#!/usr/bin/env python3

"""
script for analyzing rRNA insertions predicted using either 16SfromHMM.py or 23SfromHMM.py
- predict ORFs using Prodigal
- identifying introns usiing cmscan against Rfam intron database
"""

import sys
import os
import subprocess
from fasta import iterate_fasta as parse_fasta
from itertools import cycle as cycle

def insertions_from_masked(seq):
    """
    get coordinates of insertions from insertion-masked sequence
    """
    insertions = []
    prev = True
    for i, base in enumerate(seq):
        if base.isupper() and prev is True:
            insertions.append([])
            prev = False
        elif base.islower():
            insertions[-1].append(i)
            prev = True
    return [[min(i), max(i)] for i in insertions if i != []]

def analyze_fa(fa):
    """
    analyze fa (names, insertions) and convert fasta to prodigal/cmscan safe file
    - find insertions (masked sequence)
    - make upper case
    - assign names to id number
    """
    if fa.name == '<stdin>':
        safe = 'temp.id'
    else:
        safe = '%s.id' % (fa.name)
    safe = open(safe, 'w')
    sequences = {} # sequences[id] = sequence
    insertions = {} # insertions[id] = [[start, stop], [start, stop], ...]
    count = 0
    id2name = {}
    names = []
    for seq in parse_fasta(fa):
        id = '%010d' % (count,)
        name = seq[0].split('>', 1)[1]
        id2name[id] = name
        id2name[name] = id
        names.append(name)
        insertions[id] = insertions_from_masked(seq[1])
        sequences[id] = seq
        print('\n'.join(['>%s' % (id), seq[1].upper()]), file=safe)
        count += 1
    safe.close()
    lookup = open('%s.id.lookup' % (fa.name), 'w')
    for i in list(id2name.items()):
        print('\t'.join(i), file=lookup)
    lookup.close()
    return safe.name, sequences, id2name, names, insertions

def seq_info(names, id2names, insertions, sequences):
    """
    get insertion information from header
    """
    seqs = {} # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns]], ...]]
    for name in names:
        id = id2names[name]
        gene = name.split('fromHMM::', 1)[0].rsplit(' ', 1)[1]
        model = name.split('fromHMM::', 1)[1].split('=', 1)[1].split()[0]
        i_gene_pos = insertions[id] # coordinates of each insertion wrt gene
        i_model_pos = name.split('fromHMM::', 1)[1].split('model-pos(ins-len)=')[1].split()[0].split(';') # model overlap
        i_info = []
        for i, ins in enumerate(i_gene_pos):
            model_pos = i_model_pos[i].split('-')[1].split('(')[0] 
            length = i_model_pos[i].split('(')[1].split(')')[0]
            iheader = '>%s_%s insertion::seq=%s type=insertion strand=n/a gene-pos=%s-%s model-pos=%s'\
                    % (id, (i + 1), (i + 1), ins[0], ins[1], model_pos)
            iseq = sequences[id][1][ins[0]:(ins[1] + 1)]
            iseq = [iheader, iseq]
            info = [ins, model_pos, length, iseq, [], []]
            i_info.append(info)
        seqs[id] = [gene, model, i_info]
    return seqs

def overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))

def check_overlap(pos, ins, thresh):
    """
    make sure thresh % feature is contained within insertion
    """
    ins_pos = ins[0]
    ins_len = ins[2]
    ol = overlap(ins_pos, pos)
    feat_len = pos[1] - pos[0] + 1
#    print float(ol) / float(feat_len)
    if float(ol) / float(feat_len) >= thresh:
        return True
    return False

def find_orfs(fa, seqs):
    """
    find orfs and see if they overlap with insertions
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns]], ...]]
    """
    faa = '%s.prodigal.faa' % (fa)
    fna = '%s.prodigal.fna' % (fa)
    gbk = '%s.prodigal.gbk' % (fa)
    if os.path.exists(faa) is False:
        p = subprocess.Popen('prodigal -q -i %s -a %s -d %s -c -f gbk -m -n -o %s -p meta' \
                % (fa, faa, fna, gbk), shell = True)
        p.communicate()
    for orf in parse_fasta(faa):
        if orf[0] == []:
            continue
        id = orf[0].split('>')[1].split('_', 1)[0]
        pos = sorted([int(i) for i in orf[0].split()[2:5] if i != '#'])
        if id not in seqs:
            continue
        for i, ins in enumerate(seqs[id][2]):
            if check_overlap(pos, ins, 0.90) is True:
                seqs[id][2][i][4].append(orf)
    return seqs, faa

def find_introns(fa, seqs, sequences, threads):
    """
    find introns by searching Rfam intron databse using cmscan
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns]], ...]]
    """
    db = '%s/rfam/Rfam.cm.1_1.introns' % (os.environ['databases'])
    out = '%s.rfam-introns.cmscan' % (fa)
    tblout = '%s.rfam-introns.cmscan-tblout' % (fa)
    if os.path.exists(tblout) is False:
        p = subprocess.Popen('cmscan --cpu %s --tblout %s %s %s > %s'\
                % (threads, tblout, db, fa, out), shell = True)
        p.communicate()
    for line in open(tblout):
        if line.startswith('#'):
            continue
        line = line.strip().split()
        if line[16] == '?': # does not pass inclusion threshold
            continue
        id = line[2]
        type, start, stop, strand = line[0], int(line[7]), int(line[8]), line[9]
        if 'intron' not in type.lower():
            continue
        pos = sorted([start, stop])
        if id not in seqs:
            continue
        for i, ins in enumerate(seqs[id][2]):
            if check_overlap(pos, ins, 0.25) is True:
                seqs[id][2][i][5].append(['>%s_%s %s %s %s-%s' % (id, (i + 1), type, strand, start, stop), sequences[id][1][pos[0]-1:pos[1]]])
    return seqs

def annotate_orfs(fa, seqs, threads, threshold = 1e-6):
    """
    search orfs against pfam
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns], orfs?, introns?, {orf: annotations}], ...]]
    """
    db = '%s/pfam/Pfam-A.hmm' % os.environ['databases']
    base = fa.rsplit('.', 1)[0]
    out = '%s.pfam.hmmscan' % (base)
    tblout = '%s.pfam.hmmscan-tblout' % (base) 
    orf2pfam = {} # orf2pfam[scaffold id][orf id] = pfam
    if os.path.exists(tblout) is False:
        p = subprocess.Popen('hmmscan --cpu %s --tblout %s -o %s %s %s'\
                % (threads, tblout, out, db, fa), shell = True)
        p.communicate()
    for line in open(tblout):
        if line.startswith('#'):
            continue
        line = line.strip().split()
        hit, id, e, inc = line[0], line[2], float(line[4]), int(line[17])
        contig = id.split('_', 1)[0]
        if e > threshold:
            continue
        if contig not in orf2pfam:
            orf2pfam[contig] = {}
        if id not in orf2pfam[contig]:
            orf2pfam[contig][id] = hit
    for seq, info in list(seqs.items()):
        for insertion in info[2]:
            annotations = {}
            for orf in insertion[4]:
                orf  = orf[0].split('>')[1].split()[0]
                contig = orf.split('_', 1)[0]
                if contig in orf2pfam and orf in orf2pfam[contig]:
                    annotations[orf] = orf2pfam[contig][orf]
                else:
                    annotations[orf] = 'n/a'
            insertion.append(annotations)
    return seqs

def analyze_insertions(fa, threads = 6):
    """
    - find ORFs using Prodigal
    - find introns using cmscan (vs Rfam intron database)
    - check that ORFs and intron overlap with insertion region
    - plot insertion length versus model position for each insertion (based on insertion type)
    """
    safe, sequences, id2name, names, insertions = analyze_fa(fa)
    seqs = seq_info(names, id2name, insertions, sequences)
    seqs, orfs = find_orfs(safe, seqs)
    seqs = find_introns(safe, seqs, sequences, threads)
    seqs = seqs2bool(seqs)
    seqs = annotate_orfs(orfs, seqs, threads)
    return seqs, id2name

def seqs2bool(seqs):
    """
    convert orf and intron information to boolean
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns]], ...]]
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns], orfs?, introns?], ...]]
    """
    for seq in seqs:
        for i, ins in enumerate(seqs[seq][2]):
            if len(ins[4]) > 0:
                ins.append(True)
            else:
                ins.append(False)
            if len(ins[5]) > 0:
                ins.append(True)
            else:
                ins.append(False)
            seqs[seq][2][i] = ins
    return seqs

def print_table(seqs, id2name, name):
    """
    print table of results
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns], orfs?, introns?], ...]]
    """
    itable = open('%s.itable' % (name.rsplit('.', 1)[0]), 'w')
    print('\t'.join(['#sequence', 'gene', 'model', 'insertion', 'gene position', 'model position', 'length', 'orf?', 'intron?', 'orf?intron?', 'insertion', 'orf', 'intron']), file=itable)
    for seq, info in list(seqs.items()):
        gene, model, insertions = info
        name = id2name[seq]
        for i, ins in enumerate(insertions, 1):
            gene_pos, model_pos, length, iseq, \
                    orfs, introns, orfs_b, introns_b, orf_annotations = ins
            # append annotation to orf header
            for orf in orfs:
                parts = orf[0].split()
                annotation = orf_annotations[parts[0].split('>')[1]]
                orf[0] = '%s %s %s' % (parts[0], annotation, ' '.join(parts[1:]))
            # get orf position
            gene_pos = '-'.join([str(j) for j in gene_pos])
            # check if orf, intron is present
            if orfs_b is True or introns_b is True:
                orfs_introns_b = True
            else:
                orfs_introns_b = False
            out = [name, gene, model, i, gene_pos, model_pos, length, orfs_b, introns_b, orfs_introns_b]
            out.append('|'.join(iseq))
            out.append('|'.join(['|'.join(orf) for orf in orfs]))
            out.append('|'.join(['|'.join(intron) for intron in introns]))
            print('\t'.join([str(i) for i in out]), file=itable)
    itable.close()

def print_seqs(seqs, id2name, name):
    """
    print fasta of introns and ORFs
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, iseq, [orfs], [introns], orfs?, introns?, [orf annotations]], ...]]
    """
    orfs = open('%s.orfs.faa' % (name.rsplit('.', 1)[0]), 'w')
    introns = open('%s.introns.fa' % (name.rsplit('.', 1)[0]), 'w')
    insertions = open('%s.insertions.fa' % (name.rsplit('.', 1)[0]), 'w')
    for seq in seqs:
        for i, ins in enumerate(seqs[seq][2], 1):
            model_pos = ins[1]
            if ins[6] is True: # orf(s) in ins[4]
                for orf in ins[4]:
                    orf_info = orf[0].split('>')[1].split()
                    id = orf_info[0].split('_', 1)[0]
                    name = id2name[id]
                    annotation = orf_info[1]
                    strand = orf_info[7]
                    if strand == '1':
                        strand = '+'
                    else:
                        strand = '-'
                    start, stop = sorted([int(orf_info[3]), int(orf_info[5])])
                    header = '>%s insertion::seq=%s type=%s strand=%s gene-pos=%s-%s model-pos=%s'\
                            % (name, i, annotation, strand, start, stop, model_pos)
                    print('\n'.join([header, orf[1]]), file=orfs)
            if ins[7] is True: # intron(s) in ins[5]
                for intron in ins[5]:
                    id, type, strand, pos = intron[0].split('>', 1)[1].split()
                    name = id2name[id.split('_')[0]]
                    header = '>%s insertion::seq=%s type=%s strand=%s gene-pos=%s model-pos=%s'\
                            % (name, i, type, strand, pos, model_pos)
                    print('\n'.join([header, intron[1]]), file=introns)
            insertion = ins[3]
            id, info = insertion[0].split('>')[1].split(' ', 1)
            name = id2name[id.split('_')[0]]
            header = '>%s %s' % (name, info)
            print('\n'.join([header, insertion[1]]), file=insertions)
    insertions.close()        
    orfs.close()
    introns.close()

def max_insertion(seqs, gene, domain):
    """
    length of largest insertion
    """
    seqs = [i[2] for i in list(seqs.values()) if i[2] != [] and i[0] == gene and i[1] == domain]
    lengths = []
    for seq in seqs:
        for ins in seq:
            lengths.append(int(ins[2]))
    if lengths == []:
        return 100 
    return max(lengths)

def model_length(gene, domain):
    """
    get length of model
    """
    if gene == '16S':
        domain2max = {'E_coli_K12': int(1538), 'bacteria': int(1689), 'archaea': int(1563), 'eukarya': int(2652)}
        return domain2max[domain]
    elif gene == '23S':
        domain2max = {'E_coli_K12': int(2903), 'bacteria': int(3146), 'archaea': int(3774), 'eukarya': int(9079)}
        return domain2max[domain]
    else:
        print(sys.stderr, '# length unknown for gene: %s, domain: %s' % (gene, domain))
        exit()

def setup_markers(seqs):
    """
    setup unique marker for every orf annotation
    - change size if necessary
    """
    family2marker = {} # family2marker[family] = [marker, size]
    markers = cycle(['^', 'p', '*', '+', 'x', 'd', '|', 'v', '>', '<', '8'])
    size = 60
    families = []
    for seq in list(seqs.values()):
        for insertion in seq[2]:
            for family in list(insertion[-1].values()):
                if family not in families:
                    families.append(family)
    for family in families:
        marker = next(markers) 
        if marker == '^':
            size = size * 0.5
        family2marker[family] = [marker, size]
    return family2marker

def plot_insertions(fname, seqs, gene, domain, tax, id2name):
    """
    plot insertions wrt model positions
    2 panels:
        1. insertions with ORF
            - triangle if has intron
            - circle if no intron
        2. insertion does not have ORF
            - triangle if has intron
            - circle if no intron
        *3. conservation of sequence in model
        *4. e. coli positions
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, orf, intron], ...]]
    """
    # import
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.cm as cm
    import matplotlib.colors as col
    from matplotlib.font_manager import FontProperties
    plt.rcParams['pdf.fonttype'] = 42 # illustrator pdf
    # axis
    max_insert = max_insertion(seqs, gene, domain)
    height = max_insert + max_insert * 0.10
    xmax = model_length(gene, domain)
    f, axarr = plt.subplots(4, sharex = True, sharey = True)
    plt.axis([0, xmax, 0, height])
    plt.xticks(np.arange(0, xmax, 100), rotation = 45)
    plt.yticks(np.arange(0, height, 200))
    # labels
    axarr[0].set_title('encodes ORF and intron')
    axarr[1].set_title('encodes ORF, no intron')
    axarr[2].set_title('encodes intron, no ORF')
    axarr[3].set_title('no intron, no ORF')
    plt.suptitle('%s %s rRNA gene insertions' % (domain, gene))
    plt.ylabel('insertion length (bp)')
    plt.xlabel('position on %s %s rRNA gene model' % (domain, gene))
    # colors
    color2tax = {}
    if tax is False:
        taxa = ['n/a']
    else:
        taxa = sorted(set([tax[id2name[i]] for i, j in list(seqs.items()) if j[2] != [] and id2name[i] in tax]))
        if 'n/a' not in taxa:
            taxa.append('n/a')
    colors = cm.spectral(np.linspace(0, 1, len(taxa)))
    colors = cycle(colors)
    for t in taxa:
        color2tax[t] = next(colors)
    # markers
    markers = setup_markers(seqs)
    # plot
    for name, seq in list(seqs.items()):
        g, d = seq[0], seq[1]
        if g != gene or d != domain or seq[2] == []:
            continue
        if tax is False or id2name[name] not in tax:
            t = 'n/a'
        else:
            t = tax[id2name[name]]
        c = color2tax[t] 
        for ins in seq[2]:
            family = [i for i in list(ins[-1].values()) if i != 'n/a']
            if len(family) != 1:
                family = 'n/a'
            else:
                family = family[0]
            x, y = int(ins[1]), int(ins[2])
            orf, intron = ins[-3], ins[-2]
            if orf is True: # has orf
                if intron is True: # has intron
                    p = 0
                else:
                    p = 1
            else:
                if intron is True: # intron, no orf
                    p = 2
                else:
                    p = 3
            marker, size = 'o', 30
            if orf is True:
                marker, size = markers[family]
            axarr[p].scatter(x, y, \
                    edgecolors = c, marker = marker, s = size, label = family, \
                    facecolors = 'none', clip_on = False)
    # legend
    handles, labels = [], []
    for ax in axarr[0:2]:
        hs, ls = ax.get_legend_handles_labels()
        for h, l in zip(hs, ls):
            if l in labels:
                continue
            handles.append(h)
            labels.append(l)
    l1 = plt.legend(handles, labels, scatterpoints = 1, \
            prop = {'size':10}, loc = 'upper left', bbox_to_anchor = (1, 0.5))
    names = [t for t in taxa]
    boxes = [matplotlib.patches.Rectangle((0, 0), 1, 1, fc = color2tax[t]) for t in taxa]
    plt.legend(boxes, names, scatterpoints = 1, \
            prop = {'size':10}, loc = 'lower left', bbox_to_anchor = (1, 0.5))
    plt.gca().add_artist(l1) # add l1 as a separate legend
    # save
#    plt.tight_layout()
    figure = plt.gcf()
    figure.set_size_inches(12, 12)
    pdf = PdfPages('%s.%s-%srRNAgene-insertions.pdf' % (fname.rsplit('.', 1)[0], domain, gene))
    pdf.savefig()
    plt.close()
    pdf.close()
#    plt.show()


def plot_insertions_two_panels(fname, seqs, gene, domain, tax, id2name):
    """
    plot insertions wrt model positions
    2 panels:
        1. insertions with ORF
            - triangle if has intron
            - circle if no intron
        2. insertion does not have ORF
            - triangle if has intron
            - circle if no intron
        *3. conservation of sequence in model
        *4. e. coli positions
    # seqs[id] = [gene, model, [[i-gene_pos, i-model_pos, i-length, orf, intron], ...]]
    """
    # import
    import matplotlib
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.cm as cm
    import matplotlib.colors as col
    from itertools import cycle as cycle
    from matplotlib.font_manager import FontProperties
    plt.rcParams['pdf.fonttype'] = 42 # illustrator pdf
    # axis
    max_insert = max_insertion(seqs, gene, domain)
    height = max_insert + max_insert * 0.10
    xmax = model_length(gene, domain)
    f, axarr = plt.subplots(2, sharex = True, sharey = True)
    plt.axis([0, xmax, 0, height])
    plt.xticks(np.arange(0, xmax, 100), rotation = 45)
    plt.yticks(np.arange(0, height, 100))
    # labels
    axarr[0].set_title('encodes ORF')
    axarr[1].set_title('does not encode ORF')
    plt.suptitle('%s %s rRNA gene insertions' % (domain, gene))
    plt.ylabel('insertion length (bp)')
    plt.xlabel('position on %s %s rRNA gene model' % (domain, gene))
    # colors
    color2tax = {}
    if tax is False:
        taxa = ['n/a']
    else:
        taxa = sorted(set([tax[id2name[i]] for i, j in list(seqs.items()) if j[2] != [] and id2name[i] in tax]))
        if 'n/a' not in taxa:
            taxa.append('n/a')
    colors = cm.spectral(np.linspace(0, 1, len(taxa)))
    colors = cycle(colors)
    for t in taxa:
        color2tax[t] = next(colors)
    # plot
    for name, seq in list(seqs.items()):
        g, d = seq[0], seq[1]
        if g != gene or d != domain or seq[2] == []:
            continue
        if tax is False or id2name[name] not in tax:
            t = 'n/a'
        else:
            t = tax[id2name[name]]
        c = color2tax[t] 
        for ins in seq[2]:
            x, y = int(ins[1]), int(ins[2])
            if ins[4] == True: # has intron, set marker
                marker, size = '^', 30
            else:
                marker, size = 'o', 30
            if ins[3] == True: # has orf, plot separately
                axarr[0].scatter(x, y, marker = marker, s = size, facecolors = 'none', \
                        clip_on = False, edgecolors = c, label = t)
            else:
                axarr[1].scatter(x, y, marker = marker, s = size, facecolors = 'none', \
                        clip_on = False, edgecolors = c, label = t)
    # legend
    boxes = [matplotlib.patches.Rectangle((0, 0), 1, 1, fc = color2tax[t]) for t in taxa]
    names = [t for t in taxa]
    plt.legend(boxes, names, prop = {'size':10}, loc='center left', bbox_to_anchor=(1, 0.5), scatterpoints = 1)
    # save
    figure = plt.gcf()
    figure.set_size_inches(20, 12)
    pdf = PdfPages('%s.%s-%srRNAgene-insertions.pdf' % (fname.rsplit('.')[0], domain, gene))
    pdf.savefig()
    plt.close()
    pdf.close()

def plot_by_gene_and_domain(name, seqs, tax, id2name):
    """
    plot insertions for each gene and domain
    """
    for gene in set([seq[0] for seq in list(seqs.values())]):
        for domain in set([seq[1] for seq in list(seqs.values())]):
            plot_insertions(name, seqs, gene, domain, tax, id2name)

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: rRNA_insertions.py <[16,23]SfromHMM-out.masked.fa> <tax-mapping.tsv or False>')
        exit()
    fa, tax = sys.argv[1:]
    if fa == '-':
        fa = sys.stdin
    else:
        fa = open(fa)
    if tax == 'False' or tax == 'false' or tax == 'FALSE':
        tax = False
    else:
        tax = {i.split('\t')[0]: i.split('\t')[1].strip() for i in open(tax) if len(i.split('\t')) > 1}
    seqs, id2name = analyze_insertions(fa)
    print_table(seqs, id2name, fa.name)
    print_seqs(seqs, id2name, fa.name)
    plot_by_gene_and_domain(fa.name, seqs, tax, id2name)
