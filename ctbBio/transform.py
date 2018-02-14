#!/usr/bin/env python3

"""
script for normalizing/standardizing data

based on rows, column, or table:
    1. 0-1
    2. percent of total
    3. standardize (X - mean)/(standard deviation)

example: normalize rows 0-1 | standardize rows | standardize table
"""

import sys
import os
import numpy as np

def zero_to_one(table, option):
    """
    normalize from zero to one for row or table
    """
    if option == 'table':
        m = min(min(table))
        ma = max(max(table))
    t = []
    for row in table:
        t_row = []
        if option != 'table':
            m, ma = min(row), max(row)
        for i in row:
            if ma == m:
                t_row.append(0)
            else:
                t_row.append((i - m)/(ma - m))
        t.append(t_row)
    return t

def pertotal(table, option):
    """
    calculate percent of total
    """
    if option == 'table':
        total = sum([i for line in table for i in line])
    t = []
    for row in table:
        t_row = []
        if option != 'table':
            total = sum(row)
        for i in row:
            if total == 0:
                t_row.append(0)
            else:
                t_row.append(i/total*100)
        t.append(t_row)
    return t

def standardize(table, option):
    """
    standardize
    Z = (X - mean) / (standard deviation)
    """
    if option == 'table':
        mean = np.mean(table)
        std = np.std(table)
    t = []
    for row in table:
        t_row = []
        if option != 'table':
            mean = np.mean(row)
            std = np.std(row)
        for i in row:
            if std == 0:
                t_row.append(0)
            else:
                t_row.append((i - mean)/std)
        t.append(t_row)
    return t

def scale(table):
    """
    scale table based on the column with the largest sum
    """
    t = []
    columns = [[] for i in table[0]]
    for row in table:
        for i, v in enumerate(row):
            columns[i].append(v)
    sums = [float(sum(i)) for i in columns]
    scale_to = float(max(sums))
    scale_factor = [scale_to/i for i in sums if i != 0]
    for row in table:
        t.append([a * b for a,b in zip(row, scale_factor)])
    return t

def norm(table):
    """
    fit to normal distribution
    """
    print('# norm dist is broken', file=sys.stderr)
    exit()
    from matplotlib.pyplot import hist as hist
    t = []
    for i in table:
        t.append(np.ndarray.tolist(hist(i, bins = len(i), normed = True)[0]))
    return t

def log_trans(table):
    """
    log transform each value in table
    """
    t = []
    all = [item for sublist in table for item in sublist]
    if min(all) == 0:
        scale = min([i for i in all if i != 0]) * 10e-10
    else:
        scale = 0
    for i in table:
        t.append(np.ndarray.tolist(np.log10([j + scale for j in i])))
    return t

def box_cox(table):
    """
    box-cox transform table
    """
    from scipy.stats import boxcox as bc
    t = []
    for i in table:
        if min(i) == 0:
            scale = min([j for j in i if j != 0]) * 10e-10
        else:
            scale = 0
        t.append(np.ndarray.tolist(bc(np.array([j + scale for j in i]))[0]))
    return t

def inh(table):
    """
    inverse hyperbolic sine transformation
    """
    t = []
    for i in table:
        t.append(np.ndarray.tolist(np.arcsinh(i)))
    return t

def diri(table):
    """
    from SparCC - "randomly draw from the corresponding posterior
    Dirichlet distribution with a uniform prior"
    """
    t = []
    for i in table:
        a = [j + 1 for j in i]
        t.append(np.ndarray.tolist(np.random.mtrand.dirichlet(a)))
    return t

def transpose(table):
    """
    transpose matrix
    """
    t = []
    for i in range(0, len(table[0])):
        t.append([row[i] for row in table])
    return t

def transform(table, option, mode):
    """
    transform data table:
    option:
        rows
        columns
        table
    mode:
        0: 0-1
        pertotal: percent of total
        stand: standardize
        scale: scale based on largest column
    """
    if option != 'rows' and option != 'columns' and option != 'table':
        print('# specify: rows, columns, or table', file=sys.stderr)
        exit()
    if option == 'columns':
        table = transpose(table)
    if mode == '0':
        transformed = zero_to_one(table, option)
    elif mode == 'pertotal':
        transformed = pertotal(table, option)
    elif mode == 'stand':
        transformed = standardize(table, option)
    elif mode == 'scale':
        if option != 'table':
            print('# scaling is done on an entire table, based on the column with the largest sum - use table', file=sys.stderr)
            exit()
        transformed = scale(table)
    elif mode == 'log':
        if option != 'table':
            print('# scaling is done on a per-value basis - use table', file=sys.stderr)
            exit()
        transformed = log_trans(table)
    elif mode == 'box':
        transformed = box_cox(table)
    elif mode == 'inh':
        transformed = inh(table)
    elif mode == 'norm':
        transformed = norm(table)
    elif mode == 'diri':
        transformed = diri(table)
    else:
        print('# specify: 0 (for 0-1), pertotal (for percent of total), stand (for standardize), scale (scale table to largest column), norm (fit to normal distribution), log (log10(x+1)-transform), box (box-cox), inh (inverse hyperbolic sine), diri (Dirichlet distribution)', file=sys.stderr)
        exit()
    if option == 'columns':
        transformed = transpose(transformed)
    return transformed

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('usage: cat table.tsv | transform.py <rows, columns, table> <method>\n', file=sys.stderr)
        print('methods: 0 (0-1), pertotal (percent of total), stand (standardize), scale (scale table to largest column), norm (fit to normal distribution), log (log10(x+1)-transform), box (box-cox), inh (inverse hyperbolic sine transformation>, diri (draw from posterior Dirichlet distribution)\n', file=sys.stderr)
        print('# example: cat table.tsv | transform.py rows 0 > table.trans.tsv', file=sys.stderr)
        exit()
    option, mode = sys.argv[1:]
    data = [i.strip().split('\t') for i in sys.stdin]
    if data[0][0][0] == '%':
        header, table, names = [], [], []
        for line in data:
            if line[0][0] == '%':
                header.append(line)
            else:
                names.append(line[0])
                table.append([float(i) for i in line[1::]])
    else:
        header = [data[0]]
        names = [i[0] for i in data[1::]]
        table = [[float(j) for j in i[1:]] for i in data[1:]]
    names.insert(0, header[-1][0])
    transformed = transform(table, option, mode)
    for line in header:
        print('\t'.join(line))
    for i, line in enumerate(transformed, 1):
        print('%s\t%s' % (names[i], '\t'.join([str(i) for i in line])))
