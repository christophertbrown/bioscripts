#!/usr/bin/env python3

import sys
import os
import numpy

def stats(nums):
    stat = {}
    stat['lines'] = len(nums)
    stat['mean'] = numpy.average(nums)
    stat['median'] = numpy.median(nums)
    stat['variance'] = numpy.var(nums)
    stat['std dev'] = numpy.std(nums)
    stat['sum'] = sum(nums)
    stat['min'] = numpy.amin(nums)
    stat['max'] = numpy.amax(nums)
    return stat

if __name__ == '__main__':
    if len(sys.argv) == 1:
        nums = [float(i) for i in sys.stdin]
    elif len(sys.argv) == 2:
        nums = [float(i) for i in open(sys.argv[1])]
    else:
        print('specify file or use stdin')
        exit()
    if len(nums) > 0:
        stats = stats(nums)
    else:
        print('no nums')
        exit()
    for stat in ['lines', 'mean', 'median', \
            'variance', 'std dev', \
            'sum', 'min', 'max']:
        print('%s:\t%s' % (stat, format(stats[stat])))
