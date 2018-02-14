#!/usr/bin/env python3

"""
script for parallel execution of scripts
"""

import sys
import os
from multiprocessing import Pool as multithread
from subprocess import Popen

def run_process(process):
    """
    execute process
    """
    p = Popen(process, shell = True)
    return p.communicate()

def parallel(processes, threads):
    """
    execute jobs in processes using N threads
    """
    pool = multithread(threads)
    pool.map(run_process, processes)
    pool.close()
    pool.join()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('usage: <commands\\n> | parallel.py <threads>')
        exit()
    threads = int(sys.argv[1])
    processes = sys.stdin
    parallel(processes, threads)
