#!/usr/bin/env python

run_list = range(0,416)
run_depth = 2
start_prob = 0.01
from subprocess import call
from concurrent.futures import ThreadPoolExecutor as Pool
import os
import sys
from random import *
import time
import random
random.seed(1)
version = str(69)
if len(sys.argv) == 2: version = sys.argv[1]
print("Updating to version", version)

parallel = 8

os.system('mkdir -p store/version/')
os.system('mkdir -p store/tmp/')
os.system('mkdir -p output')


def outdated(i):
    fn = 'store/version/%d.txt'%i
    os.system('touch '+fn)
    f = open(fn, 'r')
    t = f.read()
    f.close()
    return t != version
def update(i):
    fn = 'store/version/%d.txt'%i
    os.system('touch '+fn)
    f = open(fn, 'w')
    f.write(version)
    f.close()

run_list = [i for i in run_list if outdated(i)]

global done
done = 0

start = time.time()

n = len(run_list)
def run(i):
    if call(['/usr/bin/time', './run', str(i), str(run_depth), str(start_prob)], stdout=open('store/tmp/%d_out.txt'%i,'w'), stderr=open('store/tmp/%d_err.txt'%i,'w')):
        print(i, "Crashed")
        return i
    os.system('cp store/tmp/%d_out.txt store/%d_out.txt'%(i,i))
    update(i)
    global done
    done += 1
    print("%d / %d     \r"%(done, n), end = "")
    sys.stdout.flush()
    return i

scores = []
with Pool(max_workers=parallel) as executor:
    for i in executor.map(run, run_list):
        pass

end = time.time()

print("Done in: ", round(end-start, 2), "s")
