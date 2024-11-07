#!/usr/bin/env python

from subprocess import *
from concurrent.futures import ThreadPoolExecutor as Pool
import os
import sys
import resource
import psutil
from random import *
import time
import math

from os import system
from glob import glob



SUCCESS, TLE, MLE, RTE, RUNNING = 0,1,2,3,-1
exit_names = ["SUCCESS", "TLE", "MLE", "RTE", "RUNNING"]

start_time = time.time()


MEMORY_LIMIT = 8*4096 * 0.95 # MB
TIME_LIMIT = 12*60*60 * 0.96 # Seconds


class Process:
    def __init__(self, cmd, timeout, maxmemory):
        fn = cmd.replace(' ','_')
        self.fout = open('store/tmp/%s.out'%fn,'w')
        self.ferr = open('store/tmp/%s.err'%fn,'w')
        # print(cmd)
        sys.stdout.flush()
        self.cmd = cmd
        self.process = Popen(cmd.split(), stdout=self.fout, stderr=self.ferr, shell=False)
        self.pid = self.process.pid
        self.mp = psutil.Process(self.pid)
        self.memused, self.timeused = 0, 0
        self.start_time = time.time()
        self.timeout = timeout
        self.maxmemory = maxmemory

    def update(self):
        self.memory = self.mp.memory_info().rss/2**20
        self.memused = max(self.memused, self.memory)
        self.timeused = time.time()-self.start_time
        if self.memory > self.maxmemory:
            return (MLE, self.timeused, self.memused)
        if self.timeused > self.timeout:
            return (TLE, self.timeused, self.memused)
        if not self.memory:
            if self.process.wait():
                return (RTE, self.timeused, self.memused)
            else:
                return (SUCCESS, self.timeused, self.memused)
        return (RUNNING, self.timeused, self.memused)

    def __del__(self):
        self.fout.close()
        self.ferr.close()


class Command:
    def __init__(self, cmd, expected_time = TIME_LIMIT, expected_memory = MEMORY_LIMIT, slack = 1.5):
        self.cmd = cmd
        self.time = expected_time
        self.mem = expected_memory
        self.slack = slack

    def __lt__(self, other):
        return self.time < other.time


def runAll(cmd_list, threads):
    THREAD_LIMIT = threads

    ret_stats = {}

    cmd_list = sorted(cmd_list)

    dt = 0.1
    running = []
    cmdi = 0

    def callback(process, status, timeused, memused):
        #assert(status != RTE)
        print(exit_names[status], process.cmd, " %.1fs"%timeused, "%.0fMB"%memused)
        sys.stdout.flush()

        ret_stats[process.cmd] = (status, timeused, memused)

    while len(running) or cmdi < len(cmd_list):
        while cmdi < len(cmd_list) and len(running) < THREAD_LIMIT:
            cmd = cmd_list[cmdi]
            process = Process(cmd.cmd, cmd.time*cmd.slack, cmd.mem*cmd.slack)
            running.append(process)
            cmdi += 1

        torem = []
        mems = []
        for r in running:
            status, timeused, memused = r.update()
            mems.append(r.memory)
            if status != RUNNING:
                callback(r, status, timeused, memused)
                torem.append(r)

        if sum(mems) > MEMORY_LIMIT:
            r = running[mems.index(max(mems))]
            r.process.kill()
            callback(r, MLE, r.timeused, r.memused)
            torem.append(r)
            #THREAD_LIMIT = 1

        for r in torem:
            running.remove(r)

        time.sleep(dt)

        if time.time()-start_time >= TIME_LIMIT:
            for r in running:
                r.process.kill()
                callback(r, TLE, r.timeused, r.memused)
            break

    return ret_stats
    
def has_two_high_score_solutions(taski, threshold=2):
    cands = []
    for fn in glob("/kaggle/working/absres-c-files/output/answer_%d_*.csv" % taski):
        if not os.path.isfile(fn):
            continue  # Skip if file does not exist
        with open(fn, 'r') as f:
            t = f.read().strip().split('\n')
        for cand in t[1:]:
            img, score = cand.rsplit(' ', 1)
            cands.append((float(score), img))

    # Filter candidates with scores above the threshold
    high_score_cands = [img for score, img in cands if score > threshold]
    # Get unique solutions
    unique_imgs = set(high_score_cands)
    # Check if there are two or more unique solutions
    return len(unique_imgs) >= 2

def filter_tasks_by_score(task_indices, threshold=2):
    high_score_tasks = []
    filtered_task_indices = []

    for i in task_indices:
        if has_two_high_score_solutions(i, threshold=threshold):
            high_score_tasks.append(i)
        else:
            filtered_task_indices.append(i)

    if len(high_score_tasks) > 0:
        print("High scores for:", high_score_tasks)

    return filtered_task_indices, high_score_tasks



system("mkdir -p output")
system("mkdir -p store/tmp")
system('make -s gen_profile')


if len(sys.argv) == 3:
    l = int(sys.argv[1])
    n = int(sys.argv[2])
    ntasks = n
    task_list = range(l, l+n)
else:
    ntasks = int(check_output('./count_tasks'))
    task_list = range(0, ntasks)
    #print("Usage: python %s <start_task> <#tasks>"%sys.argv[0])

run_list = []
task_indices = list(range(ntasks))

for i in task_indices:
    run_list.append(Command("./run %d 1"%i))
runAll(run_list, 12)

system('make -s use_profile')

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 21"%i))
runAll(run_list, 12)

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 31"%i))
runAll(run_list, 12)

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 2"%i))
runAll(run_list, 12)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 22"%i))
runAll(run_list, 12)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 32"%i))
runAll(run_list, 12)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 3"%i))
runAll(run_list, 4)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 23"%i))
runAll(run_list, 4)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 33"%i))
runAll(run_list, 4)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 4 0.05"%i))
runAll(run_list, 4)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 24 0.05"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 34 0.05"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 5 0.025"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 25 0.025"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

depth34 = []
for i in task_indices:
    run_list.append(Command("./run %d 35 0.025"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 4 0.05"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 24 0.05"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 34 0.05"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 5 0.025"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

run_list = []
for i in task_indices:
    run_list.append(Command("./run %d 25 0.025"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

depth34 = []
for i in task_indices:
    run_list.append(Command("./run %d 35 0.025"%i))
runAll(run_list, 2)

filtered_task_indices, high_score_tasks = filter_tasks_by_score(task_indices, threshold=3)
task_indices = filtered_task_indices

def read(fn):
    f = open(fn)
    t = f.read()
    f.close()
    return t

combined = ["output_id,output"]
for taski in task_list:
    ids = set()
    cands = []
    for fn in glob("output/answer_%d_*.csv"%taski):
        t = read(fn).strip().split('\n')
        ids.add(t[0])
        for cand in t[1:]:
            img, score = cand.split()
            cands.append((float(score), img))

    assert(len(ids) == 1)
    id = ids.pop()

    cands.sort(reverse=True)
    best = []
    for cand in cands:
        score, img = cand
        if not img in best:
            best.append(img)
            if len(best) == 3:
                break
    if not best: best.append('|0|')
    combined.append(id+','+' '.join(best))

outf = open('submission_part2.csv', 'w')
for line in combined:
    print(line, file=outf)
outf.close()
# TODO Make json directly here for submission
outf = open('/kaggle/working/submission_part2.csv', 'w')
for line in combined:
    print(line, file=outf)
outf.close()