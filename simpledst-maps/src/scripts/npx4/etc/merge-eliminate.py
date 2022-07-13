#!/usr/bin/env python

import glob, os

if __name__ == "__main__":

    badFile = 'cleaner/badRuns.txt'
    with open(badFile, 'r') as f:
        badJobs = f.readlines()
    badJobs = [job.strip() for job in badJobs]

    configs = ['IT81-III','IT81-IV']
    badFiles = []
    for config in configs:

        prefix = '/data/user/fmcnally/showerllh/%s_data/files' % config
        for jobID in badJobs:
            if config not in jobID.split('_'):
                continue
            date = jobID.split('_')[-1]
            testFile = '%s/DataLLH_%s_logdist.hdf5' % (prefix, date)
            if os.path.isfile(testFile):
                badFiles += [testFile]

    badFiles = sorted(list(set(badFiles)))
    yn = raw_input('%i bad files found. Eliminate? [y|n] ' % len(badFiles))
    if yn == 'y':
        for f in badFiles:
            os.remove(f)
