#!/usr/bin/env python

import os

if __name__ == "__main__":

    badFiles = 'cleaner/badRuns.txt'
    with open(badFiles, 'r') as f:
        badJobs = f.readlines()
    badJobs = [j.strip() for j in badJobs]

    existingFiles = []
    for badJob in badJobs:

        exeFile = 'npx4-execs/%s.sh' % badJob
        with open(exeFile, 'r') as f:
            exeLines = f.readlines()
        exeLines = [l.strip() for l in exeLines]

        exeLine = [l for l in exeLines if ' -o ' in l][0]
        exeLine = exeLine.split(' ')
        idx = exeLine.index('-o') + 1
        outFile = exeLine[idx]
        if os.path.isfile(outFile):
            existingFiles += [outFile]

    print 'The following bad files were found:'
    for f in existingFiles:
        print f

    if len(existingFiles) != 0:
        yn = raw_input('Would you like to remove them? [y|n]: ')
        if yn == 'y':
            for f in existingFiles:
                os.remove(f)
