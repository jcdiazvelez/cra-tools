#!/usr/bin/env python

import os

if __name__ == "__main__":

    badFile = 'cleaner/badRuns.txt'
    with open(badFile, 'r') as f:
        badJobs = f.readlines()
    badJobs = [job.strip() for job in badJobs]

    existingFiles = []
    for job in badJobs:

        exeFile = 'npx4-execs/%s.sh' % job
        with open(exeFile, 'r') as f:
            exeLines = f.readlines()
        myline = exeLines[-2].strip()

        a, origfile, dest = myline.split()
        basename = os.path.basename(origfile)
        outFile = '%s/%s' % (dest, basename)
        if os.path.isfile(outFile):
            existingFiles += [outFile]

    print 'The following existing files were found:'
    for efile in sorted(existingFiles):
        print efile

    if len(existingFiles) != 0:
        yn = raw_input('Do you want to delete them? [y|n]: ')
        if yn == 'y':
            for efile in existingFiles:
                os.remove(efile)
        
