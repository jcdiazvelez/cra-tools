#!/usr/bin/env python

import glob

from icecube import dataio

if __name__ == "__main__":

    # Build list of files to test from executables
    exeFiles = glob.glob('npx4-execs/*')
    exeFiles.sort()

    checkList = []
    for exeFile in exeFiles:
        with open(exeFile,'r') as f:
            lines = f.readlines()
        exeLine = lines[-2].strip()
        checkList += [l for l in exeLine.split() if '/data/sim' in l]

    badFiles = []
    nfiles = len(checkList)
    print 'Looking at %s files' % nfiles

    for i, file in enumerate(checkList):

        i3file = dataio.I3File(file)
        try: 
            while i3file.pop_frame():
                continue
        except RuntimeError:
            print 'Stream error encountered'
            badFiles += ['%s\n' % file]
            continue

    print 'Finished'
    print 'Stream errors were found in the following files:'
    for badFile in badFiles:
        print badFile.strip()

    # Update stream error file
    outFile = 'streamErrors.txt'
    outLines = []
    if os.path.isfile('streamErrors.txt'):
        with open('streamErrors.txt','r') as f:
            outLines = f.readlines()
    outLines += badFiles
    outLines = sorted(list(set(outLines)))
    with open(outFile, 'w') as f:
        f.writelines(outLines)


