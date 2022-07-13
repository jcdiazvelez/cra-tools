#!/usr/bin/env python

import os, argparse
from glob import glob
from npx4.cleaner import goodFile

if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Remove map files produced by (or causing) bad runs')
    p.add_argument('--npxdir', dest='npxdir',
            default='/home/fmcnally/npx4',
            help='Location of npx4 folders')
    p.add_argument('--strictError', dest='strictError',
            default=False, action='store_true',
            help='ONLY look at lines with the word "Error" in them')
    args = p.parse_args()

    eliminateList = []
    errFiles = sorted(glob('{}/npx4-error/*.error'.format(args.npxdir)))

    for errFile in errFiles:

        if goodFile(errFile, args.strictError):
            continue

        with open(errFile, 'r') as f:
            lines = f.readlines()

        # Isolate lines with a fits file name
        ext = '.fits'
        lines = [l.strip() for l in lines if ext in l]
        if lines == []:
            continue

        badOutput = [i for l in lines for i in l.split() if ext in i]
        eliminateList += [f for f in badOutput if os.path.isfile(f)]

    eliminateList = sorted(list(set(eliminateList)))
    if len(eliminateList) != 0:

        print 'The following files will be deleted:'
        for f in eliminateList:
            print '  %s' % f
        yn = raw_input('Continue [y|n]?:  ')
        if yn == 'y':
            for f in eliminateList:
                os.remove(f)
        
