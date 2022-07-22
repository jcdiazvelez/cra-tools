#!/usr/bin/env python

import os, glob

if __name__ == "__main__":

    eliminateList = []
    badFiles = glob.glob('/home/fmcnally/npx4/npx4-error/*')
    badFiles.sort()
    outDir =  '/data/ana/CosmicRay/Anisotropy/IceCube/IC86/2014/simple-dst'

    for f in badFiles:

        basename = os.path.basename(f)
        basename = basename.replace('.error', '.root')
        basename = 'ic86_' + basename
        badOutput = '%s/%s' % (outDir, basename)

        if os.path.isfile(badOutput):
            eliminateList += [badOutput]

    if len(eliminateList) != 0:

        print 'The following files will be deleted:'
        for f in eliminateList:
            print '  %s' % f
        yn = raw_input('Continue [y|n]?:  ')
        if yn == 'y':
            for f in eliminateList:
                os.remove(f)
        
