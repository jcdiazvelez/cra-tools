#!/usr/bin/env python

import glob
import argparse

if __name__ == "__main__":

    p = argparse.ArgumentParser()
    p.add_argument('-d', '--dir', dest='npxdir',
            default='/home/fmcnally/npx4',
            help='Submitter directory to check')
    args = p.parse_args()

    errLines = []
    alldates = []

    files = sorted(glob.glob(f'{args.npxdir}/npx4-error/timescramble*.error'))
    badFiles = []

    safeErrors = ['Error in <TSystem::ExpandFileName>: input: $HOME/.root.mimes, output: $HOME/.root.mimes']

    for filename in files:

        with open(filename, 'r') as f:
            lines = f.readlines()

        for l in lines:
            if '2018' in l:
                ldata = l.split(' ')
                for l_i in ldata:
                    if '2018' in l_i:
                        l_i = l_i.replace(',','')
                        alldates += [l_i.strip()]

            if 'Error' in l and l.strip() not in safeErrors:
                badFiles += [filename]
                if l not in errLines:
                    #print(filename)
                    #print(l)
                    #print()
                    errLines += [l]

    alldates = sorted(set(alldates))
    for d in alldates:
        print(d)

    errList1, errList2 = [],[]
    for l in errLines:
        if '.root' in l and 'badread' in l:
            f = l.split('.root')[0]
            f = f.split('File: ')[1]
            errList1 += [f+'.root']
        if '.root' in l and 'TDSTTrigger' in l:
            f = l.split('.root')[0]
            f = f.split('file ')[1]
            errList2 += [f+'.root']

    print('\nBadread files')
    for f in sorted(set(errList1)):
        print(f.split('/')[-1])

    print('\nMissing TDSTTrigger')
    for f in sorted(set(errList2)):
        print(f.split('/')[-1])

    errList = errList1 + errList2
    badFiles = sorted(set(badFiles))
    for filename in badFiles:
        with open(filename, 'r') as f:
            lines = f.readlines()
        if not any([root in l for root in errList for l in lines]):
            print(filename)
            print('Uh oh')

