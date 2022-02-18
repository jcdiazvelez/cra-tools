#!/usr/bin/env python

import os, re

from glob import glob

if __name__ == "__main__":

    # Collect all daily map files
    #prefix = '/data/user/fmcnally/anisotropy/maps'
    prefix = '/data/user/fmcnally/anisotropy/maps/test'
    masterList = []
    for cfg in ['IC','IT']:
        masterList += glob('{}/{}*/*_????-??-??.fits'.format(prefix, cfg))
    masterList = [os.path.basename(f) for f in masterList]
    masterList.sort()

    omit = ['IT73_24H_sid_STA8_6.2_{}'.format(c) for c in ['p','h','o','f']]

    configs = sorted(list(set([f.split('_')[0] for f in masterList])))
    for cfg in configs:
        print('Working on {}...'.format(cfg))
        cfgList = [f for f in masterList if cfg in f]
        dates = list(set([re.split('_|\.',f)[-2] for f in cfgList]))
        params = list(set(['_'.join(f.split('_')[:-1]) for f in cfgList]))
        params.sort()
        for p in params:
            if p in omit:
                continue
            for d in sorted(dates):
                testFile = '{}_{}.fits'.format(p, d)
                if testFile not in masterList:
                    print('{} not found!'.format(testFile))
