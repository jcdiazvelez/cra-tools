#!/usr/bin/env python

import os

if __name__ == "__main__":

    configs = ['IC86-{}'.format(yy) for yy in range(2011,2021)]
    methods = ['sid', 'anti', 'ext', 'solar']

    ex = './public_maker.py -c {} -m {}'

    for c in configs:
        for m in methods:
            print(f'\nWorking on {c} {m}')
            os.system(ex.format(c, m))

        print(f'\nWorking on {c} ebins')
        os.system(ex.format(c, 'sid') + ' --ebins')
