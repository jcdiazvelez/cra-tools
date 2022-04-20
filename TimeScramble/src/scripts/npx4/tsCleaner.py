#!/usr/bin/env python

import glob
import os

if __name__ == "__main__":

    files = glob.glob('npx4-execs/*.sh')
    files = sorted([f for f in files if 'timescramble' in f])
    
    dates = [f.split('_')[-1] for f in files]
    dates = [d.split('.')[0] for d in dates]

    print(dates)
