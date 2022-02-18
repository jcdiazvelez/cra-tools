#!/usr/bin/env python

import os

if __name__ == "__main__":

    with open('badDays.txt', 'r') as f:
        badFiles = f.readlines()

    badFiles = [f.strip() for f in badFiles if '.root' in f]
    badDates = sorted(set([f.split('_')[2] for f in badFiles]))
    print(badDates)

