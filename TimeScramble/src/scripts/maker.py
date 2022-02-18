#!/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_7_x86_64/bin/python

import os, glob, re, argparse
import datetime as dt
import random

# npx4 submission from /home/fmcnally/npx4/pysubmit.py
from npx4.pysubmit import pysubmit

def getDSTfiles(config):

    it = 'IceCube' if config[:2]=='IC' else 'IceTop'
    fpath = '/data/ana/CosmicRay/Anisotropy/{}/{}'.format(it, config)
    fileList = glob.glob('{}/*.root'.format(fpath))
    # Newer years have part files nested in daily folders
    if fileList == []:
        fileList = glob.glob('{}/*/*.root'.format(fpath))
    return sorted(fileList)

## Special cases for binning in energy, s125, and composition. Calling the flag
## without values will default to the listed values.
class defaultValues(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        bins = {}
        bins['ebins'] = [4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 6, 6.5, 100]
        bins['sbins'] = [-2, -1, -0.75, -0.5, -0.25, 1, 4]
        bins['cbins'] = ['p','h','o','f']
        if not values:
            setattr(namespace, self.dest, bins[self.dest])
        else:
            setattr(namespace, self.dest, values)


if __name__ == "__main__":

    p = argparse.ArgumentParser(
            description='Makes daily healpix maps using time-scrambling code')

    # General parameters
    p.add_argument('-c', '--config', dest='config',
            help='Detector configuration [IC59,IT73,IC86-2013...]')
    p.add_argument('-d', '--dates', dest='dates', type=str,
            default=None, nargs='*',
            help='Dates to process [yyyymmdd] (optional, inclusive)')
    p.add_argument('-i', '--inclusive', dest='inclusive',
            default=False, action='store_true',
            help='Treat two dates given as inclusive. Otherwise treat as list')
    p.add_argument('-t', '--nInt', dest='nInt', type=int,
            default=24,
            help='Integration time (in hours)')
    p.add_argument('-m', '--method', dest='method',
            default='sid',
            help='Time transformation method [sid|anti|ext|solar]')
    p.add_argument('--sd', dest='sd',
            default=False, action='store_true',
            help='Correct for solar dipole on event-by-event basis')

    # IceTop specific
    p.add_argument('-f', '--filter', dest='filter',
            help='IceTop filter to apply [STA3|STA8|NotSTA8]')
    p.add_argument('--sbins', dest='sbins',
            default=False, action=defaultValues, nargs='*',
            help='Option to do all s125 bins')
    #p.add_argument('--comp', dest='comp',
    #        default=None, action=defaultValues, nargs='*',
    #        help='Option to do all composition maps')
    p.add_argument('--emin', dest='emin',
            help='Optional minimum reconstructed energy value')

    # Energy cut options
    p.add_argument('--ebins', dest='ebins',
            default=False, action=defaultValues, nargs='*',
            help='Option to do all energy bins')

    # Additional options
    p.add_argument('-o', '--outDir', dest='outDir',
            default='/data/user/fmcnally/anisotropy/maps',
            help='Destination directory for output files')
    p.add_argument('--test', dest='test',
            default=False, action='store_true',
            help='Option for running off cluster to test')
    p.add_argument('--batchFile', dest='batchFile',
            help='Name for output file containing file lists')
    p.add_argument('--overwrite', dest='overwrite',
            default=False, action='store_true',
            help='Option to overwrite existing map files')

    args = p.parse_args()
    c_opts = vars(args).copy()
    if not args.config:
        p.error('Detector configuration not given')

    # Defaults for submission batch size
    ## NOTE: currently disabled, supports submission of one day at a time
    #defaultBatch = 1 if args.config[:2]=='IC' else 10
    #if args.test and not args.nBatch:
    #    args.nBatch = 1
    #nBatch = defaultBatch if not args.nBatch else args.nBatch

    # Default spline files
    prefix = '/data/user/fmcnally/anisotropy/sim'
    # Removed IC59/IC79 functionality, only works with IC86 for 2020 paper
    # IC59: 4046, IC79: 6451, IC86 (old): 10649
    if args.ebins and args.config[:2]=='IC':
        sp_file = '{}/IC86_20904_hist_spline.fits'.format(prefix)
        c_opts['spline'] = sp_file

    # Base name for outfiles
    outBase = '{outDir}/{config}/{config}_{nInt:02d}H_{method}'.format(**c_opts)
    if args.sd:
        outBase += '_sd'
    if args.filter:
        outBase += '_{}'.format(args.filter)
    if args.emin:
        outBase += '_emin'
    c_opts['outBase'] = outBase

    # Base name for batch file
    batchPrefix = '/home/fmcnally/anisotropy/timeScramble_10yr/tempFiles'
    if not args.batchFile:
        c_opts['batchFile'] = '%s/%s.txt' % \
                (batchPrefix, os.path.basename(outBase))

    # Create output directories if they don't already exist
    if not os.path.isdir('{outDir}/{config}'.format(**c_opts)):
        print('Creating output directory {outDir}/{config}'.format(**c_opts))
        os.makedirs('{outDir}/{config}'.format(**c_opts))

    # Create list of desired dates
    runDates = []
    if args.dates != None and not args.inclusive:
        runDates = ['{}-{}-{}'.format(d[:4],d[4:6],d[6:8]) for d in args.dates]
    if args.dates != None and args.inclusive:
        s = args.dates[0]
        e = args.dates[-1]
        sdate = dt.date(int(s[:4]), int(s[4:6]), int(s[6:8]))
        edate = dt.date(int(e[:4]), int(e[4:6]), int(e[6:8]))
        delta = edate - sdate
        runDates = [sdate + dt.timedelta(days=i) for i in range(delta.days+1)]
        runDates = [date.strftime('%Y-%m-%d') for date in runDates]

    # Collect input files
    fileList = []
    masterList = getDSTfiles(args.config)
    for rootFile in masterList:

        # Filter by desired dates
        date = re.findall('\d{4}-\d{2}-\d{2}', rootFile)[-1]
        if (args.dates != None) and (date not in runDates):
            continue

        # Apply naming conventions to find existing files
        testFiles = [outBase]
        if args.ebins:
            testFiles = ['%s_%s-%sGeV' % (f, args.ebins[i], args.ebins[i+1]) \
                    for i in range(len(args.ebins)-1) for f in testFiles]
        if args.sbins:
            testFiles = ['%s_%sto%ss125' % (f, args.sbins[i], args.sbins[i+1]) \
                    for i in range(len(args.sbins)-1) for f in testFiles]
        #if args.comp:
        #    testFiles = ['%s_%s' % (f, comp) \
        #            for comp in args.comp for f in testFiles]
        testFiles = ['%s_%s.fits' % (f, date) for f in testFiles]

        # Overwrite or omit existing files
        # Files for bins are all created at once - only omit if all exist
        if all([os.path.isfile(f) for f in testFiles]) and not args.overwrite:
            continue
        for testFile in testFiles:
            if os.path.isfile(testFile):
                os.remove(testFile)

        fileList.append(rootFile)

    # Split into days for submission
    dateList = [re.findall('\d{4}-\d{2}-\d{2}', f)[-1] for f in fileList]
    dateList = sorted(list(set(dateList)))  # Limit to unique values
    prev_dates, next_dates = [], []
    for date_str in dateList:
        yyyy, mm, dd = [int(i) for i in date_str.split('-')]
        date = dt.datetime(yyyy, mm, dd)
        prev_dates += [(date + dt.timedelta(days=-1)).strftime('%Y-%m-%d')]
        next_dates += [(date + dt.timedelta(days=1)).strftime('%Y-%m-%d')]

    # Collect surrounding files to allow runs from midnight to midnight
    # ("Daily" files grouped by run, not actually 24 hours)
    dayFiles = []
    for day, prev, post in zip(dateList, prev_dates, next_dates):
        dayFiles += [[f for f in masterList if (day in f)
                                              or (prev in f)
                                              or (post in f)]]

    # Limit files to only runs immediately flanking desired day
    ## NOTE: this should be tested to make sure we still get 24 hours
    for i, (day, prev, post) in enumerate(zip(dayFiles,prev_dates,next_dates)):
        prev_idx = 0
        post_idx = len(day)
        if any([prev in f for f in day]):
            prev_idx = max(j for j, f in enumerate(day) if prev in f)
        if any([post in f for f in day]):
            post_idx = min(j for j, f in enumerate(day) if post in f) + 1
        dayFiles[i] = day[prev_idx:post_idx]

    dayFiles = [' '.join(oneDay)+'\n' for oneDay in dayFiles]
    if args.test:
        dayFiles = dayFiles[:1]

    with open(c_opts['batchFile'], 'w') as f:
        f.writelines(dayFiles)

    # Set up parameters to feed C++ script
    print('Parameters for submission:')
    validArgs  = ['config','nInt','method','filter','outBase','emin']
    validArgs += ['spline','batchFile','sd']
    c_opts = {k:v for k,v in c_opts.items() if v!=None and k in validArgs}
    for key in sorted(c_opts.keys()):
        print('  --%s %s' % (key, c_opts[key]))
    c_opts = [['--'+key, c_opts[key]] for key in sorted(c_opts.keys())]
    # Python magic: flatten list of arbitrary depth into list of strings
    flatList = lambda *n: (str(e) for a in n \
            for e in (flatList(*a) if isinstance(a, (tuple, list)) else (a,)))
    c_opts = ' '.join(list(flatList(c_opts)))
    # Ebins/Sbins need to go last (multitoken parameters)
    if args.ebins:
        print('  --ebins ' + ' '.join([str(i) for i in args.ebins]))
        c_opts += ' --ebins %s' % ' '.join([str(i) for i in args.ebins])
    if args.sbins:
        print('  --sbins ' + ' '.join(args.sbins))
        c_opts += ' --sbins %s' % ' '.join(args.sbins)

    # Increase requested memory for jobs that produce multiple maps
    ## NOTE: exceeded memory limit (1000 MB) with single day.
    ## Don't remember this being a problem, changing from double to float
    ## does *not* fix it :(
    ## Appears to cap out around 2500 MB
    sublines = ["request_memory = 4000"]
    #sublines = ["getenv = True"]
    if args.ebins or args.sbins: #or args.comp:
        sublines = ["request_memory = 5000"]

    # Environment for script
    cvmfs = '/cvmfs/icecube.opensciencegrid.org/users/juancarlos/tools/' \
            + 'py2-v1/icetray-start'
    metaproject = '/data/user/fmcnally/offline/V04-08-00/build'
    header = ['#!/bin/sh '+cvmfs, '#METAPROJECT '+metaproject]

    # Submit files
    print('Submitting %i files...' % len(dayFiles))
    #for i, oneDay in enumerate(dayFiles):
    for i in range(len(dayFiles)):

        jobID = 'timescramble_%s_%05i' % (args.config, random.randint(0,99999))
        cmd  = '%s/TimeScramble' % os.getcwd()
        ex = [cmd, '--batch_idx', str(i), '--yyyymmdd', dateList[i], c_opts]
        #ex = ex + ['--inFiles', oneDay]
        ex = ' '.join(ex)
        print(ex)
        print('')
        pysubmit(ex, sublines=sublines, test=args.test, jobID=jobID,
                header=header)

