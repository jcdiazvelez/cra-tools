#!/usr/bin/env python

#######################################################
# Runs through error files, returns a list of files where the process ran
# correctly, and can remove files associated with good runs
#######################################################

import os, sys, re, argparse
from glob import glob

from getTime import getTime


# List of any safe errors you want to ignore
def getErrorList():
    safeErrors = []
    # General safe errors
    safeErrors.append('Error in <TSystem::ExpandFileName>: input: $HOME/.root.mimes, output: $HOME/.root.mimes')
    # IceTop time-scrambling errors
    safeErrors.append('Error in <TTree::SetBranchStatus>: unknown branch -> ShowerLLH_proton.energy')
    safeErrors.append('Error in <TTree::SetBranchStatus>: unknown branch -> ShowerLLH_iron.energy')
    safeErrors.append('Error in <TTree::SetBranchStatus>: unknown branch -> maxLLH_proton.value')
    safeErrors.append('Error in <TTree::SetBranchStatus>: unknown branch -> maxLLH_iron.value')
    safeErrors.append('Error in <TTree::SetBranchStatus>: unknown branch -> LaputopStandardParams.s125')
    safeErrors.append('Error in <TTree::SetBranchStatus>: unknown branch -> LaputopSmallShowerParams.s125')
    # Time-scrambling errors assoc. w/ file already exists?
    safeErrors.append('FITS error')
    safeErrors.append("terminate called after throwing an instance of 'Message_error'")
    safeErrors.append("/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/RHEL_7_x86_64/lib/python3.7/site-packages/astropy/config/configuration.py:532: ConfigurationMissingWarning: Configuration defaults will be used due to OSError:Could not find unix home directory to search for astropy config dir on None")
    return safeErrors


def goodFile(filename, strictError):

    # Check the error file for any lines that aren't harmless
    with open(filename, 'r') as f:
        err = f.readlines()
        err = [l.strip() for l in err]

    # Remove non-error-related I3Tray output
    oks = ['NOTICE','INFO','WARN']
    err = [l for l in err if l.split(' ')[0] not in oks]
    # Option: ignore all lines that don't explicity list as "Error"

    if strictError:
        err = [l for l in err if any([e in l for e in ['error','Error']])]

    safeErrors = getErrorList()
    if any([l not in safeErrors for l in err]):
        return False

    return True
    

def resubmit(npxdir, jobID):

    ## NOTE: Assumes the last jobs you submitted involved the same options
    ## (example: request_memory), and that the sub script is '2sub.sub'

    # Condor submission script
    d = {}
    d['executable'] = '%s/npx4-execs/%s.sh\n' % (npxdir, jobID)
    d['log'] = '%s/npx4-log/%s.log\n' % (npxdir, jobID)
    d['output'] = '%s/npx4-out/%s.out\n' % (npxdir, jobID)
    d['error'] = '%s/npx4-error/%s.error\n' % (npxdir, jobID)

    with open('2sub.sub','r') as f:
        lines = f.readlines()

    for key in d.keys():
        for i, l in enumerate(lines):
            if l.split(' ')[0].lower() == key:
                lines[i] = '%s = %s' % (key, d[key])

    with open('2sub.sub', 'w') as f:
        f.writelines(lines)

    os.system('condor_submit %s' % '2sub.sub')



if __name__ == "__main__":

    parser = argparse.ArgumentParser(
            description='Scan npx4 files and remove good runs')
    parser.add_argument('--npxdir', dest='npxdir',
            default='/home/fmcnally/npx4',
            help='Location of npx4 folders')
    parser.add_argument('--purge', dest='purge',
            default=False, action='store_true',
            help='Remove all npx4 files')
    parser.add_argument('--strictError', dest='strictError',
            default=False, action='store_true',
            help='ONLY look at lines with the word "Error" in them')
    parser.add_argument('--badInfo', dest='badInfo',
            default=False, action='store_true',
            help='Print bad output lines')
    parser.add_argument('--rerun', dest='rerun',
            default=False, action='store_true',
            help='Rerun bad run files')
    args = parser.parse_args()

    # List of all npx4 files and jobIDs
    files = glob('%s/npx4-*/*' % args.npxdir)
    jobIDs = sorted(set([re.split('/|\.', f)[-2] for f in files]))

    # Input function dependent on python version
    inFunc = input if sys.version_info.major == 3 else raw_input

    # Rerun option
    if args.rerun:
        yn = inFunc('Rerun all bad executables? [y/n]: ')
        if yn != 'y':
            print('rerun canceled. aborting...')
            sys.exit(0)

    # Purge option
    if args.purge:
        yn = inFunc('Purge all %i npx4 files? [y/n]: ' % len(files))
        if yn == 'y':
            for f in files:
                os.remove(f)
        else:
            print('purge canceled')
        sys.exit(0)

    goodRuns, badRuns, tList = [],[],[]
    t = 0.

    for jobID in jobIDs:

        #badRuns.append(jobID)
        #continue

        # Skip if not started/finished
        outFile = '{}/npx4-out/{}.out'.format(args.npxdir, jobID)
        if not os.path.isfile(outFile):
            continue
        with open('%s/npx4-out/%s.out' % (args.npxdir, jobID), 'r') as f:
            lines = f.readlines()
        if lines == []:
            print('{} is empty!'.format(outFile))
            continue
        #if lines[-1].strip() != 'Fin':
        #    continue

        # If it doesn't have an error file, remove orphan files
        errFile = '{}/npx4-error/{}.error'.format(args.npxdir, jobID)
        if not os.path.isfile(errFile):
            os.system('rm -rf %s/npx4-*/%s.*' % (args.npxdir, jobID))
            continue

        # Check the error file for any lines that aren't harmless
        isGood = goodFile(errFile, args.strictError)

        # Append run number to good or bad run list
        t0 = getTime('%s/npx4-logs/%s.log' % (args.npxdir, jobID))
        if isGood and t0:
            goodRuns.append(jobID)
            t += t0
            tList += [t0]
        #if isGood and not t0:       # Job still running?
        #    continue
        else:
            badRuns.append(jobID)

    # Print information on state of npx4 files
    nFinished = len(goodRuns) + len(badRuns)
    nRunning = len(jobIDs) - nFinished
    if nRunning != 0:
        print('Running (%i file(s)):' % nRunning)
        for jobID in jobIDs:
            if jobID not in goodRuns and jobID not in badRuns:
                print('  - %s' % jobID)
        print('')

    # Good runs
    if len(goodRuns) != 0 and not args.badInfo:
        print('Good runs (%i files(s)):' % len(goodRuns))
        for jobID in goodRuns:
            print(' --- %s ---' % jobID)
        print('  Average time per job: %.01f seconds' % (t/len(goodRuns)))
        tList.sort()
        print('    Min time: %.01f' % tList[0])
        print('    Max time: %.01f' % tList[-1])
        print('')

    # Bad runs
    if len(badRuns) != 0:
        print('Bad runs (%i file(s)):' % len(badRuns))
        for jobID in badRuns:
            print(' --- %s ---' % jobID)

            # Option to print detailed information for each file
            if args.badInfo:
                errFile = '%s/npx4-error/%s.error' % (args.npxdir, jobID)
                with open(errFile, 'r') as f:
                    lines = f.readlines()
                for l in lines:
                    try:
                        if l.split()[0] not in ['NOTICE','INFO','WARN','\n']:
                            print('    %s' % l.strip())
                    except IndexError:
                        continue
                print()

            # Option to resubmit failed executables
            if args.rerun:
                resubmit(args.npxdir, jobID)

        print('')

    # Option to remove good runs
    if len(goodRuns) != 0 and not args.badInfo and not args.rerun:
        yn = inFunc('Do you want to remove the good runs? [y/n]: ')
        if yn == 'y':
            for jobID in goodRuns:
                rmFiles = glob('%s/npx4-*/%s.*' % (args.npxdir, jobID))
                for f in rmFiles:
                    os.remove(f)


##===========================================================================##
## Relic code
##

## For reading information from condor_q output
##  user = getpass.getuser()
##  bashCommand = 'condor_q %s -wide' % user
##  process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
##  output = process.communicate()[0]
##  output = [i.strip() for i in output.split('\n')]
##
##  print(output)
##
##  # Make sure call to condor went correctly
##  # Fails with new update. Check later
##  for phrase in ['jobs','completed','removed','idle','running','held']:
##      if phrase not in output[-2]:
##          print('Fetch from condor failed')
##          return []
##  #0 jobs; 0 completed, 0 removed, 0 idle, 0 running, 0 held, 0 suspended

##===========================================================================##

