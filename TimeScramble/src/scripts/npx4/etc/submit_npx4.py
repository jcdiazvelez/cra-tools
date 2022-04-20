#!/usr/bin/env python

# An attempt at rewriting submit_npx4.sh in python
import os, stat, subprocess
import numpy as np

import myGlobals as my

def py_submit(exelines, jobID=None, sublines=None, test=False):

    # Option for testing off cluster
    if test:
        #if len(exelines) > 1:
        #    raise SystemExit('Multiple exelines not available in test')
        #os.system(exelines[0])
        os.system(exelines)
        return

    # Setup global variables
    my.setupGlobals(verbose=False)

    if jobID == None:
        randint = np.random.uniform(100000)
        jobID = 'npx4-%05d' % randint

    outexe = '%s/npx4-execs/%s.sh' % (my.npx4, jobID)
    condor_script = '%s/2sub.sub' % my.npx4

    # Run eval statement as it doesn't run by default when fed to script
    #setPath = "echo eval `/cvmfs/icecube.opensciencegrid.org/setup.sh`"
    setPath = "echo eval `/cvmfs/icecube.opensciencegrid.org/py2-v1/setup.sh`"
    p = subprocess.Popen(setPath, stdout=subprocess.PIPE, shell=True)
    path, err = p.communicate()
    path = path.strip()

    # Setup execution script
    lines = [
        "#!/bin/bash",
        "date",
        "hostname",
        "",
        #"cd %s" % os.getcwd(),
        #"eval `/cvmfs/icecube.opensciencegrid.org/setup.sh`",
        "%s" % path,
        ""
    ]
    lines += exelines
    lines += ["date"]
    lines = [l + '\n' for l in lines]

    with open(outexe, 'w') as f:
        f.writelines(lines)

    # Make file executable
    st = os.stat(outexe)
    os.chmod(outexe, st.st_mode | stat.S_IEXEC)

    # Condor submission script
    lines = [
        "Universe = vanilla\n",
        "Executable = %s/npx4-execs/%s.sh\n" % (my.npx4, jobID),
        "Log = %s/npx4-logs/%s.log\n" % (my.npx4, jobID),
        "Output = %s/npx4-out/%s.out\n" % (my.npx4, jobID),
        "Error = %s/npx4-error/%s.error\n" % (my.npx4, jobID),
        "Notification = NEVER\n",
        "Queue\n"
    ]

    # Option for additional lines to submission script
    if sublines != None:
        for l in sublines:
            lines.insert(-1, '%s\n' % l)

    with open(condor_script, 'w') as f:
        f.writelines(lines)

    os.system('condor_submit %s' % condor_script)

