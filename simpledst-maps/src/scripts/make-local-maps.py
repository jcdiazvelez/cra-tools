#!/usr/bin/env python

import os
import sys
import glob


def call(args, cwd=None, stdout=None, stderr=None, env=None):
    """
    Wrap the call function of the subprocess module for older pythons.
    """

    if not cwd:
        cwd = os.getcwd()
    if not stdout:
        stdout = sys.stdout
    if not stderr:
        stderr = sys.stderr
    if not env:
        env = os.environ

    try:
        import subprocess
        return subprocess.call(args,cwd=cwd,stdout=stdout,stderr=stderr,env=env)
    except:
        import popen2
        if cwd:
            os.chdir(cwd)
        cmd = " ".join(args)
        proc = popen2.Popen4(cmd)
        for line in proc.fromchild:
            stdout.write(line)
        rc = proc.wait()
        return rc


if __name__ == "__main__":
    """
    Main function for the script.
    """

    
    outdir = "/data/user/juancarlos/anisotropy/IceCube/localmaps/all-years/newcuts-v2/N64/"
    call(["mkdir", "-p", outdir])

    files = glob.glob("/data/ana/CosmicRay/Anisotropy/IceCube/IC86-2015/simple-dst/ic86_simpleDST_*.root")
    files += glob.glob("/data/ana/CosmicRay/Anisotropy/IceCube/IC86-2014/simple-dst/ic86_simpleDST_*.root")
    files += glob.glob("/data/ana/CosmicRay/Anisotropy/IceCube/IC86-2013/simple-dst/ic86_simpleDST_*.root")
    files += glob.glob("/data/ana/CosmicRay/Anisotropy/IceCube/IC86-2012/simple-dst/ic86_simpleDST_*.root")
    #files += glob.glob("/data/ana/CosmicRay/Anisotropy/IceCube/IC86-2011/simple-dst/ic86_simpleDST_*.root")
    files += glob.glob("/data/exp/IceCube/2012/filtered/DST_IC86/simple_dst/ic86_simpleDST_*.root") #2011

    for filepath in files:
        
        outpath = outdir
        subdir = filepath.replace(".root","").replace("ic86_simpleDST_","").split('/')[len(filepath.split('/')) - 1]
        
        outpath += subdir
        
        print "Submiting file: ", filepath
        print "Output path: ", outpath
        
        env = "/data/user/juancarlos/simpledst-maps/build/scripts/env-shell.sh"
        # Don't recreate it if it already exists
        call(["/data/user/juancarlos/simpledst-maps/build/scripts/submit.sh",
                   env, "/data/user/juancarlos/simpledst-maps/build/make-local-maps",
                   "--outdir="+outpath,  "--input="+filepath, "--nsideout=64",
                   "--spline=/data/user/fmcnally/anisotropy/sim/IC86-III_10649_spline.fits",
                   "--elogmax=4.5",
                   "--rloglmax=15",
                   "--ndirc_min=9",
                   "--ldirc_min=200",
                   ])
            
     
