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
    nmaps = 100
    max= 100

    outdir = "/data/user/juancarlos/anisotropy/isotropic/poisson.v8/"
    icecube_maps="/data/user/juancarlos/anisotropy/IceCube/localmaps/all-years/newcuts/N64/combined/"
    hawc_maps="/data/user/juancarlos/anisotropy/HAWC/local/v8/combined/" 

    call(["mkdir", "-p", outdir])

    for bin in range(0,nmaps):
        outpath = outdir+"%u"%bin
        
        print "Submiting bin: ", bin
        print "Output path: ", outpath
        
        submit = "/data/user/juancarlos/anisotropy/isotropic/scripts/submit.sh"
        env = "/data/user/juancarlos/anisotropy/isotropic/scripts/env-shell.sh"
        exe = "/data/user/juancarlos/anisotropy/isotropic/poisson.py"
        # Don't recreate it if it already exists
        call([submit, env, 'python', exe, "--nside=64",
             "--icecube="+icecube_maps, "-j %u" % bin,"--hawc="+hawc_maps, "--dest="+outpath])
            
     
