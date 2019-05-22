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

    
    outdir = "/data/user/juancarlos/anisotropy/IceCube/localmaps/all-years/newcuts/N64/combined/"
    call(["mkdir", "-p", outdir])

    for bin in range(0,360):
        #glob.glob("/data/user/juancarlos/anisotropy/IceCube/localmaps/2015/cuts/n256/*/CR_ICECUBE_LOCAL_NSIDE256_degbin-%02d.fits.gz" % bin)
        binfiles = glob.glob("/data/user/juancarlos/anisotropy/IceCube/localmaps/all-years/newcuts/N64/*/CR_ICECUBE_LOCAL_0.00-4.50GeV_NSIDE64_degbin-%02d.fits.gz" % bin)
        outfile = "CR_ICECUBE_LOCAL_NSIDE64_degbin-%03d.fits.gz" % bin
        outpath = outdir
        outpath += outfile
        
        print "Submiting bin: ", bin
        print "Output path: ", outpath
        
        submit = "/data/user/juancarlos/simpledst-maps/build/scripts/submit.sh"
        env = "/data/user/juancarlos/simpledst-maps/build/scripts/env-shell.sh"
        exe = "/data/user/juancarlos/simpledst-maps/build/combine-local-maps"
        # Don't recreate it if it already exists
        call([submit, env, exe, "--outfile="+outpath,  "--nside=64"]+binfiles)
        #call([env, exe, "--outfile="+outpath,  "--nside=64"]+binfiles)
            
     
