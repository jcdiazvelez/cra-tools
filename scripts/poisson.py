#!/usr/bin/env python
import os, os.path
import healpy as hp
import numpy as np
from numpy import random
from optparse import OptionParser
import glob
import multiLH


def generate_local_variation(localmap):
  alpha = 0.1;
  
  npix = localmap.size
  maps = []
  randmap = np.zeros(npix)
  for i in range(npix):
          # Throw a Poisson random number using the data as a mean
          Nb = localmap[i]
          Nd = random.poisson(Nb) # need to seed RNG
          if Nd > 0 and Nb > 0:
             randmap[i] = Nd
  return randmap


if __name__ == "__main__":
    # Set up command line options

    usage  = "usage: %prog [option] MAPS.fits"
    parser = OptionParser(usage)
    parser.add_option("-j", "--job", dest="job", help="job ID", type="int")
    parser.add_option("-n", "--nside", default=64, dest="nside", help="HEALPix NSide", type="int")
    parser.add_option("--icecube", dest="icecube", help="IceCube Maps", type="str")
    parser.add_option("--hawc", dest="hawc", help="HAWC Maps", type="str")
    parser.add_option("--dest", dest="dest", default='.', help="Detination", type="str")
    parser.add_option("--scratch", dest="scratch", default='.', help="Scratch", type="str")
    parser.add_option("--iterations", dest="niter", default=20, help="Number of iterations", type="int")

    (options, args) = parser.parse_args()

    random.seed(options.job)
    icecube_local_map = np.zeros(hp.nside2npix(options.nside))
    hawc_local_map = np.zeros(hp.nside2npix(options.nside))

    icecube_maps = glob.glob(options.icecube+"/*.fits.gz")
    hawc_maps = glob.glob(options.hawc+"/*.fits.gz")
    print "icecube maps:", ", ".join(icecube_maps[0:4]), ",..."
    print "hawc maps:", ", ".join(hawc_maps[0:4]), ",..."

    hawc_poisson_maps = []
    i3_poisson_maps = []
    for lm in icecube_maps:
        icecube_local_map += hp.read_map(lm,0)
        icecube_poisson_local_map = generate_local_variation(icecube_local_map)
        i3_poisson_maps.append(icecube_poisson_local_map)

    for lm in hawc_maps:
        hawc_local_map = hp.read_map(lm,0)
        hawc_poisson_local_map = generate_local_variation(hawc_local_map)
        hawc_poisson_maps.append(hawc_poisson_local_map)
 
    poissondir = options.dest
    os.makedirs(poissondir)

    multiLH.multiLH(hawc_poisson_maps, i3_poisson_maps, dest=poissondir,
        #HAWC position
        HAWClon1 = -97.0/180.*np.pi,
        HAWClat1 =  19.0/180.*np.pi,
        THETAMAX1 = 57.,

        #IceCube location
        HAWClon2 =  270.0/180.*np.pi,
        HAWClat2 =  -90.0/180.*np.pi,
        THETAMAX2 = 65.,
        Niteration = options.niter )
