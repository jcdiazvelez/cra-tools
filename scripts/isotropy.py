#!/usr/bin/env python
import os, os.path
import healpy as hp
import numpy as np
from numpy import random
from optparse import OptionParser
import glob
import multiLH


def generate_iso_local_maps(localmap,nbins=360):
  alpha = 0.1;
  
  npix = localmap.size
  maps = []
  for ibin in range(nbins):
      isomap = np.zeros(npix)
      for i in range(npix):
          # Throw a Poisson random number using the data as a mean
          Nb = localmap[i]/nbins
          Nd = random.poisson(Nb) # need to seed RNG
          if Nd > 0 and Nb > 0:
             isomap[i] = Nd
      maps.append(isomap)
  return maps


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

    for lm in icecube_maps:
        icecube_local_map += hp.read_map(lm,0)

    for lm in hawc_maps:
        hawc_local_map += hp.read_map(lm,0)
 
    hawc_iso = generate_iso_local_maps(hawc_local_map, nbins=360)
    i3_iso = generate_iso_local_maps(icecube_local_map, nbins=360)

    isodir = options.dest
    os.makedirs(isodir)

    multiLH.multiLH(hawc_iso, i3_iso, dest=isodir,
        #HAWC position
        HAWClon1 = -97.0/180.*np.pi,
        HAWClat1 =  19.0/180.*np.pi,
        THETAMAX1 = 57.,

        #IceCube location
        HAWClon2 =  270.0/180.*np.pi,
        HAWClat2 =  -90.0/180.*np.pi,
        THETAMAX2 = 65.,
        Niteration = options.niter )
