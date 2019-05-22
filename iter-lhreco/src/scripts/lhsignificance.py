#!/usr/bin/python

import healpy as H
import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import pyfits
import pickle
import os
degree = np.pi / 180

def TopHatSmooth(nside,origMap,radius,average=False):
  pixelList = []
  dsum = 0.0
  newmap = np.zeros(len(origMap))

  # Loop over all pixels
  for ipix in range(0, len(origMap)):
    theta,phi = H.pix2ang(nside, ipix)
    center = H.dir2vec(theta,phi)

    # Grab pixels with given angular radius of ith pixel
    pixelList = H.query_disc(nside, center, radius)

    #Sum up pixels in the angular cap; set value of ith pixel to this sum
    dsum = 0.
    nsum = 0.
    for jpix in pixelList:
      if origMap[jpix] != H.UNSEEN:
         dsum += origMap[jpix]
         nsum += 1.0

    if average and nsum>0:
        newmap[ipix] = dsum/nsum # set to average
    else:
        newmap[ipix] = dsum
  return newmap

if __name__ == "__main__":

	import optparse

	usage = "usage: %prog [options]"
	parser = optparse.OptionParser(usage)
	parser.add_option("-t", "--timesteps", dest="timesteps",default=360,
                      help="time steps")
	parser.add_option("-e", "--exposure", dest="exposure",
                      help="exposure file (FITS)")
	parser.add_option("-E", "--initexposure", dest="initexposure",
                      help="initial exposure file (FITS)")
	parser.add_option("-n", "--norm", dest="norm",
                      help="normalization file (.dat)")
	parser.add_option("-N", "--initnorm", dest="initnorm",
                      help="initial normalization file (.dat)")
	parser.add_option("-p", "--prefix", dest="prefix",
                      help="prefix to local maps")
	parser.add_option("-r", "--relintfile", dest="relintfile",
                      help="relative intensity file (FITS)")
	parser.add_option("-o", "--output", dest="output",
                      help="output file (FITS)")
	parser.add_option("-l", "--localmaps", dest="localmaps", default="",
                      help="location of local maps")
	parser.add_option("-s", "--smoothrad", dest="smoothrad", default=10,type=float,
                      help="smoothing radius")
	opts, args = parser.parse_args()
	print opts

	smooth_rad = opts.smoothrad*degree

	# time resolution
	Ntimestep = opts.timesteps
	minstep = 0
	maxstep = Ntimestep
	
	# number of sectors
	Nsector = 1
	THETAMAX = 60.
	isovalue = 1.0
	
	Emap = H.read_map(opts.exposure)
	Emap0 = H.read_map(opts.initexposure)

		
	isovalue = 1.0
	

	#CR anisotropy
	diffCRmap = H.read_map(opts.relintfile)

	npix = diffCRmap.size
	nside = H.npix2nside(npix)

	#Emap0 = np.zeros(npix, dtype=np.double)
	Nmap = np.zeros((Ntimestep,npix), dtype=np.double)

	#CRmap = H.read_map(opts.datafile)
	CRmap = np.ones(npix,dtype=np.double)*isovalue

	norm0 = np.zeros(maxstep, dtype=np.double)
	norm = np.zeros(maxstep, dtype=np.double)

	norm0file = open(opts.initnorm,'r')
	normfile = open(opts.norm,'r')
	sourcefolder = opts.localmaps

	FOV = np.zeros(npix, dtype=np.double)

	for i in range(0,npix) :
		FOV[i] = 1
		theta,phi = H.pix2ang(nside,i)
		
		if theta > THETAMAX/180.*np.pi :
			FOV[i] = 0
		

	# initialize with previous results :

	#initial normalization :
	for timeidx in range(0,Ntimestep) :
		norm0[timeidx] = float(norm0file.readline().strip())
		norm[timeidx] = float(normfile.readline().strip())

		namefile = sourcefolder + '/hawclocal-v202-local-cutFile-lsa-10-03-2017-allbins-bin0_'+ '{:03d}'.format(timeidx) + '_of_360.fits.gz'
		namefile = opts.prefix + '{:03d}'.format(timeidx) + '_of_360.fits.gz'
		print timeidx
		tempmap = H.read_map(namefile)
		npixinput = tempmap.size
		Nmap[timeidx] = H.ud_grade(tempmap,nside)*(npixinput/(1.*npix))

			
	#HAWC position
	HAWClon = -97.0/180.*np.pi
	HAWClat =  19.0/180.*np.pi

	#IceCube location
	#HAWClon = 0.0/180.*np.pi
	#HAWClat = -90.0/180.*np.pi

	clat = np.cos(HAWClat)
	slat = np.sin(HAWClat)
	
		
	# calculate significance
	variancemap = np.zeros(npix, dtype=np.double)
	significancemap = np.zeros(npix, dtype=np.double)
	significancemap2 = np.zeros(npix, dtype=np.double)
		
	for timeidx in range(minstep,maxstep) :
		
		#print timeidx
		
		beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon
		cb = np.cos(beta)
		sb = np.sin(beta)
		
		# rotation matrix
		rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  			
		for i in range(0,npix) :
	
			vx,vy,vz = H.pix2vec(nside,i)
		
			#rotation from local frame to Equatorial (ra,dec)
			vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
			vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
			vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
			
			j = H.vec2pix(nside,vpx,vpy,vpz)
			
			# global significance
			if Emap0[j]*norm0[timeidx] > 0.0 and Emap[j]*norm[timeidx] > 0.0 :
				variancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap[j]*norm[timeidx] 
				variancemap[i] += +2.0*CRmap[i]*Emap0[j]*norm0[timeidx]
				temp1 = Emap[j]/Emap0[j]*norm[timeidx]/norm0[timeidx]
				variancemap[i] += 2.0*Nmap[timeidx][j]*np.log(temp1*(1.0+diffCRmap[i]/CRmap[i]))


	print "cleaning out nonsensical values"
	for i in range(npix):
		if variancemap[i] < 0:
			variancemap[i] = 0.

	print "smoothing relint"
	smoothed_relint = TopHatSmooth(nside,diffCRmap,radius=smooth_rad,average=True)
	print "smoothing abs"
	smoothed_abs = TopHatSmooth(nside,np.abs(diffCRmap),radius=smooth_rad,average=True)
	sign = smoothed_relint/smoothed_abs
	print "smoothing variance"
	smoothed_variance = TopHatSmooth(nside, variancemap,smooth_rad, average=False)
	significancemap = np.sqrt(smoothed_variance)*sign

	if os.path.exists(opts.output): 
	   os.remove(opts.output)
	H.write_map(opts.output,significancemap)
