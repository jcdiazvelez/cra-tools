#!/usr/bin/python

import healpy as H
import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import pyfits
import pickle

if __name__ == "__main__":
	
	# input parameters
	nside = 64
	npix = H.nside2npix(nside)
	Ntimestep = 360
	lmax = 30
	totevents = 5e10
	sourcefolder = './sample1_nside64_360steps/'
	HAWClon = -97.0/180.*np.pi
	HAWClat =  19.0/180.*np.pi

	clat = np.cos(HAWClat)
	slat = np.sin(HAWClat)
	
	lmap = np.zeros((lmax+1,npix), dtype=np.double)	
	Cltrue = np.zeros(lmax+1, dtype=np.double)
	truediffCRmap = np.zeros(npix, dtype=np.double)
	trueEmap =  np.zeros(npix, dtype=np.double)

	Cltrue[0] = 0.0 # no monopole
	for i in range(1,lmax+1) :
		Cltrue[i] = 1e-7*18./(2.*i+1.)/(i+1.)/(i+2.) # or whatever you like
	
	
	for i in range(0,lmax+1) :

		Cl = np.zeros(lmax+1, dtype=np.double)
		Cl[i] = 1.
		
		tempmap = H.sphtfunc.synfast(Cl,nside,lmax=lmax)
		Cltemp = H.anafast(tempmap)
			
		for j in range(0,npix) :
			truediffCRmap[j] += tempmap[j]*np.sqrt(Cltrue[i]/Cltemp[i])
		
	H.write_map(sourcefolder + "true_relintensity_map.fits", truediffCRmap)	
	
	temp = 0.0
	for i in range(0,npix) :
		theta,phi = H.pix2ang(nside,i) 
		
		if theta/np.pi*180. < 60.0 : # hard 60deg zenith cut
			trueEmap[i] = np.cos(theta)*(1.0 + 0.2*np.sin(theta)*np.sin(phi+10./180.*np.pi)**2) # or whatever you like
		temp +=trueEmap[i] 
	
	for i in range(0,npix) :
		trueEmap[i] = trueEmap[i]/temp
		
	H.write_map(sourcefolder + "true_acceptance_map.fits", trueEmap)	
	
	truenorm = []
	devents = totevents/(1.*Ntimestep)
	for timeidx in range(0,Ntimestep) :
		truenorm.append(devents*(1.0 + 0.1*np.sin((70+timeidx)/(1.*Ntimestep)*np.pi*2))) # or whatever you like
	
	pickle.dump(truenorm, open(sourcefolder + "true_norm.dat", "w" ) )
	
	#exit(3)
	
	for timeidx in range(0,Ntimestep) :
	
		tempmap = np.zeros(npix, dtype=np.double)
		
		namefile = sourcefolder + 'local_timeindex' + '{:03d}'.format(timeidx) + '.fits'

		beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon
		cb = np.cos(beta)
		sb = np.sin(beta)
		
		# rotation matrix
		rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
			
		for i in range(0,npix) :
	
			vx,vy,vz = H.pix2vec(nside,i)

			#rotation from Equatorial (ra,dec) to local frame
			vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
			vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
			vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
			j = H.vec2pix(nside,vpx,vpy,vpz)
			
			mu = (1.0 + truediffCRmap[j])*trueEmap[i]*truenorm[timeidx]
			
			tempmap[i] = np.random.poisson(lam=mu)

		H.write_map(namefile,tempmap)	
		close()
