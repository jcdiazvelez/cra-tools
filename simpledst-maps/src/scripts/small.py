#!/usr/bin/python

import healpy as H
import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import pyfits
import pickle


degree = np.pi/180.0

def maskIceCubeLocal(skymap):
    """Set pixels outside the IceCube mask region to UNSEEN"""
    npix = skymap.size
    nside = H.npix2nside(npix)

    minPix = H.ang2pix(nside, 0 * degree, 0.)
    maxPix = H.ang2pix(nside, 90 * degree, 0.)
    skymap[0:minPix] = 0.0
    skymap[maxPix:npix] = 0.0

    return skymap

if __name__ == "__main__":

	# resolution of input maps
	nsideinput = 64
	npixinput = H.nside2npix(nsideinput)
	sourcefolder1 = '/data/user/juancarlos/anisotropy/HAWC/local/v7.1_local/'
	sourcefolder2 = '/data/user/juancarlos/anisotropy/IceCube/localmaps/all-years/newcuts/N64/combined/'
	
	# time resolution
	Ntimestep = 360
	minstep = 0
	maxstep = 360
	
	# resolution of output maps
	nside = 64
	npix = H.nside2npix(nside)
	
	# output paramters
	Nstart = 20        # if different from 0 indicates the continuation of a previous run
	foldername = './HAWCv7.1_and_IceCube_v2_all_IC86.1_small/' 
	
	FOV1 = np.zeros(npix, dtype=np.double)
	FOV2 = np.zeros(npix, dtype=np.double)
	
	Nmap1 = np.zeros((Ntimestep,npix), dtype=np.double)
	Nmap2 = np.zeros((Ntimestep,npix), dtype=np.double)

	mumap1 = np.zeros(npix, dtype=np.double)
	mu0map1 = np.zeros(npix, dtype=np.double)
	Nmapbin1 = np.zeros(npix, dtype=np.double)

	mumap2 = np.zeros(npix, dtype=np.double)
	mu0map2 = np.zeros(npix, dtype=np.double)
	Nmapbin2 = np.zeros(npix, dtype=np.double)

	
	for timeidx in range(0,Ntimestep) :
		namefile = sourcefolder1 + 'hawclocal-v7_bins1-2_N64_'+'{:03d}'.format(timeidx)+'_of_360.fits.gz'
		print namefile
		tempmap = H.read_map(namefile)
		Nmap1[timeidx] = H.ud_grade(tempmap,nside)*(npixinput/(1.*npix))
		
		namefile = sourcefolder2 + 'CR_ICECUBE_LOCAL_NSIDE64_degbin-'+'{:03d}'.format(timeidx)+'.fits.gz'
		print namefile
		tempmap = H.read_map(namefile)
		Nmap2[timeidx] = H.ud_grade(tempmap,nside)*(npixinput/(1.*npix))

	Emap01 = np.zeros(npix, dtype=np.double)
	norm01 = np.zeros(maxstep, dtype=np.double)
	
	Emap02 = np.zeros(npix, dtype=np.double)
	norm02 = np.zeros(maxstep, dtype=np.double)

	Emap1 = np.zeros(npix, dtype=np.double)
	norm1 = np.zeros(maxstep, dtype=np.double)
	newnorm1 = np.zeros(maxstep, dtype=np.double)
	
	Emap2 = np.zeros(npix, dtype=np.double)
	norm2 = np.zeros(maxstep, dtype=np.double)
	newnorm2 = np.zeros(maxstep, dtype=np.double)
	
	diffCRmap = np.zeros(npix, dtype=np.double)
	CRmap = np.zeros(npix, dtype=np.double)
	
	#HAWC position
	HAWClon1 = -97.0/180.*np.pi
	HAWClat1 =  19.0/180.*np.pi
	THETAMAX1 = 57.
	
	#IceCube location
	HAWClon2 =  270.0/180.*np.pi
	HAWClat2 =  -90.0/180.*np.pi
	THETAMAX2 = 65.
	SMOOTHRAD = 10.0/180.*np.pi
	
	
	for i in range(0,npix) :
		FOV1[i] = 1
		FOV2[i] = 1
		theta,phi = H.pix2ang(nside,i)
		
		if theta > THETAMAX1/180.*np.pi :
			FOV1[i] = 0
			
		if theta > THETAMAX2/180.*np.pi :
			FOV2[i] = 0	
	
	isovalue = 1.0
	
	#initial CR distribution isotropic:
	for i in range(0,npix) :
	    CRmap[i] = isovalue
	
	
	# initialize with previous results :
	
	norm1 = np.array(pickle.load(open(foldername + "norm1_32_360_iteration" + str(Nstart) + ".dat", "r" )))
	Emap1 = np.array(pickle.load(open(foldername + "exposure1_32_360_iteration" + str(Nstart) + ".dat", "r" )))
		
	norm2 = np.array(pickle.load(open(foldername + "norm2_32_360_iteration" + str(Nstart) + ".dat", "r" )))
	Emap2 = np.array(pickle.load(open(foldername + "exposure2_32_360_iteration" + str(Nstart) + ".dat", "r" )))
	namefits = foldername + 'CR_32_360_iteration' + str(Nstart) + '.fits'
	diffCRmap = H.read_map(namefits)
		
	
	LMAX=180
	alms3 = H.map2alm(diffCRmap, lmax=3 )
	#fitmap3 = H.alm2map(alms3, nside,lmax=LMAX)
	fitmap3 = H.alm2map(alms3, nside)
	diffCRmap = CRmap - fitmap3

	# substract l<=3
	subLMAX = 3
	redLMAX = 3
	out = H.anafast(tempmap,alm=True,lmax=subLMAX)
	outr = H.anafast(tempmap,alm=True,lmax=redLMAX)
	submap = H.alm2map(out[1],nside,lmax=subLMAX)

	diffCRmap = diffCRmap - submap
		
	# calculate new norm	
		
	#caluculate new acceptance
		
	# calculate signifiance
		
	significancemap = np.zeros(npix, dtype=np.double)
	significancemap2 = np.zeros(npix, dtype=np.double)
	smoothsignificancemap= np.zeros(npix, dtype=np.double)
		
	for timeidx in range(minstep,maxstep) :
		print maxstep-timeidx
		
		beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon1
		cb = np.cos(beta)
		sb = np.sin(beta)
			
		clat = np.cos(HAWClat1)
		slat = np.sin(HAWClat1)
		
		# rotation matrix
		rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  			
		for i in range(0,npix) :
	
				vx,vy,vz = H.pix2vec(nside,i)
		
				#rotation from local frame to Equatorial (ra,dec)
				vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
				vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
				vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
			
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if FOV1[j] == 0 :
					continue

				# global significance
				if Emap1[j]*norm1[timeidx] > 0.0 and CRmap[i] > 0:

					mumap1[i] += (diffCRmap[i]+submap[i]+CRmap[i])*Emap1[j]*norm1[timeidx]
					mu0map1[i] += (CRmap[i]+submap[i])*Emap1[j]*norm1[timeidx]
					Nmapbin1[i] += Nmap1[timeidx][j]

					
		beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon2
		cb = np.cos(beta)
		sb = np.sin(beta)
			
		clat = np.cos(HAWClat2)
		slat = np.sin(HAWClat2)
		
		# rotation matrix
		rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  			
		for i in range(0,npix) :
	
				vx,vy,vz = H.pix2vec(nside,i)
		
				#rotation from local frame to Equatorial (ra,dec)
				vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
				vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
				vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
			
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if FOV2[j] == 0 :
					continue
					
				# global significance
				if Emap2[j]*norm2[timeidx] > 0.0 and CRmap[i] > 0:

					mumap2[i]  += (diffCRmap[i]+submap[i]+CRmap[i])*Emap2[j]*norm2[timeidx]
					mu0map2[i] += (CRmap[i]+submap[i])*Emap2[j]*norm2[timeidx]
					Nmapbin2[i] += Nmap2[timeidx][j]
					

	for j in range(0,npix) :

	    vec = H.pix2vec(nside,j)
	    set = H.query_disc(nside,vec,SMOOTHRAD)
	    mutemp = 0.0
	    mu0temp = 0.0
	    Nmaptemp = 0.0

	    for k in set :

	       if mu0map1[k] > 0.0:
	            mutemp += mumap1[k]
	            mu0temp += mu0map1[k]
	            Nmaptemp += Nmapbin1[k]

#	       if mu0map2[k] > 0.0 :
#	            mutemp += mumap2[k]
#	            mu0temp += mu0map2[k]
#	            Nmaptemp += Nmapbin2[k]


	    if mutemp > 0.0:
	       smoothsignificancemap[j] = -2.0*mutemp+2.0*mu0temp+2.0*Nmaptemp*np.log(mutemp/mu0temp)
	       print -2.0*mutemp+2.0*mu0temp+2.0*Nmaptemp*np.log(mutemp/mu0temp)

	nameSIGfits = foldername + 'small_scale_significance_32_360_iteration' + str(Nstart) + '.fits.gz'
	H.write_map(nameSIGfits,smoothsignificancemap)
