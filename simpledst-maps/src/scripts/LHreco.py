#!/usr/bin/python

import healpy as H
import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import pyfits
import pickle

if __name__ == "__main__":

	# resolution of input maps
	nsideinput = 64
	nsideinput = 256
	npixinput = H.nside2npix(nsideinput)
	#sourcefolder = './paper3_newfreq_nside64_360steps/'
	sourcefolder = '/data/user/juancarlos/anisotropy/IceCube/localmaps/2015/combined/'
	
	# time resolution
	Ntimestep = 360
	minstep = 0
	maxstep = 360
	
	# number of sectors
	Nsector = 1
	THETAMAX = 60.
	
	# resolution of output maps
	nside = 64
	npix = H.nside2npix(nside)
	
	# output paramters
	Nstart = 0        # if different from 0 indicates the continuation of a previous run
	Niteration = 20  #
	foldername = './icecube_newfreq_180/'
	
	sectormap = np.zeros(npix, dtype=np.int)
	Nmap = np.zeros((Ntimestep,npix), dtype=np.double)
	
	for timeidx in range(0,Ntimestep) :
		#use this for real data with nsideinput=256
		#namefile = '/home/fiorino/hawc100-local-maps/hawclocal-v126-nhit30-run000630-001089_N256_' + '{:03d}'.format(timeidx) + '_of_360.fits.gz'
		
		# fake data:
		#namefile = sourcefolder + 'local_timeindex' + '{:03d}'.format(timeidx) + '.fits'
		namefile = sourcefolder + 'CR_ICECUBE_LOCAL_NSIDE256_degbin-' + '{:03d}'.format(timeidx) + '.fits.gz'
		tempmap = H.read_map(namefile)
		Nmap[timeidx] = H.ud_grade(tempmap,nside)*(npixinput/(1.*npix))
			
		# additional zenith cut? :
		
		#for i in range(0,npix) :
		#	theta,phi = H.pix2ang(nside,i)
		#	if theta > 60.0/180.*np.pi :
		#		Nmap[timeidx][i] = 0.0

	Emap0 = np.zeros((Nsector,npix), dtype=np.double)
	norm0 = np.zeros((Nsector,maxstep), dtype=np.double)
	
	totsector = np.zeros(Nsector, dtype=np.double)
	
	Emap = np.zeros((Nsector,npix), dtype=np.double)
	norm = np.zeros((Nsector,maxstep), dtype=np.double)
	newnorm = np.zeros((Nsector,maxstep), dtype=np.double)
	
	CRmap = np.zeros(npix, dtype=np.double)
	diffCRmap = np.zeros(npix, dtype=np.double)
	
	#HAWC position
	#HAWClon = -97.0/180.*np.pi
	#HAWClat =  19.0/180.*np.pi

	#IceCube location
	HAWClon = 0.0/180.*np.pi
	HAWClat = -90.0/180.*np.pi

	clat = np.cos(HAWClat)
	slat = np.sin(HAWClat)
	
	#create weight function for sector
	
	# trivial weight :
	
	for i in range(0,npix) :
		theta,phi = H.pix2ang(nside,i)
		if theta > THETAMAX/180.*np.pi :
			sectormap[i] = -1
	
	#Equi-Zenith Angle method :
	
	#for i in range(0,npix) :
	#	theta,phi = H.pix2ang(nside,i)
	#	if theta > THETAMAX/180.*np.pi :
	#		sectormap[i] = -1
	#	else :
	#		sectormap[i] = min(Nsector-1,np.floor((theta)*Nsector/((1.*THETAMAX/180.*np.pi))))

	#H.mollview(sectormap)
	#show()
	#exit(3)
	
	#normalization of isotropic flux
	isovalue = 1.0
	
	#initial CR distribution isotropic:
	for i in range(0,npix) :
		CRmap[i] = isovalue
		
	#initial exposure :	
	
	for i in range(0,npix) :

		idx = sectormap[i]
		if idx < 0 :
			continue
			
		for timeidx in range(minstep,maxstep) :
			Emap0[idx][i] += Nmap[timeidx][i]
			totsector[idx] += Nmap[timeidx][i]	
			
	for i in range(0,npix) :
		idx = sectormap[i]
		if idx < 0 :
			continue
		Emap0[idx][i] = Emap0[idx][i]/totsector[idx]
		Emap[idx][i] = Emap0[idx][i]
	
	#initial normalization :
	
	for timeidx in range(minstep,maxstep) :
		
		for i in range(0,npix) :
			idx = sectormap[i]
			if idx < 0 :
				continue
			norm[idx][timeidx] += Nmap[timeidx][i]/isovalue
			norm0[idx][timeidx] += Nmap[timeidx][i]/isovalue
		
	#initial CR anisotropy
		
	diffCRmap = np.zeros(npix, dtype=np.double)
	FOVmap = np.zeros(npix, dtype=np.int)
	
	for i in range(0,npix) :
				
		temp1 = temp2 = 0.0	
		vx,vy,vz = H.pix2vec(nside,i)
			
		for timeidx in range(minstep,maxstep) :
			
			beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon
			cb = np.cos(beta)
			sb = np.sin(beta)
		
			# rotation matrix
			rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
			
			#rotation from local frame to Equatorial (ra,dec)
			vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
			vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
			vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
		
			j = H.vec2pix(nside,vpx,vpy,vpz)
			
			idx = sectormap[j]
			if idx < 0 :
				continue
				
			if Emap0[idx][j] > 0.0 :
				temp2 += norm[idx][timeidx]*Emap0[idx][j]
					
		if temp2 > 0.0 :
			FOVmap[i] = 1
	
	npixFOV = 1.*sum(FOVmap)
	diffCRmap = np.zeros(npix, dtype=np.double)
	
	# write data or initialize with previous results :
	
	if Nstart > 0 :
		norm = np.array(pickle.load(open(foldername + "norm_32_360_iteration" + str(Nstart) + ".dat", "r" )))
		Emap = np.array(pickle.load(open(foldername + "exposure_32_360_iteration" + str(Nstart) + ".dat", "r" )))
	else :
		pickle.dump(norm, open(foldername + "norm_32_360_iteration0.dat", "w" ) )
		pickle.dump(Emap0,open(foldername + "exposure_32_360_iteration0.dat","w") )		
	
	for iteration in range(Nstart,Niteration) :
	
		namefits = foldername + 'CR_32_360_iteration' + str(iteration+1) + '.fits'
		nameEfits = foldername + 'exposure_32_360_iteration' + str(iteration+1) + '.dat'
		nameNfits = foldername + 'norm_32_360_iteration' + str(iteration+1) + '.dat'
		nameSIGfits = foldername + 'significance_32_360_iteration' + str(iteration+1) + '.fits'
		nameSIG2fits = foldername + 'significance2_32_360_iteration' + str(iteration+1) + '.fits'

		# calculate new CR anisotropy
		
		diffCRmap = np.zeros(npix, dtype=np.double)
	
		for i in range(0,npix) :
				
			temp1 = temp2 = 0.0	
			vx,vy,vz = H.pix2vec(nside,i)
			
			for timeidx in range(minstep,maxstep) :
		
				beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon
				cb = np.cos(beta)
				sb = np.sin(beta)
		
				#rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  				
				#rotation from local frame to Equatorial (ra,dec)
				vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
				vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
				vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
				idx = sectormap[j]
				if idx < 0 :
					continue
			
				if Emap0[idx][j] > 0.0 :
					temp1 += Nmap[timeidx][j]
					temp2 += norm[idx][timeidx]*Emap[idx][j]
					
			if temp2 > 0.0 :
				diffCRmap[i] = temp1/temp2-CRmap[i]
			
			
		# remove m=0 multipole moments :
		
		LMAX=180
		out  = H.anafast(diffCRmap,alm=True,lmax=LMAX)	
	
		for i in range(0,LMAX+1) :
	
			index = H.sphtfunc.Alm.getidx(LMAX,i,0)
			out[1][index] = 0.0
	
		diffCRmap = H.alm2map(out[1],nside,lmax=LMAX)
		
		# calculate new norm	
		
		for timeidx in range(minstep,maxstep) :
		
			temp1 = np.zeros(Nsector, dtype=np.double)
			temp2 = np.zeros(Nsector, dtype=np.double)
			
			beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon
			cb = np.cos(beta)
			sb = np.sin(beta)
			
			for i in range(0,npix) :
		
				idx = sectormap[i]
				if idx < 0 :
					continue
					
				vx,vy,vz = H.pix2vec(nside,i)
				
				# rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  				
				#rotation from Equatorial (ra,dec) to local frame
				vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
				vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
				vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if Emap0[idx][i] > 0.0 :
					temp1[idx] += Nmap[timeidx][i]
					temp2[idx] += Emap[idx][i]*(CRmap[j]+diffCRmap[j])
			
			for idx in range(0,Nsector) :
			
				if temp2[idx] > 0.0 :
					newnorm[idx][timeidx] = temp1[idx]/temp2[idx]
		
		norm = np.zeros((Nsector,maxstep), dtype=np.double)
		
		for timeidx in range(minstep,maxstep) : 
			for idx in range(0,Nsector) :
				norm[idx][timeidx] = newnorm[idx][timeidx]


		#caluculate new acceptance
		
		newEmap = np.zeros((Nsector,npix), dtype=np.double)
		
		for i in range(0,npix) :
		
			temp1 = temp2 = 0.0
			vx,vy,vz = H.pix2vec(nside,i)
		
			idx = sectormap[i]
			if idx < 0 :
				continue
					
			for timeidx in range(minstep,maxstep) :
		
				beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon
				cb = np.cos(beta)
				sb = np.sin(beta)
		
				# rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
				
				#rotation from Equatorial (ra,dec) to local frame
				vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
				vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
				vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if Emap0[idx][i] > 0.0 :
					temp1 += Nmap[timeidx][i]
					temp2 += norm[idx][timeidx]*(CRmap[j]+diffCRmap[j])
			
			if temp2 > 0.0 :
				newEmap[idx][i] = temp1/temp2
			else : 
				newEmap[idx][i] = 0.0
		
		for i in range(0,npix) :
			idx = sectormap[i]
			Emap[idx][i] = newEmap[idx][i]
			
		# cleaning
		
		for idx in range(0,Nsector) :
			temp = sum(Emap[idx])
			
			for timeidx in range(minstep,maxstep) :
				norm[idx][timeidx] = temp*norm[idx][timeidx]
			
			Emap[idx] = Emap[idx]/temp	
			
		# write iteration results
		
		H.write_map(namefits,diffCRmap/isovalue)	
		pickle.dump(Emap,open(nameEfits, "w" ))	
		pickle.dump(norm, open(nameNfits, "w" ))
		
		# calculate signifiance
		
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
			
				idx = sectormap[j] 
				
				if idx < 0 :
					continue
					
				# global significance
				if Emap0[idx][j]*norm0[idx][timeidx] > 0.0 and Emap[idx][j]*norm[idx][timeidx] > 0.0 :
					significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap[idx][j]*norm[idx][timeidx] 
					significancemap[i] += +2.0*CRmap[i]*Emap0[idx][j]*norm0[idx][timeidx]
					temp1 = Emap[idx][j]/Emap0[idx][j]*norm[idx][timeidx]/norm0[idx][timeidx]
					significancemap[i] += 2.0*Nmap[timeidx][j]*np.log(temp1*(1.0+diffCRmap[i]/CRmap[i]))

		H.write_map(nameSIGfits,significancemap)
