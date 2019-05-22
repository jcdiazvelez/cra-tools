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

#if __name__ == "__main__":
def multiLH(Nmap1, Nmap2, dest='.',
	#HAWC position
	HAWClon1 = -97.0/180.*np.pi,
	HAWClat1 =  19.0/180.*np.pi,
	THETAMAX1 = 57.,
	
	#IceCube location
	HAWClon2 =  270.0/180.*np.pi,
	HAWClat2 =  -90.0/180.*np.pi,
	THETAMAX2 = 65.,
	Niteration = 20  #
	):

	# resolution of input maps
	nsideinput = 64
	npixinput = H.nside2npix(nsideinput)
	
	# time resolution
	Ntimestep = 360
	minstep = 0
	maxstep = 360
	
	# resolution of output maps
	nside = 64
	npix = H.nside2npix(nside)
	
	# output paramters
	Nstart = 0        # if different from 0 indicates the continuation of a previous run
	#Niteration = 10  #
	#foldername = './HAWCv7.1_and_IceCube_v2_all_IC86.1/'
	foldername = dest + '/'
	
	FOV1 = np.zeros(npix, dtype=np.double)
	FOV2 = np.zeros(npix, dtype=np.double)
	
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
	
	CRmap = np.zeros(npix, dtype=np.double)
	diffCRmap = np.zeros(npix, dtype=np.double)
	
	
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
		
	#initial exposure :	
	
	totsector1 = 0
	totsector2 = 0
	
	for i in range(0,npix) :
	
		if FOV1[i] == 0 :
			continue
			
		for timeidx in range(minstep,maxstep) :
			Emap01[i] += Nmap1[timeidx][i]
			totsector1 += Nmap1[timeidx][i]	
			
	for i in range(0,npix) :
		
		if FOV1[i] == 0 :
			continue
			
		Emap01[i] = Emap01[i]/totsector1
		Emap1[i] = Emap01[i]		
	
	for i in range(0,npix) :
	
		if FOV2[i] == 0 :
			continue
			
		for timeidx in range(minstep,maxstep) :
			Emap02[i] += Nmap2[timeidx][i]
			totsector2 += Nmap2[timeidx][i]	

	for i in range(0,npix) :
		
		if FOV2[i] == 0 :
			continue
			
		Emap02[i] = Emap02[i]/totsector2
		Emap2[i] = Emap02[i]		
		
	#initial normalization :
	
	for timeidx in range(minstep,maxstep) :
		
		for i in range(0,npix) :
		
			if FOV1[i] == 0 :
				continue
				
			norm1[timeidx] += Nmap1[timeidx][i]/isovalue
			norm01[timeidx] += Nmap1[timeidx][i]/isovalue
			
		for i in range(0,npix) :
		
			if FOV2[i] == 0 :
				continue
				
			norm2[timeidx] += Nmap2[timeidx][i]/isovalue
			norm02[timeidx] += Nmap2[timeidx][i]/isovalue	
		
	#initial CR anisotropy
		
	diffCRmap = np.zeros(npix, dtype=np.double)
	
	# write data or initialize with previous results :
	
	if Nstart > 0 :
		norm1 = np.array(pickle.load(open(foldername + "norm1_32_360_iteration" + str(Nstart) + ".dat", "r" )))
		Emap1 = np.array(pickle.load(open(foldername + "exposure1_32_360_iteration" + str(Nstart) + ".dat", "r" )))
		
		norm2 = np.array(pickle.load(open(foldername + "norm2_32_360_iteration" + str(Nstart) + ".dat", "r" )))
		Emap2 = np.array(pickle.load(open(foldername + "exposure2_32_360_iteration" + str(Nstart) + ".dat", "r" )))
		
	else :
		pickle.dump(norm1, open(foldername + "norm1_32_360_iteration0.dat", "w" ) )
		pickle.dump(Emap01,open(foldername + "exposure1_32_360_iteration0.dat","w") )	
		
		pickle.dump(norm2, open(foldername + "norm2_32_360_iteration0.dat", "w" ) )
		pickle.dump(Emap02,open(foldername + "exposure2_32_360_iteration0.dat","w") )	
	
	for iteration in range(Nstart,Niteration) :
	
		namefits = foldername + 'CR_32_360_iteration' + str(iteration+1) + '.fits'
		
		nameE1fits = foldername + 'exposure1_32_360_iteration' + str(iteration+1) + '.dat'
		nameN1fits = foldername + 'norm1_32_360_iteration' + str(iteration+1) + '.dat'
		
		nameE2fits = foldername + 'exposure2_32_360_iteration' + str(iteration+1) + '.dat'
		nameN2fits = foldername + 'norm2_32_360_iteration' + str(iteration+1) + '.dat'
		
		nameSIGfits = foldername + 'significance_32_360_iteration' + str(iteration+1) + '.fits'
		nameSIG2fits = foldername + 'significance2_32_360_iteration' + str(iteration+1) + '.fits'

		# calculate new CR anisotropy
		
		diffCRmap = np.zeros(npix, dtype=np.double)
	
		for i in range(0,npix) :
				
			temp1 = temp2 = 0.0	
			vx,vy,vz = H.pix2vec(nside,i)
			
			for timeidx in range(minstep,maxstep) :
		
				beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon1
				cb = np.cos(beta)
				sb = np.sin(beta)
				
				clat = np.cos(HAWClat1)
				slat = np.sin(HAWClat1)
		
				#rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  				
				#rotation from local frame to Equatorial (ra,dec)
				vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
				vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
				vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
				
				if FOV1[j] == 0 :
					continue	
			
				if Emap01[j] > 0.0 :
					temp1 += Nmap1[timeidx][j]
					temp2 += norm1[timeidx]*Emap1[j]
			
			for timeidx in range(minstep,maxstep) :
		
				beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon2
				cb = np.cos(beta)
				sb = np.sin(beta)
		
				clat = np.cos(HAWClat2)
				slat = np.sin(HAWClat2)
				
				#rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  				
				#rotation from local frame to Equatorial (ra,dec)
				vpx = rot[0][0]*vx+rot[0][1]*vy+rot[0][2]*vz
				vpy = rot[1][0]*vx+rot[1][1]*vy+rot[1][2]*vz
				vpz = rot[2][0]*vx+rot[2][1]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
				
				if FOV2[j] == 0 :
					continue	
			
				if Emap02[j] > 0.0 :
					temp1 += Nmap2[timeidx][j]
					temp2 += norm2[timeidx]*Emap2[j]
					
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
		
			temp1 = 0.0
			temp2 = 0.0
			
			beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon1
			
			cb = np.cos(beta)
			sb = np.sin(beta)
			
			clat = np.cos(HAWClat1)
			slat = np.sin(HAWClat1)
				
			for i in range(0,npix) :
		
				if FOV1[i] == 0 :
					continue	
					
				vx,vy,vz = H.pix2vec(nside,i)
				
				# rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  				
				#rotation from Equatorial (ra,dec) to local frame
				vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
				vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
				vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if Emap01[i] > 0.0 :
					temp1 += Nmap1[timeidx][i]
					temp2 += Emap1[i]*(CRmap[j]+diffCRmap[j])
			
			if temp2 > 0.0 :
				newnorm1[timeidx] = temp1/temp2
		
			temp1 = 0.0
			temp2 = 0.0
			
			beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon2
				
			cb = np.cos(beta)
			sb = np.sin(beta)
			
			clat = np.cos(HAWClat2)
			slat = np.sin(HAWClat2)
			
			for i in range(0,npix) :
		
				if FOV2[i] == 0 :
					continue	
					
				vx,vy,vz = H.pix2vec(nside,i)
				
				# rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
  				
				#rotation from Equatorial (ra,dec) to local frame
				vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
				vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
				vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if Emap02[i] > 0.0 :
					temp1 += Nmap2[timeidx][i]
					temp2 += Emap2[i]*(CRmap[j]+diffCRmap[j])
			
			if temp2 > 0.0 :
				newnorm2[timeidx] = temp1/temp2
				
		norm1 = np.zeros(maxstep, dtype=np.double)
		norm2 = np.zeros(maxstep, dtype=np.double)
		
		for timeidx in range(minstep,maxstep) : 
			norm1[timeidx] = newnorm1[timeidx]
			norm2[timeidx] = newnorm2[timeidx]

		#caluculate new acceptance
		
		newEmap = np.zeros(npix, dtype=np.double)
		
		for i in range(0,npix) :
		
			temp1 = temp2 = 0.0
			vx,vy,vz = H.pix2vec(nside,i)
		
			if FOV1[i] == 0 :
				continue
					
			for timeidx in range(minstep,maxstep) :
		
				beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon1
				cb = np.cos(beta)
				sb = np.sin(beta)
				
				clat = np.cos(HAWClat1)
				slat = np.sin(HAWClat1)
		
				# rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
				
				#rotation from Equatorial (ra,dec) to local frame
				vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
				vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
				vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if Emap01[i] > 0.0 :
					temp1 += Nmap1[timeidx][i]
					temp2 += norm1[timeidx]*(CRmap[j]+diffCRmap[j])
			
			if temp2 > 0.0 :
				newEmap[i] = temp1/temp2
			else : 
				newEmap[i] = 0.0
		
		for i in range(0,npix) :
			Emap1[i] = newEmap[i]
			
		newEmap = np.zeros(npix, dtype=np.double)
		
		for i in range(0,npix) :
		
			temp1 = temp2 = 0.0
			vx,vy,vz = H.pix2vec(nside,i)
		
			if FOV2[i] == 0 :
				continue
					
			for timeidx in range(minstep,maxstep) :
		
				beta = timeidx/(1.*Ntimestep)*np.pi*2 + HAWClon2
				cb = np.cos(beta)
				sb = np.sin(beta)
				
				clat = np.cos(HAWClat2)
				slat = np.sin(HAWClat2)
		
				# rotation matrix
				rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
				
				#rotation from Equatorial (ra,dec) to local frame
				vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
				vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
				vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
		
				j = H.vec2pix(nside,vpx,vpy,vpz)
			
				if Emap02[i] > 0.0 :
					temp1 += Nmap2[timeidx][i]
					temp2 += norm2[timeidx]*(CRmap[j]+diffCRmap[j])
			
			if temp2 > 0.0 :
				newEmap[i] = temp1/temp2
			else : 
				newEmap[i] = 0.0
		
		for i in range(0,npix) :
			Emap2[i] = newEmap[i]
			
		# cleaning
		
		temp = sum(Emap1)
			
		for timeidx in range(minstep,maxstep) :
			norm1[timeidx] = temp*norm1[timeidx]
			
		Emap1 = Emap1/temp	
		
		temp = sum(Emap2)
			
		for timeidx in range(minstep,maxstep) :
			norm2[timeidx] = temp*norm2[timeidx]
			
		Emap2 = Emap2/temp
			
		# write iteration results
		
		#H.mollview(diffCRmap)
		#show()
		
		H.write_map(namefits,diffCRmap/isovalue)	
		
		pickle.dump(Emap1,open(nameE1fits, "w" ))	
		pickle.dump(norm1, open(nameN1fits, "w" ))
		
		pickle.dump(Emap2,open(nameE2fits, "w" ))	
		pickle.dump(norm2, open(nameN2fits, "w" ))
		
		# calculate signifiance
		
		significancemap = np.zeros(npix, dtype=np.double)
		significancemap2 = np.zeros(npix, dtype=np.double)
		
		for timeidx in range(minstep,maxstep) :
		
		#print timeidx
		
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
				if Emap01[j]*norm01[timeidx] > 0.0 and Emap1[j]*norm1[timeidx] > 0.0 :
					significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap1[j]*norm1[timeidx] 
					significancemap[i] += +2.0*CRmap[i]*Emap01[j]*norm01[timeidx]
					temp1 = Emap1[j]/Emap01[j]*norm1[timeidx]/norm01[timeidx]
					significancemap[i] += 2.0*Nmap1[timeidx][j]*np.log(temp1*(1.0+diffCRmap[i]/CRmap[i]))
					
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
				if Emap02[j]*norm02[timeidx] > 0.0 and Emap2[j]*norm2[timeidx] > 0.0 :
					significancemap[i] += -2.0*(diffCRmap[i]+CRmap[i])*Emap2[j]*norm2[timeidx] 
					significancemap[i] += +2.0*CRmap[i]*Emap02[j]*norm02[timeidx]
					temp1 = Emap2[j]/Emap02[j]*norm2[timeidx]/norm02[timeidx]
					significancemap[i] += 2.0*Nmap2[timeidx][j]*np.log(temp1*(1.0+diffCRmap[i]/CRmap[i]))			
					

		H.write_map(nameSIGfits,significancemap)
		#H.mollview(significancemap)
		#show()
