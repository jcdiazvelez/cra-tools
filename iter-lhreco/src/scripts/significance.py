#!/usr/bin/python

import healpy as H
import sys
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import pyfits
import pickle
import subprocess
from scipy.special import gammainc
from scipy.special import gammaincc
from scipy.special import gamma
from scipy.special import erfcinv
from scipy.special import erfinv
from matplotlib.colors import ListedColormap

               
def makegrid() :
	H.graticule(coord='C',dpar=30,dmer=30,color="black")
	H.projtext(np.pi/2.+0.02,np.pi/3.,r'$60^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='center',verticalalignment='top')
	H.projtext(np.pi/2.+0.02,np.pi*2./3.,r'$120^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='center',verticalalignment='top')
	H.projtext(np.pi/2.+0.02,np.pi,r'$180^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='center',verticalalignment='top')
	H.projtext(np.pi/2.+0.02,-np.pi/3.,r'$300^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='center',verticalalignment='top')
	H.projtext(np.pi/2.+0.02,-np.pi*2./3.,r'$240^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='center',verticalalignment='top')
	H.projtext(np.pi/2.*1/3,0.2,r'$60^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='right',verticalalignment='center')
	H.projtext(np.pi/2.*2/3,0.05,r'$30^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='right',verticalalignment='center')
	H.projtext(np.pi/2.*4/3,0.02,r'$-30^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='right',verticalalignment='center')
	H.projtext(np.pi/2.*5/3,0.02,r'$-60^\circ$',coord='C',color="black",fontsize=10,\
	rotation=0,horizontalalignment='right',verticalalignment='center')
	
if __name__ == "__main__":
	
	
	my_cmap = cm.jet
	my_cmap.set_under("w") # sets background to white
	
	rc('text', usetex=True)
	rc('font',**{'family':'serif','serif':['Times']})

	# resolution of input maps
	nsidehigh = 256
	npixhigh = H.nside2npix(nsidehigh)
	
	# resolution of output maps
	nside = 64
	npix = H.nside2npix(nside)

	
	isovalue = 1.0
	Ntimestep = 360
	minstep = 0
	maxstep = 360
	Niteration = 20
	#foldername = './HAWCcall_output64/'
	#foldername = '../Skymap/HAWC_output_full/'
        foldername = '/data/scratch/fiorino/hawc/local/markus/hawc-100/'
	
	CRmap = np.zeros(npix, dtype=np.double)
	
	#HAWC position
	HAWClon = -97.0/180.*np.pi
	HAWClat =  19.0/180.*np.pi

	ZENITHCUT = 45.0/180.*np.pi
	SMOOTHRAD = 10.0/180.*np.pi
	
	clat = np.cos(HAWClat)
	slat = np.sin(HAWClat)
	
	#initial CR distribution isotropic:
	for i in range(0,npix) :
		CRmap[i] = isovalue
	
	# choose the step from iteration method
	iteration  = 20
	
	namefits = foldername + 'CR_HAWC_64_360_iteration' + str(iteration) + '.fits.gz'
	nameEfits = foldername + 'exposure_HAWC_64_360_iteration0.fits.gz'
	namenorm = foldername + 'norm_HAWC_64_360_iteration' + str(iteration-1) + '.dat'
	
	diffCRmap = H.read_map(namefits) 
	Emap0 = H.read_map(nameEfits)

	norm = pickle.load( open(namenorm, "r" ) )
	
	Nmap = np.zeros((Ntimestep,npix), dtype=np.double)
	
	for timeidx in range(0,Ntimestep) :
		#namefile = '../Skymap/HAWCoriginal/hawclocal-v126-nhit30-run000630-001089_N256_' + '{:03d}'.format(timeidx) + '_of_360.fits.gz'
		namefile = '/data/scratch/fiorino/hawc/local/v126-nhit30/combined/hawclocal-v126-nhit30-run000630-001089_N256_' + '{:03d}'.format(timeidx) + '_of_360.fits.gz'
		
		#namefile = './HAWCtalk_nside64_360steps/local_timeindex' + '{:03d}'.format(timeidx) + '.fits'
		tempmap = H.read_map(namefile)
		Nmap[timeidx] = H.ud_grade(tempmap,nside)*(npixhigh/(1.*npix))
	
	
	tempmap = np.zeros(npix, dtype=np.double)
	for j in range(0,npix) :
		theta,phi = H.pix2ang(nside,j)

		if theta>np.pi/2.-HAWClat-ZENITHCUT and theta<np.pi/2.-HAWClat+ZENITHCUT :
			tempmap[j] = diffCRmap[j]	
		else :
			tempmap[j] = H.UNSEEN
	
	# substract l<=3
	out = H.anafast(tempmap,alm=True,lmax=3)
	submap = H.alm2map(out[1],nside,lmax=3)
		
	for j in range(0,npix) :
		diffCRmap[j] = diffCRmap[j] - submap[j]
		
	mumap = np.zeros(npix, dtype=np.double)
	mu0map = np.zeros(npix, dtype=np.double)
	Nmapbin = np.zeros(npix, dtype=np.double)
	significancemap = np.zeros(npix, dtype=np.double)
	
	for timeidx in range(minstep,maxstep) :

		print timeidx
		
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
			
			mu0map[i] += (CRmap[i]+submap[i])*Emap0[j]*norm[timeidx]
			mumap[i]  += (diffCRmap[i]+submap[i]+CRmap[i])*Emap0[j]*norm[timeidx] 
			Nmapbin[i] += Nmap[timeidx][j]
	
	smoothsignificancemap = np.zeros(npix, dtype=np.double)
	smoothdiffCRmap = np.zeros(npix, dtype=np.double)
	
	for j in range(0,npix) :
	
		print j
		
		vec = H.pix2vec(nside,j)
		set = H.query_disc(nside,vec,SMOOTHRAD)
		mutemp = 0.0
		mu0temp = 0.0
		Nmaptemp = 0.0
		for k in set :
			mutemp += mumap[k]
			mu0temp += mu0map[k]
			Nmaptemp += Nmapbin[k]
			smoothdiffCRmap[j] += diffCRmap[k]/len(set)
			
		if mutemp > 0.0 and mu0temp > 0.0 :	
			smoothsignificancemap[j] = -2.0*mutemp+2.0*mu0temp+2.0*Nmaptemp*np.log(mutemp/mu0temp)
			
	for j in range(0,npix) :
				
		if smoothsignificancemap[j] > 0.0 :
			smoothsignificancemap[j] = np.sign(smoothdiffCRmap[j])*np.sqrt(smoothsignificancemap[j])
		else :
			smoothsignificancemap[j] = 0.0
		
	mymax=max(np.amax(smoothsignificancemap),-np.amin(smoothsignificancemap))	
	H.mollview(smoothsignificancemap,title=r'significance (units of $\sigma$), iteration ' + str(iteration) ,coord='C',cbar=True,rot=180,notext=True,max=mymax,min=-mymax)
		
	makegrid()
	show()

	tempmap = np.zeros(npix, dtype=np.double)
	for j in range(0,npix) :
		theta,phi = H.pix2ang(nside,j)

		if theta>np.pi/2.-HAWClat-ZENITHCUT and theta<np.pi/2.-HAWClat+ZENITHCUT :
			tempmap[j] = smoothdiffCRmap[j]
		else :
			tempmap[j] = H.UNSEEN
			
	H.mollview(tempmap,title=r'relative intensity, iteration ' + str(iteration) ,coord='C',cbar=True,rot=180,notext=True)
		
	makegrid()
	show()
		

