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
#import colormaps as cmaps
               
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
	
	
	#my_cmap = cm.Spectral_r
	my_cmap = cm.RdBu_r
	my_cmap.set_under("w") # sets background to white
	
	rc('text', usetex=True)
	rc('font',**{'family':'serif','serif':['Times']})

	# resolution of input maps
	nsidehigh = 64
	npixhigh = H.nside2npix(nsidehigh)
	
	# resolution of output maps
	nside = 64
	npix = H.nside2npix(nside)

	
	isovalue = 1.0
	Ntimestep = 360
	minstep = 0
	maxstep = 360
	Niteration = 20
	foldername = './paper1_newfreq_180/'
	#foldername = '../Skymap/HAWC_output_full/'
	foldername = 'HAWCv7.1_and_IceCube_v2_all_IC86.1_small/'
	
	CRmap = np.zeros(npix, dtype=np.double)
	
	#HAWC position
	HAWClon = -97.0/180.*np.pi
	HAWClat =  19.0/180.*np.pi

	ZENITHCUT = 60.0/180.*np.pi
	SMOOTHRAD = 10.0/180.*np.pi
	
	clat = np.cos(HAWClat)
	slat = np.sin(HAWClat)
	
	#initial CR distribution isotropic:
	for i in range(0,npix) :
		CRmap[i] = isovalue
	
	# choose the step from iteration method
	iteration  = 20
	
	namefits = foldername + 'CR_64_360_iteration' + str(iteration) + '.fits'
	nameEfits = foldername + 'exposure_64_360_iteration' + str(iteration) +'.dat'
	namenorm = foldername + 'norm_64_360_iteration' + str(iteration) + '.dat'
	
	diffCRmap = H.read_map(namefits) 
	norm = np.array(pickle.load(open(namenorm, "r" )))[0]
	Emap0 = np.array(pickle.load(open(nameEfits, "r" )))[0]
	
	Nmap = np.zeros((Ntimestep,npix), dtype=np.double)
	
	for timeidx in range(0,Ntimestep) :
		#namefile = '../Skymap/HAWCoriginal/hawclocal-v126-nhit30-run000630-001089_N256_' + '{:03d}'.format(timeidx) + '_of_360.fits.gz'
		
		namefile = './paper1_newfreq_nside64_360steps/local_timeindex' + '{:03d}'.format(timeidx) + '.fits'
		tempmap = H.read_map(namefile)
		Nmap[timeidx] = H.ud_grade(tempmap,nside)*(npixhigh/(1.*npix))
	
	
	tempmap = np.zeros(npix, dtype=np.double)
	for j in range(0,npix) :
		theta,phi = H.pix2ang(nside,j)

		if theta>np.pi/2.-HAWClat-ZENITHCUT and theta<np.pi/2.-HAWClat+ZENITHCUT :
			tempmap[j] = diffCRmap[j]	
		else :
			tempmap[j] = 0.0 #H.UNSEEN
	
	# substract l<=3
	subLMAX = 3
	redLMAX = 3
	out = H.anafast(tempmap,alm=True,lmax=subLMAX)
	outr = H.anafast(tempmap,alm=True,lmax=redLMAX)	
	submap = H.alm2map(out[1],nside,lmax=subLMAX)
	
	#print out[1]
	
	"""
	#ALMmatrix = pickle.load( open("../inverse/ALMmatrix_paper1_nside64_LMAX10_60deg_tophat.dat", "r" ) )	
	ALMmatrix = pickle.load( open("../inverse/ALMmatrix_final_LMAX30_60deg_tophat.dat", "r" ) )	
	ALMmatrixred = ALMmatrix[:(redLMAX+1)*redLMAX,:(redLMAX+1)*redLMAX]
	
	
	inverse = np.linalg.inv(ALMmatrixred)
	
	mymax=max(np.amax(ALMmatrixred),-np.amin(ALMmatrixred))
	imshow(ALMmatrixred,interpolation='nearest',vmin=-mymax,vmax=mymax,cmap= my_cmap)
	ticks = []
	emp = []
	start = -0.5
	for l in range(1,redLMAX+1) :
		start += 2*l
		ticks.append(start)
		emp.append('')
	plt.xticks(ticks,emp)
	plt.yticks(ticks,emp)
	grid(True)
	colorbar()
	show()
	
	mymax=max(np.amax(inverse),-np.amin(inverse))
	imshow(inverse,interpolation='nearest',vmin=-mymax,vmax=mymax,cmap= my_cmap)
	ticks = []
	emp = []
	start = -0.5
	for l in range(1,redLMAX+1) :
		start += 2*l
		ticks.append(start)
		emp.append('')
	plt.xticks(ticks,emp)
	plt.yticks(ticks,emp)
	grid(True)
	colorbar()
	show()
	
	#exit(3)
	
	newalm = []
	for i in range(0,len(out[1])) :
		newalm.append(0+0j)
		
	index1 = 0
	for l1 in range(1,subLMAX+1) :
		for m1 in range(-l1,l1+1) :
			if m1 in [0] :
				continue
			
			index2 = 0
			for l2 in range(1,redLMAX+1) :
				for m2 in range(-l2,l2+1) :
					if m2 in [0] :
						continue
				
					if m2 > 0  and m1 > 0 :
						hi1 = H.sphtfunc.Alm.getidx(subLMAX,l1,m1)
						hi2 = H.sphtfunc.Alm.getidx(redLMAX,l2,m2)
						alm = outr[1][hi2]
						#print l1,m1,l2,m2,index1,index2,inverse[index1,index2],alm
						newalm[hi1] += inverse[index1,index2]*alm
					index2 += 1
			index1 += 1
		
	print out[1]	
	print np.array(newalm)
	submap2 = H.alm2map(np.array(newalm),nside,lmax=subLMAX)
	"""
	
	H.mollview(submap,title=r'large-scale anisotropy (pseudo, $\ell \leq 3$) / iteration ' + str(iteration) ,coord='C',cbar=True,rot=180,notext=True,min=-0.0005,max=0.0005,cmap=my_cmap)
		
	makegrid()
	show()
	
	"""
	H.mollview(submap2,title=r'large-scale anisotropy (truncated, $\ell \leq 3$) / iteration ' + str(iteration) ,coord='C',cbar=True,rot=180,notext=True,min=-0.0005,max=0.0005,cmap=my_cmap)
		
	"""
	makegrid()
	show()
	
	#exit(3)
	#submap = submap2
	
	for j in range(0,npix) :
		#diffCRmap[j] = diffCRmap[j] - submap2[j]
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
		totset = 0.0
		for k in set :
		
			if mu0map[k] > 0.0 : 
				mutemp += mumap[k]
				mu0temp += mu0map[k]
				Nmaptemp += Nmapbin[k]
				totset += 1.0
				smoothdiffCRmap[j] += diffCRmap[k] # /len(set)
			
		smoothdiffCRmap[j] = smoothdiffCRmap[j]/totset
			
		if mutemp > 0.0 and mu0temp > 0.0 :	
			smoothsignificancemap[j] = -2.0*mutemp+2.0*mu0temp+2.0*Nmaptemp*np.log(mutemp/mu0temp)
			
	for j in range(0,npix) :
		theta,phi = H.pix2ang(nside,j)
		
		#if theta>np.pi/2.-HAWClat-50./180.*np.pi and theta<np.pi/2.-HAWClat+50./180.*np.pi :
		if  theta<np.pi/2.-HAWClat+60./180.*np.pi :
			if smoothsignificancemap[j] > 0.0 : pass

			else :
				smoothsignificancemap[j] = 0.0
		else :	
			smoothsignificancemap[j] = H.UNSEEN
		
	mymax=max(np.amax(smoothsignificancemap),-np.amin(smoothsignificancemap))	
	mymax=9.0
	H.mollview(smoothsignificancemap,title=r'significance (units of $\sigma$, $10^\circ$ smoothed) / iteration ' + str(iteration) ,coord='C',cbar=True,rot=180,notext=False,max=mymax,min=-mymax,cmap=my_cmap)
		
	makegrid()
	show()

	tempmap = np.zeros(npix, dtype=np.double)
	for j in range(0,npix) :
		theta,phi = H.pix2ang(nside,j)

		if theta>np.pi/2.-HAWClat-50./180.*np.pi and theta<np.pi/2.-HAWClat+50./180.*np.pi :
			tempmap[j] = smoothdiffCRmap[j]
		else :
			tempmap[j] = H.UNSEEN
			
	H.mollview(tempmap,title=r'small-scale anisotropy ($\ell \leq 3$ removed, $10^\circ$ smoothed) / iteration ' + str(iteration) ,coord='C',cbar=True,rot=180,notext=False,min=-0.0004,max=0.0004,cmap=my_cmap)
		
	makegrid()
	show()
		
