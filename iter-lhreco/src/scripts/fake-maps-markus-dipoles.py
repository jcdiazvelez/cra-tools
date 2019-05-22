#!/usr/bin/python

import healpy as h
import sys
import numpy as n
from pylab import *
import matplotlib.pyplot as plt
import pyfits
import pickle

if __name__ == "__main__":
    
    # input parameters
    nside = 64
    npix = h.nside2npix(nside)
    Ntimestep = 360
    lmax = 40
    totevents = 5e10
    sourcefolder = '/data/scratch/fiorino/hawc/local/powspec-dipoles-wholesky/'
    HAWClon = -97.0/180.*n.pi
    HAWClat =  19.0/180.*n.pi
    clat = n.cos(HAWClat)
    slat = n.sin(HAWClat)

    lmap = n.zeros((lmax+1,npix), dtype=n.double)    
    truediffCRmap = n.zeros(npix, dtype=n.double)
    trueEmap =  n.zeros(npix, dtype=n.double)

    temp = 0.0
    for i in range(0,npix) :
        theta,phi = h.pix2ang(nside,i) 

        #if theta/n.pi*180. < 60.0 : # hard 60deg zenith cut
        #        trueEmap[i] = n.cos(theta)*(1.0 + 0.2*n.sin(theta)*n.sin(phi+10./180.*n.pi)**2) # or whatever you like

        # flat whole sky        
        trueEmap[i] = 1.
        temp +=trueEmap[i] 

    for i in range(0,npix) :
        trueEmap[i] = trueEmap[i]/temp
        
    h.write_map(sourcefolder + "true_acceptance_map.fits.gz", trueEmap)    
    
    truenorm = []
    devents = totevents/(1.*Ntimestep)
    for timeidx in range(0,Ntimestep) :
        truenorm.append(devents*(1.0 + 0.1*n.sin((70+timeidx)/(1.*Ntimestep)*n.pi*2))) # or whatever you like
    
    pickle.dump(norm, open(sourcefolder + "true_norm.dat", "w" ) )
    
    #exit(3)

    for g in range(0,11):

        almtemp = h.map2alm(n.zeros(len(truediffCRmap)),lmax=lmax)
        almtemp[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)]=g
        almtemp[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=1)]=10-g
        truediffCRmap=h.alm2map(almtemp,nside,lmax=lmax)
        truediffCRmap*= 1e-2 / max(truediffCRmap)
        almtemp=h.map2alm(truediffCRmap)

        Cldipole = h.anafast(truediffCRmap, lmax=lmax)
        dipolepower = Cldipole[1]
        Cltrue = n.zeros(len(Cldipole))
        Cltrue[0] = 0.0 # no monopole
        for i in range(1,lmax+1) :
                Cltrue[i] = 18*dipolepower/(2.*i+1.)/(i+1.)/(i+2.)
        psmap = h.synfast(Cltrue,64,lmax=lmax) 
        psalm = h.map2alm(psmap, lmax=lmax)
        
        psalm[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)]=almtemp[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)]
        psalm[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=1)]=almtemp[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=1)]
        truediffCRmap = h.alm2map(psalm,64,lmax=lmax)
         
        h.write_map(sourcefolder + "true_relint_dipole%02d.fits.gz" % g, truediffCRmap)    
        
        for timeidx in range(0,Ntimestep) :
        
            tempmap = n.zeros(npix, dtype=n.double)
            
            namefile = sourcefolder + 'local_timeindex' + '{:03d}'.format(timeidx) + 'dipole{:02d}'.format(g) + '.fits.gz'
    
            beta = timeidx/(1.*Ntimestep)*n.pi*2 + HAWClon
            cb = n.cos(beta)
            sb = n.sin(beta)
            
            # rotation matrix
            rot = [[cb*slat,sb*slat,-clat],[-sb,cb,0],[cb*clat,clat*sb,slat]]
            
            for i in range(0,npix) :
        
                vx,vy,vz = h.pix2vec(nside,i)
    
                #rotation from Equatorial (ra,dec) to local frame
                vpx = rot[0][0]*vx+rot[1][0]*vy+rot[2][0]*vz
                vpy = rot[0][1]*vx+rot[1][1]*vy+rot[2][1]*vz
                vpz = rot[0][2]*vx+rot[1][2]*vy+rot[2][2]*vz
            
                j = h.vec2pix(nside,vpx,vpy,vpz)
                
                mu = (1.0 + truediffCRmap[j])*trueEmap[i]*truenorm[timeidx]
                
                tempmap[i] = n.random.poisson(lam=mu)
    
            h.write_map(namefile,tempmap)    
            close()
