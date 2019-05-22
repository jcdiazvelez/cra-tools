#!/usr/bin/python
from pylab import *
import matplotlib.pyplot as plt

import healpy as H
import sys
import os
import numpy as np
import pyfits
import pickle
import optparse
    
if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option("-e", "--events", dest="totevents",type=float,default=5e10)
    parser.add_option("-o", "--outdir", dest="outdir", default="./foo/")
    parser.add_option("--nside", dest="nside",type=int,default=64)
    parser.add_option("--lmax", dest="lmax",type=int,default=40)
    parser.add_option("--timesteps", dest="Ntimestep",type=int,default=360)
    parser.add_option("--lon", dest="lon",type=float,default=-97.0)
    parser.add_option("--lat", dest="lat",type=float,default=19.0)
    parser.add_option("--thetamax", dest="thetamax",type=float,default=60.0)
    parser.add_option("--ngaussians", dest="ngaussians",type=int, default=10)
    parser.add_option("--amplitude", dest="amplitude", type=float, default=1e-2)
    parser.add_option("--width", dest="width", type=float, default=90.)

    parser.add_option("--flat", action="store_true", dest="flat", default=False)
    opts, args = parser.parse_args()
    print opts

    # input parameters
    lmax = opts.lmax
    nside = opts.nside
    npix = H.nside2npix(nside)
    Ntimestep = opts.Ntimestep
    totevents = opts.totevents
    lon = opts.lon/180.*np.pi
    lat = opts.lat/180.*np.pi
    clat = np.cos(lat)
    slat = np.sin(lat)
    thetamax = opts.thetamax
    flat = opts.flat
    sourcefolder = opts.outdir #'/data/scratch/fiorino/hawc/local/gaussians-wholesky/'
    ngaussians = opts.ngaussians
    amplitude = opts.amplitude
    width = opts.width

    if not os.path.isdir(sourcefolder):
        os.makedirs(sourcefolder)
    lmap = np.zeros((lmax+1,npix), dtype=np.double)    
    truediffCRmap = np.zeros(npix, dtype=np.double)
    trueEmap =  np.zeros(npix, dtype=np.double)

    temp = 0.0
    for i in range(0,npix) :
        theta,phi = H.pix2ang(nside,i) 

        if theta/np.pi*180. < thetamax : # hard zenith cut
            if flat:
                trueEmap[i] = 1.
            else:
                trueEmap[i] = np.cos(theta)*(1.0 + 0.2*np.sin(theta)*np.sin(phi+10./180.*np.pi)**2) # typical zenith distribution
        temp +=trueEmap[i] 

    for i in range(0,npix) :
        trueEmap[i] = trueEmap[i]/temp
        
    H.write_map(sourcefolder + "/true_acceptance_map.fits.gz", trueEmap)    
    
    truenorm = []
    devents = totevents/(1.*Ntimestep)

    for timeidx in range(0,Ntimestep) :
        truenorm.append(devents*(1.0 + 0.1*np.sin((70+timeidx)/(1.*Ntimestep)*np.pi*2))) # sinusoidal rate in sidereal time
    
    pickle.dump(truenorm, open(sourcefolder + "/true_norm.dat", "w" ) )
    
    #exit(3)

    for g in range(0,ngaussians+1):

        gaussiandir = sourcefolder + '/gaussian{:02d}'.format(g)
        if not os.path.isdir(gaussiandir):
            os.makedirs(gaussiandir)
        # create signal
        truediffCRmap = np.zeros(npix)
        truediffCRmap[ H.ang2pix( nside, .5*np.pi*(1.- 1.*g/ngaussians), 0.) ] = 10.
        truediffCRmap = H.smoothing( truediffCRmap, sigma=width*np.pi/180. )
        truediffCRmap*= amplitude / max(truediffCRmap)

        H.write_map(sourcefolder + '/gaussian{:02d}'.format(g) +  "/true_relint_gaussian%02d.fits.gz" % g, truediffCRmap)    
        
        for timeidx in range(0,Ntimestep) :
        
            tempmap = np.zeros(npix, dtype=np.double)
            
            namefile = sourcefolder + '/gaussian{:02d}'.format(g) + '/local_timeindex' + '{:03d}'.format(timeidx) + '_gaussian{:02d}'.format(g) + '.fits.gz'
    
            beta = timeidx/(1.*Ntimestep)*np.pi*2 + lon
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
