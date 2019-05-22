# usage: python dipole-strength2.py
import healpy as h
import numpy as np
import matplotlib.pyplot as p
import sys

degree=np.pi/180.

# Presets
nside    = 64
npix     = h.nside2npix(nside)
lmax     = 40
strength = 1e-3
nDipoles = 1
#nDipoles = 11
fovs     = [60., 90., 120., 150., 180.]
lats     = [0.,15.,30.,45.,60.,75.,90.]
cols     = ['k', 'b', 'r',]

n = []
m = []

def GetMaxMin(map):
    """return min/max value and angle for input map"""  
    mini=0.
    maxi=0.
    maxang=0.
    minang=0.
    for n in range(0,len(map)):
        if map[n]>maxi:
            maxi   = map[n]
            maxang = h.pix2ang(nside,n)
        if map[n]<mini:
            mini   = map[n]
            minang = h.pix2ang(nside,n)
    return maxi,maxang,mini,minang

iterations =[]
# Initialize data arrays
D     = np.zeros((1+len(iterations),nDipoles))
raAvg = np.zeros((len(lats), len(fovs), nDipoles))
covs  = np.zeros((len(lats), len(fovs), nDipoles))

# Loop over dipoles orientations
for dipole in range(0,nDipoles):
               
    almtemp = h.map2alm(np.zeros(12*nside*nside),lmax=lmax)
    almtemp[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)] = dipole
    almtemp[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=1)] = 10-dipole
    truemap = h.alm2map(almtemp,nside,lmax=lmax)
    truemap *= strength / max(truemap)

    print "========="
    print "Dipole %d" % dipole
    maxi,maxang,mini,minang = GetMaxMin(truemap)
    D[0,dipole] = 90.-maxang[0]/degree                                                 

    # Strength for RA average subtracted
    # Loop over possible detector latitudes
    for ilat, lat in enumerate(lats):
      for ifov, fov in enumerate(fovs):
        raAvgmap = truemap.copy() 
        raAvgFOVmap = truemap.copy() 

        decHi = lat+fov*.5 
        decLo = lat-fov*.5
        if decHi>=90.:
            decHi = 90.
        if decLo<=-90.:
            decLo = -90.
        pixLo = h.ang2pix(nside, np.pi/2.- decHi*degree,0.)
        pixHi = h.ang2pix(nside, np.pi/2.- decLo*degree,0.)

        out  = h.anafast(truemap,alm=True,lmax=lmax)       
        for i in range(0,lmax+1) :
           index = h.sphtfunc.Alm.getidx(lmax,i,0)
           out[1][index] = 0.0
        reducedmap = h.alm2map(out[1],nside,lmax=lmax)
        reducedmap[0:pixLo] = 0.
        reducedmap[pixHi:npix] = 0.
        fovPix = pixHi-pixLo

         
        cov = 1.*fovPix/len(truemap)
        covs[ilat][ifov][dipole] = cov

        # 'Fit' the dipole
        almtemp = h.map2alm(reducedmap,lmax=1)
        #almtemp = h.map2alm(raAvgFOVmap,lmax=1)
        tempmap = h.alm2map(almtemp,nside,lmax=1)
        maxi,maxang,mini,minang = GetMaxMin(tempmap)

        raAvg[ilat][ifov][dipole] = maxi

        corr = 16./(9*(np.sin(decHi*degree)- np.sin(decLo*degree)) + np.sin(3.*decHi*degree) - np.sin(3.*decLo*degree))
        print " RA average subtracted (lat=%d, fov=%d) %.04f, %.04f" % (lat, fov, maxi / np.cos(degree*D[0,dipole]) / strength, 1./corr)

        n.append(maxi / np.cos(degree*D[0,dipole]) / strength)
        m.append(1./corr)

# Options
p.xlabel("Fraction of Recoverable Signal in Fit")
p.ylabel("K1111")
p.plot(n,m,'ko')
p.plot(np.arange(0,5),np.arange(0,5),'r--')
p.xlim(0,1.2)
p.ylim(0,1.2)

p.show()
