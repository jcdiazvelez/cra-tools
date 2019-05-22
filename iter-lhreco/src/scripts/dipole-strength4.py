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
d1s = [0.,-30.,-60.-89.]
d2s = [60.]
cols     = ['k', 'b', 'r',]

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
ks    = np.zeros((len(d1s), len(d2s), nDipoles))
fracs = np.zeros((len(d1s), len(d2s), nDipoles))

# Loop over dipoles orientations
for dipole in range(0,nDipoles):
               
    almtrue = h.map2alm(np.zeros(12*nside*nside),lmax=lmax)
    almtrue[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)] = dipole
    almtrue[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=1)] = 10-dipole
    truemap = h.alm2map(almtrue,nside,lmax=lmax)
    truemap *= strength / max(truemap)

    print "========="
    print "Dipole %d" % dipole
    maxi,maxang,mini,minang = GetMaxMin(truemap)
    D[0,dipole] = 90.-maxang[0]/degree                                                 

    # Strength for RA average subtracted
    # Loop over possible detector latitudes
    for id1, d1 in enumerate(d1s):
      for id2, d2 in enumerate(d2s):

        print np.pi/2. - d2*degree
        print np.pi/2. - d1*degree
        pixLo = h.ang2pix(nside, np.pi/2.- d2*degree, 0.)
        pixHi = h.ang2pix(nside, np.pi/2.- d1*degree, 0.)

        out  = h.anafast(truemap,alm=True,lmax=lmax)       
        for i in range(0,lmax+1) :
           index = h.sphtfunc.Alm.getidx(lmax,i,0)
           out[1][index] = 0.0
        reducedmap = h.alm2map(out[1],nside,lmax=lmax)
        reducedmap[0:pixLo] = 0.
        reducedmap[pixHi:npix] = 0.
        fovPix = pixHi-pixLo

        # 'Fit' the dipole
        almreducedfit = h.map2alm(reducedmap,lmax=1)
        reducedfitmap = h.alm2map(almreducedfit,nside,lmax=1)
        maxi,maxang,mini,minang = GetMaxMin(reducedfitmap)

        corr = 16./(9*(np.sin(d2*degree)- np.sin(d1*degree)) + np.sin(3.*d2*degree) - np.sin(3.*d1*degree))
        print " Reduced Map (d1=%d, d2=%d) %.04f, %.04f" % (d1, d2, maxi / np.cos(degree*D[0,dipole]) / strength, 1./corr)

        ks[id1][id2][dipole] = 1./corr
        fracs[id1][id2][dipole] = maxi / np.cos(degree*D[0,dipole]) / strength
       
for id1, d1 in enumerate(d1s):
    p.plot( ks[id1][:][:], fracs[id1][:][:], label="d1=%d" % d1)
        

# Options
p.xlabel("Fraction of Recoverable Signal in Fit")
p.ylabel("K1111")
p.plot(np.arange(0,5),np.arange(0,5),'r--')
p.xlim(0,1.2)
p.ylim(0,1.2)
p.legend()

p.show()
