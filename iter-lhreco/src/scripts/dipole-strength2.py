# usage: python dipole-strength2.py
import healpy as h
import numpy as np
import matplotlib.pyplot as p
import sys

degree=np.pi/180.

# Presets
nside    = 64
lmax     = 40
strength = 1e-3
nDipoles = 11
#nDipoles = 11
fovs     = [90.,120.]#[60., 90., 120., 150., 180.]
lats     = [19.]#[0.,45.,90.]
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
#    try:
#        print "%.05f (%.02f, %.02f), %.05f (%.02f, %.02f)" % \
#              (maxi,90.-maxang[0]/degree,maxang[1]/degree,mini,90.-minang[0]/degree,minang[1]/degree)
#    except TypeError:
#        print "..."
    return maxi,maxang,mini,minang

iterations =[]
# Initialize data arrays
D     = np.zeros((1+len(iterations),nDipoles))
raAvg = np.zeros((len(lats), len(fovs), nDipoles))
covs  = np.zeros((len(lats), len(fovs), nDipoles))
noM   = np.zeros(nDipoles)

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

#    # Strength for m=0 terms removed
#    print "  m=0 terms removed"
#    almtrue = h.map2alm(truemap,lmax=lmax)
#    almtrue[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)]=0
#    maptemp = h.alm2map(almtrue, nside, lmax=lmax)
#    maxi,maxang,mini,minang = GetMaxMin( maptemp )
#    noM[dipole] = maxi 
#    p.plot(D[0], noM, 'r.', label='m=0 terms removed')
#
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

        fovPix =0.
        ringLo = int( 4*nside*(90. - decHi)/180. )
        ringHi = int( 4*nside*(90. - decLo)/180. )
        for ring in range(1,4*nside+1):
            ringInfo = h.ringinfo(nside,np.array([ring]))
            startpix = ringInfo[0]
            ringpix  = ringInfo[1]
            avg = sum(truemap[startpix:startpix+ringpix] / ringpix)
            for pix in range(startpix,startpix+ringpix):
                raAvgmap[pix]-=avg
                # Limit FOV
                if ring in range(ringLo, ringHi+1):
                    raAvgFOVmap[pix] = raAvgmap[pix]
                    fovPix+=1
                else:
                    raAvgFOVmap[pix] = 0.

         
        cov = 1.*fovPix/len(raAvgmap)
        covs[ilat][ifov][dipole] = cov

        almtemp = h.map2alm(raAvgFOVmap,lmax=1)
        tempmap = h.alm2map(almtemp,nside,lmax=1)
        maxi,maxang,mini,minang = GetMaxMin(tempmap)

        raAvg[ilat][ifov][dipole] = maxi

        a = degree*(90.-decHi)
        b = degree*(90.-decLo)
        corr = 16./(9.*(np.cos(a) - np.cos(b)) - np.cos(3.*a) + np.cos(3.*b))
        print " RA average subtracted (lat=%d, fov=%d, coverage=%.04f, frac=%.04f), %.04f, %.04f" % (lat, fov, cov, maxi / np.cos(degree*D[0,dipole]) / strength, maxi, corr)

# Plot theoretical lines
p.plot(D[0], np.ones(nDipoles)*strength, 'k-', label='True')
p.plot(D[0], np.cos(degree*D[0])*strength, 'b-', label='RA Projection')

# Options
p.xlabel("Dipole Orientation [deg]")
p.ylabel("Dipole Strength")
p.legend(loc="upper right", numpoints=1)
p.xlim( min(D[0])*.9, max(D[0])*1.1 )
p.ylim( 0., strength*2. )

for ilat, lat in enumerate(lats):
  for ifov, fov in enumerate(fovs):

    # Correction
    thetaHi = 90. - (lat+fov/2.)
    if thetaHi < 0.:
        thetaHi = 0.
    thetaLo = 90. - (lat-fov/2.)
    if thetaLo > 180.:
        thetaLo = 180.
    corr = 16./(9*(np.cos(thetaHi*degree)- np.cos(thetaLo*degree)) - np.cos(3.*thetaHi*degree) + np.cos(3.*thetaHi*degree))
    print corr, lat, fov
    #p.plot(covs[ilat][ifov], raAvg[ilat][ifov] / np.cos(degree*D[0])*strength , '%ss' % cols[ilat], label=' (lat=%d, fov=%d)' % (lat,fov)) 
    #p.plot(covs[ilat][ifov], raAvg[ilat][ifov]*corr, '%ss' % cols[ilat], label=' (lat=%d, fov=%d)' % (lat,fov)) 
    p.plot(D[0], raAvg[ilat][ifov]*corr, '%ss' % cols[ilat], label=' (lat=%d, fov=%d)' % (lat,fov)) 

p.legend(loc="upper right", numpoints=1)

p.show()
