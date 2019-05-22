# usage: python gaussian-strength.py [list of iteration steps]
import healpy as h
import numpy as np
import matplotlib.pyplot as p
import sys

degree=np.pi/180.

# Presets
nside    = 64
lmax     = 40
strength = 1e-2
nGaussians = 11
fovs     = [60., 90., 120., 150., 180.]
lats     = [0.,45.,90.]
cols     = ['k', 'b', 'r', 'g']

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
    try:
        print "%.05f (%.02f, %.02f), %.05f (%.02f, %.02f)" % \
              (maxi,90.-maxang[0]/degree,maxang[1]/degree,mini,90.-minang[0]/degree,minang[1]/degree)
    except TypeError:
        print "..."
    return maxi,maxang,mini,minang

# Parse input
iterations = []
for it in sys.argv[1:]:
    iterations.append(int(it))
print iterations

# Initialize data arrays
D     = np.zeros((1+len(iterations),nGaussians))
raAvg = np.zeros((len(lats), len(fovs), nGaussians))
covs  = np.zeros((len(lats), len(fovs), nGaussians))
noM   = np.zeros(nGaussians)

# Loop over gaussians orientations
for iwidth, width in enumerate([10, 20, 30, 40]):
    for gaussian in range(1,nGaussians):
                   
        #truemap = h.read_map("/data/scratch/fiorino/hawc/local/gaussians-width%d-wholesky/gaussian%02d/true_relint_gaussian%02d.fits.gz" % (width,gaussian,gaussian) )
        truemap = h.read_map("/data/scratch/fiorino/hawc/local/gaussians-width%d-hawcsky/gaussian%02d/true_relint_gaussian%02d.fits.gz" % (width,gaussian,gaussian) )
    
        print "========="
        print "Gaussian %d Width %d" % (gaussian, width)
        maxi,maxang,mini,minang = GetMaxMin(truemap)
        D[0,gaussian] = 90.*gaussian/(nGaussians-1.)
    
        # Strength observed with method after i iterations
        # Loop over Niterations
        for j, iteration in enumerate(iterations):
            print "  Observed"
            #r=h.read_map("/data/scratch/fiorino/hawc/local/gaussian%02d-width%d-wholesky/CR_HAWC_64_360_iteration%d.fits" % (gaussian,width,iteration))
            r=h.read_map("/data/scratch/fiorino/hawc/local/gaussian%02d-width%d-wholesky_A10_relint.fits.gz" % (gaussian,width))
            maxi,maxang,mini,minang = GetMaxMin( r )
            D[j+1,gaussian] = maxi 
    
    #    # Strength for m=0 terms removed
    #    print "  m=0 terms removed"
    #    almtrue = h.map2alm(truemap,lmax=lmax)
    #    almtrue[h.sphtfunc.Alm.getidx(lmax=lmax,l=1,m=0)]=0
    #    maptemp = h.alm2map(almtrue, nside, lmax=lmax)
    #    maxi,maxang,mini,minang = GetMaxMin( maptemp )
    #    noM[gaussian] = maxi 
    
    #    # Strength for RA average subtracted
    #    # Loop over possible detector latitudes
    #    for ilat, lat in enumerate(lats):
    #      for ifov, fov in enumerate(fovs):
    #        raAvgmap = truemap.copy() 
    #        raAvgFOVmap = truemap.copy() 
    #
    #        decHi = lat+fov*.5 
    #        decLo = lat-fov*.5
    #        if decHi>=90.:
    #            decHi = 90.
    #        if decLo<=-90.:
    #            decLo = -90.
    #
    #        fovPix =0.
    #        ringLo = int( 4*nside*(90. - decHi)/180. )
    #        ringHi = int( 4*nside*(90. - decLo)/180. )
    #        for ring in range(1,4*nside+1):
    #            ringInfo = h.ringinfo(nside,np.array([ring]))
    #            startpix = ringInfo[0]
    #            ringpix  = ringInfo[1]
    #            avg = sum(truemap[startpix:startpix+ringpix] / ringpix)
    #            for pix in range(startpix,startpix+ringpix):
    #                raAvgmap[pix]-=avg
    #                # Limit FOV
    #                if ring in range(ringLo, ringHi+1):
    #                    raAvgFOVmap[pix] = raAvgmap[pix]
    #                    fovPix+=1
    #                else:
    #                    raAvgFOVmap[pix] = 0.
    #
    #         
    #        cov = 1.*fovPix/len(raAvgmap)
    #        covs[ilat][ifov][gaussian] = cov
    #        almtemp = h.map2alm(raAvgFOVmap,lmax=1)
    #        tempmap = h.alm2map(almtemp,nside,lmax=1)
    #        maxi,maxang,mini,minang = GetMaxMin(tempmap)
    #        raAvg[ilat][ifov][gaussian] = maxi
    #        print " RA average subtracted (lat=%d, fov=%d, coverage=%.02f, frac=%.02f)" % (lat, fov, cov, maxi / np.cos(degree*D[0,gaussian]) / strength)
    # Plot projections
    #for j, lat in enumerate(lats):
    #    p.plot(D[0], raAvg[j], 'bs', alpha=(.5 + .5*j/len(lats)), label='RA average subtracted (lat=%d)' % lat)
    #p.plot(D[0], noM, 'r.', label='m=0 terms removed')
    
    # Plot gaussian strengths
    for i, iteration in enumerate(iterations):
        #p.plot(D[0], D[i+1], '%so' % cols[iwidth], alpha=(.5 + .5*iteration/max(iterations)), label='Fit (it. %d, width %d)' % (iteration,width) )
        p.plot(D[0], D[i+1], '%so' % cols[iwidth], alpha=(.5 + .5*iteration/max(iterations)), label='Maximum (width %d)' % (width) )

# Plot theoretical lines
p.plot(D[0], np.ones(nGaussians)*strength, 'k-', label='True')

# Options
p.xlabel("Gaussian Orientation [deg]")
p.ylabel("Gaussian Strength")
p.legend(loc="upper right", numpoints=1)
p.xlim( min(D[0])*.9, max(D[0])*1.1 )
p.ylim( 0., strength*2. )

#for ilat, lat in enumerate(lats):
#  for ifov, fov in enumerate(fovs):
#    p.plot(covs[ilat][ifov], raAvg[ilat][ifov] / np.cos(degree*D[0])*strength , '%ss' % cols[ilat], label=' (lat=%d, fov=%d)' % (lat,fov)) 
#
#p.legend(loc="upper left", numpoints=1)

p.show()
