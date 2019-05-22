# usage: python dipole-strength.py [list of iteration steps]
import healpy as h
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys

nside=64
npix=12*nside*nside
lmax=40
degree=np.pi/180.

def GetMaxMin(map):
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
        print "%.05f (%.02f, %.02f), %.05f (%.02f, %.02f)" % \
              (maxi,90.-maxang[0]/degree,maxang[1]/degree,mini,90.-minang[0]/degree,minang[1]/degree)
        return maxi,maxang,mini,minang

iterations = int(sys.argv[1])

Cltrue = np.zeros(lmax+1, dtype=np.double)
Cltrue[0] = 0.0 # no monopole
for i in range(1,lmax+1) :
        Cltrue[i] = 1e-7*18./(2.*i+1.)/(i+1.)/(i+2.) # or whatever you like
print Cltrue
mpl.rc("font", family="serif", size=14)
fig = plt.figure(1, figsize=(12,6))
ax = fig.add_subplot(111)

orientations =[]
Clavgs = []
for i in range(0,iterations):
               
    truediffCRmap =np.zeros(nside*nside*12,dtype=np.double)
    for j in range(0,lmax+1) :

            Cl = np.zeros(lmax+1, dtype=np.double)
            Cl[j] = 1.

            tempmap = h.sphtfunc.synfast(Cl,nside,lmax=lmax,new=True)
            Cltemp = h.anafast(tempmap,lmax=lmax)

            for k in range(0,npix) :
                    truediffCRmap[k] += tempmap[k]*np.sqrt(Cltrue[j]/Cltemp[j])

    almdipole = h.map2alm(truediffCRmap,lmax=1)
    crdipole  = h.alm2map(almdipole, nside, lmax=1)
    maxi, maxang,mini,minang = GetMaxMin(crdipole)
    orientations.append(90.-maxang[0]/degree)
 
    raAvgCRmap = truediffCRmap.copy() 
#    for ring in range(1,4*nside+1):
    for ring in range(int(4*nside*((90.-float(sys.argv[2]))/180.)),int((4*nside+1)*((90.-float(sys.argv[3]))/180.))):
        ringInfo = h.ringinfo(nside,np.array([ring]))
        startpix = ringInfo[0]
        ringpix = ringInfo[1]
        avg = sum(truediffCRmap[startpix:startpix+ringpix] / ringpix)
        for pix in range(startpix,startpix+ringpix):
            raAvgCRmap[pix]-=avg
    Clavg = h.anafast(raAvgCRmap,lmax=lmax)
    Clavgs.append(Clavg)
    ax.plot(range(0,41),Clavg,'ro',alpha=.25)
    Clavgs.append(Clavg)

Clavgs = np.array(Clavgs)

ax.plot(range(0,41),Cltrue,'k-',label='True Power')
ax.set_xlabel(r"multipole $\ell$")
ax.set_ylabel(r"$C_\ell$")
ax.set_xlim([0,40])
ax.set_ylim([1e-11,1e-5])
ax.set_yscale("log")
ax.legend(frameon=False,numpoints=1,loc="upper right")
ticks = [1., 4., 9., 18., 36.]
cl2deg = lambda cl: [(r"%g$^\circ$" % (180./c)) for c in cl]
at = ax.twiny()
at.set_xlim([0,40])
at.set_xticks(ticks)
at.set_xticklabels(cl2deg(ticks))
fig.tight_layout()
#plt.show()
plt.savefig("PowSpecsTheory_FOV_%s_to_%s.png"% (sys.argv[2],sys.argv[3]))

print orientations
fig2 = plt.figure(2, figsize=(12,6))
ax2 = fig2.add_subplot(111)
for l in [1,2,3,10]:
    Clplot = np.zeros(iterations)
    for j in range(0,iterations):
        tempCl = Clavgs[j]
        Clplot[j]=tempCl[l]/Cltrue[l]
    
    print Clplot
    ax2.plot((np.array(orientations)), Clplot, "*", label="Fraction of ell %d" %l)

ax2.set_xlabel("Dipole Orientation [deg]")
ax2.set_ylabel(r"$C_\ell$ observed / $C_\ell$ true")
ax2.set_xlim([-90,90])
ax2.set_ylim([0.,1.2])
ax2.legend(frameon=False,numpoints=1,loc="lower right")
plt.savefig("PowSpecsTheory_FOV_%s_to_%s_orientations.png"% (sys.argv[2],sys.argv[3]))


