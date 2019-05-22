import healpy
import healpy as hp
import math
import numpy
import numpy as np
import os
import os.path
import ephem
import pylab
import copy
import astropy.units as u
from astropy.coordinates import Galactic
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from numpy import log, sqrt, cos, pi


degree = math.pi / 180
default_binsize = 20


def kcorrection(bottom=-90*degree,top=90*degree):
    d1=bottom
    d2=top
    K1111 =.25*(3*(np.sin(d2) - np.sin(d1)) + pow(np.sin(d1),3) - pow(np.sin(d2),3))
    return K1111


def maskIceCube(skymap):
    """Set pixels outside the IceCube mask region to UNSEEN"""
    npix = skymap.size
    nside = healpy.npix2nside(npix)

    minDec = -89 * degree
    maxDec = -26 * degree
    maxDec = 0 * degree
    minPix = healpy.ang2pix(nside, 90 * degree - maxDec, 0.)
    maxPix = healpy.ang2pix(nside, 90 * degree - minDec, 0.)

    skymap[0:minPix] = healpy.UNSEEN
    skymap[maxPix:npix] = healpy.UNSEEN

    return skymap

def maskCircle(skymap,center0=59*degree,center1=-9*degree,radius=15*degree):
    """Set pixels outside the IceCube mask region to UNSEEN"""
    npix = skymap.size
    nside = healpy.npix2nside(npix)
    masked_skymap = copy.deepcopy(skymap)
    for ipix in xrange(npix):
        theta, phi = healpy.pix2ang(nside, ipix)
        dtheta = 90*degree - theta - center1
        dphi = phi - center0
        rad = dtheta*dtheta + dphi*dphi
        if rad < radius*radius:
           masked_skymap[ipix] = healpy.UNSEEN

    return masked_skymap




def mkMask(skymap, minDec=-85 * degree, maxDec=70 * degree, mask_value=healpy.UNSEEN):
    """Set pixels outside the IceCube mask region to UNSEEN"""
    npix = skymap.size
    nside = healpy.npix2nside(npix)
    masked_skymap = copy.deepcopy(skymap)

    minPix = healpy.ang2pix(nside, 90 * degree - maxDec, 0.)
    maxPix = healpy.ang2pix(nside, 90 * degree - minDec, 0.)
    masked_skymap[minPix:maxPix] = mask_value

    return masked_skymap




def icecube_mask(imap, top=-20 * degree, bottom=-80 * degree, mask_value=healpy.UNSEEN):
    top_mask = mkMask(imap, top, 90 * degree, mask_value=mask_value)
    return mkMask(top_mask, -90 * degree, bottom, mask_value=mask_value)


def hawc_mask(imap, top=70 * degree, bottom=-30 * degree, mask_value=healpy.UNSEEN):
    top_mask = mkMask(imap, top, 90 * degree, mask_value=mask_value)
    return mkMask(top_mask, -90 * degree, bottom, mask_value=mask_value)


def mask(imap, top=80 * degree, bottom=-80 * degree, mask_value=healpy.UNSEEN):
    top_mask = mkMask(imap, top, 90 * degree, mask_value=mask_value)
    return mkMask(top_mask, -90 * degree, bottom, mask_value=mask_value)


def plot_chi2(signals, backgrounds, legends, title='1D Projection Sidereal-time (all-sky)', binsize=10):
    from scipy.stats import chi2
    plt.figure(figsize=(15, 8))
    ax = plt.subplot(211)
    xchi2 = []

    x0, y0, yerr0 = get_1d_bin(signal[0], background[1])
    x1, y1, yerr1 = get_1d_bin(signal[1], background[1])
    xchi2 = map(chi2, zip(x0, x1))
    p = plt.plot(x * binsize, xchi2, fmt='o', label="$\chi^2$")

    plt.legend(loc='lower right', numpoints=1)
    plt.xlim(plt.xlim()[::-1])
    ax.set_title(title)
    ax.set_xlabel('RA (deg)')
    plt.show()





def generate_map(bgMap, dataMap):
  import numpy.random
  alpha = 0.1;

  npix = dataMap.size
  signalMap = np.zeros(npix)
  for i in range(npix):
    # Throw a Poisson random number using the data as a mean
    Nd = numpy.random.poisson(bgMap[i])
    Nb = bgMap[i]
    if Nd > 0 and Nb > 0:
      #signalMap[i] = (Nd - alpha*Nb) / numpy.sqrt(Nd + alpha*Nb)
      #signalMap[i] = (Nd - alpha*Nb) / (alpha*Nb+1e-16)
      signalMap[i] = Nd 
      #signalMap[i] = (Nd - alpha*Nb) / (alpha*Nb+1e-16)
  return signalMap



def isotropic_band(dataMap, bgMap=None, n=100, lmax=40, percentile=90, 
                       isomap_dir=None,lhiter=10, 
                       mask_top=70 * degree, mask_bottom=-80 * degree):
    from tqdm import tqdm
    import contextlib
    print "generating isotropic band"
    @contextlib.contextmanager
    def capture():
        import sys
        from cStringIO import StringIO
        oldout,olderr = sys.stdout, sys.stderr
        try:
            out=[StringIO(), StringIO()]
            sys.stdout,sys.stderr = out
            yield out
        finally:
            sys.stdout,sys.stderr = oldout, olderr
            out[0] = out[0].getvalue()
            out[1] = out[1].getvalue()


    from tqdm import tqdm
    with capture() as out:
       print 'hi'
       cls = []
       p0 = (100 - percentile)/2
       for c in xrange(lmax+1):
          cls.append([])
       clmins = numpy.ones(lmax+1)*1e6
       clmaxs = numpy.zeros(lmax+1)
       for it in tqdm(range(n)):
           if isomap_dir:
              iso_relint = healpy.read_map(os.path.join(isomap_dir,str(it),"CR_32_360_iteration%u.fits" % lhiter),0)
           elif bgMap is not None:
              iso = generate_map(bgMap, dataMap)
              iso_relint = (iso-bgMap)/(bgMap+1e-16) 
           else:
              raise "No background map provided!"

           iso_ccl = hp.anafast(mask(iso_relint, top=mask_top, bottom=mask_bottom))
           for c in xrange(lmax+1):
              cls[c].append(iso_ccl[c])

    for c in xrange(lmax+1):
      pmin, pmax = numpy.percentile(cls[c], [p0,100-p0])
      clmins[c]=pmin
      clmaxs[c]=pmax

    return clmins,clmaxs



def error_bars(dataMap, bgMap=None, n=100, lmin=0, lmax=40, percentile=90, poissonmap_dir=None,lhiter=10, do_chisq=False):
    import scipy.stats
    from tqdm import tqdm
    import contextlib
    print "error bars"
    @contextlib.contextmanager
    def capture():
        import sys
        from cStringIO import StringIO
        oldout,olderr = sys.stdout, sys.stderr
        try:
            out=[StringIO(), StringIO()]
            sys.stdout,sys.stderr = out
            yield out
        finally:
            sys.stdout,sys.stderr = oldout, olderr
            out[0] = out[0].getvalue()
            out[1] = out[1].getvalue()


    from tqdm import tqdm
    with capture() as out:
       print 'hi'
       cls = []
       p0 = (100 - percentile)/2
       for c in xrange(lmax+1):
          cls.append([])
       yerr = numpy.zeros(lmax+1)
       for it in tqdm(range(n)):
           if poissonmap_dir:
              poisson_relint = healpy.read_map(os.path.join(poissonmap_dir,str(it),"CR_32_360_iteration%u.fits" % lhiter),0)
           elif bgMap:
              poisson = generate_map(bgMap, dataMap)
              poisson_relint = (poisson-bgMap)/(bgMap+1e-16) 
           else:
              raise "No background map provided!"
           if lmin > 1:
              alms = healpy.map2alm(poisson_relint, lmax=lmin-1 )
              npix = len(poisson_relint)
              nside = healpy.npix2nside(npix)
              poisson_relint -= healpy.alm2map(alms, nside)
              

           #poisson_ccl = hp.anafast(poisson_relint, lmax=lmax)
           poisson_ccl = hp.anafast(mask(poisson_relint, top=70 * degree, bottom=-80 * degree))
           for c in xrange(lmax+1):
              cls[c].append(poisson_ccl[c])

    chisq = []
    for c in xrange(lmax+1):
      sigma = numpy.std(cls[c])
      mu = numpy.mean(cls[c])
      yerr[c] = sigma
      if do_chisq:
         hist, bin_edges = np.histogram(np.array(cls))
         pdf = 1./(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bin_edges[1:] - mu)**2 / (2 * sigma**2) )
         chisq.append(scipy.stats.chisquare(hist, f_exp=pdf, ddof=0, axis=0))
      


    if lmin > 1:
       print "Using subtracted multipole l=",lmin-1

    return yerr, chisq




def TopHatSmooth(nside,origMap,radius,average=False):
  pixelList = []
  dsum = 0.0
  newmap = numpy.zeros(len(origMap))

  # Loop over all pixels
  for ipix in range(0, len(origMap)):
    theta,phi = healpy.pix2ang(nside, ipix)
    center = healpy.dir2vec(theta,phi)

    # Grab pixels with given angular radius of ith pixel
    pixelList = healpy.query_disc(nside, center, radius)

    #Sum up pixels in the angular cap; set value of ith pixel to this sum
    dsum = 0.
    nsum = 0.
    for jpix in pixelList:
      if origMap[jpix] != healpy.UNSEEN:
         dsum += origMap[jpix]
         nsum += 1.0

    if average and nsum>0:
        newmap[ipix] = dsum/nsum # set to average
    else:
        newmap[ipix] = dsum
  return newmap


# Maps with celestial features

def SetupAbsThresholdColormap(amin, amax, threshold):
        """ Create a color map for "two-sided" thresholds.  Below the threshold,
            the map is a cool green-blue palette.  Between the lower and upper
            threshold, the map is gray-white-gray.  Above the upper threshold,
            the map is a warm red-yellow palette.
        """
        x1 = (-threshold - amin) / (amax - amin)
        x3 = (amax - threshold) / (amax - amin)
        x2 = 1. - x1 - x3
        gvl = 0.5
        thrDict = {
            "red"   : ((0.0, 1.0, 0.5), (x1, 0.0, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                        (x1 + x2, gvl, 0.7), (1.0, 1.0, 1.0)),
            "green" : ((0.0, 1.0, 1.0), (x1, 0.0, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                        (x1 + x2, gvl, 0.0), (1.0, 1.0, 1.0)),
            "blue"  : ((0.0, 1.0, 1.0), (x1, 0.7, gvl), (x1 + 0.5*x2, 1.0, 1.0),
                        (x1 + x2, gvl, 0.0), (1.0, 0.5, 1.0)) }
        return LinearSegmentedColormap("thresholdColormap", thrDict, 256)

def set_color_bar(title=None, ticks=None, coord="C"):
    """Create the color bar for a HEALPix figure
    """
    fig = pylab.figure(1)
    for ax in fig.get_axes():
        if type(ax) is hp.projaxes.HpxMollweideAxes:
            tickLabels = []
            if ticks:
                tickLabels = [TickToLabel(t) for t in ticks]

            cb = fig.colorbar(ax.get_images()[0], ax=ax,
                              orientation="horizontal",
                              shrink=0.8, aspect=35,
                              pad=0.05, fraction=0.1,
                              ticks=ticks)
            cb.ax.set_xticklabels(tickLabels)
            cb.set_label(title)

    

def plot_map_with_features(skymap,title,proj='C',gc=True,label='', show_ibex=False,
                          fitvalue=37.8421671285,CG_corr=None,lang='eng',filename=None,watermark="", 
                           center_on_fit=False,thresh=None,dMin=None,dMax=None):
    
    limf= ephem.Ecliptic(227*degree,34.62*degree)
    limferr = [5.7,4.4]
    limferr_vec = [[],[]]
    for th in range(0,360,30):
        limferr_vec[0].append(227+limferr[0]*np.cos(th*degree))
        limferr_vec[1].append(34.62+limferr[1]*np.sin(th*degree))
    
    colormap = pylab.get_cmap("bwr")
    if thresh != None:
        colormap = SetupAbsThresholdColormap(dMin, dMax, thresh)
        colormap.set_under("w")    
            
        
    ibex = ephem.Ecliptic(219.2*degree,39.9*degree)

    glimf = ephem.Galactic(limf)
    gibex = ephem.Galactic(ibex)
    climf = ephem.Equatorial(limf)
    
    cibex = ephem.Equatorial(ibex)
    print 'LIMF',climf.ra/degree, climf.dec/degree
    print 'LIMF',climf.ra,climf.dec

    
    
    

    notext=False
    if proj=='C0':
        rotation = (0,0,0)
    elif proj=='C':
        rotation = (-180,0,0)
    elif proj=='LIMF':
        rotation = (climf.ra/degree,climf.dec/degree)
        notext=True
    elif proj=='ANTILIMF':
        rotation = (climf.ra/degree-180,climf.dec/degree-34.62)
        rotation = (climf.ra/degree-180,-climf.dec/degree) 
        notext=True
    else:
        rotation = (0,0,0)
    hp.mollview(mask(skymap,top=67*degree,bottom=-90*degree), 
            fig=1,
            title=title, 
            rot=rotation,coord=['C'],
            unit=label,
            notext=notext, cmap = colormap,min=dMin,max=dMax)
    #set_color_bar(title="Neutral Hydrogen Flux [1/(cm$^2$ s sr keV)]", ticks=None, coord="C")
    fig = pylab.figure(1)

    x = [227.28]
    y = [34.62]
    x =[0]
    y=[-89]
 
    labelsmap = {}
    labelsmap['eng'] = {'c-g': 'G-C', 'gc':'GC', 'ibex_center':'IBEX Ribbon','v_ism':r'$v_{ISM}$',
                        'v_lsr':r'$v_{LSR}$','limf':r'$\vec{B}_\mathrm{LIMF}$'}
    labelsmap['esp'] = {'c-g': 'C-G', 'gc':'Cent. Galact.', 'ibex_center':'Cinta IBEX','v_ism':r'$v_{MIS}$',
                        'v_lsr':r'$v_{LSR}$', 'limf':r'$\vec{B}_\mathrm{LIMF}$'}
    labels = labelsmap[lang]


    print glimf.lat, glimf.lon

    hp.projscatter(x, y, lonlat=True,color='black', coord='G',rot=(glimf.lon/degree,glimf.lat/degree))
    #hp.projscatter(limferr_vec[0], limferr_vec[1], lonlat=True,color='grey', coord='E')

    hp.projscatter([0,263.5520,195.1338,071.3350], [0, -02.7873,04.2658,03.0668], lonlat=True, coord='G')
    hp.projtext(195.1338, 04.2658, 'Geminga', lonlat=True, coord='G')
    hp.projtext(263.5520, -02.7873, 'Vela', lonlat=True, coord='G')
    hp.projtext(071.3350, 03.0668, 'Cygnus-X',lonlat=True, coord='G')
    hp.projtext(0, 0, 'GC', lonlat=True, coord='G',fontsize=15)


    for i in range(0,360,30):

        lmflat = np.arange(-90,90)
        lmflon = np.zeros(180)+i
        marker = '--'
        #f i == 210: marker = 'b-'
        hp.projplot(lmflon, lmflat, marker,lonlat=True,coord='G',rot=(glimf.lon/degree,glimf.lat/degree) ,color='gray')

    #hp.projplot(equateur_lon, equateur_lat, 'r-', lonlat=True, coord='C')
    lmfeqlon = np.arange(0,360)
    lmfeqlat = np.zeros(360)
    hp.projplot(lmfeqlon, lmfeqlat, 'b-',lonlat=True,coord='G',rot=(glimf.lon/degree,glimf.lat/degree) )
    
    
    #ibex ribbon
    if show_ibex:
        hp.projplot(lmfeqlon, lmfeqlat+90-73.26, '-',lonlat=True,coord='G', color='grey', rot=(gibex.lon/degree,gibex.lat/degree) )
        hp.projplot(lmfeqlon, lmfeqlat+90-79.18, '-',lonlat=True,coord='G', color='grey', rot=(gibex.lon/degree,gibex.lat/degree) )
        hp.projscatter([gibex.lon/degree], [gibex.lat/degree], lonlat=True,color='black', coord='G')
        hp.projtext(gibex.lon/degree,gibex.lat/degree, labels['ibex_center'], lonlat=True, coord='G')
        
    #2D fit
    if fitvalue>0:

        
        fitlat = np.arange(-90,90)
        fitlon = np.zeros(180)
        
        hp.projplot(fitlon+fitvalue, fitlat, '-+',lonlat=True, color='red' )
        hp.projplot(fitlon+180+fitvalue, fitlat, '-+',lonlat=True, color='blue' )
        if CG_corr is not None:
            hp.projplot(fitlon+fitvalue+CG_corr, fitlat, '-+',lonlat=True, color='grey' )
            hp.projplot(fitlon+180+CG_corr+fitvalue, fitlat, '-+',lonlat=True, color='grey' )
    

    #37.8421671285


    gplon = np.arange(0,360)
    gplat = np.zeros(360)
    if gc:
        hp.projplot(gplon, gplat, 'r-',lonlat=True, coord='G')

    #B-V plane
    lmflat = np.arange(-90,90)
    lmflon = np.zeros(180)-31.25
    marker = 'b-'
    hp.projplot(lmflon, lmflat, marker,lonlat=True,coord='G',rot=(glimf.lon/degree,glimf.lat/degree) )
    lmflon = (lmflon+180)%360
    hp.projplot(lmflon, lmflat, marker,lonlat=True,coord='G',rot=(glimf.lon/degree,glimf.lat/degree) )
    
    x = [227.28]
    y = [34.62]
    x =[glimf.lon/degree]
    y=[glimf.lat/degree]

    print glimf.lat, glimf.lon
    #p.projscatter(x, y, lonlat=True,color='black', coord='G',rot=(glimf.lon/degree,glimf.lat/degree))
    hp.projscatter(x, y, lonlat=True,color='black', coord='G')
    hp.projtext(glimf.lon/degree,glimf.lat/degree, labels['limf'], lonlat=True, coord='G',fontsize=18)
    
 
    xv=[5.25]
    yv=[12]
    solapex = [[47.9],[23.8]]
   
    hp.projscatter(solapex[0],solapex[1],lonlat=True,color='black', coord='G')
    hp.projtext(solapex[0],solapex[1], labels['v_lsr'], lonlat=True, coord='G',fontsize=18)
    
    vism = [[5.25],[12]]
    hp.projscatter(vism[0],vism[1],lonlat=True,color='black', coord='G')
    hp.projtext(vism[0],vism[1], labels['v_ism'], lonlat=True, coord='G',fontsize=18)
    
    hp.projtext([np.array(lmflon).mean()],[180+45], '$B-V$', lonlat=True, coord='G', 
                rot=(glimf.lon/degree,glimf.lat/degree),fontsize=15)

    r = hp.Rotator(coord=['G'],rot=(glimf.lon/degree,glimf.lat/degree)) 
    ri = r.get_inverse()
    vcr_lon, vcr_lat = r(5.25*degree, 12*degree)

    hp.projtext(95*degree+(rotation[0]-180)*degree, 280*degree-rotation[1]*degree,
                    watermark,
                    #coord=coords,
                    color="grey",
                    alpha=0.5,
                    rotation=0,
                    fontdict={"family":"sans-serif", "weight":"bold", "size":42})
    print ri(vcr_lon,0)
    lmflonrc = np.arange(-90,90)
    lmflatrc = np.zeros(180)
    lxrc = [];lyrc=[]
    if filename:
        fig.savefig(filename, dpi=100)
    plt.show()

