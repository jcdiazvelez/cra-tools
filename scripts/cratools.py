import setupMaps
import healpy
import numpy.random
import healpy as hp
import math
import numpy
import numpy as np
import os
import os.path
import ephem
import pylab

import copy
from scipy.stats import ks_2samp
import astropy.units as u
from astropy.coordinates import Galactic
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from numpy import log, sqrt, cos, pi


degree = math.pi / 180
default_binsize = 20

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
    for ipix in range(npix):
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


def get_1d_bin(data, bg, minDec=-60 * degree, maxDec=40 * degree, binsize=10):

    npix = data.size
    nside = healpy.npix2nside(npix)

    minPix = healpy.ang2pix(nside, 90 * degree - maxDec, 0.)
    maxPix = healpy.ang2pix(nside, 90 * degree - minDec, 0.)

    y0 = np.zeros(360 / binsize)
    y1 = np.zeros(360 / binsize)
    xaxis = np.array(range(0, 360 / binsize))
    err_sqr = np.zeros(360 / binsize)
    sig = (data - bg) / (bg + 1e-16)

    for ipix in range(minPix, maxPix):
        theta, phi = healpy.pix2ang(nside, ipix)
        # print theta/degree,phi/degree, (data[ipix]-bg[ipix])/(bg[ipix]+1e-16)
        y0[int(np.floor(phi / degree)) / binsize] += data[ipix]
        y1[int(np.floor(phi / degree)) / binsize] += bg[ipix]
    for ibin in range(0, 360 / binsize):
        err_sqr[ibin] += 1.0 / y1[ibin]

    return (xaxis, (y0 - y1) / y1, np.sqrt(err_sqr))


"""
Bin Healpix relative intensity pixels in bins along a declination band
"""
def get_1d_bin_relint(relintdata, bg, minDec=-30 * degree, 
		maxDec=-20 * degree, xbinsize=10, ybinsize=1, name=None):


    npix = relintdata.size
    nside = healpy.npix2nside(npix)
    bnpix = bg.size
    bnside = healpy.npix2nside(bnpix)

    y0 = np.zeros(360 / xbinsize)
    y1 = np.zeros(360 / xbinsize)
    err_sqr = np.zeros(360 / xbinsize)
    avg = np.zeros(360 / xbinsize)
    inbinpix = np.zeros(360 / xbinsize)
    xaxis = np.array(range(0, 360 / xbinsize))*xbinsize + .5*xbinsize

    dec_strip = healpy.query_strip(nside, 90*degree - maxDec, 90*degree - minDec, inclusive=True)
    for ipix in dec_strip:
            theta, phi = healpy.pix2ang(nside, ipix)
            idx = int(phi / degree) / xbinsize
            y0[idx] += relintdata[ipix]
            inbinpix[idx] += 1.0
    dec_strip = healpy.query_strip(bnside, 90*degree - maxDec, 90*degree - minDec, inclusive=True)
    for ipix in dec_strip:
            theta, phi = healpy.pix2ang(bnside, ipix)
            idx = int(phi / degree) / xbinsize
            y1[idx] += bg[ipix]

    avg = y0/inbinpix
    err_sqr = 1/y1
    return (xaxis, avg, np.sqrt(err_sqr))


def get_1d_bin_relint_new(
		relintdata, bg, minDec=-30 * degree, 
		maxDec=-20 * degree, xbinsize=10, ybinsize=1, name=None):

    npix = relintdata.size
    nside = healpy.npix2nside(npix)

    minPix = healpy.ang2pix(nside, 90 * degree - maxDec, 0.)
    maxPix = healpy.ang2pix(nside, 90 * degree - minDec, 0.)

    minTheta = 90 - maxDec/degree
    maxTheta = 90 - minDec/degree
    ybins = int((maxTheta-minTheta)/ybinsize)

    y0 = np.zeros(360 / xbinsize)
    y1 = np.zeros(360 / xbinsize)
    nbinpix = np.zeros(360 / xbinsize)
    y2 = [[]]*(360 / xbinsize)
    xaxis = np.array(range(0, 360 / xbinsize))*xbinsize
    err_sqr = np.zeros(360 / xbinsize)

    #print "ybins",ybins
    #print "minT",minTheta, minDec
    #print "maxT",maxTheta, maxDec

    for iy in range(ybins):
        theta1 = (minTheta+iy)*degree
        theta2 = (minTheta+iy+ybinsize)*degree
        phi = 0
        dec_strip = healpy.query_strip(nside, theta1, theta2, inclusive=True)
        iy0 = np.zeros(360 / xbinsize)
        iy1 = np.zeros(360 / xbinsize)
        inbinpix = np.zeros(360 / xbinsize)
        for ipix in dec_strip:
            theta, phi = healpy.pix2ang(nside, ipix)
            #print theta/degree, phi/degree
            idx = int(np.floor(phi / degree)) / xbinsize
            iy0[idx] += relintdata[ipix]
            iy1[idx] += bg[ipix]
            inbinpix[idx] += 1.0
        avg = np.mean(inbinpix)
        for ibin in range(0, 360 / xbinsize):
            y0[ibin] +=  iy0[ibin]/inbinpix[idx]
            if iy1[ibin] > 0:
               err_sqr[ibin] += 1.0/iy1[ibin]/(inbinpix[idx]*inbinpix[idx])
            nbinpix[ibin] += 1.0

    for ibin in range(0, 360 / xbinsize):
        y0[ibin] /=  nbinpix[ibin]
    bsum = 0
    dec_strip = healpy.query_strip(nside, minTheta, maxTheta, inclusive=True)
    for ipix in dec_strip:
        bsum += bg[ipix]

    return (xaxis, y0, np.sqrt(err_sqr)/ybins)





def plot_projection(signals,
                    backgrounds,
                    legends,
                    title='1D Projection Sidereal-time',
                    binsize=10,
                    fmts='-o',
                    minmaxes = [(-40*degree,60*degree)],
                    split=True,
                    ymin=-2,
                    ymax=2,
                    log_scale=False):


    colors = ['red','blue','green','magenta']
    if split:
        # Two subplots, the axes array is 1-d
        f, axarr = plt.subplots(len(signals), 1, sharex=True, figsize=(10,7))
        ax = axarr[0]

    else:
        plt.figure(figsize=(20, 10))
        ax = plt.subplot(211)

    if type(fmts) != list:
        fmts = [fmts] * len(signals)

    iplot = 0
    ax.set_title(title, fontsize=20)
    ax.invert_xaxis()
    if log_scale:
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e6)
    for signal, background, legend, fmt, minmax in zip(signals, backgrounds, legends, fmts, minmaxes):
        x, y, yerr = get_1d_bin(signal, background,minDec=minmax[0],maxDec=minmax[1],binsize=binsize)
        if split:
           ax = axarr[iplot]
        x = np.array(x)*binsize
        p = ax.errorbar(x , y*1e3, yerr=np.array(yerr)*1e3, fmt=fmt, label=legend,color=colors[iplot])
        ax.legend(loc='upper left', numpoints=1, fontsize=20)
        iplot += 1

        for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)

        ax.set_ylim([ymin, ymax])
        ax.set_ylabel('Relint [$10^{-3}$]', fontsize=20)

    #ax.xlim(ax.xlim()[::-1])
    ax.set_xlabel('RA ($^\circ$)', fontsize=20)
    plt.show()


from numpy import log, sqrt, cos, pi


def LMSig(Non, Noff, alpha):
    Noff /= alpha
    if (Non < alpha * Noff):
        sign = -1.
    else:
        sign = 1.
    # if Non & Noff are close FPE may occur to make sqrt(neg num)
    if((Non * log(((1. + alpha) * Non) / (alpha * (Non + Noff))) + Noff * log(((1. + alpha) * Noff) / (Non + Noff))) < 0.):
        return 0.

    sigma = sign * sqrt(2. * (Non * log(((1. + alpha) * Non) / (alpha *
                                                                (Non + Noff))) + Noff * log(((1. + alpha) * Noff) / (Non + Noff))))
    return sigma


def alpha(r, dec, int_hrs=24.):
    return pi * r / (2. * int_hrs * 15 * degree * cos(dec))


def CalcLMSignificance(md, mb, boundary=-20 * degree, smoothing=5. * degree):
    npix = md.size
    nside = healpy.npix2nside(npix)
    pixrad = healpy.max_pixrad(nside)
    SignificanceSkymap = numpy.zeros(healpy.nside2npix(nside))
    for i in range(npix):
        Nd = float(md[i])
        Nb = float(mb[i])
        theta, phi = healpy.pix2ang(nside, i)
        if Nd > 0 and Nb > 0:
            if theta > boundary:
                SignificanceSkymap[i] = LMSig(
                    Nd, Nb, alpha=pi * smoothing / (2. * 24. * 15. * degree * cos(90. * degree - theta)))
            else:
                SignificanceSkymap[i] = LMSig(Nd, Nb, alpha=1. / 20.)
        else:
            SignificanceSkymap[i] = 0.
    return SignificanceSkymap


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


def plot_proj_dec(
        signals,
        backgrounds,
        legends,
        lang='eng',
        colors = ['red','blue','magenta','green'],
        title='1D Projection $\delta$',
        binsize=10,
        mindec=-90 * degree,
        maxdec=90 * degree,
        fmts='-o',
        log_scale=False,
        xmin=None,
        xmax=None,
        ymin=-1,
        ymax=-1,
        xspacing=20,
        output=None,
        figsize=(15,8)):

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)

    if type(fmts) != list:
        fmts = [fmts] * len(signals)

    icolor=0
    patterns = ['-', '+', 'x', '\\', '*', 'o', 'O', '.']
    for signal, background, legend, fmt in zip(signals, backgrounds, legends, fmts):
        x, y, yerr = get_1d_bin_dec(
            signal, background, binsize=binsize, minRA=mindec, maxRA=maxdec)
        yerr = np.array(yerr)
        p = ax.errorbar(x, y, yerr=yerr, fmt=fmt, label=legend,color=colors[icolor])
        y0 = np.zeros(y.size)
        ax.fill_between(x, y0, y, where=y >= y0, facecolor=colors[icolor], interpolate=True, alpha=0.3)
        #p = plt.plot(x*binsize,y,'o',label=legend)
        icolor = (icolor+1)%len(colors)

    plt.legend(loc='upper right', numpoints=1, fontsize=20)
    # plt.xlim(plt.xlim()[::-1])
    if log_scale:
        ax.set_yscale('log')
    if ymin >= 0:
        ax.set_ylim(bottom=ymin)
    if ymax >= 0:
        ax.set_ylim(top=ymax)

    if xmin is None : xmin = x[0]
    if xmax is None : xmax = x[-1]
    ax.set_xlim((xmin,xmax))
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    ax.set_title(title, fontsize=20)
    plt.xticks(np.arange(xmin, xmax+1, xspacing))
    ax.set_xlabel('$\delta$ [$^\circ$]', fontsize=25)
    if lang == 'esp':
       ax.set_ylabel(r'$\rm{N\acute{u}mero}$ de sucesos', fontsize=25)
    else:
       ax.set_ylabel('Number of events', fontsize=25)
    if output != None:
        if output.find("jpg") > 0:
            fig.savefig(output, dpi=100)
        else:
            #savefig(output)
            jpgOutput = output.replace("png", "jpg")
            fig.savefig(jpgOutput, dpi=100)

    plt.show()




def get_1d_bin_dec(data, bg, minRA=-80 * degree, maxRA=80 * degree, binsize=1):
    from math import sin, asin
    import random
    npix = data.size
    nside = healpy.npix2nside(npix)

    minbin = int(minRA / degree)
    maxbin = int(maxRA / degree)
    array_size = (maxbin - minbin) / binsize
    y0 = np.zeros(array_size)
    y1 = np.zeros(array_size)
    err_sqr = np.zeros(array_size)

    print (minbin, maxbin, binsize)
    xaxis = np.array(range(minbin, maxbin, binsize)) + binsize*.5

    print(minbin + 90, maxbin + 90, binsize)
    for ipix in range(data.size):
        theta, phi = healpy.pix2ang(nside, ipix)
        #theta += (0.45 - random.random()) * sqrt(healpy.nside2pixarea(nside))
        thetaprime = 90 - (theta) / degree
        ibin = int(thetaprime - minbin) / binsize
        if ibin >= 0 and ibin < y0.size:
            y0[ibin] += data[ipix]
            y1[ibin] += bg[ipix]

    for ibin in range(y1.size):
        if y1[ibin] > 0.0:
            err_sqr[ibin] += sqrt(1.0 / y1[ibin])

    norm = np.sum(y0)
    return (xaxis, y0, err_sqr/y1.size)


def significance(OrigData, OrigBG, smooth_radius=5, boundary=-20 * degree):
    npix = OrigData.size
    nside = hp.npix2nside(npix)
    SmoothData = numpy.zeros(hp.nside2npix(nside))
    SmoothBG = numpy.zeros(hp.nside2npix(nside))
    SmoothSignal = numpy.zeros(hp.nside2npix(nside))
    hp_boundary = 90 * degree - boundary

    for ipix in range(npix):
        ListOfPixels = []
        theta, phi = hp.pix2ang(nside, ipix)
        vec = hp.ang2vec(theta, phi)
        ListOfPixels = hp.query_disc(nside, vec, smooth_radius * degree)

        # Sum up pixels in the angular cap; set value of ith pixel to this sum
        sumData = 0.
        sumBG = 0.

        for jpix in range(ListOfPixels.size):
            sumData += OrigData[ListOfPixels[jpix]]
            sumBG += OrigBG[ListOfPixels[jpix]]

        SmoothData[ipix] = sumData
        SmoothBG[ipix] = sumBG

    print (" Computing  significance... ")
    for ipix in range(SmoothData.size):

        Nd = SmoothData[ipix]
        Nb = SmoothBG[ipix]

        if Nd > 0 and Nb > 0:
            theta, phi = healpy.pix2ang(nside, ipix)
            if Nd > 0 and Nb > 0:
                if theta > hp_boundary:
                    alpha = pi * smooth_radius * degree / \
                        (2. * 24. * 15. * degree * cos(90. * degree - theta))
                else:
                    alpha = 1. / 20.
                SmoothSignal[ipix] = LMSignificance(Nd, Nb, alpha=alpha)
        else:
            SmoothSignal[ipix] = 0
        if SmoothSignal[ipix] != SmoothSignal[ipix]:
            SmoothSignal[ipix] = 0

    return SmoothSignal


def LMSignificance(Non, Nback, alpha=1. / 20.):

    Noff = Nback / alpha
    sign = 1.
    if (Non < alpha * Noff):
        sign = -1.

    # Li&Ma (eq 5)
    #  return (Non - alpha*Noff)/sqrt(Non +alpha*alpha*Noff);

    # Li&Ma (eq 9)
    # return (Non - alpha*Noff)/sqrt(alpha*(Non +Noff));

    # Li&Ma (eq 17)
    return sign * numpy.sqrt(2 * (Non * numpy.log(((1 + alpha) * Non) / (alpha * (Non + Noff))) +
                                  Noff * numpy.log(((1 + alpha) * Noff) / (Non + Noff))))



def generate_map(bgMap, dataMap):
  alpha = 1.0;

  npix = dataMap.size
  signalMap = np.zeros(npix) + hp.UNSEEN
  for i in range(npix):
    # Throw a Poisson random number using the data as a mean
    Nd = numpy.random.poisson(dataMap[i])
    Nb = bgMap[i]
    if Nd > 0 and Nb > 0:
      #signalMap[i] = (Nd - alpha*Nb) / numpy.sqrt(Nd + alpha*Nb)
      #signalMap[i] = (Nd - alpha*Nb) / (alpha*Nb+1e-16)
      #signalMap[i] = Nd 
      signalMap[i] = (Nd - alpha*Nb) / (alpha*Nb)
  return signalMap


def galactic(equaMap):
     from astropy.coordinates import SkyCoord
     from astropy import units as u

     npix = equaMap.size
     nside = healpy.npix2nside(npix)
     galacticMap = np.zeros(npix)
     for i in range(npix):
         theta, phi = healpy.pix2ang(nside, i)
         c_icrs = SkyCoord(ra=phi*u.radian, dec=theta*u.radian, frame='icrs')
         g_theta = c_icrs.galactic.b.rad
         g_phi = c_icrs.galactic.l.rad
         try:
            g_i = healpy.ang2pix(nside, g_theta+90*degree, g_phi)
         except Exception:
            print (g_phi/degree,g_theta/degree )
            raise Exception
         galacticMap[g_i] = equaMap[i]
     return galacticlMap

def isotropic_band(dataMap, bgMap=None, n=100, lmax=40, percentile=90, 
                       isomap_dir=None,lhiter=10, 
                       mask_top=70 * degree, mask_bottom=-80 * degree):
    from tqdm import tqdm
    import contextlib
    @contextlib.contextmanager
    def capture():
        import sys
        from io import StringIO
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
       cls = []
       p0 = (100 - percentile)/2
       for c in range(lmax+1):
          cls.append([])
       clmins = numpy.ones(lmax+1)*1e6
       clmaxs = numpy.zeros(lmax+1)
       for it in tqdm(range(n)):
           if isomap_dir:
              iso_relint = healpy.read_map(os.path.join(isomap_dir,str(it),"CR_32_360_iteration%u.fits" % lhiter),0)
           elif bgMap is not None:
              iso_relint = generate_map(bgMap, bgMap)
           else:
              raise BaseException("No background map provided!")

           #iso_ccl = hp.anafast(iso_relint, lmax=lmax)
           iso_ccl = hp.anafast(mask(iso_relint, top=mask_top, bottom=mask_bottom))
           for c in range(lmax+1):
              cls[c].append(iso_ccl[c])

    for c in range(lmax+1):
      pmin, pmax = numpy.percentile(cls[c], [p0,100-p0])
      clmins[c]=pmin
      clmaxs[c]=pmax

    return clmins,clmaxs



def error_bars(dataMap, bgMap=None, n=100, lmin=0, lmax=40, percentile=90, 
              isomap_dir=None,lhiter=10, do_chisq=False, thetamin=-80*degree,thetamax=70*degree):
    import scipy.stats
    from tqdm import tqdm
    import contextlib
    @contextlib.contextmanager
    def capture():
        import sys
        from io import StringIO
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
       cls = []
       p0 = (100 - percentile)/2
       for c in range(lmax+1):
          cls.append([])
       yerr = numpy.zeros(lmax+1)
       yerr1 = numpy.zeros(lmax+1)
       yerr2 = numpy.zeros(lmax+1)
       for it in tqdm(range(n)):
           if isomap_dir:
              iso_relint = healpy.read_map(os.path.join(isomap_dir,str(it),"CR_32_360_iteration%u.fits" % lhiter),0)
           elif bgMap is not None:
              iso_relint = generate_map(bgMap, dataMap)
           else:
              raise BaseException("No background map provided!")
           if lmin > 1:
              alms = healpy.map2alm(iso_relint, lmax=lmin-1 )
              npix = len(iso_relint)
              nside = healpy.npix2nside(npix)
              iso_relint -= healpy.alm2map(alms, nside)
              

           #iso_ccl = hp.anafast(iso_relint, lmax=lmax)
           iso_ccl = hp.anafast(mask(iso_relint, top=thetamax, bottom=thetamin))
           for l in range(lmax+1):
              cls[l].append(iso_ccl[l])

    chisq = []
    means = np.zeros(lmax+1)
    for l in range(lmax+1):
      sigma = numpy.std(cls[l])
      mu = numpy.mean(cls[l])
      means[l] = mu
      yerr[l] = sigma
      #q = [.16, .5, .84]
      q = [.1, .5, .9]
      #quatiles = np.quantile(cls[c], q)
      quantiles = np.percentile(cls[l], q)
      yerr1[l] = quantiles[0]
      yerr2[l] = quantiles[2]
      #print "quantiles", quantiles

      if do_chisq:
         hist, bin_edges = np.histogram(np.array(cls))
         pdf = 1./(sigma * np.sqrt(2 * np.pi)) * np.exp( - (bin_edges[1:] - mu)**2 / (2 * sigma**2) )
         chisq.append(scipy.stats.chisquare(hist, f_exp=pdf, ddof=0, axis=0))
      

    if lmin > 1:
       print ("Using subtracted multipole l=",lmin-1)

    return means, yerr, chisq, yerr1, yerr2, np.array(cls)




####
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


def plot_projection_relint(relintdata,
                    backgrounds,
                    legends,
                    unsmoothed_relintdata=[None,None,None],
                    title='1D Projection Sidereal-time',
                    xbinsize=10,
                    ybinsize=1,
                    ybands= [(-20*degree,-17*degree)],
                    fmts='-o',
                    minmaxes = [(-19*degree,-18*degree)],
                    colors = ['blue','red','green','magenta'],
                    alpha = [1.0],
                    watermark= "",
                    log_scale=False, step=False,
                    scale=1.0,
                    errors=True,
                    ylim=None,
                    ylabel=None,
                    output=None):

    fig = plt.figure(figsize=(16, 8))
    ax = plt.subplot(111)
    rets =[]
    from scipy.stats import ks_2samp


    if type(fmts) != list:
        fmts = [fmts] * len(relintdata)

    c = 0
    ks_sample1 = None
    ks_sample2 = None

    for signal, unsmoothed_signal, background, legend, fmt, minmax, yband, in zip(relintdata,unsmoothed_relintdata, backgrounds, legends, fmts, minmaxes, ybands):
        
        nside = healpy.npix2nside(signal.size)
        #signal = TopHatSmooth(nside,unsmoothed_signal,radius=5.0*degree,average=True)
        if unsmoothed_signal is None:
            unsmoothed_signal = signal
        #x, y, yerr = get_1d_bin_relint_lima(
        #x, y, yerr = get_1d_bin_relint_new(
        x, y, yerr = get_1d_bin_relint(
			signal, background, minDec=minmax[0], 
			maxDec=minmax[1], xbinsize=xbinsize, 
			ybinsize=ybinsize,name=legend)
        yerr = np.array(yerr)

        ks_sample2 = ks_sample1
        ks_sample1 = cdf(x,y)( numpy.random.uniform(size=1000) )
      
        if ks_sample2 is not None:
           print ('ks-test',ks_2samp(ks_sample2,ks_sample1))

        syserr = np.zeros(len(y))
        y = np.array(y)
 
        myplot = {'x':x,'y':y }
        if legend:
           if errors:
              myplot = { 'x':x,'y':y,'yerr':yerr, 'fmt':fmt,'legend':legend, 'color':colors[c] }
              p = plt.errorbar(x, y*scale, yerr=yerr*scale, fmt=fmt, label=legend, color=colors[c], alpha=alpha[c%len(alpha)])
           else:
              myplot = { 'x':x,'y':y,'yerr':None, 'fmt':fmt,'legend':legend, 'color':colors[c] }
              p = plt.plot(x, y*scale, fmt, label=legend, color=colors[c], alpha=alpha[c%len(alpha)])
        else:
           p = plt.errorbar(x, y*scale, yerr=yerr*scale, fmt=fmt, color=colors[c], alpha=alpha[c%len(alpha)])
           myplot = { 'x':x,'y':y,'yerr':yerr, 'fmt':fmt,'legend':None, 'color':colors[c] }
        if step:
           plt.step(x, y, where='mid', linewidth=2, color=colors[c])
        rets.append(myplot)
        ymin = 1.0e12*np.ones(y.size)
        ymax = -1.0e12*np.ones(y.size)
        dy = minmax[1]-minmax[0]
        if yband > 0:
           for yib in np.arange(yband[0], yband[1],1*degree):
               #xb, yb, yerrb = get_1d_bin_relint_new(
               xb, yb, yerrb = get_1d_bin_relint(
			unsmoothed_signal, background, minDec=yib-0.5*dy, 
			maxDec=yib+0.5*dy, xbinsize=xbinsize, 
			ybinsize=ybinsize)
               for i in range(y.size):
                   ymin[i] = min(yb[i]-yerrb[i],ymin[i])
                   ymax[i] = max(yb[i]+yerrb[i],ymax[i])
           for i in range(y.size):
               syserr[i] = max(abs(y[i]-ymin[i]), abs(y[i]-ymax[i]))
           #ax.fill_between(x, (y+syserr)*scale, (y+syserr)*scale, where=ymax >= ymin, facecolor=colors[c], interpolate=True, alpha=0.05)
           ax.fill_between(x, (y-syserr)*scale, (y+syserr)*scale,  facecolor=colors[c], interpolate=True, alpha=0.1)

        c += 1

    ax.text(180, 1.5, watermark, alpha=0.3,ha= 'left', fontsize=30)

    plt.legend(loc='lower left', numpoints=1, fontsize=15)
    plt.xlim([360,0])
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    ax.set_title(title, fontsize=20)
    if ylabel: 
       ax.set_ylabel(ylabel, fontsize=20)
    ax.set_xlabel('RA ($^\circ$)', fontsize=20)
    if ylim: 
        ax.set_ylim(ylim)
    elif log_scale:
        ax.set_yscale('log')
        ax.set_ylim(bottom=1e6)

    if output != None:
        if output.find("jpg") > 0:
            fig.savefig(output, dpi=100)
        else:
            #savefig(output)
            jpgOutput = output.replace("png", "jpg")
            fig.savefig(jpgOutput, dpi=100)

    plt.show()
    return rets


def hist_proj_dec(
        signals,
        backgrounds,
        legends,
        lang='eng',
        sine=False,
        colors = ['red','blue','magenta','green'],
        title='1D Projection $\delta$',
        binsize=10,
        mindec=-90 * degree,
        maxdec=90 * degree,
        fmts='-o',
        log_scale=False,
        xmin=None,
        xmax=None,
        ymin=-1,
        ymax=-1,
        alpha=0.3,
        xspacing=20,
        output=None,
        patterns = ['', '', '/', '\\', '*', 'o', 'O', '.'],
        bars= [False, False, True, True, True, True, True, True],
        watermark="",
        figsize=(15,8)):

    fig = plt.figure(figsize=figsize)
    ax = plt.subplot(111)

    if type(fmts) != list:
        fmts = [fmts] * len(signals)

    icolor=0
    nbins = (maxdec-mindec)/binsize/degree
    if sine:
        nbins = (np.sin(maxdec)-np.sin(mindec))/binsize
        print ("nbins",maxdec,mindec,binsize,nbins )
        nbins = int(nbins)
    print (maxdec/degree, mindec/degree, nbins, binsize/degree)
    for signal, background, legend, fmt, pattern, bar in zip(signals, backgrounds, legends, fmts,patterns, bars):
        y1, weights = get_hist_dec(
            signal, background, binsize=binsize, minRA=mindec, maxRA=maxdec,sine=sine)
        #n, bins, patches = plt.hist( 
        #       np.array(y1), nbins, weights=np.array(weights), 
        #       facecolor=colors[icolor],alpha=0.3, histtype='stepfilled')
        n, bins = np.histogram(
               np.array(y1), nbins, weights=np.array(weights)) 
        yerr, errbins = np.histogram(
               np.array(y1), nbins, weights=np.square(np.array(weights)))
        #print len(n), len(bins)
        xaxis = np.array(bins[0:-1]) #+ binsize*.5
        #print legend,"min",min(weights)
        #print legend,"max",max(weights)

        #p = ax.plot(xaxis, n, fmt, label=legend,color=colors[icolor])
        #p = ax.errorbar(xaxis , n, yerr=np.sqrt(yerr), fmt=fmt, label=legend,color=colors[icolor])
        #y0 = np.zeros(xaxis.size)
        #ax.fill_between(xaxis, y0, n, where=n >= y0, facecolor=colors[icolor], interpolate=True, alpha=0.3)
        if bar:
        	ax.bar(xaxis, n, width=binsize, color=colors[icolor], yerr=1.0/np.sqrt(yerr), 
			alpha=alpha,label=legend,hatch=pattern) 
        else:
        	ax.step(xaxis,n, where='post', linewidth=2, color=colors[icolor],label=legend)
        icolor = (icolor+1)%len(colors)

    plt.legend(loc='upper right', numpoints=1, fontsize=20)
    # plt.xlim(plt.xlim()[::-1])
    if log_scale:
        ax.set_yscale('log')
    if ymin >= 0:
        ax.set_ylim(bottom=ymin)
    if ymax >= 0:
        ax.set_ylim(top=ymax)
    if xmin is None : xmin = bins[0]
    if xmax is None : xmax = bins[-1]

    x0 = [-16,-16]
    y0 = [1e6,1e12]
    if ymax > 0:
       y0 = [ymin,ymax]
    if sine:
       x0 = [np.sin(-16*degree),np.sin(-16*degree)]
    #ax.plot(x0, y0, '-', color='red')



    ax.set_xlim((xmin,xmax))
    for item in (ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(20)
    ax.set_title(title, fontsize=20)
    plt.xticks(np.arange(xmin, xmax+binsize, xspacing))
    if sine:
       ax.set_xlabel('$\sin(\delta)$', fontsize=25)
    else:
       ax.set_xlabel('$\delta$ [$^\circ$]', fontsize=25)
    if lang == 'esp':
       ax.set_ylabel(r'$\rm{N\acute{u}mero}$ de sucesos', fontsize=25)
    else:
       ax.set_ylabel('Number of events', fontsize=25)

    #print "loc", 0.25*(xmin+xmax), 0.3*(ymin+ymax), watermark
    ax.text(xmin+10, 0.3*(ymin+ymax), watermark, alpha=0.3,ha= 'left', fontsize=26)
    if output != None:
        if output.find("jpg") > 0:
            fig.savefig(output, dpi=100)
        else:
            #savefig(output)
            jpgOutput = output.replace("png", "jpg")
            fig.savefig(jpgOutput, dpi=100)

    plt.show()




def get_hist_dec(data, bg, minRA=-80 * degree, maxRA=80 * degree, binsize=1,sine=False):
    y0 = []
    weights = []
    errors = []

    minbin = int(minRA / degree)
    maxbin = int(maxRA / degree)
    if sine:
       maxbin = np.sin(maxRA)
       minbin = np.sin(minRA)
    array_size = int((maxbin - minbin) / binsize)
    npix = data.size
    nside = healpy.npix2nside(npix)

    print (minbin, maxbin, binsize)
    for ipix in range(data.size):
        theta, phi = healpy.pix2ang(nside, ipix)
        thetaprime = 90 - (theta) / degree
        if sine:
           thetaprime = np.sin(0.5*np.pi - theta)
        y0.append(thetaprime)
        weights.append(data[ipix])

    return (y0, weights)

def get_1d_bin_relint_lima(
		relintdata, bg, minDec=-30 * degree, 
		maxDec=-20 * degree, xbinsize=10, ybinsize=1, name=None):

    npix = relintdata.size
    nside = healpy.npix2nside(npix)

    minTheta = 90 - maxDec/degree
    maxTheta = 90 - minDec/degree

    y0 = np.zeros(360 / xbinsize)
    y1 = np.zeros(360 / xbinsize)
    counts = np.zeros(360 / xbinsize)
    nbinpix = np.zeros(360 / xbinsize)
    xaxis = np.array(range(0, 360 / xbinsize))*xbinsize
    err_sqr = np.zeros(360 / xbinsize)

    #print "minT",minTheta, minDec
    #print "maxT",maxTheta, maxDec

    dec_strip = healpy.query_strip(nside, minTheta*degree, maxTheta*degree, inclusive=True)
    bsum = 0
    for ipix in dec_strip:
        bsum += bg[ipix]
        theta, phi = healpy.pix2ang(nside, ipix)
        idx = int(np.floor(phi / degree)) / xbinsize
        y0[idx] += relintdata[ipix]
        y1[idx] += bg[ipix]
        counts[idx] += 1

    #err_sqr = y0*np.sqrt(y1)/counts
    err_sqr = 1/np.sqrt(y1)
    if name:
       print (name, "%0.2g" % bsum, "events in between ", minDec/degree, " and ", maxDec/degree )

    return (xaxis, y0/counts, err_sqr)


#how galactic plane
import astropy.units as u
from astropy.coordinates import Galactic
import ephem
import pylab
from matplotlib.colors import LinearSegmentedColormap

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

#            ax.annotate(0, xy=(-1.19,0.03), color="white", weight="bold")
#            ax.annotate(90, xy=(1.77,0.03), color="white", weight="bold")
#            ax.annotate(180, xy=(0.75,0.03), color="white", weight="bold")
#            ax.annotate(270, xy=(-0.24,0.03), color="white", weight="bold")
    


def MaskDisk(nside,origMap,radius=0,dec=0,ra=0):
  pixelList = []
  dsum = 0.0
  newmap = copy.deepcopy(origMap)

  # Loop over all pixels
  theta = 90*degree - dec
  phi = ra
  center = healpy.dir2vec(theta,phi)

  # Grab pixels with given angular radius of ith pixel
  pixelList = healpy.query_disc(nside, center, radius)

  #print len(pixelList)
  for jpix in pixelList:
      newmap[jpix] = healpy.UNSEEN
  return newmap


def plot_map_with_features(skymap,title,proj='C',gc=True,label='', show_ibex=False,fit=None,
                           fitvalue=37.8421671285,CG_corr=None,lang='eng',filename=None,watermark="", 
                           features=True,moll=True, graticule=True,
                           center_on_fit=False,thresh=None,dMin=None,dMax=None,fieldlines=True,
                           limf_l=227.28*degree, limf_b=34.62*degree,fitcolor="darkgreen",
                           whiteb=False,
                           notext=False,
                           exposure_hawc=[],exposure_icecube=[]):
    
    
    limf = ephem.Ecliptic(limf_l, limf_b)
    bbox=(dict(facecolor='white', alpha=0.5, edgecolor='black', boxstyle="round,pad=0.1"))
        
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


    
    gibex = ephem.Galactic(ibex)
    
    glimf = ephem.Galactic(limf)
    climf = ephem.Equatorial(limf)
    
    cibex = ephem.Equatorial(ibex)
    #print 'LIMF',climf.ra/degree, climf.dec/degree
    #print 'ANTILIMF',climf.ra/degree-180, -climf.dec/degree
    #print 'LIMF',climf.ra,climf.dec
    
    tail = ephem.Ecliptic(21*degree, -6.*degree)
    ctail = ephem.Equatorial(tail)
    gtail = ephem.Galactic(tail)
    #vism = [[5.25],[12]]
    evism = ephem.Ecliptic(79.00*degree, -4.98*degree)
    gvism = ephem.Galactic(evism)
    cvism = ephem.Equatorial(evism)
    vism = [[cvism.ra/degree],[cvism.dec/degree]]

    

    if proj=='C0':
        rotation = (0,0,0)
    elif proj=='C':
        rotation = (-180,0,0)
    elif proj=='VISM':
        rotation = (180+cvism.ra/degree,-cvism.dec/degree)
    elif proj=='LIMF':
        rotation = (climf.ra/degree,climf.dec/degree)
        notext=True
    elif proj=='ANTILIMF':
        rotation = (climf.ra/degree-180,climf.dec/degree-34.62)
        rotation = (climf.ra/degree-180,-climf.dec/degree) 
        notext=True
    else:
        rotation = (0,0,0)
    if moll:
        hp.mollview(mask(skymap,top=67*degree,bottom=-90*degree), 
            fig=1,
            title=title, 
            rot=rotation,coord=['C'],
            unit=label,
            notext=notext, cmap = colormap,min=dMin,max=dMax)

    else:
        healpy.cartview(mask(skymap,top=67*degree,bottom=-90*degree), 
            fig=1,
            title=title, 
            rot=rotation,coord=['C'],
            unit=label,
            notext=notext, cmap = colormap,min=dMin,max=dMax)

 
    #set_color_bar(title="Neutral Hydrogen Flux [1/(cm$^2$ s sr keV)]", ticks=None, coord="C")
    fig = pylab.figure(1)
    for ax in fig.get_axes():
            if proj=='C0':
                ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="x-large")
                ax.annotate("360$^\circ$",xy=(-1.95, 0.625), size="x-large")
            elif proj=='C':
                ax.annotate("0$^\circ$", xy=(1.8, 0.625), size="x-large")
                ax.annotate("360$^\circ$",xy=(-1.95, 0.625), size="x-large")

    x = [227.28]
    y = [34.62]
    x =[0]
    y=[-89]
 
    labelsmap = {}
    labelsmap['eng'] = {'c-g': 'G-C', 'gc':'GC', 'ibex_center':'IBEX Ribbon','v_ism':r'$v_{ISM}$',
                        'v_lsr':r'$v_{LSR}$','limf':r'$\vec{B}_\mathrm{LIMF}$','antilimf':r'$-\vec{B}_\mathrm{LIMF}$',
                        'cygx':'Cygnus X-1','tail':'Heliotail'}
    labelsmap['esp'] = {'c-g': 'C-G', 'gc':'Cent. Galact.', 'ibex_center':'Cinta IBEX','v_ism':r'$v_{MIS}$',
                        'v_lsr':r'$v_{LSR}$', 'limf':r'$\vec{B}_\mathrm{LIMF}$','antilimf':r'$-\vec{B}_\mathrm{LIMF}$',
                        'cygx':'Cygnus X-1','tail':'Cola'}
    labels = labelsmap[lang]


    #print glimf.lat, glimf.lon

    
    #hp.projscatter(limferr_vec[0], limferr_vec[1], lonlat=True,color='grey', coord='E')
    if features:
       hp.projscatter([0,263.5520,195.1338,071.3350], [0, -02.7873,04.2658,03.0668], lonlat=True, coord='G')
       hp.projtext(195.1338+10, 04.2658+5, 'Geminga', lonlat=True, coord='G')
       hp.projtext(263.5520, -02.7873, 'Vela', lonlat=True, coord='G')
       hp.projtext(071.3350+1, 03.0668+1, labels['cygx'],lonlat=True, coord='G')
       hp.projtext(0, 0, 'GC', lonlat=True, coord='G',fontsize=15)


    #Exposure
    for hc in exposure_hawc:
        hawc_lat = np.zeros(360) + hc/degree
        hawc_lon = np.arange(360)
        hp.projplot(hawc_lon, hawc_lat, '--',lonlat=True, color='green' )
        #hp.projtext(90, hc/degree, 'HAWC', lonlat=True, coord='C',color='green',fontsize=12)
    for ic in exposure_icecube:
        ic_lat = np.zeros(360) + ic/degree
        ic_lon = np.arange(360)
        hp.projplot(ic_lon, ic_lat, '--',lonlat=True, color='black' )
        #hp.projtext(90, ic/degree, 'IceCube', lonlat=True, coord='C',color='black', fontsize=12)
        
    #2D fit
    if fitvalue>0:

        
        fitlat = np.arange(-90,90)
        fitlon = np.zeros(180)
        
        hp.projplot(fitlon+fitvalue, fitlat, '-',lonlat=True, color='red', linewidth=3 )
        hp.projplot(fitlon+180+fitvalue, fitlat, '-',lonlat=True, color='blue', linewidth=3 )
        if CG_corr is not None:
            hp.projplot(fitlon+fitvalue+CG_corr, fitlat, '-',lonlat=True, color='grey' )
            hp.projplot(fitlon+180+CG_corr+fitvalue, fitlat, '-',lonlat=True, color='grey' )
     
    
    #Show fit to circle
    if fit:
            
        cfit = ephem.Equatorial(fit[0],fit[1])
        gfit = ephem.Galactic(cfit)
        efit = ephem.Ecliptic(cfit)
        #print "fit(ecliptic)=", efit.lon/degree, efit.lat/degree
        
        #fitlat = np.arange(-90,90)
        #fitlon = np.zeros(180)
        
        #hp.projplot(fitlon+fit[0], fitlat, '-+',lonlat=True, color='red' )
        #hp.projplot(fitlon+180+fit[0], fitlat, '-+',lonlat=True, color='blue' )
            
        fiteqlat = np.zeros(360)
        fiteqlon = np.arange(360)
        hp.projplot(fiteqlon, fiteqlat, '-+',lonlat=True, color=fitcolor, rot=(cfit.ra/degree,cfit.dec/degree) )
        
        hp.projscatter([cfit.ra/degree], [cfit.dec/degree], lonlat=True,color=fitcolor, coord='C')
        hp.projtext(cfit.ra/degree,cfit.dec/degree-10, 'Fit',color=fitcolor, lonlat=True, coord='C',fontsize=18)#,bbox=bbox) 
    

    #37.8421671285


    gplon = np.arange(0,360)
    gplat = np.zeros(360)
    if gc:
        hp.projplot(gplon, gplat, 'r-',lonlat=True, coord='G')

    marker = 'b-'
    
    # BField
    bcolor = "black"
    x =[glimf.lon/degree]
    y=[glimf.lat/degree]


    #BFiled lines
    if fieldlines:
        for i in range(0,360,30):

            lmflat = np.arange(-90,90)
            lmflon = np.zeros(180)+i
            marker = '--'
            #f i == 210: marker = 'b-'
            hp.projplot(lmflon, lmflat, marker,lonlat=True,coord='G',rot=(glimf.lon/degree,glimf.lat/degree) ,color='gray')

            
        # BField-Equator
        lmfeqlon = np.arange(0,360)
        lmfeqlat = np.zeros(360)
        hp.projplot(lmfeqlon, lmfeqlat, 'b-',lonlat=True,coord='G',color=bcolor,rot=(glimf.lon/degree,glimf.lat/degree) )
    
        #B-V plane
        lmflat = np.arange(-90,90)
        lmflon = np.zeros(180)-31.25
        hp.projplot(lmflon, lmflat, 'b-',lonlat=True,coord='G',color=bcolor, rot=(glimf.lon/degree,glimf.lat/degree) )
        lmflon = (lmflon+180)%360
        hp.projplot(lmflon, lmflat, 'b-',lonlat=True,coord='G',color=bcolor, rot=(glimf.lon/degree,glimf.lat/degree) )
        hp.projtext([np.array(lmflon).mean()],[180+45], '$B-V$', lonlat=True,color=bcolor, coord='G', 
                rot=(glimf.lon/degree,glimf.lat/degree),fontsize=15)
    
    
    elif graticule:
        hp.graticule()
        

    #print glimf.lat, glimf.lon
    if features:
    
       hp.projscatter(x, y, lonlat=True,color=bcolor, coord='G')
       hp.projtext(climf.ra/degree-5,climf.dec/degree, labels['limf'],color=bcolor, fontweight='bold',lonlat=True, coord='C',fontsize=19)
    
       hp.projscatter(x, [glimf.lat/degree+180], lonlat=True,color=bcolor, coord='G')
       hp.projtext(glimf.lon/degree,glimf.lat/degree+180, labels['antilimf'],color=bcolor, fontweight='bold', lonlat=True, coord='G',fontsize=19)
    #hp.projscatter(x, y, lonlat=True,color='black', coord='G',rot=(glimf.lon/degree,glimf.lat/degree))        
        
       #ibex ribbon
       if show_ibex:
           hp.projplot(lmfeqlon, lmfeqlat+90-73.26, '-*',lonlat=True,coord='G', color='grey', rot=(gibex.lon/degree,gibex.lat/degree) )
           hp.projplot(lmfeqlon, lmfeqlat+90-79.18, '-*',lonlat=True,coord='G', color='grey', rot=(gibex.lon/degree,gibex.lat/degree) )
           hp.projscatter([gibex.lon/degree], [gibex.lat/degree], lonlat=True,color='black', coord='G')
           hp.projtext(gibex.lon/degree,gibex.lat/degree, labels['ibex_center'], lonlat=True, coord='G')
    

       #V_ism
       hp.projscatter(180+vism[0][0],-vism[1][0], lonlat=True, coord='C')
       hp.projtext(180+vism[0][0]-5,-vism[1][0] , labels['v_ism'], lonlat=True, coord='C',fontsize=18)
       hp.projscatter(vism[0],vism[1],lonlat=True,color='black', coord='C')
       hp.projtext(vism[0][0]-5,vism[1][0], "-"+labels['v_ism'], lonlat=True, coord='C',fontsize=18)
    
    
    
    #hp.projscatter(gtail.lon/degree,gtail.lat/degree, lonlat=True, coord='G')
    #hp.projtext(gtail.lon/degree,gtail.lat/degree , labels['tail'], lonlat=True, coord='G',fontsize=11)
    
    #print 'tail', gtail.lon/degree,gtail.lat/degree
    #print 'tail', gtail.lon/degree,gtail.lat/degree
    #print 'tail', ctail.ra/degree,ctail.dec/degree
    #print 'vism',vism[0][0],vism[1][0]
    #print '-vism', 180+vism[0][0],-vism[1][0]
    #print 'vism-tail',vism[0][0]-gtail.lon/degree,vism[1][0]-gtail.lat/degree
    
    xv=[5.25]
    yv=[12]
    solapex = [[47.9],[23.8]]
   
    if features:
       hp.projscatter(solapex[0],solapex[1],lonlat=True,color='black', coord='G')
       hp.projtext([solapex[0][0]+1],[solapex[1][0]+1], labels['v_lsr'], lonlat=True, coord='G',fontsize=18)
    

    


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
    #print ri(vcr_lon,0)
    lmflonrc = np.arange(-90,90)
    lmflatrc = np.zeros(180)
    lxrc = [];lyrc=[]
    if filename:
        fig.savefig(filename, dpi=100)
    plt.show()
    #return fig

"""
Solid angle
"""
def solid_angle(thetamin,thetamax):
    return 2*np.pi*(np.cos(thetamin)-np.cos(thetamax))

"""
Determine noise level
"""
def noise(m, dmin,dmax,nside=64):
    #nsum = np.sum(m)
    
    thetamin = -dmax+np.pi*.5
    thetamax = -dmin+np.pi*.5
    #print thetamin,thetamax
    nsum = 0.
    o = solid_angle(thetamin,thetamax)
    osqr = o*o*.25/np.pi
    
    npix = hp.nside2npix(nside)
    for ipix in range(npix):
        if m[ipix] == hp.UNSEEN or m[ipix] <= 0:
            continue
        theta, phi = hp.pix2ang(nside, ipix)
        if theta > thetamin and theta < thetamax:
            nsum+=1./m[ipix]*(o/(1.*npix))**2/4./np.pi
        
    return nsum
    
def find_max(m):

    npix = m.size
    nside = healpy.npix2nside(npix)
    
    imax = np.argmax(m)
    
    theta, phi = H.pix2ang(nside,imax)

 
    return 90-theta/degree, phi/degree              


def cdf(x,y):
    yy = []
    s = 0
    for yi in y:
       s += yi
       yy.append(s)
    yy = np.array(yy)*1.0/yy[-1]
    return np.vectorize(lambda xx: np.interp(xx,yy,x))

def get_1d_bin_ks_data(relintdata,bg, minDec=-30 * degree, 
		maxDec=-20 * degree, xbinsize=10):

    npix = relintdata.size
    nside = healpy.npix2nside(npix)
    bnpix = bg.size
    bnside = healpy.npix2nside(bnpix)
    nevents = 0
    dat = []
    weights = []

    res = hp.nside2resol(nside)

    dec_strip = healpy.query_strip(nside, 90*degree - maxDec, 90*degree - minDec, inclusive=True)
    for ipix in dec_strip:
            theta, phi = healpy.pix2ang(nside, ipix)
            phip = np.random.normal(phi,res)
            dat.append(phip)
            weights.append((1.0+relintdata[ipix])*bg[ipix])
            nevents += bg[ipix]

    return (np.array(dat), np.array(weights), nevents)


