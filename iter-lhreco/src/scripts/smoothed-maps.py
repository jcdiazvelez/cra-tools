import healpy as h
import numpy as n 
import pylab as p
import sys

degree = n.pi/180.

r  =h.read_map(sys.argv[1])
#d  =h.read_map(sys.argv[2])
prefix =       sys.argv[3]
radius = float(sys.argv[4])

nside = int(n.sqrt(len(r)/12))

ravg  = n.zeros(len(r))

#b     = d / (r+1.)
asigsm = n.zeros(len(r))
#dsm   = n.zeros(len(r))
#bsm   = n.zeros(len(r))
#rsm   = n.zeros(len(r))
#alpha = n.zeros(len(r))


#multipole subtracted
#alm = h.map2alm(r, lmax=3)
#fit = h.alm2map(alm, nside, lmax=3)
#rsub = n.zeros(len(r))
#sigsub = n.zeros(len(r))
#


for i in range(0,len(r)):
    pixels = h.query_disc(nside, h.pix2vec(nside,i), radius*degree)

    ravg[i] = sum( r[pixels]) / len(pixels)
   # rsub[i] = (sum(r[pixels]) - sum(fit[pixels])) / len(pixels)
   # dsm[i] = sum(d[pixels])
   # bsm[i] = sum(b[pixels])
   # rsm[i] = dsm[i]/bsm[i] -1.

   # ang=h.pix2ang(nside,i)
   # dec = n.pi/2. - ang[0]
   # 
   # onexposure  = 2.*n.pi*(1-n.cos(degree*radius))
   # offexposure = 15.*degree*24.*(n.cos(n.pi/2.-dec-radius*degree) - n.cos(n.pi/2.-dec+radius*degree))
   # alpha[i]= onexposure/offexposure

   # offcorrected =   bsm[i] - alpha[i]*(dsm[i]-bsm[i])
   # if offcorrected<0.:
   #   offcorrected=0.;

   # if (offcorrected > 0.0):
   #   xon   =  dsm[i]
   #   xoff  =  offcorrected/alpha[i]
   #   Z = 0.
   #   if xon > 0. and xoff> 0.:
   #     logterm1 = ((1.+alpha[i])/alpha[i])*(xon/(xon+xoff))
   #     logterm2 = (1.+alpha[i])*(xoff/(xon+xoff))
   #     sqarg    = xon*n.log(logterm1) + xoff*n.log(logterm2)
   #     if sqarg < 0.:
   #       sqarg = 0.;
   #     Z = n.sqrt(2.)*n.sqrt(sqarg);
   #     if offcorrected > xon:
   #       Z = -Z
   #      
   #   sigsm[i]= Z;
   # 


#h.write_map("%s_S%02d_sig.fits.gz"   %(prefix,int(radius)),sigsm)
#h.write_map("%s_S%02d_alpha.fits.gz" %(prefix,int(radius)),alpha)
#h.write_map("%s_S%03d_relint.fits.gz"%(prefix,int(radius)),rsm)
h.write_map("%s_A%02d_relint.fits.gz"%(prefix,int(radius)),ravg)
#h.write_map("%s_A%02d_relint_multsub.fits.gz"%(prefix,int(radius)),rsub)
