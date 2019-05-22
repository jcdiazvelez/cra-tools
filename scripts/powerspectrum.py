import cratools
from cratools import *
import pylab

watermark = "PRELIMINARY"
marker_size=9

# Read Iter data
datdir = "/data/maps/new_hawcv8.sdpatm_icecube_all/"
combined_llh_iter20 = healpy.read_map(datdir+"//CR_HAWC-IceCube_64_360_iteration20.fits.gz")

# Smooth
combined_llh_smooth_iter20 = TopHatSmooth(64,combined_llh_iter20,radius=5.0*degree,average=True)

# Generate isotropic band
lh_ymin,lh_ymax = isotropic_band(combined_llh_iter20, n=100, lmax=40, percentile=90, isomap_dir="/data/maps/isomaps.v8",lhiter=20)
lh_ymin,lh_ymax = isotropic_band(combined_llh_iter20, n=100, lmax=40, percentile=90, isomap_dir="/data/maps/isomaps.v8",lhiter=20)

# Generate error bars
yerr,chisq = error_bars(combined_llh_iter20, n=100, lmax=40, poissonmap_dir="/data/maps/poisson.v8",lhiter=20)
yerr_l3,chisq_l3 = error_bars(combined_llh_iter20, n=100, lmin=4, lmax=40, poissonmap_dir="/data/maps/poisson.v8",lhiter=10)

# Generate powerspectrum
combined_lhreco_ccl = hp.anafast(mask(combined_llh_iter20, top=70 * degree, bottom=-85 * degree), lmax=40,iter=10)
combined_lhreco_ell = np.arange(len(combined_lhreco_ccl))

alm3 = hp.map2alm(combined_llh_iter20, lmax = 3) 
fitmap3 = healpy.alm2map(alm3, 64)

combined_lhreco_l3_ccl = hp.anafast(mask(combined_llh_iter20-fitmap3, top=70 * degree, bottom=-80 * degree), lmax=40,iter=10)
combined_lhreco_l3_ell = np.arange(len(combined_lhreco_l3_ccl))


# Make plot
fig = plt.figure(figsize=(12, 7))
ax = plt.subplot(111)

# errors
ax.errorbar(combined_lhreco_ell[1:], 
            combined_lhreco_ccl[1:]/kcorrection(top=70 * degree, bottom=-90*degree),
            yerr=yerr[1:],fmt='o',
            label="LHReco",color='black',markersize=marker_size)

# ell=3 subtracted errors
ax.errorbar(combined_lhreco_l3_ell[1:], 
            combined_lhreco_l3_ccl[1:]/kcorrection(top=70 * degree, bottom=-90*degree),
            yerr=yerr_l3[1:],fmt='p', 
            label="LHReco sans $\ell$=1,2,3", color='red',markersize=marker_size)

# Isotropic band
ax.fill_between(
       combined_lhreco_ell[1:],
       lh_ymin[1:]/kcorrection(top=70 * degree, bottom=-90*degree), 
       lh_ymax[1:]/kcorrection(top=70 * degree, bottom=-90*degree),
       facecolor='gray', alpha=.3)

# Add legened
plt.legend(loc='upper right', numpoints=1, fontsize=18)

# Set labels
for item in (ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(20)
                    
plt.xlabel('$\ell$',fontsize=30);
plt.ylabel('$\~c_\ell$',fontsize=32); #plt.grid()


# Add labels
x0 = [0,40]
hawc_noise = 2.94556668757e-10
icecube_noise = 1.85329616887e-10
y0 = [hawc_noise,hawc_noise]
y1 = [icecube_noise,icecube_noise]

ax.plot(x0, y0, '-', color='gray')
ax.plot(x0, y1, '-', color='gray')

ax.text(.5, hawc_noise, 'HAWC300 noise level', ha= 'left', fontsize=14)
ax.text(.5, icecube_noise*.6, 'IC86 noise level', ha= 'left', fontsize=14)

# Watermark
ax.text(0.25*(0+40), 0.3*(7e-11+1e-5), watermark, alpha=0.3,ha= 'left', fontsize=26)

ax.set_yscale("log")

# Add second x axis on top with angular scale
ax2 = ax.twiny()  # ax2 is responsible for "top" axis and "right" axis
angle_ticks = [180,45,20,10,5]
angle_l_ticks = map(lambda x:180/x,angle_ticks )
ax2.set_xticks([1]+angle_l_ticks+[40])
angle_ticks_labels = map(lambda x:"%u$^\circ$"%x,angle_ticks)
ax2.set_xticklabels([""]+angle_ticks_labels+[], fontsize=17)
ax.set_xlim(0, 40)
ax.set_ylim(7e-11, 1e-5)
ax.grid(True)
fig.savefig("CRAI3HW2yr_ang_power_spectrum.pdf", dpi=100)
plt.show()
