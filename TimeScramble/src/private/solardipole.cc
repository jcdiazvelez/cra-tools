#include <solardipole.h>



double solar_dipole(double mjd, double ra, double dec) { 
        double dvh[3],dph[3]; // heliocentric coords 
        double dvb[3],dpb[3]; //place holders barycentric coords 
        float deqx = 2000.0; 
        double km = 1.496e+8; // conversion factor to from AU to km
        double gamma = 2.659; // # differential spectral index  \pm 0.061

        palEvp(mjd, deqx, dvb, dpb, dvh, dph);
        double vmag = sqrt(dvh[0]*dvh[0]+ dvh[1]*dvh[1]+dvh[2]*dvh[2])*km;

        double dr,dd;      // J2000.0 mean RA,Dec (radians)
        double sundr,sundd,diam;      // J2000.0 mean RA,Dec (radians)
        palRdplan(mjd, 0, 0, 0, &sundr, &sundd, &diam);

        double dl = 0.;   //  ecliptic longitude 
        double db = 0.;   //  ecliptic latitude

        palEqecl(sundr, sundd, mjd, &dl, &db);
        dl += 270.*degrees;   //  ecliptic longitude 

        //slaEqecl(dr, dd, mjd, &dl, &db);
        palEcleq(dl, db, mjd, &dr, &dd);
        //slaEcleq(0, 0, mjd, &sundr, &sundd);
        // calculate angular distance between 
        double theta = acos(sin(dec)*sin(dd)+cos(dec)*cos(dd)*cos(ra-dr));
        //double suntheta = acos(sin(sundd)*sin(dd)+cos(sundd)*cos(dd)*cos(sundr-dr));
        double speedratio = vmag/2.99792458e5; //# v/c

        double dipole = (gamma+2.0)*speedratio*cos(theta);
        //printf("%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n", 
        //                ra,dec,dr,dd,sundr,sundd,theta*rad2deg,suntheta*rad2deg, 1-dipole);
        return 1-dipole;
}
