#include <solardipole.h>
#include <iostream>
#include <units.h>


namespace polyfit { 
    double curve1(double x) { 
        /*
         * Coefficients from fit done by Thomas Yan
         */
        double coeff[4] = { 
            -0.000000014952662,
            0.000005520043609,
            -0.000679788839714,
            0.027953661402772
        };
        return x*(x*(coeff[0]*x +coeff[1])+coeff[2])+coeff[3];
    } 
    
    double curve2(double x) { 
        /*
         * Coefficients from fit done by Thomas Yan
         */
        double coeff[4] = { 
            -0.000000009321591,
             0.000004212543836,
            -0.000633849115097,
            0.031730170853338
        };
        return x*(x*(coeff[0]*x +coeff[1])+coeff[2])+coeff[3];
    } 
    
    double line(double x) { 
        /*
         * Coefficients from fit done by Thomas Yan
         */
        return 0.000000835883*x -0.000152000265487;
    } 

    double polyfit(double x) { 
        if (x > 154)
            return line(x);
        else if (x > 134) 
            return curve2(x);
        return curve1(x);
    }
}


double solar_2nd_order(double dec) { 
        return polyfit::polyfit(dec/units::degree);
}

double solar_dipole(double mjd, double ra, double dec, bool soc) { 
        double dvh[3],dph[3]; // heliocentric coords 
        double dvb[3],dpb[3]; //place holders barycentric coords 
        float deqx = 2000.0; 
        double km = 1.496e+8; // conversion factor to from AU to km
        double gamma = 2.659; // # differential spectral index  \pm 0.061

        palEvp(mjd, deqx, dvb, dpb, dvh, dph);
        double vmag = sqrt(dvh[0]*dvh[0]+ dvh[1]*dvh[1]+dvh[2]*dvh[2])*km;

        double dr,dd;      // J2000.0 mean RA,Dec (radians)
        double sundr,sundd,diam;      // J2000.0 mean RA,Dec (radians)
        //palRdplan(mjd, 0, 0, 0, &sundr, &sundd, &diam);
        palRdplan(mjd, 0, 0, -0.5*constants::pi, &sundr, &sundd, &diam);
        std::cout << "mjd "<< mjd << ", Sun: lon " << sundr/units::degree << ", lat " << sundd/units::degree;

        double dl = 0.;   //  ecliptic longitude 
        double db = 0.;   //  ecliptic latitude

        palEqecl(sundr, sundd, mjd, &dl, &db);
        std::cout << ", SunEc: lon " << dl/units::degree << ", lat " << db/units::degree ;
        dl += 270.*units::degree;   //  ecliptic longitude 

        //slaEqecl(dr, dd, mjd, &dl, &db);
        palEcleq(dl, db, mjd, &dr, &dd);
        std::cout << ", SD_eq: lon " << dr/units::degree << ", lat " << dd/units::degree << std::endl;
        //slaEcleq(0, 0, mjd, &sundr, &sundd);
        // calculate angular distance between 
        double theta = acos(sin(dec)*sin(dd)+cos(dec)*cos(dd)*cos(ra-dr));
        //double suntheta = acos(sin(sundd)*sin(dd)+cos(sundd)*cos(dd)*cos(sundr-dr));
        double speedratio = vmag/2.99792458e5; //# v/c

        double dipole = (gamma+2.0)*speedratio*cos(theta);

        if (soc) { // apply second-order correction
            dipole *= (1-solar_2nd_order(dec));
        }
        //printf("%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n", 
        //                ra,dec,dr,dd,sundr,sundd,theta*rad2deg,suntheta*rad2deg, 1-dipole);
        return 1-dipole;
}




