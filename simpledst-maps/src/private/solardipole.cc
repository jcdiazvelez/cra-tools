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
        static double speedratio = 0.0;
        static double dr_ = 0; // diple RA
        static double dd_ = 0; // dipole Dec

        static int imjd_  = 0; 

        // update sun velocity & position (Eq) every 0.1 days (2.4h)
        if (imjd_ == 0) std::cout << "initializing solar dipole" << std::endl;

        int imjd  = static_cast<int>(floor(mjd*100));
        if (imjd_ != imjd) {

            imjd_  = imjd;
            palEvp(mjd, deqx, dvb, dpb, dvh, dph);
            double vmag = sqrt(dvh[0]*dvh[0]+ dvh[1]*dvh[1]+dvh[2]*dvh[2])*km;
            speedratio = vmag/2.99792458e5; //# v/c

            double dr,dd;      // J2000.0 mean RA,Dec (radians)
            double sundr,sundd,diam;      // J2000.0 mean RA,Dec (radians)
            //palRdplan(mjd, 0, 0, 0, &sundr, &sundd, &diam);
            palRdplan(mjd, 0, 0, -0.5*constants::pi, &sundr, &sundd, &diam);

            double dl = 0.;   //  ecliptic longitude 
            double db = 0.;   //  ecliptic latitude

            palEqecl(sundr, sundd, mjd, &dl, &db);
            dl += 270.*units::degree;   //  ecliptic longitude 

            //slaEqecl(dr, dd, mjd, &dl, &db);
            palEcleq(dl, db, mjd, &dr, &dd);
            dr_ = dr;
            dd_ = dd;

        }
        //slaEcleq(0, 0, mjd, &sundr, &sundd);
        // calculate angular distance between 
        double theta = acos(sin(dec)*sin(dd_)+cos(dec)*cos(dd_)*cos(ra-dr_));
        //double suntheta = acos(sin(sundd)*sin(dd)+cos(sundd)*cos(dd)*cos(sundr-dr));

        double dipole = (gamma+2.0)*speedratio*cos(theta);

        if (soc) { // apply second-order correction
            return 1-dipole * (1-solar_2nd_order(dec));
        }
        //printf("%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n", 
        //                ra,dec,dr,dd,sundr,sundd,theta*rad2deg,suntheta*rad2deg, 1-dipole);
        return 1-dipole;
}




