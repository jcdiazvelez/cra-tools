#ifndef SOLARDIPOLE_H
#define SOLARDIPOLE_H


#include <star/pal.h>
#include <math.h>
#define degrees M_PI/180.


/*
 * solar dipole correction from theory 
 * \frac{\Delta I}{I} = (\gamma + 2) \frac{v}{c} cos(\theta)
 */
double solar_dipole(double mjd, double ra, double dec, bool soc=true);


/*
 * Second-order solar dipole correction (taken empirically from data) accounts
 * for angular resolution as a function of declination (zenith angle)
 */
double solar_2nd_order(double dec);




#endif
