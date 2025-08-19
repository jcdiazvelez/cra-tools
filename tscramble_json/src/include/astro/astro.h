#ifndef ASTRO_H_INCLUDED
#define ASTRO_H_INCLUDED

#include <units.h>

/**
 * @brief High-level functions for astronomical conversions
 * for IceCube datatypes 
 *
 * @copyright (C) 2015 The IceCube Collaboration
 * 
 * $Id$
 *
 * @file Astro.h
 * @author Kevin Meagher
 * @date August 2015
 * $Rev$
 * $LastChangedBy$
 * $LastChangedDate$
 *
 */

class Direction;

///Nominal latitude of the icecube detector in radians
const double ICECUBE_LATITUDE  = -89.9944*units::degree;
///Nominal longitude of the icecube detector in radians
const double ICECUBE_LONGITUDE = -62.6081*units::degree;


/**
 * @brief Container class for storing equatorial (J2000) coordinates
 */
struct Equatorial 
{
 Equatorial() : ra(NAN),dec(NAN) {};
 Equatorial(double r,double d) : ra(r),dec(d) {}; 
  
  /// Right Ascension in radians
  double ra;
  /// Declination in radians in radians 
  double dec;
};

/**
 * @brief Container class for storing galactic coordinates (IAU 1958)
 */
struct Galactic 
{
 Galactic() : l(NAN),b(NAN) {};
 Galactic(double gl,double gb) : l(gl),b(gb) {};
  
  /// Galactic Longitude in radians
  double l;
  /// Galactic Latitude in radians
  double b;
};

/**
 * @brief Container class for storing supergalactic (de Vaucouleurs) coordinates
 */
struct SuperGalactic 
{
 SuperGalactic() : l(NAN),b(NAN) {};
 SuperGalactic(double SGL,double SGB) : l(SGL),b(SGB) {};

  /// Supergalactic Longitude in radians 
  double l;
  /// Supergalactic Latitude in radians
  double b;
};

/**
 * @brief Container class for storing ecliptic coordinates 
 */
struct Ecliptic
{
 Ecliptic() : l(NAN),b(NAN) {};
 Ecliptic(double el,double eb) : l(el),b(eb) {};
  
  /// Galactic Longitude in radians
  double l;
  /// Galactic Latitude in radians
  double b;
};



/**
 * @brief Convert from double to sidereal mjd (GMST)
 *
 * @param eventdouble - the mjd of the observation in double
 * @returns GMST in decimal hours
 *
 */
double GetGMST(const double eventdouble);


/**
 * @brief Convert from double to antisidereal mjd (GMAST)
 *
 * @param eventdouble - the mjd of the observation in double
 * @returns GMAST in decimal hours
 *
 */
double GetGMAST(const double eventdouble);

/**
 * @brief Convert from double to extended-sidereal mjd (GMEST)
 *
 * @param eventdouble - the mjd of the observation in double
 * @returns GMEST in decimal hours
 *
 */

double GetGMEST(const double eventdouble);

/**
 * @brief Gets the direction of the Moon in local IceCube coordinates 
 * at mjd
 *
 * @param mjd - the mjd of the observation in double
 * @returns Direction object
 */
Direction GetMoonDirection(const double mjd);


/**
 * @brief Gets the direction of the Sun in local IceCube coordinates 
 * at mjd
 *
 * @param mjd - double of the astronomical observation
 * @returns Direction object
 */
Direction GetSunDirection(const double mjd);

/**
 * @brief Convert an astronomical position in J2000 equatorial 
 * coordinates to local IceCube coordnates at a given mjd
 * 
 * @param equatorial - Equatorial position of astronomical observation 
 *        in J2000 coordinate system
 * @param mjd - double of the astronomical observation
 * @returns Direction containing local IceCube direction
 */
Direction GetDirectionFromEquatorial(const Equatorial& equatorial,const double mjd);

/**
 * @brief Convert an IceCube direction to 
 * J2000 equatorial coordinates system at a given mjd
 *
 * @param direction - IceCube local coordinate to convert
 * @param mjd - double of the astronomical observation
 * @returns Equatorial astronomical position of the 
 *          event in J2000 coordinate system
 */
Equatorial GetEquatorialFromDirection(const Direction& direction,const double mjd);

/**
 * @brief Convert from Equatorial (J2000) to Ecliptic coordinate system
 * 
 * @param equatorial - Equatorial position of astronomical observation 
 *        in J2000 coordinate system
 * @param mjd - double of the astronomical observation
 * @returns Ecliptic position in ecliptic coordinate system
 */
Ecliptic GetEclipticFromEquatorial(const Equatorial& equatorial, const double mjd);

/**
 * @brief Convert from Ecliptic to Equatorial (J2000) coordinate system
 * 
 * @param ecliptic - Ecliptic position in galactic coordinate system
 * @param mjd - double of the astronomical observation
 * @returns Equatorial position of astronomical observation 
 *          in J2000 coordinate system
 */
Equatorial GetEquatorialFromEcliptic(const Ecliptic& ecliptic, const double mjd);



/**
 * @brief Convert from Equatorial (J2000) to Galactic (IAU 1958) coordinate system
 * 
 * @param equatorial - Equatorial position of astronomical observation 
 *        in J2000 coordinate system
 * @returns Galactic position in galactic coordinate system
 */
Galactic GetGalacticFromEquatorial(const Equatorial& equatorial);

/**
 * @brief Convert from Galactic (IAU 1958) to Equatorial (J2000) coordinate system
 * 
 * @param galactic - Galactic position in galactic coordinate system
 * @returns Equatorial position of astronomical observation 
 *          in J2000 coordinate system
 */
Equatorial GetEquatorialFromGalactic(const Galactic& galactic);

/**
 * @brief Convert from Galactic (IAU 1958) to SuperGalactic (de Vaucouleurs) coordinate system
 *
 * @param galactic - Galactic position in galactic coordinate system
 * @returns SuperGalactic position of astronomical observation
 *
 */
SuperGalactic GetSuperGalacticFromGalactic(const Galactic& galactic);

/**
 * @brief Convert from SuperGalactic (de Vaucouleurs) to Galactic (IAU 1958) coordinate system
 *
 * @param supergalactic - SuperGalactic position in SuperGalactic coordinate system
 * @returns Galactic position of astronomical observation in IAU 1958 coordinate system
 *
 */
Galactic GetGalacticFromSuperGalactic(const SuperGalactic& supergalactic);

/**
 * @brief Convert from Equatorial (J2000) to SuperGalactic (de Vaucouleurs) coordinate system
 *
 * @param equatorial - Equatorial position in Equatorial (J200) coordinate system
 * @returns SuperGalactic position of astronomical observation in de Vaucouleurs coordinate system
 *
 */
SuperGalactic GetSuperGalacticFromEquatorial(const Equatorial& equatorial);

/**
 * @brief Convert from SuperGalactic (de Vaucouleurs) to Equatorial (J2000) coordinate system
 *
 * @param supergalactic - SuperGalactic position in Supergalactic (de Vaucouleurs) coordinate system
 * @returns Equatorial position of astronomical observation in J2000 coordinate system
 *
 */
Equatorial GetEquatorialFromSuperGalactic(const SuperGalactic& supergalactic);




#endif //ASTRO_H_INCLUDED

