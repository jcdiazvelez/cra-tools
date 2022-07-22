/**
 * @brief 
 *
 * @copyright (C) 2015 The Icecube Collaboration
 * 
 * $Id$
 *
 * @file Astro.cxx
 * @author Kevin Meagher
 * @date August 2015
 * $Rev$
 * $LastChangedBy$
 * $LastChangedDate$
 * 
 */

#include <star/pal.h>
#include <Direction.h>
#include <astro/astro.h>
#include <units.h>
#include <cassert>

namespace {

  /**
   * Floored modular division rounding toward negative infinity
   *
   * It is totally ridiculous that c++ doesn't have this function
   * but the relativitly useless fmod and remainder are part of
   * the specification
   */
  double floored_division_remainder(double x,double y)
  {
    return fmod(fmod(x,y)+y,y);
  }

  /** 
   * Convert apparent equatorial cooridnates to local coordinates
   * This function is hidden in an anonymous namespace because 
   * apparent coordinates are too similar to J2000 and it confuses 
   * people.
   */
  Direction direction_from_apparent(double ra_rad,
					double dec_rad,
					double mjd)
  {
    double PalAz=NAN;
    double PalZen=NAN;
    double PalHA=NAN;
    double PalDec=NAN;
    double PalRA=NAN;

    // SLA_AOP : Apparent to Observed
    // ACTION  : Apparent to observed place, for sources distant from the solar system.
    palAop(ra_rad,dec_rad,   //geocentric apparent [α,δ] (radians)
	   mjd,              // date/time (Modified Julian Date, JD−2400000.5)
	   0,                //UT1−UTC (UTC seconds)
	   ICECUBE_LONGITUDE,//observer’s mean longitude (radians, east +ve)
	   ICECUBE_LATITUDE, //observer’s mean geodetic latitude (radians)
	   0,                //observer’s height above sea level (metres)
	   0,0,              //polar motion [x,y] coordinates (radians)
	   0,                //local ambient temperature (K; std=273.15D0)
	   0,                //local atmospheric pressure (mb; std=1013.25D0)
	   0,                //local relative humidity (in the range 0D0 – 1D0)
	   0,                //effective wavelength (μm, e.g. 0.55D0)
	   0,                //tropospheric lapse rate (K per metre, e.g. 0.0065D0)
	   //RETURNED:
	   &PalAz,           //observed azimuth (radians: N=0, E=90∘)
	   &PalZen,          //observed zenith distance (radians)
	   &PalHA,           //observed Hour Angle (radians)
	   &PalDec,          //observed δ (radians)
	   &PalRA            //observed α (radians)
	   );
	   
    
    //convert to IceCube Coordinates
    double ic_az  = M_PI/2. - PalAz - ICECUBE_LONGITUDE;
    
    //normaize azimuth to [0, 360)
    ic_az = floored_division_remainder(ic_az,2*M_PI);
    
    assert ( ic_az  <= 2*M_PI);
    assert ( ic_az  >= 0 );
    assert ( PalZen <= M_PI);
    assert ( PalZen >= 0);
    
    return Direction(PalZen,ic_az);
  }

  Direction direction_of_planet(int planet_number, const double mjd)
  {
    //PAL_DTT : TT minus UTC PAL_DTT
    //ACTION  : Compute ∆TT, the increment to be applied to Coordinated
    //          Universal Time UTC to give Terrestrial Time TT.
    double tt = mjd+palDtt(mjd)/86400.;

    double raRad=NAN;
    double decRad=NAN;
    double diameterRad=NAN;
    //PAL_RDPLAN : Apparent [ α, δ ] of Planet PAL_RDPLAN
    //ACTION     : Approximate topocentric apparent [ α, δ ]
    //           : and angular size of a planet.
    palRdplan( tt,               //MJD of observation (JD−2400000.5) TT can be used 
	       planet_number,    // 0 = Sun, 3 = Moon
	       ICECUBE_LONGITUDE,// observer’s longitude (east +ve)
	       ICECUBE_LATITUDE, // and latitude (radians)
	       //RETURNED
	       &raRad, &decRad,  //topocentric apparent [α,δ] (radians)
	       &diameterRad      //angular diameter (equatorial, radians)
	       );
    
    return direction_from_apparent(raRad,decRad,mjd);
  }
  
}//anonymous namespace

double GetGMST(const double mjd)
{
  // NOTE: MJD is expected with respect to UT
  double gmst = palGmst(mjd)*units::radian/units::degree;
  gmst = fmod(gmst, 360);
  gmst = (gmst < 0.) ? gmst + 360: gmst;

  return gmst*24/360; // in hours
}

double GetGMAST(const double mjd)
{
  // NOTE: MJD is expected with respect to UT
  const double gmst = GetGMST(mjd);
  double dt = mjd - gmst;
  double ast = mjd + dt;
  ast = fmod(ast, 24);
  ast = (ast < 0.) ? ast + 24: ast;

  return ast; // in hours
}

double GetGMEST(const double mjd)
{
  // NOTE: MJD is expected with respect to UT
  const double gmst = GetGMST(mjd);
  double dt = mjd - gmst;
  double ast = mjd - 2*dt;
  ast = fmod(ast, 24);
  ast = (ast < 0.) ? ast + 24: ast;

  return ast; // in hours
}



Direction GetMoonDirection(const double mjd)
{
  //Moon = 3
  return direction_of_planet(3,mjd);
}

Direction GetSunDirection(const double mjd)
{
  //Sun = 0
  return direction_of_planet(0,mjd);
}

Direction GetDirectionFromEquatorial(const Equatorial &eq, const double mjd)
{

  //PAL_DTT : TT minus UTC PAL_DTT
  //ACTION  : Compute ∆TT, the increment to be applied to Coordinated
  //          Universal Time UTC to give Terrestrial Time TT.
  double tt = mjd+palDtt(mjd)/86400.;
  
  double raRad  = NAN;
  double decRad = NAN;
  //SLA_MAP : Mean to Apparent SLA_MAP
  //ACTION  : Transform star [ α, δ ] from mean place to geocentric apparent. The reference
  //        : frames and timescales used are post IAU 1976.
  palMap( eq.ra,    //mean RA 
  	  eq.dec,   //mean dec
  	  0.0,      //proper motions RA
  	  0.0,      //proper motions Dec
  	  0.0,      //parallax			       
  	  0.0,      // radial velocity
  	  2000.0,   //epoch and eqinox of star data (Julian)
  	  tt,       //TDB for apparent place (JD−2400000.5)
	  // RETURNS
  	  &raRad,  //apparent RA
  	  &decRad  //apparent Dec			       
  	  );

  return direction_from_apparent(raRad,decRad,mjd);
}

Equatorial GetEquatorialFromDirection(const Direction & dir, const double mjd)
{

  // Convert local coordinates from IceCube convention to Astronomy convention
  double PalAz = M_PI/2. - dir.GetAzimuth() - ICECUBE_LONGITUDE;

  double raRad  = NAN;
  double decRad = NAN;

  // SLA_OAP  : Observed to Apparent
  // ACTION   : Observed to apparent place.
  char coordtype[] ={'A'};
  palOap (coordtype,              //type of coordinates – ‘R’, ‘H’ or ‘A’ 
	  PalAz,            //observed Az, HA or RA (radians; Az is N=0, E=90∘)
	  dir.GetZenith(),  //observed zenith distance or δ (radians)
	  mjd,              //date/time (Modified Julian Date, JD−2400000.5)
	  0,                //UT1−UTC (UTC seconds)
	  ICECUBE_LONGITUDE,//observer’s mean longitude (radians, east +ve)
	  ICECUBE_LATITUDE, //observer’s mean geodetic latitude (radians)
	  0,                //observer’s height above sea level (metres)
	  0,0,              //polar motion [x,y] coordinates (radians)
	  0,                //local ambient temperature (K; std=273.15D0)
	  0,                //local atmospheric pressure (mb; std=1013.25D0)
	  0,                //local relative humidity (in the range 0D0 – 1D0)
	  0,                //effective wavelength (μm, e.g. 0.55D0)
	  0,                //tropospheric lapse rate (K per metre, e.g. 0.0065D0)
	  //RETURNED
	  &raRad,&decRad    //Returned geocentric apparent [α,δ]
	  );

  //SLA_AMP  : Apparent to Mean SLA_AMP
  //ACTION   : Convert star [ α, δ ] from geocentric apparent to mean place (post IAU 1976).
  double mra  = NAN;
  double mdec = NAN;
  palAmp(raRad,decRad, // apparent [ α, δ ] (radians)
	 mjd,          // TDB for apparent place (JD−2400000.5)
	 2000.0,       // equinox: Julian epoch of mean place
	 //RETURNED 
	 &mra, &mdec   // mean [ α, δ ] (radians)
	 );
  
  Equatorial eq(mra,mdec);

  assert(eq.ra  >=  0);
  assert(eq.ra  <=  2*M_PI);
  assert(eq.dec >= -M_PI/2);
  assert(eq.dec <= +M_PI/2);
  
  return eq;
}

Ecliptic GetEclipticFromEquatorial(const Equatorial& eq, const double mjd)
{
  double ecl_longitude(NAN);
  double ecl_latitude(NAN);

  //PAL_EQECL : J2000 α, δ to Ecliptic
  //ACTION    : Transform from J2000.0 equatorial coordinates 
  //          : to ecliptic coordinates 
  palEqecl( eq.ra, eq.dec, mjd, &ecl_longitude, &ecl_latitude);

  Ecliptic ecl(ecl_longitude, ecl_latitude);

  assert(ecl.l >=  0);
  assert(ecl.l <=  2*M_PI);
  assert(ecl.b >= -M_PI/2);
  assert(ecl.b <= +M_PI/2);
  
  return ecl;
}

Equatorial GetEquatorialFromEcliptic(const Ecliptic& ecl, const double mjd)
{
  double eq_ra(NAN);
  double eq_dec(NAN);

  //PAL_ECLEQ : Ecliptic to J2000 α, δ 
  //ACTION    : Transform from ecliptic coordinates 
  //          : to J2000.0 FK5 equatorial coordinates.
  palEcleq(ecl.l,ecl.b,mjd,&eq_ra,&eq_dec);
  
  Equatorial eq(eq_ra, eq_dec);

  assert(eq.ra  >=  0);
  assert(eq.ra  <=  2*M_PI);
  assert(eq.dec >= -M_PI/2);
  assert(eq.dec <= +M_PI/2);

  
  return eq;

}

Galactic GetGalacticFromEquatorial(const Equatorial& eq)
{
  double gal_longitude(NAN);
  double gal_latitude(NAN);
  
  //PAL_EQGAL : J2000 α, δ to Galactic
  //ACTION    : Transformation from J2000.0 FK5 equatorial coordinates
  //            to IAU 1958 galactic coordinates
  palEqgal( eq.ra, eq.dec, &gal_longitude, &gal_latitude);

  Galactic gal(gal_longitude,gal_latitude);

  assert(gal.l >=  0);
  assert(gal.l <=  2*M_PI);
  assert(gal.b >= -M_PI/2);
  assert(gal.b <= +M_PI/2);
  
  
  return gal;
}

Equatorial  GetEquatorialFromGalactic(const Galactic& gal)
{
  double eq_ra(NAN);
  double eq_dec(NAN);

  //PAL_GALEQ : Galactic to J2000 α, δ 
  //ACTION    : Transformation from IAU 1958 galactic coordinates
  //          : to J2000.0 FK5 equatorial coordinates.
  palGaleq(gal.l,gal.b,&eq_ra,&eq_dec);
  
  Equatorial eq(eq_ra, eq_dec);

  assert(eq.ra  >=  0);
  assert(eq.ra  <=  2*M_PI);
  assert(eq.dec >= -M_PI/2);
  assert(eq.dec <= +M_PI/2);
  
  return eq;
}

SuperGalactic GetSuperGalacticFromGalactic(const Galactic& gal)
{
  double sg_longitude(NAN);
  double sg_latitude(NAN);

  palGalsup( gal.l, gal.b, &sg_longitude, &sg_latitude );

  SuperGalactic SG(sg_longitude, sg_latitude);

  assert(SG.l >= 0);
  assert(SG.l <= 2*M_PI);
  assert(SG.b >= -M_PI/2);
  assert(SG.b <= M_PI/2);

  return SG;
}

Galactic GetGalacticFromSuperGalactic(const SuperGalactic& SG)
{
  double gal_longitude(NAN);
  double gal_latitude(NAN);

  palSupgal( SG.l, SG.b, &gal_longitude, &gal_latitude );

  Galactic gal(gal_longitude, gal_latitude);

  assert(gal.l >= 0);
  assert(gal.l <= 2*M_PI);
  assert(gal.b >= -M_PI/2);
  assert(gal.b <= M_PI/2);

  return gal;
}

SuperGalactic GetSuperGalacticFromEquatorial(const Equatorial& eq)
{
  Galactic gal = GetGalacticFromEquatorial(eq);
  SuperGalactic SG = GetSuperGalacticFromGalactic(gal);
  return SG;
}

Equatorial GetEquatorialFromSuperGalactic(const SuperGalactic& SG)
{
  Galactic gal = GetGalacticFromSuperGalactic(SG);
  Equatorial eq = GetEquatorialFromGalactic(gal);
  return eq;
}
