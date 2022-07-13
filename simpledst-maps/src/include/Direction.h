/**
 * copyright  (C) 2004
 * the icecube collaboration
 * @version $Id$
 * @file Direction.h
 * @date $Date$
 */

//***********************************************************
//-- Created: Dusan Turcan, UMD, Sep 2, 2004
//***********************************************************
#include <units.h>

// $Id$

#ifndef DIRECTION_H_INCLUDED
#define DIRECTION_H_INCLUDED

#define EPSILON 1e-15

bool inline approximatelyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool inline essentiallyEqual(float a, float b, float epsilon)
{
    return fabs(a - b) <= ( (fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool inline definitelyGreaterThan(float a, float b, float epsilon)
{
    return (a - b) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool inline definitelyLessThan(float a, float b, float epsilon)
{
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}



class Direction 
{
 public:

  Direction():
  zenith_(NAN),
  azimuth_(NAN),
  xDir_(NAN),
  yDir_(NAN),
  zDir_(NAN),
  isCalculated_(false)
  {}

  /**
   * Additional constructor: 2 arguments mean construct dir. with zen,azi
   */
  Direction(double zen, double azi):
  zenith_(zen),
  azimuth_(azi),
  xDir_(NAN),
  yDir_(NAN),
  zDir_(NAN),
  isCalculated_(false)
  {}

  /**
   * Additional constructor: 3 arguments mean construct dir. with x,y,z
   */
  Direction(double x, double y, double z):
  zenith_(NAN),
  azimuth_(NAN),  
  xDir_(x),
  yDir_(y),
  zDir_(z),
  isCalculated_(false)    
  {
    CalcSphFromCar();
  }
  
  

  //--------------

  /**
   * Store direction with theta, phi
   */
  void SetThetaPhi(double theta, double phi);

  //--------------

  /**
   * Provide Zenith of direction
   */
  inline double GetZenith() const {return zenith_;}

  /**
   * Provide Azimuth of direction
   */
  inline double GetAzimuth() const {return azimuth_;}

  /**
   * Provide X of direction in cartesian ref frame
   */
  inline double GetX() const {
    if (!isCalculated_) CalcCarFromSph();
    return xDir_;
  }

  /**
   * Provide Y of direction in cartesian ref frame
   */
  inline double GetY() const {
    if (!isCalculated_) CalcCarFromSph();
    return yDir_;
  }

  /**
   * Provide Z of direction in cartesian ref frame
   */
  inline double GetZ() const {
    if (!isCalculated_) CalcCarFromSph();
    return zDir_;
  }

  /**
   * Calculate Theta of direction
   */
  inline double CalcTheta() const {
    return constants::pi - static_cast<double>(zenith_);
  }

  /**
   * Calculate Phi of direction
   */
  inline double CalcPhi() const {
    double phi = constants::pi + static_cast<double>(azimuth_);
    if (phi >= 2.*constants::pi) phi -= 2.*constants::pi;
    return phi;
  }


  //--------------

  /**
   * Rotate direction around X axis by angle
   */
  void RotateX(double angle);

  /**
   * Rotate direction around Y axis by angle
   */
  void RotateY(double angle);

  /**
   * Rotate direction around Z axis by angle
   */
  void RotateZ(double angle);
  
  /**
   * Vector inversion (makes the vector point in the opposite direction)
   */
  Direction operator-() const{
    Direction d;
    d.zenith_=constants::pi-zenith_;
    d.azimuth_=azimuth_+constants::pi;
    if(d.azimuth_>=2*constants::pi)
      d.azimuth_-=2*constants::pi;
    d.isCalculated_=false;
    return d;
  }

  /**
   * Cross product of this x d
   */
  Direction Cross(const Direction& d) const{
    if (!isCalculated_) CalcCarFromSph();
    return Direction (yDir_*d.GetZ() - zDir_*d.GetY(),
                        zDir_*d.GetX() - xDir_*d.GetZ(),
                        xDir_*d.GetY() - yDir_*d.GetX());
  }
  

  /**
   * Scalar (dot) product
   */
  double operator*(const Direction& other) const{
    if(isCalculated_ && other.isCalculated_)
      return xDir_*other.xDir_ + yDir_*other.yDir_ + zDir_*other.zDir_;
    //if one of the directions doesn't already have cartesian
    //coordinates calculated, calculating them would be less efficient
    double cad=cos(azimuth_ - other.azimuth_);
    double czd=cos(zenith_ - other.zenith_);
    double czs=cos(zenith_ + other.zenith_);
    return 0.5*(cad*(czd-czs)+czd+czs);
  }
  
  

  bool operator==(const Direction& rhs) const {
    return approximatelyEqual(zenith_, rhs.zenith_, EPSILON);
  }
  bool operator!=(const Direction& rhs) const {
    return !(*this == rhs);
  }

  //---
  
  /** returns the angle between in two Directions
   * @return the angle in rad (native units)
   */
  inline double Angle(const Direction& rhs) const
  {
    return acos(*this*rhs);
  }

 protected:
  /**
   * direction coordinates -- spherical
   */ 
  double zenith_;
  double azimuth_;

  /**
   * direction coordinates -- cartesian (direction cosines)
   */ 
  mutable double xDir_; //!
  mutable double yDir_; //!
  mutable double zDir_; //!

  /**
   * Did we calculate the directions before?
   */
  mutable bool isCalculated_; //!

 private:
  /**
   * Change zenith,azimuth coordinates into x,y,z directional coordinates.
   * The three numbers x,y,z are calculated to add up (in quadrature) to 1.
   * theta=pi-zenith and phi=azimuth-pi in these IceCube coordinates.
   */
  void CalcCarFromSph() const;

  /**
   * Change x,y,z directional coordinates to zenith,azimuth coordinates.
   * zenith=pi-theta and azimuth=phi+pi in these IceCube coordinates.
   * The three numbers DO NOT have to add up (in quadrature) to 1.
   * Even if they don't, the direction that they define is stored.
   * SO BE CAREFUL IF YOU SPECIFY THE DIRECTION IN THIS WAY!
   */
  void CalcSphFromCar();


};




#endif //DIRECTION_H_INCLUDED
