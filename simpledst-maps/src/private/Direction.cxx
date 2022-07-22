
// $Id$

#include <Direction.h>
#include <cmath>
#include <units.h>



//-----------------------------------------------------------
void Direction::SetThetaPhi(double theta, double phi)
{
  if (theta>constants::pi) theta = 2.*constants::pi-theta;
  zenith_ = constants::pi-theta;
  double azimuth = constants::pi+phi;
  if (azimuth>=2.*constants::pi) azimuth -= 2*constants::pi;
  azimuth_=azimuth;
  isCalculated_=false;
}

//-----------------------------------------------------------
void Direction::RotateX(double angle)
{
// Rotate around x-axis by angle
  if (!isCalculated_) CalcCarFromSph();
  const double s=std::sin(angle);
  const double c=std::cos(angle);
  const double y=yDir_;
  yDir_=c*y-s*zDir_;
  zDir_=s*y+c*zDir_;
  CalcSphFromCar();
}

//-----------------------------------------------------------
void Direction::RotateY(double angle)
{
// Rotate around y-axis by angle
  if (!isCalculated_) CalcCarFromSph();
  const double s=std::sin(angle);
  const double c=std::cos(angle);
  const double z=zDir_;
  zDir_=c*z-s*xDir_;
  xDir_=s*z+c*xDir_;
  CalcSphFromCar();
}

//-----------------------------------------------------------
void Direction::RotateZ(double angle)
{
// Rotate around z-axis by angle
  if (!isCalculated_) CalcCarFromSph();
  const double s=std::sin(angle);
  const double c=std::cos(angle);
  const double x=xDir_;
  xDir_=c*x-s*yDir_;
  yDir_=s*x+c*yDir_;
  CalcSphFromCar();
}


//-----------------------------------------------------------
void Direction::CalcCarFromSph() const
{
  // Calculate Cartesian coordinates from Spherical
  // Direction is stored on disk in Spherical coordinates only.
  // theta=pi-zenith and phi=azimuth-pi in these IceCube coordinates.
  const double theta = constants::pi-zenith_;
  const double phi = azimuth_-constants::pi;
  const double rho = std::sin(theta);
  xDir_ = rho*std::cos(phi);
  yDir_ = rho*std::sin(phi);
  zDir_ = std::cos(theta);
  isCalculated_=true;
}

//-----------------------------------------------------------
void Direction::CalcSphFromCar()
{
  // Calculate Spherical coordinates from Cartesian
  // Direction is stored on disk in Spherical coordinates only
  // zenith=pi-theta and azimuth=phi+pi in these IceCube coordinates.
  const double r = std::sqrt(xDir_*xDir_+yDir_*yDir_+zDir_*zDir_);
  double theta = 0.;
  if (r && std::abs(zDir_/r)<=1.) {
    theta=std::acos(zDir_/r);
  } else {
    if (zDir_<0.) theta=constants::pi;
  }
  if (theta<0.) theta+=2.*constants::pi;
  double phi=0;
  if ((xDir_!=0.) || (yDir_!=0.)) phi=std::atan2(yDir_,xDir_);
  if (phi<0.) phi+=2.*constants::pi;
  double zenith = constants::pi-theta;
  double azimuth = phi+constants::pi;
  if (zenith>constants::pi) zenith -= constants::pi-(zenith-constants::pi);
  azimuth -= (int)(azimuth/(2.*constants::pi))*(2.*constants::pi);

  zenith_=zenith;
  azimuth_=azimuth;

  isCalculated_=false;
}


