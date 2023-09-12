
// $Id: Position.cxx 122541 2014-08-15 16:34:33Z nega $


#include <iostream>
#include "Position.h"
#include "Direction.h"
#include "units.h"
#include <iostream>



//-----------------------------------------------------------
Position::Position(const Direction& d):
x_(d.GetX()),
y_(d.GetY()),
z_(d.GetZ()),
r_(1),
theta_(d.CalcTheta()),
phi_(d.CalcPhi()),
rho_(sin(theta_)),
isCalculated_(true)
{}

//-----------------------------------------------------------
Position::Position(double r1, double r2, double r3, RefFrame frame){
  // Store position according to reference frame f
  isCalculated_=false;

  switch (frame) {
    case car: // Input given in Cartesian coordinates
      x_=r1;
      y_=r2;
      z_=r3;
      break;

    case sph: // Input given in Spherical coordinates
      r_=r1;
      theta_=r2;
      phi_=r3;
      CalcCarCylFromSph();
      break;

    case cyl: // Input given in Cylindrical coordinates
      rho_=r1;
      phi_=r2;
      z_=r3;
      CalcCarSphFromCyl();
      break;

    default: // Unsupported reference frame
      std::cerr << "Unsupported reference frame passed to Position::SetPosition: " <<frame <<std::endl;
      exit(1);
      break;
  }
}

//-----------------------------------------------------------
void Position::RotateX(double angle)
{
// Rotate around x-axis by angle
  const double s=std::sin(angle);
  const double c=std::cos(angle);
  const double y=y_;
  const double z=z_;
  y_=c*y-s*z;
  z_=s*y+c*z;
  isCalculated_ = false;
}

//-----------------------------------------------------------
void Position::RotateY(double angle)
{
// Rotate around y-axis by angle
  const double s=std::sin(angle);
  const double c=std::cos(angle);
  const double x=x_;
  const double z=z_;
  z_=c*z-s*x;
  x_=s*z+c*x;
  isCalculated_ = false;
}

//-----------------------------------------------------------
void Position::RotateZ(double angle)
{
// Rotate around z-axis by angle
  const double s=std::sin(angle);
  const double c=std::cos(angle);
  const double x=x_;
  const double y=y_;
  x_=c*x-s*y;
  y_=s*x+c*y;
  isCalculated_ = false;
}

//-----------------------------------------------------------
double Position::operator*(const Direction& rhs) const{
  return x_*rhs.GetX() + y_*rhs.GetY() + z_*rhs.GetZ();
}

//-----------------------------------------------------------
Position Position::Cross(const Direction& d) const{
  return Position (y_*d.GetZ() - z_*d.GetY(),
                     z_*d.GetX() - x_*d.GetZ(),
                     x_*d.GetY() - y_*d.GetX());
}

//-----------------------------------------------------------
void Position::CalcSphCylFromCar() const
{
  // Calculate Spherical and Cylindrical coordinates from Cartesian
  // Position is stored on disk in Cartesian coordinates only
  r_=std::sqrt(x_*x_+y_*y_+z_*z_);
  theta_=0;
  if ((r_!=0.) && std::abs(z_/r_)<=1.) {
    theta_=std::acos(z_/r_);
  } else {
    if (z_<0.) theta_=constants::pi;
  }
  if (theta_<0.) theta_+=2.*constants::pi;
  phi_=0;
  if ((x_!=0.f) || (y_!=0.f)) phi_=std::atan2(y_,x_);
  if (phi_<0.) phi_+=2.*constants::pi;
  rho_=r_*std::sin(theta_);
  isCalculated_=true;
}

//-----------------------------------------------------------
void Position::CalcCarCylFromSph()
{
  // Calculate Cartesian and Cylindrical coordinates from Spherical
  rho_=r_*std::sin(theta_);
  x_=rho_*std::cos(phi_);
  y_=rho_*std::sin(phi_);
  z_=r_*std::cos(theta_);
  isCalculated_=true;
}

//-----------------------------------------------------------
void Position::CalcCarSphFromCyl()
{
  // Calculate Cartesian and Spherical coordinates from Cylindrical
  r_=std::sqrt(rho_*rho_+z_*z_);
  if (phi_<0.) phi_+=2.*constants::pi;
  theta_=0;
  if ((r_!=0.) && std::abs(z_/r_)<=1.) {
    theta_=std::acos(z_/r_);
  } else {
    if (z_<0.) theta_=constants::pi;
  }
  if (theta_<0.) theta_+=2.*constants::pi;
  x_=rho_*std::cos(phi_);
  y_=rho_*std::sin(phi_);
  isCalculated_=true;
}

//-----------------------------------------------------------
Position operator*(double a, const Position& p){
  return(p*a);
}

//-----------------------------------------------------------
std::ostream& operator<<(std::ostream& oss, const Position& p){
  oss << "Position(" << p.GetX() << "," <<p.GetY() << "," << p.GetZ() << ")";
  return oss;
}

//-----------------------------------------------------------
double abs(const Position& p){
	return(p.Magnitude());
}
