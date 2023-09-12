/**
 * copyright  (C) 2004
 * the icecube collaboration
 * @version $Id: Position.h 118949 2014-04-16 19:38:12Z nega $
 * @file Position.h
 * @date $Date: 2014-04-16 14:38:12 -0500 (Wed, 16 Apr 2014) $
 */

//***********************************************************
//-- Created: Dusan Turcan UMD 26-05-2004
//   Taken from: Nick van Eijndhoven 06-feb-1999 UU-SAP Utrecht
//***********************************************************

// $Id: Position.h 118949 2014-04-16 19:38:12Z nega $

#ifndef I3POSITION_H_INCLUDED
#define I3POSITION_H_INCLUDED

#include <cmath>
#include <sstream>

/**
 * @brief The basic position class for IceCube.
 *
 * All positions in IceCube should be written with this class.
 * Positions can be given in cartesian, spherical, or cylindrical coordinates.
 *
 * @todo implement "print out" of all information in a uniform way...
 * @todo insure that the temporary data isn't written to disk.
 */

//Forward declaration
class Direction;

class Position 
{
 public:

  /**
   * Possible reference frames.
   */
  enum RefFrame { car = 0, sph = 1, cyl = 2 };

  //--------------

  /**
   * Default constructor
   */
  Position():
  x_(NAN),
  y_(NAN),
  z_(NAN),
  isCalculated_(false)
  {}

  /**
   * Constructor for different coordinate systems
   *
   * The meaning of this constructor depends on the value of RefFrame.
   * If it is Position::car the three coordinates are treated as
   *   cartesian x, y, and z
   * If it is Position::sph the three coordinates are treated as
   *   spherical r, theta, and phi
   * If it is Position::cyl the three coordinates are treated as
   *   cylindrical rho, phi, and z
   */
	Position(double x, double y, double z, RefFrame f);

  /**
   * Constructor from cartesian coordinates
   */
  Position(double x, double y, double z):
  x_(x),
  y_(y),
  z_(z),
  isCalculated_(false)
  {}

  /**
   * Copy constructor
   */
  Position(const Position& p):
  x_(p.x_),
  y_(p.y_),
  z_(p.z_),
  isCalculated_(false)
  {}

  explicit Position(const Direction& d);

  //--------------

  /**
   * Provide X of position in cartesian ref frame
   */
  inline double GetX() const {return x_;}

  /**
   * Provide Y of position in cartesian ref frame
   */
  inline double GetY() const {return y_;}

  /**
   * Provide Z of position in cartesian ref frame
   */
  inline double GetZ() const {return z_;}

  /**
   * Provide R of position in spherical ref frame
   * If non-cartesian have not been calculated, then calculate them first
   */
  inline double GetR() const {
    if (!isCalculated_) CalcSphCylFromCar();
    return r_;
  }

  /**
   * Provide Theta of position in spherical ref frame
   * If non-cartesian have not been calculated, then calculate them first
   */
  inline double GetTheta() const {
    if (!isCalculated_) CalcSphCylFromCar();
    return theta_;
  }

  /**
   * Provide Phi of position in spherical or cylindrical ref frame
   * If non-cartesian have not been calculated, then calculate them first
   */
  inline double GetPhi() const {
    if (!isCalculated_) CalcSphCylFromCar();
    return phi_;
  }

  /**
   * Provide Rho of position in cylindrical ref frame
   * If non-cartesian have not been calculated, then calculate them first
   */
  inline double GetRho() const {
    if (!isCalculated_) CalcSphCylFromCar();
    return rho_;
  }

  //--------------

  /**
   * Set X position while keeping Y,Z constant.  Recalculate SPH and CYL.
   */
  inline void SetX(double x) {
    x_=x;
    isCalculated_=false; // when accessing CYL/SPH, they will be recalculated
  }

  /**
   * Set Y position while keeping X,Z constant.  Recalculate SPH and CYL.
   */
  inline void SetY(double y) {
    y_=y;
    isCalculated_=false; // when accessing CYL/SPH, they will be recalculated
  }

  /**
   * Set Z position while keeping X,Y constant.  Recalculate SPH and CYL.
   */
  inline void SetZ(double z) {
    z_=z;
    isCalculated_=false; // when accessing CYL/SPH, they will be recalculated
  }

  //--------------

  /**
   * Rotate position around X axis by angle
   */
  void RotateX(double angle);

  /**
   * Rotate position around Y axis by angle
   */
  void RotateY(double angle);

  /**
   * Rotate position around Z axis by angle
   */
  void RotateZ(double angle);

  /**
   * Computes the distance from this position to the origin of the
   * coordinate system (it's magnitude as a vector)
   */
  double Magnitude() const{
    if(isCalculated_)
      return r_;
    //otherwise use self dot-product
    return sqrt(*this * *this);
  }

  /**
   * Computes the square of the vector magnitude of the position
   */
  double Mag2() const{
    if(isCalculated_)
      return r_*r_;
    //otherwise use self dot-product
    return *this * *this;
  }

  /**
   * Vector inversion (makes the vector point in the opposite direction)
   */
  Position operator-() const{
    Position p;
    p.x_=-x_;
    p.y_=-y_;
    p.z_=-z_;
    p.isCalculated_=false;
    return p;
  }

  /**
   * Vector addition
   */
  Position& operator+=(const Position& rhs){
    x_+=rhs.x_;
    y_+=rhs.y_;
    z_+=rhs.z_;
    isCalculated_=false;
    return *this;
  }

  /**
   * Vector subtraction
   */
  Position& operator-=(const Position& rhs){
    x_-=rhs.x_;
    y_-=rhs.y_;
    z_-=rhs.z_;
    isCalculated_=false;
    return *this;
  }

  /**
   * Vector addition
   */
  Position operator+(const Position& rhs) const{
    return Position(*this)+=rhs;
  }

  /**
   * Vector subtraction
   */
  Position operator-(const Position& rhs) const{
    return Position(*this)-=rhs;
  }

  /**
   * Scalar (dot) product
   */
  double operator*(const Position& rhs) const{
    return x_*rhs.x_ + y_*rhs.y_ + z_*rhs.z_;
  }

  /**
   * Scalar (dot) product
   */
  double operator*(const Direction&) const;

  /**
   * Multiplication by a scalar
   */
  Position& operator*=(double a){
    x_*=a;
    y_*=a;
    z_*=a;
    isCalculated_=false;
    return *this;
  }

  /**
   * Divison by a scalar
   */
  Position& operator/=(double a){
    x_/=a;
    y_/=a;
    z_/=a;
    isCalculated_=false;
    return *this;
  }

  /**
   * Multiplication by a scalar
   */
  Position operator*(double a) const{
    return Position(*this)*=a;
  }

  /**
   * Division by a scalar
   */
  Position operator/(double a) const{
    return Position(*this)/=a;
  }

  /**
   * Vector (cross) product
   */
  Position Cross(const Position& d) const{
    return Position (y_*d.z_ - z_*d.y_,
                       z_*d.x_ - x_*d.z_,
                       x_*d.y_ - y_*d.x_);
  }

  /**
   * Vector (cross) product
   */
  Position Cross(const Direction&) const;

 protected:
  /**
   * cartesian (car)
   */
  double x_;
  double y_;
  double z_;

  /**
   * spherical (sph)
   */
  mutable double r_;
  mutable double theta_;
  mutable double phi_;

  /**
   * cylindrical (cyl) - Z and Phi are same.
   */
  mutable double rho_;

  /**
   * Whether the coordinates in secondary coordinates systems
   * (sph and cyl) are already computed
   */
  mutable bool isCalculated_;

 private:

  void CalcSphCylFromCar() const;
  void CalcCarCylFromSph();
  void CalcCarSphFromCyl();

};

class I3Position: public Position {

  public:
    explicit I3Position() : Position() {};
    explicit I3Position(double x, double y, double z): Position(x, y, z) {};
};

inline bool operator==(const Position& lhs, const Position& rhs) {
  return ((lhs.GetX() == rhs.GetX()) &&
          (lhs.GetY() == rhs.GetY()) &&
          (lhs.GetZ() == rhs.GetZ()));
}

inline bool operator!=(const Position& lhs, const Position& rhs) {
  return ((lhs.GetX() != rhs.GetX()) ||
          (lhs.GetY() != rhs.GetY()) ||
          (lhs.GetZ() != rhs.GetZ()));
}

Position operator*(double, const Position&);

std::ostream& operator<<(std::ostream& oss, const Position& p);

double abs(const Position& p);





#endif //I3POSITION_H_INCLUDED
