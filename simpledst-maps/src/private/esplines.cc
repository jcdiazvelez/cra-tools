#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>


bool ICenergyCut(unsigned NChannels, photospline::splinetable<> &spline, double zenith, double emin, double emax) {

  // Setup basic parameters
  double x = cos(zenith);
  double y = log10(NChannels);

  // Boundary check (energy cut tables go to 0.3 in cos(zenith))
  if (x < 0.3)
    return -1;

  // Catch additional outliers
  double coords[2] = {x, y};
  int centers[spline.get_ndim()];
  if (spline.searchcenters(coords, centers) != 0) {
    std::cout << "Variables outside of table boundaries" << std::endl;
    std::cout << "x: " << x << " y: " << y << std::endl;
    return false;
  }

  // Calculate reconstructed energy
  double median = spline.ndsplineeval(coords, centers, 0);
  // Make sure we're in the energy bin range
  if (emax <= 0)
     return (median >= emin); 
  return ((median >= emin) && (median < emax));
}


