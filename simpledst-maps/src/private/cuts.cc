#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "cuts.h"


bool ICenergyCut(SimpleDST dst, photospline::splinetable<> &spline, double zenith, double emin, double emax) {

  // Setup basic parameters
  double x = cos(zenith);
  double y = log10(dst.NChannels);

  // Boundary check (energy cut tables go to 0.3 in cos(zenith))
  if (x < 0.3)
    return false;

  // Catch additional outliers
  double coords[2] = {x, y};
  int centers[spline.get_ndim()];
  if (!spline.searchcenters(coords, centers)) {
    std::cout << "Variables outside of table boundaries" << std::endl;
    std::cout << "x: " << x << " y: " << y << std::endl;
    return false;
  }

  // Calculate reconstructed energy
  double median = spline.ndsplineeval(coords, centers, 0);
  // Make sure we're in the energy bin range
  if (emax <= 0.0)
     return (median >= emin); 
  return ((median >= emin) && (median < emax));
}


int ITenergyCut(SimpleDST dst, double emin, double emax) {

  // Get most likely energy value
  double llhEnergy = (dst.pLLH >= dst.fLLH) ? dst.pEnergy : dst.fEnergy;
  double logEnergy = log10(llhEnergy);

  if (emax <= 0.0)
    return (logEnergy >= emin);
  // Get energy bin
  return ((logEnergy >= emin) && (logEnergy < emax));
}

int ITs125Cut(SimpleDST dst, double smin, double smax) {

  // Get desired s125 value
  double s125 = (dst.nStations >= 5) ? dst.s125 : dst.ss125;
  double logS125 = log10(s125);

  // Make sure we're in the bin range
  if (smax <= 0.0)
    return (logS125 >= smin);
  return ((logS125 >= smin) && (logS125 < smax));

}

int ITstats(Simble DST dst, double loen, double hien) {
  //low energy is 3<numStat<8
  //high energy is 8<numStat
  //how to separate these into bins???
  double stations = (dst.numStat >= 8) ? dst.s124 : dst.ss125;
  double logstations = log10(stations);
  
  // correct?



