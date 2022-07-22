#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <units.h>
#include <photospline/splinetable.h>
#include <photospline/bspline.h>


int main(int argc, char* argv[]) 
{
  photospline::splinetable<> spline; 
  std::string splineFile("/data/user/fmcnally/anisotropy/sim/IC86_20904_hist_spline.fits");
  spline.read_fits(splineFile.c_str()); 

  double emin = 0;
  double emax = 4.5;
  double radians = 180/constants::pi;

  for (unsigned NChannels = 1; NChannels < 5000; NChannels *= 2)
  {
    for (int zenith_bin = 0; zenith_bin < 8; zenith_bin++)
    { 
        // Setup basic parameters 
        double zenith = zenith_bin*0.5*constants::pi/10;
        double x = cos(zenith); 
        double y = log10(NChannels);

        // Boundary check (energy cut tables go to 0.3 in cos(zenith))
        if (x < 0.3) 
        {
            std::cerr << "Variables outside of table boundaries" << std::endl; 
            std::cerr << "x: " << x << " y: " << y << std::endl; 
            std::cerr << "zenith: " << zenith*radians << "deg, NChan: " <<  NChannels << std::endl; 
            continue;
        }

        // Catch additional outliers
        double coords[2] = {x, y};
        int centers[spline.get_ndim()];
        if (!spline.searchcenters(coords, centers)) { 
            std::cout << "Variables outside of table boundaries" << std::endl; 
            std::cout << "x: " << x << " y: " << y << std::endl; 
            //return false; 
        }

        // Calculate reconstructed energy
        double median = spline.ndsplineeval(coords, centers, 0);
        std::cout << "NChan:" << NChannels << ", zenith: " << zenith << ", emedian: " << median << std::endl;
        // Make sure we're in the energy bin range
        //if (emax <= 0) 
        //    return (median >= emin); 
        //return ((median >= emin) && (median < emax));
    }
  }
}


