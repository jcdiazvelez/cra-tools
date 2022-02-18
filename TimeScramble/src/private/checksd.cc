#include <SimpleDST.h>

#include <TChain.h>
#include <TH1D.h>
#include <TMath.h>
#include <TRandom.h>
#include <TStopwatch.h>

#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/pointing.h>

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

#include <boost/program_options.hpp>


#include <astro/astro.h>
#include <astro/time.h>

#include <solardipole.h>
#include <Direction.h>


using namespace std;
const double deg2rad = TMath::DegToRad();
const double rad2deg = TMath::RadToDeg();


int main(int argc, char* argv[]) {

  int NSide = 16;

  Healpix_Map<float> testmap;
  testmap.SetNside(NSide, RING);
  testmap.fill(0.);

  pointing sphereDir;
  int pixelID;


  // Get info from previous- and next-day files
  int yy = 2019;
  int mm = 5;
  int dd = 16;
  astro::Time t(yy, mm, dd, 0, 0, 0);
  astro::Time t1(yy, mm, dd, 0, 0, 0);

  Double_t mjd1 = t.GetMJD();

  // Integration time



  //*********************************************************************//
  // Begin iterating through events
  //*********************************************************************//
  int npix = testmap.Npix();
  double sdweight;

  for (double j=0; j<365; j++) {
    t.SetTime(mjd1+j);
    t1.SetTime(mjd1+j-.25);

    double maxsd=0;
    double minsd=2;
    double maxra=-1;
    double minra=-1;
    double maxdec=-1;
    double mindec=-1;
    double maxra1=-1;
    double maxdec1=-1;

    double dr,dd;

    for (unsigned ipix=0; ipix<npix; ++ipix) { 
            pointing p = testmap.pix2ang(ipix);
            Direction dir(p.theta, p.phi);
            Equatorial eq = GetEquatorialFromDirection(dir, t.GetMJD());

            sdweight = solar_dipole(mjd1+j, eq.ra, eq.dec);

            Equatorial eq1 = GetEquatorialFromDirection(dir, t1.GetMJD());
            if (sdweight > maxsd) 
            {
                    maxsd = sdweight;
                    maxra = eq.ra;
                    maxdec = eq.dec;
            }
            if (sdweight < minsd) 
            {
                    minsd = sdweight;
                    minra = eq.ra;
                    mindec = eq.dec;
                    maxra1 = eq1.ra;
                    maxdec1 = eq1.dec;
            }
            testmap[ipix] += sdweight;
    } 
    double sundr,sundd,diam;      // J2000.0 mean RA,Dec (radians)
    palRdplan(mjd1+j, 0, 0, 0, &sundr, &sundd, &diam);
    double theta = acos(sin(sundd)*sin(maxdec)+cos(sundd)*cos(maxdec)*cos(sundr-maxra));
    printf("%2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f, %2.4f\n", 
                            mjd1+j,maxra*rad2deg,maxdec*rad2deg,
                            minra*rad2deg,mindec*rad2deg, sundr*rad2deg, sundd*rad2deg, theta*rad2deg,
                            maxra1*rad2deg,maxdec1*rad2deg );


    
  }
  
    fitshandle fitsOut;
    arr<std::string> colname(1);
    fitsOut.create("testmap.fits");
    fitsOut.add_comment("Maps: sd");
    colname[0] = "test map";
    prepare_Healpix_fitsmap(fitsOut, testmap, PLANCK_FLOAT64, colname);

    fitsOut.write_column(1, testmap.Map());
    fitsOut.close();

  return 0;
}


