#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/vec3.h>

#include <TMinuit.h>

#include <cmath>
#include <iostream>
#include <iomanip>

typedef Healpix_Map<double> HMap;

#define TDOUBLE PLANCK_FLOAT64

using std::cout;
using std::cerr;
using std::endl;
using std::setw;
using std::sqrt;

const double degree = 4*atan(1.) / 180.;

static HMap map;
static HMap mapVar;

static void
chi2(int& npar, double* gin, double& f, double* par, int iflag)
{
  double m  = par[0]; // Monopole
  double p1 = par[1]; // Dipole
  double p2 = par[2];
  double p3 = par[3];
  double Q1 = par[4]; // Quadrupole
  double Q2 = par[5];
  double Q3 = par[6];
  double Q4 = par[7];
  double Q5 = par[8];

  double sum = 0;
  double df;
  double sigma2;
  vec3 v;

  const double maxZ = sin(70*degree);
  const double minZ = sin(-90*degree);

  for (int i = 0; i < map.Npix(); ++i) {
    v = map.pix2vec(i);
    // Mask out the northern hemisphere
    if (v.z >= minZ && v.z <= maxZ  && map[i] > -1.e30) {
      df = map[i] - m -(p1*v.x + p2*v.y + p3*v.z + 
                         Q1 * 0.5*(3*v.z*v.z-1.) +
                         Q2 * 2*v.x*v.z +
                         Q3 * 2*v.y*v.z +
                         Q4 * (v.x*v.x - v.y*v.y) +
                         Q5 * 2*v.x*v.y);

      sigma2 = mapVar[i];
      sum += df*df / sigma2;
      //sum += fabs(df);
    }
  }
  f = sum;
}

int main(int argc, char* argv[])
{
  if (argc != 3) {
    cerr << "\nUsage: " << argv[0] << " [bkg.FITS] [dat.FITS]\n\n";
    return 1;
  }

  HMap bkg;
  HMap dat;

  // Read in the relative intensity map
  read_Healpix_map_from_fits(argv[1], bkg);
  read_Healpix_map_from_fits(argv[2], dat);

  map.SetNside(bkg.Nside(), bkg.Scheme());
  mapVar.SetNside(bkg.Nside(), bkg.Scheme());

  HMap mapRes;
  mapRes.SetNside(bkg.Nside(), bkg.Scheme());
  
  HMap dipmap;
  dipmap.SetNside(bkg.Nside(), bkg.Scheme());
  
  HMap dqmap;
  dqmap.SetNside(bkg.Nside(), bkg.Scheme());

  int npix = 0;
  double Nb;
  double Nd;
  double minZ = sin(-90*degree);
  double maxZ = sin(70*degree);
  //const double alpha = 1./20.;
  const double alpha = 1.;
  vec3 v;

  for (int i = 0; i < bkg.Npix(); ++i) {
    v = map.pix2vec(i);
    if ( v.z >= minZ && v.z <= maxZ) {

      Nb = bkg[i];
      map[i] = dat[i];
      Nd = Nb*(1.-map[i]);
      npix++;
      if (Nb > 0){
         //mapVar[i] = Nd*(Nb + alpha*Nd) / (Nb*Nb*Nb);
         mapVar[i] = 1.0/Nb;
      }
      if (map[i] != map[i])
        map[i] = 0.;
      if (mapVar[i] != mapVar[i])
        mapVar[i] = 1e99;
    }
    else {
      map[i] = 0;
      mapVar[i] = 0;
    }
  }

  // Initialize the MINUIT fit
  TMinuit minuit(9);
  minuit.SetFCN(chi2);

  double argList[10];
  argList[0] = 1;
  int iErrFlag = 0;

  // Set parameter start values, mimizer step size, and min/max limits
  double vstart[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  double step[9] = { 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6 };
  double bmin[9] = { -1e6, -1e6, -1e6, -1e6, -1e6, -1e6, -1e6, -1e6, -1e6 };
  double bmax[9] = {  1e6,  1e6,  1e6, 1e6,  1e6,  1e6,  1e6,  1e6,  1e6 };
  const char* name[9] = { "$m$", "$p_1$", "$p_2$", "$p_3$", "$Q_1$", "$Q_2$", "$Q_3$", "$Q_4$", "$Q_5$" };
  for (int i = 0; i < 9; ++i)
    minuit.mnparm(i, name[i], vstart[i], step[i], bmin[i], bmax[i], iErrFlag);

  // Fix a parameter for the fit
  // argList[0] = 0; // Change to parameter number;
  // minuit.mnexcm("FIX", argList, 1, iErrFlag);
  // argList[0] = 8;
  // minuit.mnexcm("FIX", argList, 1, iErrFlag);
  // argList[0] = 4;
  // minuit.mnexcm("FIX", argList, 1, iErrFlag);

  // Set minimization strategy (1=standard; 2=improved, but slow)
  argList[0] = 1;
  minuit.mnexcm("SET STR", argList, 2, iErrFlag);

  // Set up MIGRAD to iterate up to 1000 times
  argList[0] = 1000;
  minuit.mnexcm("MIGRAD", argList, 1, iErrFlag);

  // Scan about the global minimum
  // argList[0] = 0;
  // minuit.mnexcm("SCAN", argList, 1, iErrFlag);

  double x, dx;
  const double scale = 1e4;
  cout << "\nFit Results (x" << 1./scale << "):\n"; 
  cout.precision(3);
  for (int i = 0; i < 9; ++i) {
    minuit.GetParameter(i, x, dx);
    cout << std::fixed
         << setw(2) << name[i] << " & $"
         << setw(6) << x * scale << " \\pm "
         << setw(5) << dx * scale << "$"
         << endl;
  }

  double m, p1, p2, p3, Q1, Q2, Q3, Q4, Q5;
  double dp1, dp2, dp3; 
  minuit.GetParameter(0, m, dx);
  minuit.GetParameter(1, p1, dp1);
  minuit.GetParameter(2, p2, dp2);
  minuit.GetParameter(3, p3, dp3);
  minuit.GetParameter(4, Q1, dx);
  minuit.GetParameter(5, Q2, dx);
  minuit.GetParameter(6, Q3, dx);
  minuit.GetParameter(7, Q4, dx);
  minuit.GetParameter(8, Q5, dx);
  
  minZ = sin(-90*degree);
  maxZ = sin(70*degree);

  double chi2sum = 0.;
  for (int i = 0; i < mapRes.Npix(); ++i) {
    v = mapRes.pix2vec(i);
    // Mask out the northern hemisphere
    mapRes[i] = 0.;
    if (v.z >= minZ && v.z <= maxZ  && map[i] > -1.e30) {
      double dipole = p1*v.x + p2*v.y + p3*v.z;
      double quad = Q1 * 0.5*(3*v.z*v.z-1.) +
                    Q2 * 2*v.x*v.z +
                    Q3 * 2*v.y*v.z +
                    Q4 * (v.x*v.x - v.y*v.y) +
                    Q5 * 2*v.x*v.y;

      mapRes[i] = map[i] - (m + dipole + quad);
      dipmap[i] = m + dipole;
      dqmap[i] = m + dipole + quad;

      map[i];
      mapRes[i];

///

      double df = map[i] - m -(p1*v.x + p2*v.y + p3*v.z + 
                         Q1 * 0.5*(3*v.z*v.z-1.) +
                         Q2 * 2*v.x*v.z +
                         Q3 * 2*v.y*v.z +
                         Q4 * (v.x*v.x - v.y*v.y) +
                         Q5 * 2*v.x*v.y);

      double sigma2 = mapVar[i];
      chi2sum += df*df / sigma2;
      //sum += fabs(df);
//
    }
  }

  fitshandle out1;
  out1.create("!relInt.fits");
  write_Healpix_map_to_fits(out1, map, TDOUBLE);

  fitshandle out2;
  out2.create("!residuals.fits");
  write_Healpix_map_to_fits(out2, mapRes, TDOUBLE);
  
  fitshandle out3;
  out3.create("!dipolefit.fits");
  write_Healpix_map_to_fits(out3, dipmap, TDOUBLE);
  
  fitshandle out4;
  out4.create("!dqfit.fits");
  write_Healpix_map_to_fits(out4, dqmap, TDOUBLE);
  cout << "npix " << npix << endl;
  cout << "*********************"  << endl;
  cout << "npix (all) " << map.Npix() << endl;
  cout << "Chi2/ndf " << chi2sum << "/" << npix-9 << endl;

  double p = sqrt(p1*p1 + p2*p2);
  double dp = sqrt((p1*p1*dp1*dp1 + p2*p2*dp2*dp2 )/(p*p));

  cout << "dipole: " << p * 1e4 << " +/- " << dp * 1e4 << endl;
  cout << "alpha: " << atan(p2/p1)/degree << endl;


  return 0;
}

