#include <healpix_cxx/alm.h>
#include <healpix_cxx/xcomplex.h>
#include <healpix_cxx/alm_healpix_tools.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <math.h>
#include <numeric>

#include <boost/shared_ptr.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/program_options.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <memory>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <iterator>
#include <vector>

#include <pickle.h>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <TChain.h>
#include <cuts.h>
#include <config.h>
#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <solardipole.h>

namespace po = boost::program_options;
#define log_info(args...) cout << args << endl

namespace fs = boost::filesystem;
using namespace std;
using namespace boost::numeric::ublas;

using boost::format;

typedef Healpix_Map<double> SkyMap;
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer


double
Sum(const SkyMap& map, unsigned int npix)
{
    double sumval = 0.;
    for (unsigned int i=0; i < npix; i++)
    {
         sumval += map[i];
    }
    return sumval;
}


void progress_bar(float progress, std::string msg)
{
    int barWidth = 70;
    std::cout << msg << " [";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
            if (i < pos) std::cout << "=";
            else if (i == pos) std::cout << ">";
            else std::cout << " ";
    }
    std::cout.precision(2);
    std::cout << "] " << progress * 100.0 << " %\r";
    std::cout.flush();
}


int main(int argc, char* argv[])
{

//*****************************************************************************
////// Input //////////////////////////////////////////////////////////////////
//*****************************************************************************
    gROOT->SetBatch( 1 );
    std::string foldername;
    std::string outfile;
    std::string method("sid");
    std::string dirreco;
    std::string config;
    int nsideOut;
    bool sundp;
    bool show_progress;
    double elogmin, elogmax, rloglmax, ldirmin, ndirmin;
    double slogmin, slogmax, smin, smax;
    std::vector< std::string > input;
    int nEvents = 0;

    po::options_description desc("Options");
    po::positional_options_description p;
    p.add("input", -1);
    po::variables_map vm;

    try {
       // Define and parse the program options
       desc.add_options()
             ("help", "produce help message")
             ("input", po::value<std::vector <std::string> >(&input)->multitoken(), "Input files")
             ("outdir", po::value<std::string>(&foldername)->default_value("./sample/"), "Directory of output")
             ("outfile", po::value<std::string>(&outfile)->default_value("CR_ITpass2"), "Base name for output healpix maps")
             ("spline", po::value<string>(), "File containing spline tables")
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum log energy (GeV)")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum log energy (GeV) (0: no cut)")
             ("slogmin", po::value<double>(&slogmin)->default_value(0), "Minimum logS125")
             ("slogmax", po::value<double>(&slogmax)->default_value(0), "Maximum logS125 (0: no cut)")
             ("smin", po::value<double>(&smin)->default_value(0), "Minimum nStations")
             ("smax", po::value<double>(&smax)->default_value(0), "Maximum nStations")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map")
             ("sundp", po::bool_switch(&sundp)->default_value(false), "subtract solar dipole")
             ("config", po::value<string>(&config), "Detector configuration")
             ("method", po::value<string>(&method), "Sidereal, Anti, Solar, Extended")
             ("dirreco", po::value<string>(&dirreco), "Directional reconstruction (ShowerPlane, Laputop, LaputopSmall) (IceTop-only)")
             ("progress", po::bool_switch(&show_progress)->default_value(false), "show progress bar")
             ;

        //po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);

        /// --help option
        if ( vm.count("help")  ) {
            std::cout << "Basic Command Line Parameter App"
                      << std::endl << desc << std::endl;
            return 0;
        }
        po::notify(vm); // throws on error, so do after help in case there are any problems

        if (vm.count("input"))
        {
            cout << "Input files: " << "\n";
        }

    } catch(po::error& e) {
            std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
            std::cerr << desc << std::endl;
            return 1;

    } catch(std::exception& e) {
            std::cerr << "Unhandled Exception reached the top of main: "
                      << e.what() << ", application will now exit" << std::endl;
            return 2;
    }

//*****************************************************************************
////// Initialize /////////////////////////////////////////////////////////////
//*****************************************************************************

    const char* masterTree;
    bool sundp2 = true;
    bool usesplines = false;

    Config cfg(vm);
    if (vm.count("spline") && (cfg.detector == Config::IceCube) )
        usesplines = true;

    if (cfg.detector == Config::IceCube) {
        masterTree = "CutDST";
    }
    if (cfg.detector == Config::IceTop) {
        sundp2 = false;
        if (cfg.cfg == Config::ITv3) {
            masterTree = "MasterTree";
        } else {
            masterTree = "master_tree";
        }
    }


    // TFile f("Event.root");
    TChain cutDST(masterTree);
    for (unsigned i = 0; i < input.size(); ++i) {
        cutDST.Add(input[i].c_str());
    }
    SimpleDST dst(&cutDST, config);


    cout << "Number of chained files: " << cutDST.GetNtrees() << endl;
    Long64_t nEntries = cutDST.GetEntries();


    // Read in spline tables if provided
    photospline::splinetable<> spline;
    if (usesplines) {
        string splineFile = vm["spline"].as<string>();
        spline.read_fits(splineFile.c_str());
    }


    unsigned int npix = 12*nsideOut*nsideOut;
    unsigned int nTimeBins = 360;
    // Import data : n_tau_i
    std::vector<SkyMapPtr> Nmap;

    for (unsigned int degbin=0;degbin< 360;degbin++)
    {
        SkyMapPtr locMapDegr(new SkyMap);
        locMapDegr->SetNside(nsideOut, RING);
        locMapDegr->fill(0.);
        Nmap.push_back(locMapDegr);
    }
    log_info("Loaded " << input[0] << ", etc.");

    // output directory
    fs::path dir(foldername);
    if(!(fs::exists(dir)) ) {
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir))
                    log_info("....successfully created !");
    }

    double eventweight;
    double weightsum = 0.;
    float zenith, azimuth;
    double mjd, lst;
    bool fitPassed;
    int timeBin;
    SkyMapPtr tmp_map;

    cout << "Reading " << nEntries << " entries...\n";
    if (usesplines) {
        cout << "spline " << vm.count("spline") << " elogmin"<< elogmin << " elogmax" << elogmax << std::endl;
    }

    for (Long64_t i = 0; i < nEntries; ++i) {

      cutDST.GetEntry(i);

      if ( show_progress && (i % 5000 == 0) )
         progress_bar(i/float(nEntries), " Reading entries");

      if ( (rloglmax > 0 ) && (dst.RLogL > rloglmax) )
         continue;


      mjd = dst.ModJulDay;
      lst = GetGMST(mjd);
      if (cfg.method == Config::antisid) {
          lst = GetGMAST(mjd);
      } else if (cfg.method == Config::extsid) {
          lst = GetGMEST(mjd);
      } else if (cfg.method == Config::solar) {
          lst= ( mjd - int(mjd) )* 24.;
      }


      // Establish directional reconstruction (IceCube)
      if (cfg.detector == Config::IceCube) {
          azimuth = dst.LLHAzimuthDeg/180.*M_PI;
          zenith = dst.LLHZenithDeg/180.*M_PI;
          fitPassed = dst.isReco;
      }

      // Establish directional reconstruction (IceTop)
      else if (cfg.detector == Config::IceTop) {
          if (cfg.reco == Config::Laputop) {
              azimuth = dst.LapAzimuth;
              zenith = dst.LapZenith;
              fitPassed = (dst.LapFitStatus == 0);
          } else if (cfg.reco == Config::LaputopSmall) {
              azimuth = dst.LSAzimuth;
              zenith = dst.LSZenith;
              fitPassed = (dst.LSFitStatus == 0);
          } else {
              azimuth = dst.SPAzimuth;
              zenith = dst.SPZenith;
              fitPassed = (dst.SPFitStatus == 0);
          }
      }

      // General directional reconstruction cuts
      if (!fitPassed || std::isnan(zenith) || std::isnan(azimuth))
          continue;

      // Additional IceCube cuts
      if ( dst.NDirHits < ndirmin*cos(zenith) )
         continue;
      if ( dst.LDir < ldirmin*cos(zenith) )
          continue;
      if (usesplines && !ICenergyCut(dst, spline, zenith, elogmin, elogmax))
          continue;

      // Additional IceTop cuts
      if (cfg.detector == Config::IceTop) {
          if (cfg.cfg == Config::ITv3) {
              if (!ITNstatCut(dst, smin,smax))
                  continue;
          } else if (!ITs125Cut(dst, slogmin,slogmax)) {
                  continue;
          }
      }

      // Check if angle above pi
      pointing pt(zenith, azimuth);
      if(pt.theta > M_PI)
          continue;

      // Optional solar dipole weighting
      eventweight = 1.0;
      if (sundp) {
          Direction dir(zenith,azimuth);
          Equatorial eq = GetEquatorialFromDirection(dir, mjd);

          eventweight = solar_dipole(mjd, eq.ra, eq.dec,sundp2);
      }
      weightsum += eventweight;

      // Increment the map (https://en.wikipedia.org/wiki/Time_in_Antarctica)
      timeBin = int( nTimeBins*lst/24.0 );
      tmp_map = Nmap[timeBin % 360];// make sure tsid < 360 deg
      int j = tmp_map->ang2pix(pt);
      (*tmp_map)[j] += eventweight;
      nEvents++;

    }

    std::cout << std::endl;

    string detector = "ICECUBE";
    if (cfg.detector == Config::IceTop)
        detector = "ICETOP";

    for (unsigned int degbin=0; degbin< 360; degbin++) {

        fitshandle fitsOut;
        stringstream namefits;

        namefits << foldername << boost::format("/%s_degbin-%03d.fits.gz") % outfile % degbin;

        if (fs::exists(namefits.str()) )
                fs::remove(namefits.str() );

        SkyMap jmap = *(Nmap[degbin]);
        std::cout << "bin: " << degbin << " events: " << Sum(jmap,npix);
        std::cout << std::endl;
        //std::cout << namefits.str() << std::endl;
        fitsOut.create(namefits.str().c_str());
        fitsOut.add_comment("Local Map");

        write_Healpix_map_to_fits(fitsOut, jmap, PLANCK_FLOAT64);
        fitsOut.close();
    }

    std::cout << "passing rate: " << nEvents << "/" << nEntries << std::endl;
    std::cout << "weightsum: " << weightsum << std::endl;

}


