#include <healpix_cxx/alm.h>
#include <healpix_cxx/xcomplex.h>
#include <healpix_cxx/alm_healpix_tools.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#include <math.h>
#include <numeric>

#include <TRandom.h>
#include <TRandom1.h>
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

#include <iter-lhreco-proj/pickle.h>
#include "TFile.h"
#include "TTree.h"
#include <TChain.h>

#include "cuts.h"

#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <solardipole.h>
#include <SimpleTrigger.h>



#ifdef HAWCNEST
#include <hawcnest/Logging.h>
#include <hawcnest/CommandLineConfigurator.h>
#else
namespace po = boost::program_options; 
#define log_info(args...) cout << args << endl
#endif


namespace fs = boost::filesystem; 
using namespace std;
using namespace boost::numeric::ublas;
using boost::format;

typedef Healpix_Map<double> SkyMap; 
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer

bool newConfig(string config);

double 
Sum(const SkyMap& map,unsigned int npix) 
{
    double sumval = 0.;
    for (unsigned int i=0; i < npix; i++) 
    {
         sumval += map[i];
    }
    return sumval;
}

bool filterCut(po::variables_map vm, SimpleDST dst) {

  string filter = vm["filter"].as<string>();
  string config = vm["config"].as<string>();
  if (filter=="STA3" && dst.isSTA3)
    return true;
  if (config=="IT59" || config=="IT73" || config=="IT81") {
    if ((filter=="STA8" && dst.isSTA8) ||
        (filter=="NotSTA8" && dst.isSTA3 && !dst.isSTA8))
      return true;
  }
  if (config=="IT81-2012" || config=="IT81-2013") {
    if ((filter=="STA8" && dst.nStations>=8) ||
        (filter=="NotSTA8" && dst.isSTA3 && dst.nStations<8))
      return true;
  }
  return false;
}



int main(int argc, char* argv[])
{
//*****************************************************************************
////// Input ////////////////////////////////////////////////////////////////// 
//*****************************************************************************
    std::string foldername;
    std::string method("sid");
    std::string config;
    std::string filter;
    unsigned int nTimesteps;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideOut;
    int nmaps;
    bool sundp;
    bool show_progress;
    double lon;
    double lat;
    double elogmin, elogmax, rloglmax,ldirmin,ndirmin;
    double slogmin, slogmax;
    std::vector< std::string > input;
    std::string coords = "Azimuth/ Zenith";
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
             ("spline", po::value<string>(), "File containing spline tables")
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum energy")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum energy (0: no cut)")
             ("slogmin", po::value<double>(&slogmin)->default_value(0), "Minimum logS125")
             ("slogmax", po::value<double>(&slogmax)->default_value(0), "Maximum logS125 (0: no cut)")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map") 
             ("nmaps", po::value<int>(&nmaps)->default_value(10), "Number of random maps togenerate") 
             ("sundp", po::bool_switch(&sundp)->default_value(false), "subtract solar dipole")
             ("config", po::value<string>(&config), "Detector configuration") 
             ("method", po::value<string>(&method), "Sidereal, Anti, Solar, Extended")
             ("filter", po::value<string>(&filter), "Filter for IceTop data")
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
    const char* triggerTree;
    string detector = config.substr(0,2); 
    if (detector == "IC") { 
        masterTree = "CutDST"; 
        triggerTree = "TDSTTriggers"; 
    } 
    if (detector == "IT") { 
        masterTree = "master_tree"; 
        triggerTree = "";   // Unused? Will probably break IT functionality...  
    }


    // TFile f("Event.root");
    TChain cutDST(masterTree);
    for (unsigned i = 0; i < input.size(); ++i) { 
        cutDST.Add(input[i].c_str()); 
    } 
    SimpleDST dst(&cutDST, config);

    TChain trigDST(triggerTree);
    if (newConfig(config)) { 
        for (unsigned i = 0; i < input.size(); ++i) { 
            trigDST.Add(input[i].c_str()); 
        }
    } 
    SimpleTrigger dst_trig(&trigDST);

    cout << "Number of chained files: " << cutDST.GetNtrees() << endl; 
    Long64_t nEntries = cutDST.GetEntries();

    double rndMJD;


    // Read in spline tables if provided
    photospline::splinetable<> spline;
    if (vm.count("spline")) {
        string splineFile = vm["spline"].as<string>();
        spline.read_fits(splineFile.c_str());
    }


    for (unsigned int rand_iter=0; rand_iter<nmaps;rand_iter++)
    {
            TRandom1 rand(rand_iter+1);


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


            // output directory
            stringstream subfoldername;
            subfoldername << foldername << boost::format("/%02d") % rand_iter ;
            fs::path subdir(subfoldername.str()); 
            if(!(fs::exists(subdir)) ) { 
                    log_info("Directory " << subfoldername.str() << " doesn't exist");
                    if (fs::create_directory(subdir)) 
                            log_info("....successfully created !");
            }

            double eventweight;
            double weightsum = 0.;
            float zenith, azimuth;
            int timeBin;
            int prevTimeBin = 0 ;
            int firstTimeBin = 0 ;
            int prevMJD = 0 ;
            int startMJD = 0 ;
            int currMJD = 0 ;
            bool init=false;
            SkyMapPtr tmp_map;

            for (int i = 0; i < nEntries; ++i) {
               cutDST.GetEntry(i); 
               if (newConfig(config)) { 
                   trigDST.GetEntry(i); 
               }

               if ( (rloglmax > 0 ) && (dst.RLogL > rloglmax) )
                  continue; 

               double mjd = dst.ModJulDay; 
               double lst = GetGMST(mjd); 
               if (method == "anti") { 
                   lst = GetGMAST(mjd); 
               } else if (method == "ext") { 
                   lst = GetGMEST(mjd); 
               } else if (method == "solar") { 
                   lst= ( mjd - int(mjd) )* 24.; 
               }

               // https://en.wikipedia.org/wiki/Time_in_Antarctica
               azimuth = dst.LLHAzimuthDeg/180.*M_PI; 
               //phi = RADeg/180.*M_PI; 
               zenith = dst.LLHZenithDeg/180.*M_PI; 
               //theta = DecDeg/180.*M_PI; 


               if ( dst.NDirHits < ndirmin*cos(zenith) ) 
                  continue;

               if ( dst.LDir < ldirmin*cos(zenith) ) 
                  continue;

               double rngMST = rand.Uniform(24.0);
               
               timeBin = int( nTimeBins*rngMST/24.0 ); 

               // Reconstruction cuts
               if (!dst.isReco || zenith != zenith || azimuth != azimuth) 
                   continue;

               // IceTop filter cut
               if (vm.count("filter") && !filterCut(vm, dst)) 
                   continue;
               
               if (vm.count("spline") && !ICenergyCut(dst, spline, zenith, elogmin, elogmax))
                  continue;

               // Energy cuts for IceTop and IceCube
               if (detector == "IT" && !ITenergyCut(dst, elogmin,elogmax)) 
                   continue;

               pointing pt(zenith, azimuth);

               eventweight = 1.0;
               if (sundp) { 
                   Direction dir(zenith,azimuth); 
                   Equatorial eq = GetEquatorialFromDirection(dir, mjd); 
                   eventweight = solar_dipole(mjd, eq.ra, eq.dec,true); 
               } 
               weightsum += eventweight;

               tmp_map = Nmap[timeBin % 360];// make sure tsid < 360 deg
               int j = tmp_map->ang2pix(pt);
               (*tmp_map)[j] += 1.0; 
               nEvents++;

            }
            std::cout << std::endl;


            for (unsigned int degbin=0;degbin< 360;degbin++) 
            { 
                 fitshandle fitsOut;
                 stringstream namefits;
                 if (elogmin > 0 || elogmax > 0)
                    namefits << foldername << boost::format("/%02d/CR_ICECUBE_LOCAL_%4.2f-%4.2fGeV_NSIDE%d_degbin-%02d.fits.gz") % rand_iter % elogmin % elogmax % nsideOut % degbin;
                 else
                    namefits << foldername << boost::format("/%02d/CR_ICECUBE_LOCAL_NSIDE%d_degbin-%02d.fits.gz") % rand_iter %  nsideOut % degbin;
                 if (fs::exists(namefits.str()) ) { 
                         fs::remove(namefits.str() ); 
                 }
                 std::cout << namefits.str() << std::endl;
                 fitsOut.create(namefits.str().c_str());
                 fitsOut.add_comment("Local Map");

                 SkyMap jmap = *(Nmap[degbin]);
                 write_Healpix_map_to_fits(fitsOut, jmap, PLANCK_FLOAT64);
                 fitsOut.close();
            }
            std::cout << "passing rate: " << nEvents << "/" << nEntries << std::endl;
            std::cout << "weightsum: " << weightsum << std::endl;
            
    }
}

bool newConfig(string config) {

  if (config=="IC86-2011" || config=="IC86-2012" || config=="IC86-2013" || 
      config=="IC86-2014" || config=="IC86-2015") {
    return false;
  }
  return true;
}


