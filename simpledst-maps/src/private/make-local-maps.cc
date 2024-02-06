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

#include <iter-lhreco-proj/pickle.h>
#include <iter-lhreco-proj/illh-utils.h>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <TChain.h>

#include "cuts.h"
#include "config.h"

#include <astro/astro.h>
#include <astro/time.h>
#include <Direction.h>
#include <solardipole.h>
#include <SimpleTrigger.h>

namespace po = boost::program_options; 
#define log_info(args...) cout << args << endl

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

SkyMapPtr load_siderial_map(std::string path)
{
    SkyMapPtr sidMap(new SkyMap);

    fitshandle handle;
    handle.open(path.c_str());
    handle.goto_hdu(2);
    read_Healpix_map_from_fits(handle, *sidMap, 1);
    handle.close();
    
    return sidMap;
}


int main(int argc, char* argv[])
{
//*****************************************************************************
////// Input ////////////////////////////////////////////////////////////////// 
//*****************************************************************************
    gROOT->SetBatch( 1 );
    std::string foldername;
    std::string smethod("sid");
    std::string config;
    std::string filter;
    config::method method_;
    unsigned int nTimesteps;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideOut;
    bool sundp;
    bool show_progress;
    bool costheta;
    bool rmsid(false);
    double detlon;
    double detlat;
    double elogmin, elogmax, rloglmax,ldirmin,ndirmin;
    double slogmin, slogmax;
    std::vector< std::string > input;
    std::string coords = "Azimuth/ Zenith";
    std::string sidpath; 
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
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum log energy (GeV)")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum log energy (GeV) (0: no cut)")
             ("slogmin", po::value<double>(&slogmin)->default_value(0), "Minimum logS125")
             ("slogmax", po::value<double>(&slogmax)->default_value(0), "Maximum logS125 (0: no cut)")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("costheta", po::bool_switch(&costheta)->default_value(false), "cuts propto costheta")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map") 
             ("sundp", po::bool_switch(&sundp)->default_value(false), "subtract solar dipole")
             ("remove-sid", po::value<string>(&sidpath), "sidereal distribution to subract (for solar map)")
             ("lon", po::value<double>(&detlon)->default_value(270.0), "Longitude of detector") 
             ("lat", po::value<double>(&detlat)->default_value(-90.0), "Latitude of detector")
             ("config", po::value<string>(&config), "Detector configuration") 
             ("method", po::value<string>(&smethod), "sidereal, anti, solar, ext")
             ("filter", po::value<string>(&filter), "Filter for IceTop data")
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


    // Read in spline tables if provided
    photospline::splinetable<> spline;
    if (vm.count("spline")) {
        string splineFile = vm["spline"].as<string>();
        spline.read_fits(splineFile.c_str());
    }


    unsigned int npix = 12*nsideOut*nsideOut; 
    unsigned int nTimeBins = 360;
    // Import data : n_tau_i
    std::vector<SkyMapPtr> Nmap;
    SkyMapPtr sidereal_map; 


    if (vm.count("remove-sid"))  {
        string sidpath = vm["remove-sid"].as<string>();
        cout << "reading relint distribution from " << sidpath << "..." ;
        sidereal_map = load_siderial_map(sidpath.c_str());
        rmsid = true;
        cout << "done." << endl;
    }

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
    int timeBin;
    int prevTimeBin = 0 ;
    int firstTimeBin = 0 ;
    int prevMJD = 0 ;
    int startMJD = 0 ;
    int currMJD = 0 ;
    bool init=false;
    SkyMapPtr tmp_map;

    if (smethod == "anti") { 
        method_ = config::antisid;
    } else if (smethod == "ext") { 
        method_ = config::extsid;
    } else if (smethod == "solar") {
        method_ = config::solar;
    } else {
        method_ = config::sidereal;
    } 
 


    cout << "Reading " << nEntries << " entries...\n";
    cout << "spline " << vm.count("spline") << " elogmin"<< elogmin << " elogmax" << elogmax << std::endl;

    for (int i = 0; i < nEntries; ++i) {
       cutDST.GetEntry(i); 
       if (newConfig(config)) { 
           trigDST.GetEntry(i); 
       }
    
       if ( show_progress && (i % 5000 == 0) ) 
          progress_bar(i/float(nEntries), " Reading entries");

       if ( (rloglmax > 0 ) && (dst.RLogL > rloglmax) )
          continue;

       
       double mjd = dst.ModJulDay;
       double lst = GetGMST(mjd);
       int mstTimeBin = int( nTimeBins*lst/24.0 ); 

       if (method_ == config::antisid) { 
           lst = GetGMAST(mjd); 
       } else if (method_ == config::extsid) { 
           lst = GetGMEST(mjd);
       } else if (method_ == config::solar) {
           lst= ( mjd - int(mjd) )* 24.; 
       }
 

       // https://en.wikipedia.org/wiki/Time_in_Antarctica
       timeBin = int( nTimeBins*lst/24.0 ); 

       azimuth = dst.LLHAzimuthDeg/180.*M_PI; 
       zenith = dst.LLHZenithDeg/180.*M_PI; 

       if ( costheta && (dst.NDirHits < ndirmin*cos(zenith) ) )
           continue;
       else if ( dst.NDirHits < ndirmin ) 
           continue;

       if ( costheta && (dst.LDir < ldirmin*cos(zenith) ) )
           continue;
       else if ( dst.LDir < ldirmin ) 
           continue;

       // Reconstruction cuts
       if (!dst.isReco || zenith != zenith || azimuth != azimuth) 
           continue;

       // IceTop filter cut
       if (vm.count("filter") && !filterCut(vm, dst))
           continue;

       // Energy cuts for IceTop and IceCube
       if (detector == "IT" && !ITenergyCut(dst, elogmin,elogmax)) 
           continue;

       if (vm.count("spline") && !ICenergyCut(dst, spline, zenith, elogmin, elogmax))
           continue;

       if (detector == "IT" &&  !ITs125Cut(dst, slogmin,slogmax))
           continue;

       pointing pt(zenith, azimuth);
       eventweight = 1.0;

       if (sundp) { 
           Direction dir(zenith,azimuth);
           Equatorial eq = GetEquatorialFromDirection(dir, mjd);
           eventweight = solar_dipole(mjd, eq.ra, eq.dec,true); 
       } else if (rmsid) { 
           int l = sidereal_map->ang2pix(pt);
           int k = illh::eq2loc_idx(l, mstTimeBin, detlat, detlon, nTimeBins, *sidereal_map); 
           eventweight = (1-(*sidereal_map)[k]);
       }
       tmp_map = Nmap[timeBin % 360];// make sure tsid < 360 deg

       while (pt.phi < 0)
          pt.phi += 2.*constants::pi;

       int j = tmp_map->ang2pix(pt);
       (*tmp_map)[j] += eventweight; 
       weightsum += eventweight;

       nEvents++;

    }
    std::cout << std::endl;


    for (unsigned int degbin=0;degbin< 360;degbin++) 
    { 
         fitshandle fitsOut;
         stringstream namefits;
         if (elogmin > 0 || elogmax > 0)
            namefits << foldername << boost::format("/CR_ICECUBE_LOCAL_%4.2f-%4.2fGeV_NSIDE%d_degbin-%03d.fits.gz") % elogmin % elogmax % nsideOut % degbin;
         else
            namefits << foldername << boost::format("/CR_ICECUBE_LOCAL_NSIDE%d_degbin-%03d.fits.gz") % nsideOut % degbin;
         if (fs::exists(namefits.str()) ) { 
                 fs::remove(namefits.str() ); 
         }
         SkyMap jmap = *(Nmap[degbin]);
         std::cout << "bin: " << degbin << " events: "<< Sum(jmap,npix) << " ";
         std::cout << namefits.str() << std::endl;
         fitsOut.create(namefits.str().c_str());
         fitsOut.add_comment("Local Map");

         write_Healpix_map_to_fits(fitsOut, jmap, PLANCK_FLOAT64);
         fitsOut.close();
    }
    std::cout << "passing rate: " << nEvents << "/" << nEntries << std::endl;
    std::cout << "weightsum: " << weightsum << std::endl;

   
}

bool newConfig(string config) {

  if (config=="IC86-2011" || config=="IC86-2012" || config=="IC86-2013" || 
      config=="IC86-2014" || config=="IC86-2015") {
    return false;
  }
  return true;
}


