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
#include <iostream>  
#include <string>
#include <stdexcept>
#include <iterator>
#include <vector>

#include <iter-lhreco-proj/pickle.h>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include <TChain.h>
#include <cuts.h>
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
    std::string method("sid");
    std::string config;
    std::string filter;
    unsigned int nTimesteps;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideOut;
    int index;
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
             ("index", po::value<int>(&index)->default_value(0), "file index")
             ("spline", po::value<string>(), "File containing spline tables")
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum log energy (GeV)")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum log energy (GeV) (0: no cut)")
             ("slogmin", po::value<double>(&slogmin)->default_value(0), "Minimum logS125")
             ("slogmax", po::value<double>(&slogmax)->default_value(0), "Maximum logS125 (0: no cut)")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map") 
             ("sundp", po::bool_switch(&sundp)->default_value(false), "subtract solar dipole")
             ("config", po::value<string>(&config), "Detector configuration") 
             ("method", po::value<string>(&method), "Sidereal, Anti, Solar, Extended")
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
    bool sundp2 = true;
    bool usesplines = false;

    Config cfg(vm);
    if (vm.count("spline") && (cfg.detector == Config::IceCube) )
        usesplines = true;

    if (cfg.detector == Config::IceCube) { 
        masterTree = "CutDST"; 
        triggerTree = "TDSTTriggers"; 
    } 
    if (cfg.detector == Config::IceTop) { 
        triggerTree = "";   // Unused? Will probably break IT functionality...  
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


    TChain trigDST(triggerTree);
    if (cfg.newConfig()) { 
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
    SkyMapPtr mymap(new SkyMap); 
    mymap->SetNside(nsideOut, RING);     
    mymap->fill(0.);

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }


//*****************************************************************************
////// Input ////////////////////////////////////////////////////////////////// 
    fitshandle fitsOut;
    stringstream namefits;

    log_info("Loaded " << input[0] << ", etc.");
    log_info("output file: " << namefits.str() );

    if (elogmin > 0 || elogmax > 0)
            namefits << foldername << boost::format("/CR_ICECUBE_SIDEREAL_%4.2f-%4.2fGeV_NSIDE%d.%03d.fits.gz") % elogmin % elogmax % nsideOut % index;
    else
            namefits << foldername << boost::format("/CR_ICECUBE_SIDEREAL_NSIDE%d.%03d.fits.gz") % nsideOut % index;



    double eventweight;
    double weightsum = 0.;
    float zenith, azimuth;

    for (int i = 0; i < nEntries; ++i) {
       cutDST.GetEntry(i); 
       if (cfg.newConfig()) { 
           trigDST.GetEntry(i); 
       }

       if ( show_progress && (i % 5000 == 0) ) 
          progress_bar(i/float(nEntries), " Reading entries");

       if ( (rloglmax > 0 ) && (dst.RLogL > rloglmax) )
          continue;

       
       double mjd = dst.ModJulDay;
       double lst = GetGMST(mjd);
       if (cfg.method == Config::antisid) { 
           lst = GetGMAST(mjd); 
       } else if (cfg.method == Config::extsid) { 
           lst = GetGMEST(mjd);
       } else if (cfg.method == Config::solar) {
           lst= ( mjd - int(mjd) )* 24.; 
       }
 

       azimuth = dst.LLHAzimuthDeg/180.*M_PI; 
       zenith = dst.LLHZenithDeg/180.*M_PI; 

       if (cfg.detector == Config::IceCube) { 
           azimuth = dst.LLHAzimuthDeg/180.*M_PI; 
           zenith = dst.LLHZenithDeg/180.*M_PI; 

       } else if (cfg.detector == Config::IceTop) { 
           zenith = dst.Zenith;
           azimuth = dst.Azimuth;

       }

       if ( dst.NDirHits < ndirmin*cos(zenith) ) 
          continue;

       if ( dst.LDir < ldirmin*cos(zenith) ) 
           continue;

       // Reconstruction cuts
       if (!dst.isReco || zenith != zenith || azimuth != azimuth) 
           continue;

       // IceTop filter cut
       if (!filterCut(cfg,dst))
           continue;

      // Energy cuts for IceTop and IceCube
        if (cfg.detector == Config::IceTop) {
           if (cfg.cfg == Config::ITv3) { 
               if (!ITNstatCut(dst, slogmin,slogmax)) 
                   continue;
           } else { 
               if (!ITenergyCut(dst, elogmin,elogmax)) 
                   continue; 
               if (!ITs125Cut(dst, slogmin,slogmax)) 
                   continue;
           }
       }
       if (usesplines && !ICenergyCut(dst, spline, zenith, elogmin, elogmax))
           continue;

       if (vm.count("spline") && !ICenergyCut(dst, spline, zenith, elogmin, elogmax))
           continue;


       Direction dir(zenith,azimuth);
       Equatorial eq = GetEquatorialFromDirection(dir, mjd);

       double phi = eq.ra; 
       double theta = constants::pi-eq.dec;
       pointing pt(theta, phi);

       eventweight = 1.0;
       if (sundp) { 
           Direction dir(zenith,azimuth);
           eventweight = solar_dipole(mjd, eq.ra, eq.dec,sundp2); 
       }
       weightsum += eventweight;



       int j = mymap->ang2pix(pt);
       (*mymap)[j] += 1.0; 
       nEvents++;

    }


    if (fs::exists(namefits.str()) ) { 
                 fs::remove(namefits.str() ); 
    }
    std::cout << namefits.str() << std::endl;
    fitsOut.create(namefits.str().c_str());
    fitsOut.add_comment("Local Map");

    write_Healpix_map_to_fits(fitsOut, *mymap, PLANCK_FLOAT64);
    //write_Healpix_map_to_fits(fitsOut, mymap, FITSUTIL<double>::DTYPE);
    fitsOut.close();
    
}



