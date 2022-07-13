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
#include "TFile.h"
#include "TTree.h"

//#include <photospline/splinetable.h>
//#include <photospline/bspline.h>
#include "esplines.h"

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


//bool ICenergyCut(unsigned NChannels, splinetable t, double zenith, double emin, double emax);

int main(int argc, char* argv[])
{
//*****************************************************************************
////// Input ////////////////////////////////////////////////////////////////// 
//*****************************************************************************
    std::string foldername;
    unsigned int nTimesteps;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideOut,index;
    double lon;
    double lat;
    double elogmin, elogmax, rloglmax,ldirmin,ndirmin;
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
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum energy")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum energy (0: no cut)")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map") 
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
//
    fitshandle fitsOut;
    stringstream namefits;
    if (elogmin > 0 || elogmax > 0)
            namefits << foldername << boost::format("/CR_ICECUBE_SIDEREAL_%4.2f-%4.2fGeV_NSIDE%d.%03d.fits.gz") % elogmin % elogmax % nsideOut % index;
    else
            namefits << foldername << boost::format("/CR_ICECUBE_SIDEREAL_NSIDE%d.%03d.fits.gz") % nsideOut % index;

    // TFile f("Event.root");
    TFile file(input[0].c_str());

    TTree* tree = (TTree*) file.Get("CutDST");
    float LLHZenithDeg, LLHAzimuthDeg, RADeg, DecDeg, RLogL;
    double ModJulDay, LocalMST;
    UShort_t NChannels; 
    UInt_t NDirHits, LDir;

    tree->SetBranchAddress("ModJulDay", &ModJulDay);
    tree->SetBranchAddress("LocalMST", &LocalMST);
    tree->SetBranchAddress("LLHZenithDeg", &LLHZenithDeg);
    tree->SetBranchAddress("LLHAzimuthDeg", &LLHAzimuthDeg);
    tree->SetBranchAddress("RADeg", &RADeg);
    tree->SetBranchAddress("DecDeg", &DecDeg);
    tree->SetBranchAddress("NChannels", &NChannels);
    tree->SetBranchAddress("NDirHits", &NDirHits);
    tree->SetBranchAddress("LDir", &LDir);

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
    log_info("Loaded " << input[0] << ", etc.");
    log_info("output file: " << namefits.str() );

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }

    float theta, phi;

    for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
       tree->GetEntry(i); 

       if ( (rloglmax > 0 ) && (RLogL > rloglmax) )
          continue;

       
       //log_info("declination: " << DecDeg );
       //if ( NDirHits < ndirmin ) 
       if ( NDirHits < ndirmin*cos(LLHZenithDeg) ) 
          continue;

       //std::cout << "LDir= " << LDir << " ldirmin = " << ldirmin << std::endl;
       //if ( LDir < ldirmin ) 
       if ( LDir < ldirmin*cos(LLHZenithDeg) ) 
          continue;

       theta = LLHZenithDeg/180.*M_PI; 

       if ((vm.count("spline")) && !ICenergyCut(NChannels, spline, theta, elogmin, elogmax))
          continue;

       phi = RADeg/180.*M_PI; 
       theta = (90-DecDeg)/180.*M_PI; 
       pointing pt(theta, phi);

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
