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

#include <photospline/splinetable.h>
#include <photospline/bspline.h>

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


bool ICenergyCut(unsigned NChannels, splinetable t, double zenith, double emin, double emax);

int main(int argc, char* argv[])
{
//*****************************************************************************
////// Input ////////////////////////////////////////////////////////////////// 
//*****************************************************************************
    std::string foldername;
    unsigned int nTimesteps;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideOut;
    bool sundp;
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
             ("spline", po::value<string>(), "File containing spline tables")
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum energy")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum energy (0: no cut)")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map") 
             ("sundp", po::bool_switch(&sundp)->default_value(false), "subtract solar dipole")
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

    // TFile f("Event.root");
    TFile file(input[0].c_str());

    TTree* tree = (TTree*) file.Get("CutDST");
    float LLHZenithDeg, LLHAzimuthDeg, RADeg, DecDeg, RLogL, RaSun, DecSun;
    double ModJulDay, LocalMST;
    UShort_t NChannels; 
    UInt_t NDirHits, LDir;

    tree->SetBranchAddress("ModJulDay", &ModJulDay);
    tree->SetBranchAddress("LocalMST", &LocalMST);
    tree->SetBranchAddress("LLHZenithDeg", &LLHZenithDeg);
    tree->SetBranchAddress("LLHAzimuthDeg", &LLHAzimuthDeg);
    tree->SetBranchAddress("RADeg", &RADeg);
    tree->SetBranchAddress("DecDeg", &DecDeg);
    tree->SetBranchAddress("RaSun", &RaSun);
    tree->SetBranchAddress("DecSun", &DecSun);
    tree->SetBranchAddress("NChannels", &NChannels);
    tree->SetBranchAddress("NDirHits", &NDirHits);
    tree->SetBranchAddress("LDir", &LDir);

    // Read in spline tables if provided
    struct splinetable table;
    if (vm.count("spline")) {
       string splineFile = vm["spline"].as<string>();
       readsplinefitstable(splineFile.c_str(), &table);
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

    float theta, phi;
    int timeBin;
    int prevTimeBin = 0 ;
    int firstTimeBin = 0 ;
    int prevMJD = 0 ;
    int startMJD = 0 ;
    int currMJD = 0 ;
    bool init=false;
    SkyMapPtr tmp_map;

    for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
       tree->GetEntry(i); 

       if ( (rloglmax > 0 ) && (RLogL > rloglmax) )
          continue;

       
       //std::cout << "NDirHits=" << NDirHits << " ndirmin = " << ndirmin << std::endl;
       //if ( NDirHits < ndirmin ) 
       if ( NDirHits < ndirmin*cos(LLHZenithDeg) ) 
          continue;

       //std::cout << "LDir= " << LDir << " ldirmin = " << ldirmin << std::endl;
       //if ( LDir < ldirmin ) 
       if ( LDir < ldirmin*cos(LLHZenithDeg) ) 
          continue;

       currMJD = ModJulDay;
       // https://en.wikipedia.org/wiki/Time_in_Antarctica
       timeBin = int( nTimeBins*LocalMST/24.0 ); 

       phi = LLHAzimuthDeg/180.*M_PI; 
       //phi = RADeg/180.*M_PI; 
       theta = LLHZenithDeg/180.*M_PI; 
       //theta = DecDeg/180.*M_PI; 
       
       if ((vm.count("spline")) && !ICenergyCut(NChannels, table, theta, elogmin, elogmax))
          continue;

       pointing pt(theta, phi);

       tmp_map = Nmap[timeBin % 360];// make sure tsid < 360 deg
       int j = tmp_map->ang2pix(pt);
       (*tmp_map)[j] += 1.0; 
       nEvents++;

    }


    for (unsigned int degbin=0;degbin< 360;degbin++) 
    { 
         fitshandle fitsOut;
         stringstream namefits;
         if (elogmin > 0 || elogmax > 0)
            namefits << foldername << boost::format("/CR_ICECUBE_LOCAL_%4.2f-%4.2fGeV_NSIDE%d_degbin-%02d.fits.gz") % elogmin % elogmax % nsideOut % degbin;
         else
            namefits << foldername << boost::format("/CR_ICECUBE_LOCAL_NSIDE%d_degbin-%02d.fits.gz") % nsideOut % degbin;
         if (fs::exists(namefits.str()) ) { 
                 fs::remove(namefits.str() ); 
         }
         std::cout << namefits.str() << std::endl;
         fitsOut.create(namefits.str().c_str());
         fitsOut.add_comment("Local Map");

         SkyMap jmap = *(Nmap[degbin]);
         write_Healpix_map_to_fits(fitsOut, jmap, PLANCK_FLOAT64);
         //write_Healpix_map_to_fits(fitsOut, jmap, FITSUTIL<double>::DTYPE);
         fitsOut.close();
    }
    
}
