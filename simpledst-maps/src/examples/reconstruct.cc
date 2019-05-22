#include <healpix-cxx/healpix/alm.h>
#include <healpix-cxx/cxxsupport/xcomplex.h>
#include <healpix-cxx/healpix/alm_healpix_tools.h>
#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
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
    double lon;
    double lat;
    std::vector< std::string > input;

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
             ("nsideout", po::value<int>(&nsideOut)->default_value(64), "Healpix Nside for output map") 
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
    float LLHZenithDeg, LLHAzimuthDeg, RADeg;
    tree->SetBranchAddress("LLHZenithDeg", &LLHZenithDeg);
    tree->SetBranchAddress("LLHAzimuthDeg", &LLHAzimuthDeg);
    tree->SetBranchAddress("RADeg", &RADeg);


    unsigned int npix = 12*nsideOut*nsideOut; 
    // Import data : n_tau_i
    std::vector<SkyMapPtr> Nmap;
    for (unsigned int degbin=0;degbin< 360;degbin++) 
    {
        SkyMapPtr locMapDegr(new SkyMap); 
        locMapDegr->SetNside(nsideOut, RING);     
        locMapDegr->fill(0.);
        Nmap.push_back(locMapDegr);
    }
    //log_info("Loaded " << input[0] << ", etc.");

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }

    double theta, phi;
    int imap;
    SkyMapPtr tmp_map;

    for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
       tree->GetEntry(i); 

       imap = int(floor(RADeg)) % 360;

       phi = LLHAzimuthDeg/180.*M_PI; 
       theta = LLHZenithDeg/180.*M_PI; 
       pointing pt(theta, phi);

       tmp_map = Nmap[imap];
       int j = tmp_map->ang2pix(pt);
       (*tmp_map)[j] += 1.0;
       //std::cout << "imap:" << imap << " LLHAzimuthDeg:" << LLHAzimuthDeg << " LLHZenithDeg:" << LLHZenithDeg << std::endl;
    }


    for (unsigned int degbin=0;degbin< 360;degbin++) 
    { 
         fitshandle fitsOut;
         stringstream namefits;
         namefits << foldername << boost::format("CR_ICECUBE_LOCAL_NSIDE%d_degbin%02d.fits") % nsideOut % degbin;
         if (fs::exists(namefits.str()) ) { 
                 fs::remove(namefits.str() ); 
         }
         std::cout << namefits.str() << std::endl;
         fitsOut.create(namefits.str().c_str());
         SkyMap jmap = *(Nmap[degbin]);
         //write_Healpix_map_to_fits(fitsOut, jmap, PLANCK_FLOAT64);
         write_Healpix_map_to_fits(fitsOut, jmap, FITSUTIL<double>::DTYPE);
         fitsOut.close();
    }
    
}
