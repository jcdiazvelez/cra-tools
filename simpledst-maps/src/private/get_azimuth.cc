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

#include <iter-lhreco/pickle.h>
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

#include "cuts.h"
#include "SimpleDST.h"

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
    std::string filename;
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
             ("out", po::value<std::string>(&filename)->default_value("picodst.root"), "output file") 
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
    TFile* outfile = new TFile(filename.c_str(), "recreate");
    

    TChain cutDST("CutDST");
    for (unsigned i = 0; i < input.size(); ++i) { 
        cutDST.Add(input[i].c_str()); 
    } 
    SimpleDST dst(&cutDST, "IC86-2012");

    TTree *newtree = new TTree("picoDST", "picoDST");
    newtree->SetAutoSave(10000);

    double llhazimuth;
    newtree->Branch("LLHAzimuthDeg", &llhazimuth, "llhazimuth/D");


    // Read in spline tables if provided
    photospline::splinetable<> spline;
    if (vm.count("spline")) {
        string splineFile = vm["spline"].as<string>();
        spline.read_fits(splineFile.c_str());
    }




    unsigned int npix = 12*nsideOut*nsideOut; 
    unsigned int nTimeBins = 360;
    // Import data : n_tau_i

    float theta, phi;
    bool init=false;

    for (int i = 0, N = cutDST.GetEntries(); i < N; ++i) {
       cutDST.GetEntry(i); 

       theta = dst.LLHZenithDeg*M_PI/180.;
       phi = dst.LLHAzimuthDeg*M_PI/180.;

       if ( (rloglmax > 0 ) && (dst.RLogL > rloglmax) )
          continue;

       if ( dst.NDirHits < ndirmin*cos(theta) ) 
          continue;

       if ( dst.LDir < ldirmin*cos(phi) ) 
          continue;
       
       if ((vm.count("spline")) && !ICenergyCut(dst, spline, theta, elogmin, elogmax))
          continue;

       llhazimuth = dst.LLHAzimuthDeg;
       newtree->Fill();
       nEvents++;

    }
    outfile->cd();
    outfile->Add(newtree);
    outfile->Write();
    outfile->Close();
    
}
