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
    bool sundp;
    double elogmin,  cxpemin, nhitmin, pincmin;
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
             ("elogmin", po::value<double>(&elogmin)->default_value(3.9), "Minimum energy")
             ("cxpemin", po::value<double>(&cxpemin)->default_value(40), "Minimum cxpe40nch (0: no cut)")
             ("pincmin", po::value<double>(&pincmin)->default_value(1.8), "Minimum pinc (0: no cut)")
             ("nhitmin", po::value<double>(&nhitmin)->default_value(75), "Minimum nHit (0: no cut)")
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

    // TFile f("Event.root");
    TFile file(input[0].c_str());

    log_info("Reading file" << input[0] << " "); 
    TTree* tree = (TTree*) file.Get("XCDF");
    log_info("getting tree XCDF"); 
    float zenith, azimuth, CxPE40XnCh; 
    float protonlheEnergy, pinc, weight;
    UInt_t nHit, angleFitStatus;

    tree->SetBranchAddress("rec.zenithAngle", &zenith);
    tree->SetBranchAddress("rec.azimuthAngle", &azimuth);
    tree->SetBranchAddress("rec.nHit", &nHit);
    tree->SetBranchAddress("rec.CxPE40XnCh", &CxPE40XnCh);
    tree->SetBranchAddress("rec.protonlheEnergy", &protonlheEnergy);
    tree->SetBranchAddress("rec.protonlheEnergy", &protonlheEnergy);
    tree->SetBranchAddress("rec.PINC", &pinc);
    tree->SetBranchAddress("rec.angleFitStatus", &angleFitStatus);
    tree->SetBranchAddress("sweets.IWgt", &weight);

    unsigned int npix = 12*nsideOut*nsideOut; 
    unsigned int nTimeBins = 360;
    // Import data : n_tau_i
    log_info("Creating map " );

    SkyMapPtr locMapDegr(new SkyMap); 
    locMapDegr->SetNside(nsideOut, RING);     
    locMapDegr->fill(0.);
    log_info("Loaded " << input[0] << ", etc.");

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }

    float theta, phi;
    bool init=false;

/*
"sweets.IWgt*(rec.nHit>=75)*(rec.CxPE40XnCh>40)*(rec.protonlheEnergy-9>=3.9)*(rec.PINC > 1.6 )*(rec.angleFitStatus == 0)*(rec.azimuthAngle < 0)*(rec.zenithAngle > 0*TMath::DegToRad())*(rec.zenithAngle < 10*TMath::DegToRad())","",
20093156, 0);
*/
    log_info("Reading " << tree->GetEntries() << " entries"); 
    for (int i = 0, N = tree->GetEntries(); i < N; ++i) {
       tree->GetEntry(i); 

       if ( (nHit < nhitmin) )
          continue;

       if ( CxPE40XnCh <= cxpemin) 
          continue;

       if ( protonlheEnergy-9 < elogmin ) 
          continue;

       if ( pinc <= pincmin ) 
          continue;

       if ( angleFitStatus != 0 ) 
          continue;

       if (azimuth == 0.)
          continue;

       phi = azimuth;
       if (azimuth < 0)
           phi = azimuth+2.0*M_PI; 

       theta = zenith;

       pointing pt(theta, phi);

       int j = locMapDegr->ang2pix(pt);
       (*locMapDegr)[j] += weight; 
       nEvents++;

    }


    fitshandle fitsOut;
    stringstream namefits;
    namefits << foldername << boost::format("/MC_HAWC_LOCAL_NSIDE%d.fits.gz") % nsideOut;
    if (fs::exists(namefits.str()) ) { 
                 fs::remove(namefits.str() ); 
    }
    std::cout << namefits.str() << std::endl;
    fitsOut.create(namefits.str().c_str());
    fitsOut.add_comment("Local MC Map");

    write_Healpix_map_to_fits(fitsOut, *locMapDegr, PLANCK_FLOAT64);
    fitsOut.close();
    
}
