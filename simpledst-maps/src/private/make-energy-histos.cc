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

#include <algorithm>           // std::for_each
//#include <boost/histogram.hpp> // make_histogram, regular, weight, indexed
#include <functional>          // std::ref

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
#include "histogram.h"

#include <astro/astro.h>
#include <astro/time.h>
#include <util/cnpy.h>
#include <Direction.h>
#include <solardipole.h>
#include <SimpleTrigger.h>


namespace po = boost::program_options; 
#define log_info(args...) cout << args << endl

namespace fs = boost::filesystem; 
using namespace std;
using namespace boost::numeric::ublas;
using boost::format;
//using namespace boost::histogram; // strip the boost::histogram prefix

typedef Healpix_Map<double> SkyMap; 
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer

// use this when axis configuration is fix to get highest performance





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
    std::string outfile;
    std::string smethod("sid");
    std::string config;
    std::string filter;
    config::method method_;
    unsigned int nTimesteps;
    unsigned int nbins;
    int nsideOut;
    bool sundp;
    bool show_progress;
    bool rmsid(false);
    double detlon, detlat;
    double decmin, decmax;
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
             ("outfile", po::value<std::string>(&outfile)->default_value("ehisto.npz"), "Output file") 
             ("spline", po::value<string>(), "File containing spline tables")
             ("elogmin", po::value<double>(&elogmin)->default_value(0), "Minimum log energy (GeV)")
             ("elogmax", po::value<double>(&elogmax)->default_value(0), "Maximum log energy (GeV) (0: no cut)")
             ("decmin", po::value<double>(&decmin)->default_value(-90), "Minimum declination (degrees)")
             ("decmax", po::value<double>(&decmax)->default_value(-25), "Maximum declination (degrees)")
             ("slogmin", po::value<double>(&slogmin)->default_value(0), "Minimum logS125")
             ("slogmax", po::value<double>(&slogmax)->default_value(0), "Maximum logS125 (0: no cut)")
             ("ldirc_min", po::value<double>(&ldirmin)->default_value(0), "Minimum ldir (0: no cut)")
             ("ndirc_min", po::value<double>(&ndirmin)->default_value(0), "Minimum ndir (0: no cut)")
             ("rloglmax", po::value<double>(&rloglmax)->default_value(0), "Maximum rlogL (0: no cut)")
             ("nsideout", po::value<int>(&nsideOut)->default_value(256), "Healpix Nside for output map") 
             ("sundp", po::bool_switch(&sundp)->default_value(false), "subtract solar dipole")
             ("nbins", po::value<unsigned int>(&nbins)->default_value(20), "Number of bins in energy histogram")
             ("lon", po::value<double>(&detlon)->default_value(270.0), "Longitude of detector") 
             ("lat", po::value<double>(&detlat)->default_value(-90.0), "Latitude of detector")
             ("config", po::value<string>(&config), "Detector configuration") 
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
    decmin = decmin*units::degree;
    decmax = decmax*units::degree;
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
    std::vector<FloatHistogramPtr> Nhist;

    SkyMapPtr locMapDegr(new SkyMap); 
    locMapDegr->SetNside(nsideOut, RING);     
    locMapDegr->fill(0.);

    pointing minpt(0.5*constants::pi - decmax, 0.0);
    cout << "mintheta:" << minpt.theta << endl;
    pointing maxpt(0.5*constants::pi - decmin, 0.0);
    cout << "maxtheta:" << maxpt.theta << endl;
    unsigned minPix = locMapDegr->ang2pix(minpt);
    unsigned maxPix = locMapDegr->ang2pix(maxpt);
    unsigned visiblenpix = maxPix-minPix+1;

    //hist_t h = make_histogram(axis::regular<double>(6, -1.0, 2.0));
    for (unsigned int ipix=0;ipix< visiblenpix;ipix++) 
    {
        Nhist.push_back(FloatHistogramPtr(new FloatHistogram(nbins, elogmin, elogmax)));
    }
    log_info("Loaded "<< input[0] << ", etc.");

    // output directorrry
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }

    double eventweight;
    double weightsum = 0.;
    float zenith, azimuth;
    double recoEnergy;
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


       if ( dst.NDirHits < ndirmin*cos(zenith) ) 
           continue;

       if ( dst.LDir < ldirmin*cos(zenith) ) 
           continue;

       // Reconstruction cuts
       if (!dst.isReco || zenith != zenith || azimuth != azimuth) 
           continue;

       // IceTop filter cut
       //if (vm.count("filter") && !filterCut(vm, dst))
       //    continue;

       // Energy cuts for IceTop and IceCube
       //if (detector == "IT" && !ITenergyCut(dst, elogmin,elogmax)) 
       //    continue;

       //if (vm.count("spline") && !ICenergyCut(dst, spline, zenith, elogmin, elogmax))
       //    continue;

       //if (detector == "IT" &&  !ITs125Cut(dst, slogmin,slogmax))
       //    continue;

       //if (vm.count("spline")) { 
           recoEnergy = ICenergy(dst, spline, zenith);
           //recoEnergy = log10(dst.NChannels);

           if (recoEnergy < 0) 
               continue;
       //}

       Direction dir(zenith,azimuth);
       Equatorial eq = GetEquatorialFromDirection(dir, mjd);


       pointing pt(0.5*constants::pi - eq.dec, eq.ra);
       while (pt.phi < 0)
          pt.phi += 2.*constants::pi;

       int jpix = locMapDegr->ang2pix(pt);

       if ( jpix < minPix || jpix > maxPix )
           continue;


       FloatHistogramPtr pixhist = Nhist[jpix-minPix];
       pixhist->insert(recoEnergy);

       nEvents++;

    }
    std::cout << std::endl;

    std::cout << "passing rate: " << nEvents << "/" << nEntries << std::endl;
    std::cout << "weightsum: " << weightsum << std::endl;
    stringstream filename;
    filename << foldername << "/" << outfile;
    std::vector<float> edges = Nhist[0]->get_edges();
    std::vector<size_t> eshape = {1,Nhist[0]->size()};
    std::vector<size_t> bshape = {visiblenpix,Nhist[0]->size()};
    std::vector<double> data;

    cnpy::npz_save(filename.str(),"edges",&edges[0],eshape,"w");

    for (unsigned jpix = 0; jpix < visiblenpix;jpix++)
    {
        //std::vector<double> data;
        std::vector<double> tmp = Nhist[jpix]->get_bins();
        data.reserve(data.size() + distance(tmp.begin(),tmp.end()));
        data.insert(data.end(),tmp.begin(),tmp.end());
    }
    cnpy::npz_save(filename.str(),"bins",&data[0],bshape,"a");
   
}

bool newConfig(string config) {

  if (config=="IC86-2011" || config=="IC86-2012" || config=="IC86-2013" || 
      config=="IC86-2014" || config=="IC86-2015") {
    return false;
  }
  return true;
}


