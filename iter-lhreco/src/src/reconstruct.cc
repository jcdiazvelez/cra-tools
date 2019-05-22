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

#ifdef HAWCNEST
#include <hawcnest/HAWCNest.h>
#include <hawcnest/Logging.h>
#include <hawcnest/CommandLineConfigurator.h>
#include <detector-service/ConfigDirDetectorService.h>
#define MyDTYPE PLANCK_FLOAT64

#else
namespace po = boost::program_options; 
#define log_info(args...) cout << args << endl
#define log_fatal(args...) throw args
#define MyDTYPE FITSUTIL<double>::DTYPE 
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
    std::string detector;
    unsigned int nTimesteps;
    unsigned int timeidxMin;
    unsigned int timeidxMax;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideIn;
    int nsideOut;
    double lon;
    double lat;
    std::vector< std::string > input;

#ifdef HAWCNEST
    CommandLineConfigurator cl;
  
    cl.AddPositionalOption<std::vector <std::string> >("input", "HEALPix local maps (theta/phi) in timesteps of a sidereal day");
    cl.AddOption<std::string>("outdir","./sample/","Directory of output");
    cl.AddOption<std::string>("detector","HAWC","name of detector");
    cl.AddOption<int>("nsideout",64,"nSide of output maps");
    cl.AddOption<int>("iterations",20,"Number of iterations");
    cl.AddOption<int>("timesteps",360,"Timesteps in sidereal day");
    cl.AddOption<int>("timestepmin",0,"first time step to use");
    cl.AddOption<int>("timestepmax",0,"last time step to use (default = timesteps)");
    cl.AddOption<int>("sectors",1,"Number of FOV sectors to use");
    cl.AddOption<double>("lon",-97.0,"Longitude of detector");
    cl.AddOption<double>("lat",19.0,"Latitude of detector");
    if (!cl.ParseCommandLine(argc, argv)) {
      return 1;
    }
    

    input = cl.GetArgument<std::vector <std::string> >("input");
    foldername = cl.GetArgument<std::string>("outdir");
    detector = cl.GetArgument<std::string>("detector");
    nTimesteps = cl.GetArgument<int>("timesteps");
    timeidxMin = cl.GetArgument<int>("timestepmin");
    timeidxMax = cl.GetArgument<int>("timestepmax");
    if (timeidxMax == 0){
      timeidxMax=nTimesteps;
    }
    nIterations = cl.GetArgument<int>("iterations");
    nSectors = cl.GetArgument<int>("sectors");
    nsideOut = cl.GetArgument<int>("nsideout");
    lon = cl.GetArgument<double>("lon");
    lat = cl.GetArgument<double>("lat");
#else

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
             ("timesteps", po::value<unsigned int>(&nTimesteps)->default_value(360), "Number of time steps") 
             ("timestepmin", po::value<unsigned int>(&timeidxMin)->default_value(0), "First time step to use") 
             ("timestepmax", po::value<unsigned int>(&timeidxMax)->default_value(0), "Last time step to use") 
             ("iterations", po::value<unsigned int>(&nIterations)->default_value(20), "Number of iterations") 
             ("sectors", po::value<unsigned int>(&nSectors)->default_value(1), "Number sectors") 
             ("lon", po::value<double>(&lon)->default_value(-97.0), "Longitude of detector") 
             ("lat", po::value<double>(&lat)->default_value(19.0), "Latitude of detector")
			("detector", 
			 po::value<string>(&detector)->default_value("HAWC"), 
			 "Name of detector"); 
     
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
     
        if (timeidxMax==0){
            timeidxMax=nTimesteps;
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
#endif 
 
//*****************************************************************************
////// Initialize ///////////////////////////////////////////////////////////// 
//*****************************************************************************

    unsigned int npix = 12*nsideOut*nsideOut; 

    fitshandle handle;
    // Import data : n_tau_i
    std::vector<SkyMapPtr> Nmap;
    std::string cuts;
    std::string coords;
    std::string timesys;
    double startMJD = -1.;
    double stopMJD = -1.;
    double totDur = 0.;
    double nEventsFio = 0.;
    int nTimeBins = -1;
    sort(input.begin(), input.end());
    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    {
        //double nEvents = 0.0; 

        // Read header info and then map
        handle.open(input[timeidx].c_str());
        log_info("reading file : " << input[timeidx] );

        std::string cutsTemp = "";
        if (handle.key_present("CUTS")){
          handle.get_key("CUTS", cutsTemp);
        }
        string coordsTemp = "";  
        if (handle.key_present("COORDS")){
          handle.get_key("COORDS", coordsTemp);
        }
        double nEventsTemp = 0;
        if (handle.key_present("NEVENTS")){
          handle.get_key("NEVENTS",  nEventsTemp);
        }
        double startMJDTemp = 0;
        if (handle.key_present("STARTMJD")){
          handle.get_key("STARTMJD", startMJDTemp);
        }
        double stopMJDTemp = 0;
        if (handle.key_present("STOPMJD")){
          handle.get_key("STOPMJD", stopMJDTemp);
        }
        double totDurTemp = 0;
        if (handle.key_present("TOTDUR")){
          handle.get_key("TOTDUR",   totDurTemp);
        }
        std::string timesysTemp = "";
        if (handle.key_present("TIMESYS")){
          handle.get_key("TIMESYS", timesysTemp);
        }
        int bin = -1;
        if (handle.key_present("TIMEBIN")){
            handle.get_key("TIMEBIN",   bin);
        }
        int nTimeBinsTemp = -1;
        if (handle.key_present("NBINS")) {
            handle.get_key("NBINS",  nTimeBinsTemp);
        }

        handle.goto_hdu(2);
        SkyMapPtr locMap(new SkyMap);
        read_Healpix_map_from_fits(handle, *locMap, 1);
        handle.close();

        if (timeidx == timeidxMin) {
            nsideIn = locMap->Nside();
            log_info("Old Nside: " << nsideIn << ", New Nside: " << nsideOut);
            cuts = cutsTemp;
            coords = coordsTemp;
            timesys = timesysTemp;
            startMJD = startMJDTemp;
            stopMJD = stopMJDTemp;
            nTimeBins = nTimeBinsTemp;
        }

        totDur += totDurTemp;
        nEventsFio += nEventsTemp;

        if (cuts != cutsTemp) {
            log_fatal("Cuts are different");
        }
        if (coords != coordsTemp && !coordsTemp.empty()) {
            log_fatal("Coordinate systems are different");
        }
        if (timesys != timesysTemp && !timesysTemp.empty()) {
            log_fatal("Time systems are different");
        }
        if (nTimeBins != nTimeBinsTemp) {
            log_fatal("Number of time bins are different");
        }
        if (bin != timeidx && bin > 0) {
            log_fatal("Time Bin ID does not match Time Index of Loop. Maps may be out of order.");
        }

        if ( nsideIn != nsideOut) {
            SkyMapPtr locMapDegr(new SkyMap); 
            locMapDegr->SetNside(nsideOut, RING);     
            locMapDegr->fill(0.);
            for (int i = 0; i < 12*nsideIn*nsideIn; i++){
                pointing pt        = locMap->pix2ang(i);
                int j              = locMapDegr->ang2pix(pt);
                (*locMapDegr)[j] += (*locMap)[i];
            }
            Nmap.push_back(locMapDegr);
        } else {
            Nmap.push_back(locMap);
        }
    }
    log_info("Loaded " << input[timeidxMin] << ", etc.");

    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }
    
    // detector position
    lon = lon/180.*M_PI;
    lat = lat/180.*M_PI;
#ifdef HAWCNEST
    if (detector == "HAWC"){
      string configDir = "";
      if (getenv("CONFIG_HAWC")){
        configDir = getenv("CONFIG_HAWC");
      } else {
        log_fatal("config-dir is not defined by CONFIG_HAWC");
      }
      HAWCNest nest;
      nest.Service<ConfigDirDetectorService>("det")
        ("configDir",configDir);
      nest.Configure();
      const det::Detector& detector = GetService<DetectorService>("det").GetDetector(TimeStamp(1111111111));

      const LatLonAlt loc = detector.GetLatitudeLongitudeHeight();
      lon = -loc.GetLongitude();
      lat = loc.GetLatitude();
    }
#endif
    double clat = cos(lat);
    double slat = sin(lat);
    log_info("Latitude=" << lat << ",  Longitude=" << lon);
    
    // normalization of isotropic flux
    double isovalue = 1.0;
    // initial CR anisotropy : I_a^(0) = 1.
    SkyMap CRmap;
    CRmap.SetNside(nsideOut, RING);
    CRmap.fill( isovalue );
                
    // initial exposure : A_i^(0)
    SkyMap Emap0;
    Emap0.SetNside(nsideOut, RING);
    Emap0.fill(0.);
    // n^th exposure : A_i^(n)
    SkyMap Emap;
    Emap.SetNside(nsideOut, RING);
    Emap.fill(0.);
    double nEventsTot = 0.0;
    for (unsigned int i=0; i < npix;i++ ) 
    { 
            double nEvents = 0.0; 
            //for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    nEvents    += (*Nmap[timeidx])[i];
                    nEventsTot += (*Nmap[timeidx])[i];
            } 
            Emap0[i] = nEvents; 
    }
            
    for (unsigned int i=0; i < npix;i++ ) 
    { 
            Emap0[i] = Emap0[i]/nEventsTot; 
            Emap[i] = Emap0[i]; 
    }
    log_info("Total events " << nEventsTot);

    // initial normalization : N_tau^(0)
    std::vector<double> norm0(nTimesteps);
    // n^th normalization : N_tau^(n)
    std::vector<double> norm(nTimesteps);
    //for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
    { 
             double nEvents = 0.0;
             for (unsigned int i=0; i < npix;i++ ) 
             { 
                     nEvents += (*Nmap[timeidx])[i];
             } 
             norm[timeidx] = nEvents/isovalue; 
             norm0[timeidx] = norm[timeidx]; 
    }
        
    // data in equatorial coords : n_a
    SkyMap dataMap;
    dataMap.SetNside(nsideOut, RING); 
    dataMap.fill(0.);
    //// window function of FOV map : F_a
    //SkyMap FOVmap;
    //FOVmap.SetNside(nsideOut, RING); 
    //FOVmap.fill(0.);
    for (unsigned int i=0; i < npix;i++ ) 
    { 
            double nBkg = 0.; 
            double nEvents = 0.; 
            vec3 v =  CRmap.pix2vec(i);
            
            //for (unsigned int timeidx=0; timeidx < nTimesteps;timeidx++ ) 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon; 
                    double cb = cos(beta); 
                    double sb = sin(beta);
        
                    // rotation matrix
                    matrix<double> rot (3, 3);

                    rot(0,0) = cb*slat;
                    rot(0,1) = sb*slat;
                    rot(0,2) = -clat;

                    rot(1,0) = -sb;
                    rot(1,1) = cb;
                    rot(1,2) = 0;

                    rot(2,0) = cb*clat;
                    rot(2,1) = clat*sb;
                    rot(2,2) = slat;

                    // rotation from local frame to Equatorial (ra,dec)
                    vec3 vp;
                    vp.x = rot(0,0)*v.x+rot(0,1)*v.y+rot(0,2)*v.z;
                    vp.y = rot(1,0)*v.x+rot(1,1)*v.y+rot(1,2)*v.z;
                    vp.z = rot(2,0)*v.x+rot(2,1)*v.y+rot(2,2)*v.z;
        
                    int j = CRmap.vec2pix(vp); // local pixel
                    if (Emap0[j] > 0.0 ){ 
                            nEvents += (*Nmap[timeidx])[j];
                            nBkg    += norm[timeidx]*Emap0[j]; 
                    }
                } 

                dataMap[i]      = nEvents;

                //if (nBkg > 0.0) { 
                //        FOVmap[i] = 1; 
                //} 
    } 

    //double npixFOV = 1.*Sum(FOVmap,npix); 
    
    //  write N_tau^(0)
    stringstream normName;
    normName << foldername << boost::format("/norm_%s_%d_%d_iteration00.fits") % detector % nsideOut % nTimesteps;
    ////pickle::dump(norm, normName.str(), "\n");
    std::ofstream fileout(normName.str().c_str());
    for(std::vector<double>::iterator it = norm.begin(); it != norm.end(); ++it) {
      fileout << *it << "\n";
    }
    fileout.close();

    // write A_i^(0)
    fitshandle fitsOut; 
    stringstream exposuremapname;
    exposuremapname <<  foldername << boost::format("/exposure_%s_%d_%d_iteration00.fits") % detector % nsideOut % nTimesteps;
    if (fs::exists(exposuremapname.str()) ) { 
            fs::remove( exposuremapname.str() ); 
    }
    fitsOut.create(exposuremapname.str().c_str()); 
    //write_Healpix_map_to_fits(fitsOut, Emap0, PLANCK_FLOAT64); 
    write_Healpix_map_to_fits(fitsOut, Emap0, MyDTYPE);
    fitsOut.close(); 

    // write n_a
    stringstream datamapname;
    datamapname <<  foldername << boost::format("/data_%s_%d_%d.fits") % detector % nsideOut % nTimesteps;
    if (fs::exists(datamapname.str()) ) {
            fs::remove( datamapname.str() );
    }
    fitsOut.create(datamapname.str().c_str());
    //write_Healpix_map_to_fits(fitsOut, dataMap, PLANCK_FLOAT64);
    write_Healpix_map_to_fits(fitsOut, dataMap, MyDTYPE);
    fitsOut.close();

//*****************************************************************************
////// Iterate //////////////////////////////////////////////////////////////// 
//*****************************************************************************

    // n^th differential CR anisotropy : I_a^(n) - isovalue
    SkyMap diffCRmap;
    diffCRmap.SetNside(nsideOut, RING);
    diffCRmap.fill(0.);
    for (unsigned int iteration = 1; iteration <= nIterations; iteration++)
    { 

            // calculate new CR anisotropy
            for (unsigned int i=0; i<npix;i++) 
            { 
                    double nEvents = 0.0;
                    double nBkg = 0.0;
                    vec3 v =  CRmap.pix2vec(i); 


                    //for (unsigned int timeidx=0;timeidx < nTimesteps;timeidx++)
                    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                    { 
                            double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon;
                            double cb = cos(beta);
                            double sb = sin(beta); 

                            // rotation matrix 
                            matrix<double> rot (3, 3); 

                            rot(0,0) = cb*slat; 
                            rot(0,1) = sb*slat; 
                            rot(0,2) = -clat; 
                                
                            rot(1,0) = -sb; 
                            rot(1,1) = cb; 
                            rot(1,2) = 0; 
                                
                            rot(2,0) = cb*clat; 
                            rot(2,1) = clat*sb; 
                            rot(2,2) = slat; 
                                
                            // rotation from local frame to Equatorial (ra,dec) 
                            vec3 vp; 
                            vp.x = rot(0,0)*v.x+rot(0,1)*v.y+rot(0,2)*v.z;
                            vp.y = rot(1,0)*v.x+rot(1,1)*v.y+rot(1,2)*v.z;
                            vp.z = rot(2,0)*v.x+rot(2,1)*v.y+rot(2,2)*v.z;
                                                
                            int j = CRmap.vec2pix(vp); // local pixel

                            if (Emap0[j] > 0.0) { 
                                    nEvents += (*Nmap[timeidx])[j];
                                    nBkg += norm[timeidx]*Emap[j]; 
                            } 

                    }
                        
                    if (nBkg > 0.0) { 
                            diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
                    }
            } 

            // remove m=0 multipole moments :
            const int LMAX=180;
            // Initialize the spherical harmonic coefficients of the map
            Alm<xcomplex<double> > alm( LMAX, LMAX);
            arr<double> weight(2*nsideOut);
            weight.fill(1.);
            const int numIter = 3;
            map2alm_iter( diffCRmap, alm, numIter, weight);
            for (unsigned int l=0; l < LMAX+1; l++)
            {    
                alm(l,0) = 0.;
            }
            alm2map( alm, diffCRmap);

            // calculate new norm 
            //for (unsigned int timeidx=0; timeidx < nTimesteps; timeidx++) 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                double nEvents = 0.0;
                double nBkg = 0.0; 

                double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon; 
                double cb = cos(beta);
                double sb = sin(beta); 
                        
                for (unsigned int i=0; i<npix; i++) 
                { 
                        vec3 v = CRmap.pix2vec(i); 
            
                        // rotation matrix 
                        matrix<double> rot (3, 3); 

                        rot(0,0) = cb*slat; 
                        rot(0,1) = sb*slat; 
                        rot(0,2) = -clat; 
                                
                        rot(1,0) = -sb; 
                        rot(1,1) = cb; 
                        rot(1,2) = 0; 
                                
                        rot(2,0) = cb*clat; 
                        rot(2,1) = clat*sb; 
                        rot(2,2) = slat; 
                                
                        // rotation from local frame to Equatorial (ra,dec)
                        vec3 vp; 
                        vp.x = rot(0,0)*v.x+rot(1,0)*v.y+rot(2,0)*v.z;
                        vp.y = rot(0,1)*v.x+rot(1,1)*v.y+rot(2,1)*v.z;
                        vp.z = rot(0,2)*v.x+rot(1,2)*v.y+rot(2,2)*v.z;

                        int j = CRmap.vec2pix(vp); // Equatorial pixel 

                        if (Emap0[i] > 0.0) { 
                                nEvents += (*Nmap[timeidx])[i];
                                nBkg += Emap[i]*(CRmap[j]+diffCRmap[j]);
                        } 
                } 
                if (nBkg > 0.0) { 
                        norm[timeidx] = nEvents/nBkg; 
                } 
            } 

            // caluculate new acceptance
            for (unsigned int i=0; i<npix; i++) 
            {
                double nEvents = 0.0;
                double nBkg = 0.0; 
                vec3 v = CRmap.pix2vec(i); 
            
                //idx = sectormap[i]
                //if idx < 0 :
                //    continue
                        
                //for (unsigned int timeidx=0; timeidx < nTimesteps; timeidx++) 
                for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                {
                    double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon;
                    double cb = cos(beta);
                    double sb = sin(beta);
            
                    // rotation matrix 
                    matrix<double> rot (3, 3); 

                    rot(0,0) = cb*slat; 
                    rot(0,1) = sb*slat; 
                    rot(0,2) = -clat; 
                            
                    rot(1,0) = -sb; 
                    rot(1,1) = cb; 
                    rot(1,2) = 0; 
                            
                    rot(2,0) = cb*clat; 
                    rot(2,1) = clat*sb; 
                    rot(2,2) = slat; 
                                
                    // rotation from Equatorial (ra,dec) to local frame
                    vec3 vp;
                    vp.x = rot(0,0)*v.x+rot(1,0)*v.y+rot(2,0)*v.z;
                    vp.y = rot(0,1)*v.x+rot(1,1)*v.y+rot(2,1)*v.z;
                    vp.z = rot(0,2)*v.x+rot(1,2)*v.y+rot(2,2)*v.z;

                    int j = CRmap.vec2pix(vp); // Equatorial pixel 

                    if (Emap0[i] > 0.0) {
                            nEvents += (*Nmap[timeidx])[i];
                            nBkg += norm[timeidx]*(CRmap[j]+diffCRmap[j]);
                    }
                }
            
                if (nBkg > 0.0) 
                {
                    Emap[i] = nEvents/nBkg;
                } else
                { 
                    Emap[i] = 0.0;
                }
            }
            
            // Renormalize
            double sumEmap = Sum(Emap, npix);
            
            //for (unsigned int timeidx=0;timeidx < nTimesteps;timeidx++) 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                norm[timeidx] = norm[timeidx]*sumEmap; 
            } 
            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    Emap[i] = Emap[i]/sumEmap;
            } 
            
            SkyMap diffCRmapNormed; 
            diffCRmapNormed.SetNside(nsideOut, RING); 
            diffCRmapNormed.fill(0.);

            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    diffCRmapNormed[i] = diffCRmap[i]/isovalue; 
            } 
            
            // write I_a^(n)
            stringstream namefits;
            namefits << foldername << boost::format("/CR_%s_%d_%d_iteration%02d.fits") % detector % nsideOut % nTimesteps % iteration;
            if (fs::exists(namefits.str()) ) { 
                    fs::remove( namefits.str() ); 
            }
            fitshandle fitsOut; 
            fitsOut.create(namefits.str().c_str()); 
            //write_Healpix_map_to_fits(fitsOut, diffCRmapNormed, PLANCK_FLOAT64); 
            write_Healpix_map_to_fits(fitsOut, diffCRmapNormed, MyDTYPE);
            fitsOut.close(); 

            // write A_i^(n)
            stringstream expmapname;
            expmapname << foldername << boost::format("/exposure_%s_%d_%d_iteration%02d.fits") % detector % nsideOut % nTimesteps % iteration;
            if (fs::exists(expmapname.str()) ) {
                    fs::remove( expmapname.str() );
            }
            fitsOut.create(expmapname.str().c_str());
            //write_Healpix_map_to_fits(fitsOut, Emap, PLANCK_FLOAT64);
            write_Healpix_map_to_fits(fitsOut, Emap, MyDTYPE);
            fitsOut.close();
        
            // write N_tau^(n)
            stringstream nameNfits;
            nameNfits << foldername <<  boost::format("/norm_%s_%d_%d_iteration%02d.dat") % detector % nsideOut % nTimesteps % iteration; 
            ////pickle::dump(norm, nameNfits.str(), "\n");
            std::ofstream fileout2(nameNfits.str().c_str());
            for(std::vector<double>::iterator it = norm.begin(); it != norm.end(); ++it) {
              fileout2 << *it << "\n";
            }
            fileout2.close();

            // calculate signifiance : S_a^(n)
            SkyMap significancemap; 
            significancemap.SetNside(nsideOut, RING);
            significancemap.fill(0.);
            
            //for (unsigned int timeidx=0;timeidx < nTimesteps;timeidx++) 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    
                    double beta = timeidx/(1.*nTimesteps)*M_PI*2 + lon; 
                    double cb = cos(beta); 
                    double sb = sin(beta);
        
                    vec3 vp;
            
                    // rotation matrix 
                    matrix<double> rot (3, 3); 

                    rot(0,0) = cb*slat; 
                    rot(0,1) = sb*slat; 
                    rot(0,2) = -clat; 
            
                    rot(1,0) = -sb; 
                    rot(1,1) = cb; 
                    rot(1,2) = 0; 
                    
                    rot(2,0) = cb*clat; 
                    rot(2,1) = clat*sb; 
                    rot(2,2) = slat; 
                                
                    for (unsigned int i=0; i < npix;i++ ) 
                    { 
                            //rotation from Equatorial (ra,dec) to local
                            vec3 v =  CRmap.pix2vec(i); 
                            vp.x = rot(0,0)*v.x+rot(0,1)*v.y+rot(0,2)*v.z; 
                            vp.y = rot(1,0)*v.x+rot(1,1)*v.y+rot(1,2)*v.z; 
                            vp.z = rot(2,0)*v.x+rot(2,1)*v.y+rot(2,2)*v.z;
                           
                            unsigned int j = CRmap.vec2pix(vp); // local pixel 
                   
                            // global significance 
                            if (Emap0[j]*norm0[timeidx] > 0.0 && Emap[j]*norm[timeidx] > 0.0) { 
                                    significancemap[i] += 
                                            -2.0*(diffCRmap[i]+CRmap[i])*Emap[j]*norm[timeidx]; 

                                    significancemap[i] += 
                                            +2.0*CRmap[i]*Emap0[j]*norm0[timeidx]; 

                                    double temp1 = norm[timeidx]/norm0[timeidx]*Emap[j]/Emap0[j]; 

                                    significancemap[i] += 
                                            2.0*(*Nmap[timeidx])[j]*log(temp1*(1.0+diffCRmap[i]/CRmap[i])); 
                            }
                    } 
            }

            //write S_a^(n)
            stringstream nameSIGfits;
            nameSIGfits << foldername <<  boost::format("/significance_%s_%d_%d_iteration%02d.fits") % detector % nsideOut % nTimesteps % iteration; 
            if (fs::exists(nameSIGfits.str()) ) {
                    fs::remove(nameSIGfits.str() );
            }
            fitsOut.create(nameSIGfits.str().c_str()); 
            //write_Healpix_map_to_fits(fitsOut, significancemap, PLANCK_FLOAT64);
            write_Healpix_map_to_fits(fitsOut, significancemap, MyDTYPE);
            fitsOut.close(); 

            log_info("Finished iteration " << iteration << " of " << nIterations << "...");

    }
    return 0;
}
