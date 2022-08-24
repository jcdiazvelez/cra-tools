/**
 * multi-llh
 *
 * @version $Id: $
 *
 * @date: $Date: $
 *
 * @author Juan Carlos Diaz-Velez <juan.diazvelez@alumnos.udg.mx>
 *
*/


#include <iter-lhreco/illh-utils.h>
#include <iter-lhreco/config.h>
#include <stdexcept>

namespace fs = boost::filesystem; 
using namespace std;
using namespace boost::numeric::ublas;
using namespace boost::random; 
using boost::format;

typedef Healpix_Map<double> SkyMap; 
typedef boost::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer

void progress_bar(float progress, std::string msg)
{
    int barWidth = 70; 
    std::cout.width(16);
    std::cout << msg << std::left << " ["; 
    int pos = barWidth * progress; 
    for (int i = 0; i < barWidth; ++i) { 
            if (i < pos) std::cout << "="; 
            else if (i == pos) std::cout << ">"; 
            else std::cout << " "; 
    } 
    //std::cout.precision(2);
    std::cout << "] " << boost::format("%02.2f") % (progress * 100.0) << "%\r"; 
    std::cout.flush(); 
}


int main(int argc, char* argv[])
{
    // declare all the variables
    std::string foldername;
    std::vector<std::string> detector;
    unsigned int nTimesteps;
    unsigned int timeidxMin;
    unsigned int timeidxMax;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideIn;
    int nsideOut;

    std::vector<double> lon;
    std::vector<double> lat;

    std::vector<double> thetamax;
    bool randfluct;
    bool iso;
    bool show_progress;
    unsigned int seed;

    std::string cfgstr;

    po::options_description desc("Options"); 
    po::positional_options_description p;
    po::variables_map vm; 
    
    try { 
       // Define and parse the program options 
       desc.add_options() 
             ("help,h", "produce help message") 
             ("cfg", po::value<std::string>(&cfgstr), "detector config file") 
             ("outdir,o", po::value<std::string>(&foldername)->default_value("./sample/"), "Directory of output") 
             ("nsideout", po::value<int>(&nsideOut)->default_value(64), "Healpix Nside for output map") 
             ("timesteps", po::value<unsigned int>(&nTimesteps)->default_value(360), "Number of time steps") 
             ("timestepmin", po::value<unsigned int>(&timeidxMin)->default_value(0), "First time step to use") 
             ("timestepmax", po::value<unsigned int>(&timeidxMax)->default_value(0), "Last time step to use") 
             ("iterations", po::value<unsigned int>(&nIterations)->default_value(20), "Number of iterations") 
             ("progress", po::bool_switch(&show_progress)->default_value(false), "show progress bar")
             ("prefix", po::value<std::string>(), "(override) prefix for input files") 
             ("suffix", po::value<std::string>(), "(override) suffix for input files") 
             ("lon", po::value<double>(), "(override) Longitude of detector") 
             ("lat", po::value<double>(), "(override) Latitude of detector")
             ("thetamax", po::value<double>(), "(override) max theta detector")
			 ("detector", po::value<string>(), "(override) Name of detector")
#if __cplusplus > 199711L
             ("fluctuate,f", po::bool_switch(&randfluct)->default_value(false), "add random fluctuations")
             ("seed", po::value<unsigned int>(&seed)->default_value(123), "RNG seed")
             ("iso", po::bool_switch(&iso)->default_value(false), "make isotropic map")
#endif
             ("sectors", po::value<unsigned int>(&nSectors)->default_value(1), 
                    "Number of sectors (i.e. observatories)");
     
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
         
        /// --help option 
        if ( vm.count("help")  ) { 
            std::cout << "Basic Command Line Parameter App" 
                      << std::endl << desc << std::endl; 
            std::cout << "The (override) parameters will supersede values in cfg for single sector." << std::endl; 
            return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        
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
 
    //*****************************************************************************
    ///// Initialize ///////////////////////////////////////////////////////////// 
    //*****************************************************************************
    cout << cfgstr <<endl;
    Config cfg(cfgstr);

    unsigned int npix = 12*nsideOut*nsideOut; 

    // Import data : n_tau_i
    std::vector< std::vector<SkyMapPtr> > Nmap;
    std::vector< std::string > prefix;
    std::vector< std::string > suffix;


    for (unsigned int isect=0; isect<nSectors;isect++) 
    { 
        if (isect+1 > cfg.detectors.size() )
            throw boost::program_options::error("Could not find additional detector configurations!!!");

        detector.push_back(cfg.detectors[isect]);
        log_info("Reading config for " << detector.back() << ": ");

        lat.push_back(cfg.latitude[isect]);
        log_info("  latitude = " << lat.back());

        lon.push_back(cfg.longitude[isect]);
        log_info("  longitude = " << lon.back());

        thetamax.push_back(cfg.thetamax[isect]);
        log_info("  thetamax = " << thetamax.back());

        prefix.push_back(cfg.prefix[isect]);
        log_info("  prefix = " << prefix.back());

        suffix.push_back(cfg.suffix[isect]);
        log_info("  suffix = " << suffix.back());
    }
    if ( vm.count("detector")  )
        detector[0] = vm["detector"].as<std::string>();
    if ( vm.count("prefix")  )
        prefix[0] = vm["prefix"].as<std::string>();
    if ( vm.count("suffix")  )
        suffix[0] = vm["suffix"].as<std::string>();
    if ( vm.count("lon")  )
        lon[0] = vm["lon"].as<double>();
    if ( vm.count("lat")  )
        lat[0] = vm["lat"].as<double>();
    if ( vm.count("thetamax")  )
        thetamax[0] = vm["thetamax"].as<double>();

    //*****************************************************************************
    ////// Read input maps //////////////////////////////////////////////////////// 
    //*****************************************************************************
    for (unsigned int isect=0; isect<nSectors;isect++) 
    {
        //std::vector<SkyMapPtr> tmpNmap;
        Nmap.push_back(std::vector<SkyMapPtr>());
        illh::loadMap( Nmap.back(), timeidxMin, timeidxMax, 
                nTimesteps, nsideIn, nsideOut, prefix[isect] , suffix[isect]);
        log_info("Read localmaps for " << detector[isect]);
    }

    if (randfluct && iso)
            log_fatal("ranfluct and iso are mutually exclusive!!!");

    if (randfluct) 
    { 
#if __cplusplus > 199711L
            boost::mt19937 rng(seed);
            for (unsigned int isect=0; isect<nSectors;isect++) 
            { 
                illh::fluctuate(Nmap[isect],rng); 
            }
#else
            log_info("isotropic function disabled");
#endif
    }
    else if (iso) 
    { 
#if __cplusplus > 199711L
            boost::mt19937 rng(seed);
            for (unsigned int isect=0; isect<nSectors;isect++) 
            {
                illh::isotropic(Nmap[isect],rng);
            }
#else 
            log_info("isotropic function disabled");
#endif
    }

    log_info("Creating output");
    // output directory
    fs::path dir(foldername); 
    if(!(fs::exists(dir)) ) { 
            log_info("Directory " << foldername << " doesn't exist");
            if (fs::create_directory(dir)) 
                    log_info("....successfully created !");
    }
    
    // detector position
    
    for (unsigned int isect=0; isect<nSectors;isect++) 
    { 
        lon[isect] = lon[isect]/180.*M_PI;
        lat[isect] = lat[isect]/180.*M_PI;
        log_info("Sector "<< isect << " ("<<detector[isect]<<"): Latitude=" << lat[isect] 
                << ",  Longitude=" << lon[isect]);
        thetamax[isect] = thetamax[isect]/180.*M_PI;
    }
    
    // normalization of isotropic flux
    double isovalue = 1.0;

    // initial CR anisotropy : I_a^(0) = 1.
    SkyMap CRmap;
    CRmap.SetNside(nsideOut, RING);
    CRmap.fill( isovalue );
                
    // initial exposure : A_i^(0)
    std::vector<SkyMapPtr> Emap0;
    std::vector<SkyMapPtr> Emap;

    //// window function of FOV map : F_a
    matrix<bool> FOV (nSectors, npix);

    for (unsigned int isect=0; isect<nSectors;isect++) 
    {
        SkyMapPtr Emap0i(new SkyMap);
        Emap0i->SetNside(nsideOut, RING);
        Emap0i->fill(0.);
        Emap0.push_back(Emap0i);

        // n^th exposure : A_i^(n)
        SkyMapPtr Emapi(new SkyMap);
        Emapi->SetNside(nsideOut, RING);
        Emapi->fill(0.);
        Emap.push_back(Emapi);

        for (unsigned i=0;i<npix;++i) 
        {
            FOV(isect,i) = 1;
            pointing pt = Emap[isect]->pix2ang(i);
            if (pt.theta > thetamax[isect])
               FOV(isect,i) = 0;
        }
    }

    // Calculate initial exposure :	
    std::vector<double> totsector(nSectors,0.0);

    for (unsigned int isect=0; isect<nSectors;isect++) 
    { 
        for (unsigned int i=0; i < npix;i++ ) 
        { 
            double nEvents = 0.0; 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    nEvents    += (*Nmap[isect][timeidx])[i]; 
            } 
            if (FOV(isect,i)) {
                    (*Emap0[isect])[i] = nEvents; 
                    totsector[isect] += nEvents;
            }
        }
    }
            
    double tmptot(0);
    for (unsigned int isect=0; isect<nSectors;isect++) 
    { 
        for (unsigned int i=0; i < npix;i++ ) 
        { 
            if (FOV(isect,i)) { 
                    (*Emap0[isect])[i] = (*Emap0[isect])[i]/totsector[isect]; 
                    (*Emap[isect])[i] = (*Emap0[isect])[i]; 
            }
        }
        log_info("Total events ("<<detector[isect]<<"): " << totsector[isect]);
        tmptot += totsector[isect];
    }
    log_info("Total events: " << tmptot);

    // initial normalization : N_tau^(0)
    std::vector< std::vector<double> > norm0;
    // n^th normalization : N_tau^(n)
    std::vector< std::vector<double> > norm;

    log_info("Initializing normalization...");
    for (unsigned int isect=0; isect<nSectors;isect++) 
    { 
        norm0.push_back(std::vector<double>(nTimesteps,0.));
        norm.push_back(std::vector<double>(nTimesteps,0.));

        for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
        { 
             double nEvents = 0.0;
             for (unsigned int i=0; i < npix;i++ ) 
             { 
                 if (FOV(isect,i)) 
                        nEvents += (*Nmap[isect][timeidx])[i];
             } 
             norm[isect][timeidx] = nEvents/isovalue; 
             norm0[isect][timeidx] = norm[isect][timeidx]; 
        }
    }
        
    // data in equatorial coords : n_a
    SkyMap dataMap;
    dataMap.SetNside(nsideOut, RING); 
    dataMap.fill(0.);

    SkyMap bkgMap;
    bkgMap.SetNside(nsideOut, RING); 
    bkgMap.fill(0.);


    // rotate and combine maps
    log_info("rotate and combine maps");
    for (unsigned int i=0; i < npix;i++ ) 
    { 
        double nBkg = 0.; 
        double nEvents = 0.; 
            
        for (unsigned int isect=0; isect<nSectors;isect++) 
        { 
            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    int j = illh::loc2eq_idx(i, timeidx, lat[isect], lon[isect], nTimesteps, CRmap);
                    if ((*Emap0[isect])[j] > 0.0 ){ 
                            nEvents += (*Nmap[isect][timeidx])[j];
                            nBkg    += norm[isect][timeidx]*(*Emap0[isect])[j]; 
                    }
            } 
        } 
        dataMap[i] = nEvents;
    } 
    
    //  write N_tau^(0), write A_i^(0)
    log_info("Writting initial exposure");
    for (unsigned int isect=0; isect<nSectors;isect++) 
        illh::save_iter(foldername, norm[isect], *Emap0[isect], detector[isect], nsideOut, nTimesteps, 0);


    // write n_a
    stringstream detectors;
    detectors <<  detector[0]; 
    for (unsigned int isect=1; isect<nSectors;isect++) 
        detectors <<  "-" << detector[isect]; 

    //*****************************************************************************
    ////// Iterate //////////////////////////////////////////////////////////////// 
    //*****************************************************************************

    // n^th differential CR anisotropy : I_a^(n) - isovalue
    SkyMap diffCRmap;
    diffCRmap.SetNside(nsideOut, RING);
    diffCRmap.fill(0.);

    for (unsigned int iteration = 1; iteration <= nIterations; iteration++)
    { 
            stringstream iterstr;
            iterstr << "Iter " << iteration;
            if (!show_progress) 
                log_info(iterstr.str());

            // calculate new CR anisotropy
            for (unsigned int i=0; i<npix;i++) 
            { 
                    if ( show_progress && (i % 100 == 0) ) 
                        progress_bar(i/float(npix), iterstr.str());

                    double nEvents = 0.0;
                    double nBkg = 0.0;

                    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                    { 
                        for (unsigned int isect=0; isect<nSectors;isect++) 
                        {
                            int j = illh::loc2eq_idx(i, timeidx, lat[isect], lon[isect], nTimesteps, CRmap);
                            if (FOV(isect,j) && ((*Emap0[isect])[j] > 0.0)) { 
                                    nEvents += (*Nmap[isect][timeidx])[j];
                                    nBkg += norm[isect][timeidx]*(*Emap[isect])[j]; 
                            } 
                        }
                    }
                        
                    if (nBkg > 0.0) { 
                            diffCRmap[i] = nEvents/nBkg-CRmap[i]; 
                            bkgMap[i] = nBkg;
                    }
            } 
            if (show_progress) {
                progress_bar(1.0, iterstr.str()); 
                std::cout << std::endl;
            }

            // remove m=0 multipole moments :
            const int LMAX=180;

            // Initialize the spherical harmonic coefficients of the map
            Alm<xcomplex<double> > alm( LMAX, LMAX);
            arr<double> weight(2*nsideOut);
            weight.fill(1.);
            
	        // generate alm coefficients with terative map2alm method
            const int numIter = 3;
            map2alm_iter( diffCRmap, alm, numIter, weight);
            for (unsigned int l=0; l < LMAX+1; l++)
            {    
                alm(l,0) = 0.;
            }

            // generate diffCRmap from alm coefficients
            alm2map( alm, diffCRmap);

            stringstream acc_iterstr;
            acc_iterstr << "Acceptance " << iteration;
            if (!show_progress) 
                log_info(acc_iterstr.str());



            float total_steps = 2.0*nSectors*(timeidxMax-timeidxMin)*npix;
            int step = 0;
            for (unsigned int isect=0; isect<nSectors;isect++) 
            { 
                // calculate new normalization 
                for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                { 

                    double nEvents = 0.0;
                    double nBkg = 0.0; 

                    // Integrate over all (rotated) pixels 
                    for (unsigned int i=0; i<npix; i++) 
                    { 
                        if ( show_progress && (step % 10000 == 0) ) 
                            progress_bar(step/total_steps, acc_iterstr.str());

                        int j = illh::eq2loc_idx(i, timeidx, lat[isect], lon[isect], nTimesteps, CRmap);
                        if (FOV(isect,i) && ((*Emap0[isect])[i] > 0.0)) { 
                                nEvents += (*Nmap[isect][timeidx])[i];
                                nBkg += (*Emap[isect])[i]*(CRmap[j]+diffCRmap[j]);
                        } 
                        step++;

                    } 
                    if (nBkg > 0.0) 
                        norm[isect][timeidx] = nEvents/nBkg; 
                } 
            } 

            // caluculate new acceptance
            for (unsigned int isect=0; isect<nSectors;isect++) 
            { 
                for (unsigned int i=0; i<npix; i++) 
                {
                    double nEvents = 0.0;
                    double nBkg = 0.0; 
            
                    bool skip=true;
                    for (unsigned int jsect=0; jsect<nSectors;jsect++) 
                        skip = skip && !((*Emap0[jsect])[i] > 0.0);

                    if ( skip ) 
                        continue;

                    for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                    {
                        if ( show_progress && (step % 10000 == 0) ) 
                            progress_bar(step/total_steps, acc_iterstr.str());

                        int j = illh::eq2loc_idx(i, timeidx, lat[isect], lon[isect], nTimesteps, CRmap);
                        nEvents += (*Nmap[isect][timeidx])[i];
                        nBkg += norm[isect][timeidx]*(CRmap[j]+diffCRmap[j]);
                        step++;
                    } 
                    if (((*Emap0[isect])[i] > 0.0) && (nBkg > 0.0)) 
                    { 
                        (*Emap[isect])[i] = nEvents/nBkg; 
                    } else { 
                        (*Emap[isect])[i] = 0.0; 
                    } 
                }
            
                // Renormalize
                double sumEmap = illh::Sum(*Emap[isect], npix);
            
                for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
                { 
                    norm[isect][timeidx] = norm[isect][timeidx]*sumEmap; 
                } 
                for (unsigned int i=0; i < npix;i++ ) 
                { 
                    (*Emap[isect])[i] = (*Emap[isect])[i]/sumEmap; 
                } 
            }
            if (show_progress) {
                progress_bar(1.0, acc_iterstr.str()); 
                std::cout << std::endl;
            }
            
                
            SkyMap diffCRmapNormed; 
            diffCRmapNormed.SetNside(nsideOut, RING); 
            diffCRmapNormed.fill(0.);

            for (unsigned int i=0; i < npix;i++ ) 
            { 
                    diffCRmapNormed[i] = diffCRmap[i]/isovalue; 
            } 
            
            // write relative intensity map I_a^(n)
            illh::write_maps(foldername,"CR", diffCRmapNormed, dataMap, bkgMap,
                detectors.str(), nsideOut, nTimesteps, iteration);

       
            for (unsigned int isect=0; isect<nSectors;isect++) 
            { 
                // write N_tau^(n) normalizaton to file
                illh::save_iter(foldername, norm[isect], *Emap[isect], detector[isect], 
                        nsideOut, nTimesteps, iteration);
            }

            // calculate statistical significance : S_a^(n)
            SkyMap significancemap; 
            significancemap.SetNside(nsideOut, RING);
            significancemap.fill(0.);
            
            stringstream sig_iterstr;
            sig_iterstr << "Significance " << iteration;
            if (!show_progress) 
                log_info(iterstr.str());

            for (unsigned int timeidx=timeidxMin; timeidx < timeidxMax;timeidx++ ) 
            { 
                    if ( show_progress ) 
                        progress_bar(timeidx/float(timeidxMax), sig_iterstr.str());

                    for (unsigned int i=0; i < npix;i++ ) 
                    { 
                        //rotation from local to Equatorial (ra,dec) 
                        for (unsigned int isect=0; isect<nSectors;isect++) 
                        { 
                            int j = illh::loc2eq_idx(i, timeidx, lat[isect], lon[isect], nTimesteps, CRmap);

                            // global significance 
                            if ((*Emap0[isect])[j]*norm0[isect][timeidx] > 0.0 && (*Emap[isect])[j]*norm[isect][timeidx] > 0.0) { 
                                 significancemap[i] += 
                                     -2.0*(diffCRmap[i]+CRmap[i])*(*Emap[isect])[j]*norm[isect][timeidx]; 
                                 significancemap[i] += 
                                     +2.0*CRmap[i]*(*Emap0[isect])[j]*norm0[isect][timeidx]; 
                                 double temp1 = (*Emap[isect])[j]/(*Emap0[isect])[j]*norm[isect][timeidx]/norm0[isect][timeidx];
                                 significancemap[i] += 
                                     2.0*(*Nmap[isect][timeidx])[j]*log(temp1*(1.0+diffCRmap[i]/CRmap[i])); 
                            }
                        }
                    } 
            }
            progress_bar(1.0, sig_iterstr.str());
            std::cout << std::endl;

            //write S_a^(n)
            stringstream nameSIGfits;
            nameSIGfits << foldername;
            illh::write_maps(foldername,"significance", significancemap, dataMap, bkgMap,
                detectors.str(), nsideOut, nTimesteps, iteration);
            log_info("Finished iteration " << iteration << " of " << nIterations << "...");

    }

    return 0;
}
