#include <ctime>
#include <sys/types.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept> 
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <random>


#include <boost/filesystem.hpp>
namespace fs = boost::filesystem; 


// linalg
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "Sampling.h"
#include "units.h"
#include "solardipole.h"

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options; 

typedef Healpix_Map<double> SkyMap; 
typedef std::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer

using namespace std;

namespace sampling {

        namespace ublas = boost::numeric::ublas;

        SkyMap TopHatSmooth(SkyMap origMap, double radius, bool average)
        { 
                double theta, phi; 
                unsigned int nside = origMap.Nside();
                unsigned int npix = 12*nside*nside;
                SkyMap newmap; 
                newmap.SetNside(nside, RING);     

                // Loop over all pixels 
                for (int ipix = 0; ipix < npix; ipix++)
                { 
                        // Grab pixels with given angular radius of ith pixel
                        pointing center  = origMap.pix2ang(ipix);

                        std::vector<int> listpix;
                        origMap.query_disc_inclusive(center, radius, listpix);

                        //Sum up pixels in the angular cap; set value of ith pixel to this sum
                        double dsum = 0.;
                        int nsum = 0;
                        vector<int>::iterator jit;
                        for (jit = listpix.begin(); jit != listpix.end(); jit++)
                        {
                                if (origMap[*jit] != Healpix_undef)
                                { 
                                        dsum += origMap[*jit];
                                        ++nsum;
                                }
                        }
                        if (average && nsum>0)
                        { 
                                newmap[ipix] = dsum/nsum; // set to average
                        } else { 
                                newmap[ipix] = dsum;
                        } 
                } 
                return newmap;
        }



        /*
         * loadMap - read local maps and generate vector of maps bined in 
         * siderial time steps
         */
        void loadMap(SkyMap& map, std::string fitspath, bool normalize)
        {
                fitshandle handle;

                std::cout <<"loading healpix fits file " << fitspath << endl ;
                // Read header info and then map
                SkyMap rawmap;
                handle.open(fitspath.c_str());
                handle.goto_hdu(2);
                read_Healpix_map_from_fits(handle, rawmap, 1);
                handle.close();

                //map = TopHatSmooth(rawmap, 10.*phys::degree, true);
                map = rawmap;
                unsigned int nside = map.Nside();

                // add 1e-16 to avoid zeros values
                bool negatives = false;
                double high = Healpix_undef;
                double low = -Healpix_undef;

                if (!normalize) return;
                
                for (int i = 0; i < 12*nside*nside; i++)
                { 
                        if (map[i]  != Healpix_undef) 
                        {
                                high = max(high,map[i]);
                                low  = min(low,map[i]);
                        }
                } 
                for (int i = 0; i < 12*nside*nside; i++)
                {
                        if (map[i]  != Healpix_undef) 
                            map[i] = (map[i]-low)/(high-low);
                }
        }

 
        /*
         * Markov Chain Monte Carlo implementation with Metropolis–Hastings algoritm.
         * Taken from example at:
         *  https://stephens999.github.io/fiveMinuteStats/MH-examples1.html
         */
        void FillRandom(
                SkyMap& outmap, 
                SkyMap& weightsmap, 
                SkyMap& bgmap, 
                const SkyMap& inmap, 
                const double mjd,
                unsigned int numpoints)
        {
          std::random_device rd;
          std::mt19937 generator(rd());
          std::uniform_real_distribution<double> rdist(0,1.0);
          std::uniform_real_distribution<double> pixeldist(-1.0,1.0);

          unsigned int nside = inmap.Nside();
          unsigned int npix = 12*nside*nside;
          double pixelarea = 4.*constants::pi/npix;
          std::cout << "starting fill random " << std::endl;




          for (int ipix=0; ipix<npix;ipix++) 
          {
                  //std::cout << "." << ipix ;
                  //if (inmap[ipix] < 0) 
                  //        continue;

                  double mu = numpoints*inmap[ipix];
                  std::poisson_distribution<> d(mu);
                  //unsigned nevents =d(generator);
                  unsigned nevents = int(mu);
                  //cout<< "numpoints " <<  numpoints << ", pixel count " << inmap[ipix] <<endl;

                  pointing pt = inmap.pix2ang(ipix); 
                  //std::cout << ":"  << nevents << std::endl;
                  

                  //double theta = fmod(pt.theta + pixeldist(generator)*sqrt(pixelarea),constants::pi);
                  double theta = pt.theta;
                  if (theta < 0) theta += constants::pi;
 
                  //double phi = fmod(pt.phi + pixeldist(generator)*sqrt(pixelarea),2*constants::pi); 
                  double phi = pt.phi;
                  if (phi < 0) phi += 2*constants::pi;

                  pointing newpt(theta,phi);
                  int newpix = weightsmap.ang2pix(newpt); 
                  outmap[newpix]+=1.;
                  weightsmap[newpix] += solar_dipole(mjd, phi, 0.5*constants::pi-theta);
                  continue;


                  for (unsigned ievent=0;ievent<nevents;ievent++)
                  {

                        double theta = fmod(pt.theta + pixeldist(generator)*sqrt(pixelarea),constants::pi);
                        if (theta < 0) theta += constants::pi;

                        double phi = fmod(pt.phi + pixeldist(generator)*sqrt(pixelarea),2*constants::pi);
                        if (phi < 0) phi += 2*constants::pi;

                        pointing newpt(theta,phi);
                        int newpix = weightsmap.ang2pix(newpt); 
                        outmap[newpix]+=1.;
                        weightsmap[newpix] += solar_dipole(mjd, phi, 0.5*constants::pi-theta);

                        pointing rand_pt(theta, acos(pixeldist(generator))*2.);
                        int rand_pix = weightsmap.ang2pix(rand_pt); 
                        bgmap[rand_pix] += 1.;
                  }
          }
          std::cout << "finished fill random " << std::endl;
        }
         
        void MCMC(
                SkyMap& outmap, 
                SkyMap& weightsmap, 
                const SkyMap& inmap, 
                const double mjd,
                unsigned int numpoints)
        {
          unsigned int proposed, current;
          unsigned int nside = inmap.Nside();
          unsigned int npix = 12*nside*nside;
          double pixelarea = 4.*constants::pi/npix;

          std::random_device rd;
          std::mt19937 generator(rd());
          //std::default_random_engine generator((unsigned)time(0)*getpid());
          std::uniform_real_distribution<double> rdist(0,1.0);
          std::uniform_real_distribution<double> pixeldist(-1.0,1.0);
          std::uniform_int_distribution<int> idist(0,npix);
               
          // if an event occurs 4 times a minute on average
          //     // how often is it that it occurs n times in one
          //     minute?
          //         std::poisson_distribution<> d(4);

          // If we are in a masked pixel move until we are not
          current = idist(generator); 
          while (inmap[current] <= Healpix_undef) 
                current = idist(generator); 
          //while (!(inmap[current] > 0))
          //      current = idist(generator); 
          for(unsigned int i = 0; i < numpoints; i++)
          { 
                  // If we are in a masked pixel move until we are not
                  proposed = idist(generator); 
                  //while (map[proposed] == Healpix_undef) 
                  //while (!(inmap[proposed] > 0) )
                  //        proposed = idist(generator); 
                  proposed = idist(generator); 

                  float A = inmap[proposed]/inmap[current]; 
                  if(rdist(generator)<A){ 
                          current = proposed; // accept move with probabily min(1,A) 
                  } 
                  pointing pt = inmap.pix2ang(current); 
                  //std::cout << i << ", ";
                  //std::cout << "(" <<pt.theta/units::degree << ", " << pt.phi/units::degree << ")" << std::endl;
                  double theta = fmod(pt.theta + pixeldist(generator)*sqrt(pixelarea),constants::pi);
                  if (theta < 0) theta += constants::pi;
                  double phi = fmod(pt.phi + pixeldist(generator)*sqrt(pixelarea),2*constants::pi);
                  if (phi < 0) phi += 2*constants::pi;
                  //double weight = 1.0/(inmap[current]+1e-16);

                  pointing newpt(theta,phi);
                  int newpix = weightsmap.ang2pix(newpt); 
                  outmap[newpix] += 1.0/numpoints;
                  weightsmap[newpix] += solar_dipole(mjd, phi, 0.5*constants::pi-theta)/numpoints;
                  //std::cout << "sd: "<< solar_dipole(mjd, phi, 0.5*constants::pi-theta) << std::endl;
          }
        }

} // namespace sampling


int main(int argc, char * argv[]) { 

    // initialize random number generator 
    srand((unsigned)time(0)*getpid());


    std::string positionfile, rootname,fitspath;
    std::string maptypestr;
    unsigned int ntraj, nthread, hpdist,nside;
    float mjd;
    bool writemap;

    po::options_description desc("Options"); 
    po::variables_map vm; 
    try { 
       // Define and parse the program options 
       desc.add_options() 
             ("help,h", "produce help message") 
             ("positions,p", po::value<std::string>(&fitspath), "Initial positions file") 
             ("ntraj,n", po::value<unsigned int>(&ntraj)->default_value(1), 
                "Number of trajectories to simulate") 
             ("mjd,m", po::value<float>(&mjd)->default_value(5753), 
                "Modified Julian Day") 
             ("nside,s", po::value<unsigned int>(&nside)->default_value(64), 
                "HEALPix NSide parameter") 
             ("prefix,o", po::value<std::string>(&rootname)->default_value("ouput.fits"), 
                "Output file prefix")
             ;
        po::store(po::command_line_parser(argc, argv).options(desc).run(), vm); 
         
        /// --help option 
        if ( vm.count("help")  ) { 
           std::cout << desc << std::endl; 
           return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        
    } catch(po::error& e) { 
            std::cerr << e.what() << std::endl;
            return 1; 
        
    } catch(std::exception& e) { 
            std::cerr << "Unhandled Exception reached the top of main: " << e.what() << std::endl;
            std::cout <<"loading healpix fits file " << fitspath << std::endl ;
            return 2; 
    } 
 

    /************************************************************************************
    Start of trajectory loop
    ************************************************************************************/
    std::string filename = rootname + "." + std::to_string(ntraj) + ".";

    SkyMap spherical_dist;
    //ifstream posfile;

    std::cout << "opening HealPix file " << fitspath<< std::endl;
    sampling::loadMap(spherical_dist, fitspath,true);

    SkyMap DataMap;
    DataMap.SetNside(nside, RING);
    DataMap.fill(0.);
    SkyMap WeightsMap;
    WeightsMap.SetNside(nside, RING);
    WeightsMap.fill(0.);
    SkyMap BGMap;
    BGMap.SetNside(nside, RING);
    BGMap.fill(0.);

    // populate pool of starting positions
    // sample distribution from healpix map
    //sampling::MCMC(DataMap, WeightsMap, spherical_dist,  mjd, ntraj);
    sampling::FillRandom(DataMap, WeightsMap, BGMap,spherical_dist,  mjd, ntraj);
    
	    // read initial position from file
 
    // Write paths to file

    /*
     * if (io::ends_with(rootname,".fits") || io::ends_with(rootname,".fits.gzip")) 
    {
            io::WriteFitsMap(rootname, npydata);
    }
     */

    fs::remove( rootname); 
    fitshandle fitsOut;
    fitsOut.create(rootname);
    fitsOut.add_comment("Maps: data, weights, bg");
    arr<std::string> colname(3);
    colname[0] = "data map";
    colname[1] = "weights map";
    colname[2] = "background map";


    prepare_Healpix_fitsmap(fitsOut, DataMap, PLANCK_FLOAT64, colname);

    fitsOut.write_column(1, DataMap.Map());
    fitsOut.write_column(2, WeightsMap.Map());
    fitsOut.write_column(3, BGMap.Map());
    fitsOut.close();


    /************************************************************************************
       End of trajectory loop
     ************************************************************************************/
    return 0;
}

