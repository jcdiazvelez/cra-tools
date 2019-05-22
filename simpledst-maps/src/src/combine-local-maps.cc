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

void
Combine(SkyMap& to, const SkyMap& from) 
{
    unsigned int npix = from.Npix();
    std::cout << "combining " << npix << " pixels" << std::endl;
    for (unsigned int i=0; i < npix; i++) 
    {
         to[i] += from[i];
    }
}


double
Sum(const SkyMap& map, unsigned int npix) 
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
    std::string outputfilename;
    unsigned int nTimesteps;
    unsigned int nIterations;
    unsigned int nSectors;
    int nsideOut;
    double lon;
    double lat;
    std::vector< std::string > input;

    po::options_description desc("Options"); 
    po::positional_options_description p;
    po::positional_options_description positionalOptions;
    //p.add("input", -1);
    po::variables_map vm; 
    
    try { 
       // Define and parse the program options 
       desc.add_options() 
             ("help", "produce help message") 
             ("input", po::value<std::vector <std::string> >(&input)->multitoken(), "Input files") 
             ("outfile", po::value<std::string>(&outputfilename)->default_value("combined_map.fits"), "Combined output file") 
             ("nsideout", po::value<int>(&nsideOut)->default_value(64), "Healpix Nside for output map") 
             ;
        positionalOptions.add("input", 3000);
     
        //po::store(po::parse_command_line(argc, argv, desc),  vm); // can throw 
        po::store(po::command_line_parser(argc, argv).options(desc).positional(positionalOptions).run(),
                                                                                                                                                          vm);
         
        /// --help option 
        if ( vm.count("help")  ) { 
            std::cout << "Basic Command Line Parameter App" 
                      << std::endl << desc << std::endl; 
            return 0; 
        } 
        po::notify(vm); // throws on error, so do after help in case there are any problems
        
        //if (vm.count("input")) 
        { 
            cout << "Input files: " << input[0] << ", " << input[1] << "..." << "\n";
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

    unsigned int npix = 12*nsideOut*nsideOut; 
    SkyMap locMap;
    //locMap.SetNside(nsideOut, RING);     
    read_Healpix_map_from_fits(input[0].c_str(), locMap);
    locMap.fill(0.);

    std::vector< std::string >::iterator it;
    for (it = input.begin();it != input.end(); it++)
    {
        std::cout << *it << std::endl;
        SkyMap tempmap;
        read_Healpix_map_from_fits(it->c_str(), tempmap);
        Combine(locMap,tempmap);
    }

    fitshandle fitsOut;
    if (fs::exists(outputfilename) ) { 
                 fs::remove(outputfilename); 
    }
    std::cout << outputfilename << std::endl;
    fitsOut.create(outputfilename.c_str());
    //write_Healpix_map_to_fits(fitsOut, locMap, FITSUTIL<double>::DTYPE);
    write_Healpix_map_to_fits(fitsOut, locMap, PLANCK_FLOAT64);
    fitsOut.close();
    
}
