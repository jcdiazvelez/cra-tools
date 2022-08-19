#ifndef ESPLINES
#define ESPLINES
#include <vector>
#if __cplusplus > 199711L
#include <photospline/splinetable.h>
#include <photospline/bspline.h>
#include <boost/program_options.hpp>
#include <SimpleDST.h>
#include <string>

namespace po = boost::program_options;

class Config
{
     public:

      enum Detector { 
        IceCube = 1,
        IceTop = 2
      };

      enum Configuration { 
        ICv1 = 5, // < 2016
        ICv2 = 6, // >= 2016+
        ITv1 = 7, // IT59,IT73,IT81
        ITv2 = 8, // IT81-2012,IT81-2013
        ITv3 = 9 // ITpass2
      };

      enum Filter { 
        None = 0, 
        STA3 = 11, 
        STA8 = 12, 
        NotSTA8 = 13, 
      };

      enum Method { 
        sidereal = 20, 
        antisid = 21, 
        extsid = 22, 
        solar = 23
      };


      int32_t detector;
      int32_t cfg;
      int32_t filter;
      int32_t method;

      Config(Detector dt = IceCube, Configuration c = ICv2, Filter f = None, Method = sidereal);

      //Config(std::string config, std::string f = "");
      Config(po::variables_map vm);

      bool newConfig();
};


/**
 *	Spline-based energy cuts for In-Ice SimpleDST
 *
 */
bool ICenergyCut(SimpleDST dst, photospline::splinetable<> &spline, double zenith, double emin, double emax);
int ICenergyCut(SimpleDST dst, photospline::splinetable<> &table, double zenith, std::vector<float> ebins);
#endif // __cplusplus


/**
 *	IceTop energy cut
 *
 */
bool ITenergyCut(SimpleDST dst, double emin, double emax);
int ITenergyCut(SimpleDST dst, std::vector<float> ebins);

/**
 *	IceTop s125 cut
 *
 */
bool ITs125Cut(SimpleDST dst, double smin, double smax);
int ITs125Cut(SimpleDST dst, std::vector<float> sbins);

/**
 *  GB IceTop stations cut --> low and high
 *
 */
bool ITNstatCut(SimpleDST dst, double smin, double smax);
//int ITNstatCut(SimpleDST dst, std::vector<float> sbins) {
//

bool filterCut(Config cfg, SimpleDST dst);

#endif // ESPLINES
