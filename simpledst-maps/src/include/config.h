#include <boost/program_options.hpp>

namespace po = boost::program_options;

class Config {

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

      enum Method {
        sidereal = 10,
        antisid = 11,
        extsid = 12,
        solar = 13
      };

      enum Reco {
        ShowerPlane = 20,
        Laputop = 21,
        LaputopSmall = 22
      };

      int32_t detector;
      int32_t cfg;
      int32_t method;
      int32_t reco;

      Config(po::variables_map vm);

      bool newConfig();
};

