#include <iostream>
#include "config.h"


Config::Config(const po::variables_map& vm) {

  // Detector configuration, defaults to newer in-ice format
  std::string config;
  if (vm.count("config") )
    config = vm["config"].as<std::string>();

  if (config=="IT59" || config=="IT73" || config=="IT81") {
    std::cout << "IceTop v1" << std::endl;
    detector = IceTop;
    cfg = ITv1;
  }
  else if (config=="IT81-2012" || config=="IT81-2013") {
    std::cout << "IceTop v2" << std::endl;
    detector = IceTop;
    cfg = ITv2;
  }
  else if (config=="ITpass2") {
    std::cout << "IceTop v3" << std::endl;
    detector = IceTop;
    cfg = ITv3;
  }
  else if (config=="IC86-2011" || config=="IC86-2012" || config=="IC86-2013" ||
      config=="IC86-2014" || config=="IC86-2015")
  {
    std::cout << "IceCube v1" << std::endl;
    detector = IceCube;
    cfg = ICv1;
  }
  else {
    std::cout << "IceCube v2" << std::endl;
    detector = IceCube;
    cfg = ICv2;
  }

  // Time frame (needs language update), defaults to sidereal
  std::string m;
  if (vm.count("method") )
    m = vm["method"].as<std::string>();
  if (m == "ext") {
      method = extsid;
  } else if (m == "anti") {
      method = antisid;
  } else if (m == "solar") {
      method = solar;
  } else {
      method = sidereal;
  }

  // Directional reconstruction, Defaults to ShowerPlane
  std::string r;
  if (vm.count("reco") )
      r = vm["reco"].as<std::string>();
  if (r == "Laputop") {
      reco = Laputop;
  } else if (r == "LaputopSmall") {
      reco = LaputopSmall;
  } else {
      reco = ShowerPlane;
  }

}


bool Config::newConfig() {
  return (cfg == ICv2);
}
