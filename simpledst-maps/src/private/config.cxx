#include "iter-lhreco/config.h"
#include <iostream>
#include <sstream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

using namespace std;

namespace config = boost::property_tree;
typedef  std::shared_ptr<config::ptree> ptree_ptr; 


Config::Config(std::string path)
{
        ReadConfig(path);
}

void Config::ReadConfig(std::string path)
{
  config::ptree cfg; 
  //config::ptree::read_json(path, cfg);
  boost::property_tree::read_json(path,cfg);
  //config::ini_parser::read_ini(path, cfg); 

  std::cout << "Input config parameters: " << std::endl;
  std::cout << "-------------------------------------------------------------" 	<< std::endl;
  try { 
          //loglevel = cfg.get<std::string>("logging.level", "INFO");

          //detectors =  cfg.get_child("detectors");
          config::ptree &detlist = cfg.get_child("detectors");
          for (config::ptree::const_iterator it = detlist.begin(); it != detlist.end(); ++it)
          {
              //std::cout << it->second.get_value<std::string>() << std::endl;
              config::ptree dt = it->second;
              detectors.push_back(dt.get_child("name").get_value<std::string>());
              thetamax.push_back(dt.get_child("thetamax").get_value<double>());
              longitude.push_back(dt.get_child("longitude").get_value<double>());
              latitude.push_back(dt.get_child("latitude").get_value<double>());
              prefix.push_back(dt.get_child("prefix").get_value<std::string>());
              suffix.push_back(dt.get_child("suffix").get_value<std::string>());
          }
          
    } catch(config::ptree_error& e) { 
            std::cerr << e.what() << std::endl;
            exit(1);
    }
}
 
