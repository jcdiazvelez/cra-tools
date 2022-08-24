#ifndef ILLH_CONFIG_H
#define ILLH_CONFIG_H


#include <vector>
#include <string>


class Config{

        public:

          
             std::vector<std::string> detectors;
             std::vector<std::string> prefix;
             std::vector<std::string> suffix;
             std::vector<double> longitude;
             std::vector<double> latitude;
             std::vector<double> thetamax;

             std::string loglevel;

             Config() {};
             Config(std::string path);
             ~Config() {};
             void ReadConfig(std::string path); 

        //private:

};


#endif
