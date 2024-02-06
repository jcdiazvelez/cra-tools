#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED


/**
 * @namespace Config
 * @brief config definitions 
 *
 */

namespace config
{

  /**
   * method
   */
  enum method
    { 
        solar,
        sidereal,
        antisid,
        extsid
    };

  enum Detector
    { 
        IceCube1,
        IceCube2,
        IceTop,
        Other
    };




}
#endif
