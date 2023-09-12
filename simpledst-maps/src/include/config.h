#ifndef CONFIG_H_INCLUDED
#define CONFIG_H_INCLUDED


/**
 * @namespace Config
 * @brief config definitions 
 *
 */

namespace Config
{

  /**
   * method
   */
  enum method
    { 
        solar,
        sidereal,
        anti_sidereal,
        extende_sidereal
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
