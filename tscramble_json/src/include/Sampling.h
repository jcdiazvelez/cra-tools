#ifndef SAMPLING_H
#define SAMPLING_H

#include <memory>
#include <iostream>
#include <fstream>

#ifdef ICECUBE
#include <healpix-cxx/healpix/alm.h>
#include <healpix-cxx/cxxsupport/xcomplex.h>
#include <healpix-cxx/healpix/alm_healpix_tools.h>
#include <healpix-cxx/cxxsupport/fitshandle.h>
#include <healpix-cxx/healpix/healpix_map.h>
#include <healpix-cxx/healpix/healpix_map_fitsio.h>
#define MyDTYPE FITSUTIL<double>::DTYPE 

#else
#include <healpix_cxx/alm.h>
#include <healpix_cxx/xcomplex.h>
#include <healpix_cxx/alm_healpix_tools.h>
#include <healpix_cxx/fitshandle.h>
#include <healpix_cxx/healpix_map.h>
#include <healpix_cxx/healpix_map_fitsio.h>
#define MyDTYPE PLANCK_FLOAT64
#endif //ICECUBE


typedef Healpix_Map<double> SkyMap; 
typedef std::shared_ptr<SkyMap> SkyMapPtr; // Map shared pointer

namespace maptype {
   enum Type { position , momentum };
}



#endif
