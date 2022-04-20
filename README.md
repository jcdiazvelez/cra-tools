# cra-tools

This is a collection of tools used for analyzing the angular distribution of cosmic-ray arrival directions.


The project is divided into three subprojects:

## TimeScramble
Time-scrambling code for IceCube CRA studies. Currently only reads ROOT format SimpleDST files.

## iter-lhreco
C++ implementation of maximum-likelihood technique for reconstructing cosmic-ray anisotropy maps
http://iopscience.iop.org/article/10.3847/0004-637X/823/1/10

## simpledst-maps
C++ code for extracting data from IceCube's reduced data format (simple-dst) in ROOT/HDF5 and generating local HEALpix maps.


## scripts
Collection of scripts for driving production of extraction, reconstruction and analysis of cosmic-ray data.


## Installation


**Prerequisites**:

iter-lhreco dependencies:

* CFITSIO: depends on cfitsio: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html 
* GNU Scientific Library: depends on gsl: https://www.gnu.org/software/gsl/doc/html/index.html
* HEALpix: depends on healpix: http://healpix.sourceforge.net/
* BOOST C++ libraries: https://www.boost.org/


simpledst-maps dependencies:

* CFITSIO: depends on cfitsio: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html 
* GNU Scientific Library: depends on gsl: https://www.gnu.org/software/gsl/doc/html/index.html
* Photospline: extraction depends on photospline to make energy cuts: https://github.com/IceCubeOpenSource/photospline
* Starlink Positional Astronomy Library: extraction depends on pal: https://github.com/Starlink/pal
* CERN ROOT: extraction depends on root: https://root.cern.ch/
* HEALpix: depends on healpix: http://healpix.sourceforge.net/

You will also need:

* [CMake](https://cmake.org) >= 3.1
* A C++11-compliant compiler (e.g. [gcc](https://gcc.gnu.org) >= 4.8.1 or [clang](https://clang.llvm.org) >= 3.3)


**Installation**:

The C++ projects are built using CMake.
Each project has a CMakeList.txt that will detect dependencies and generate a MakeFile. To do this cd into the build directory and excecute the commands:

  cmake ../src;
  make



