# cra-tools

This is a collection of tools used for analyzing cosmic-ray arrival distributions.


The project is divided into three subprojects:

## iter-lhreco
C++ implementation of maximum-likelihood technique for reconstructing cosmic-ray anisotropy
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


simpledst-maps dependencies:

* CFITSIO: depends on cfitsio: https://heasarc.gsfc.nasa.gov/fitsio/fitsio.html 
* GNU Scientific Library: depends on gsl: https://www.gnu.org/software/gsl/doc/html/index.html
* Photospline: extraction depends on photospline to make energy cuts: https://github.com/cnweaver/photospline
* Starlink Positional Astronomy Library: extraction depends on pal: https://github.com/Starlink/pal
* CERN ROOT: extraction depends on root: https://root.cern.ch/
* HEALpix: depends on healpix: http://healpix.sourceforge.net/

scripts:



**Installation**:





