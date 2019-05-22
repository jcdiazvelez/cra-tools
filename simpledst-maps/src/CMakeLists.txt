
MESSAGE( STATUS "SOFTWARE_NAME: " ${SOFTWARE_NAME} )
if (SOFTWARE_NAME STREQUAL "AERIE")
   HAWC_PROJECT (iter-lhreco-proj CXX C)

   HAWC_ADD_LIBRARY (iter-lhreco-proj
     HEADERS include/iter-lhreco-proj/*.h
     SOURCES src/*.cc
     USE_PROJECTS hawcnest data-structures rng-service aerie-io
               astro-service detector-service rng-service track-fitter
     USE_PACKAGES CFITSIO HEALPix ROOT Boost pal)

   HAWC_ADD_EXECUTABLE (make-local-maps
     SOURCES src/make-local-maps.cc src/pickle.cc src/esplines.cc
     USE_PROJECTS hawcnest
     USE_PACKAGES CFITSIO HEALPix Boost pal)

   HAWC_ADD_EXECUTABLE (get_azimuth
     SOURCES src/get_azimuth.cc src/pickle.cc src/esplines.cc
     USE_PROJECTS hawcnest
     USE_PACKAGES CFITSIO HEALPix Boost pal)


   HAWC_ADD_EXECUTABLE (make-sidereal-maps
     SOURCES src/make-sidereal-maps.cc src/pickle.cc src/esplines.cc
     USE_PROJECTS hawcnest
     USE_PACKAGES CFITSIO HEALPix Boost pal)

   HAWC_ADD_EXECUTABLE (combine-local-maps
     SOURCES examples/combine-local-maps.cc src/pickle.cc
     USE_PROJECTS hawcnest
     USE_PACKAGES CFITSIO HEALPix Boost)

else (SOFTWARE_NAME STREQUAL "AERIE")
        cmake_minimum_required(VERSION 2.8)

        set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

        project(iter-lhreco-proj CXX C)

        include(CheckIncludeFile)

        FIND_PACKAGE(HEALPIX REQUIRED)
        FIND_PACKAGE(CFITSIO REQUIRED)
        #FIND_PACKAGE(Boost REQUIRED)
        FIND_PACKAGE(Boost COMPONENTS system filesystem program_options REQUIRED)
        FIND_PACKAGE(ROOT)
        FIND_PACKAGE(PHOTOSPLINE)
        #FIND_PACKAGE(PAL)
        #FIND_PACKAGE(ASTRO)

        #foreach(library ${HEALPIX_LIBRARIES})
        #  MESSAGE( STATUS "Boost libraries: " ${library} )
        #endforeach(library $HEALPIX_LIBRARIES})

        #include_directories(${ASTRO} ${PHOTOSPLINE_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${HEALPIX_INCLUDE_DIRS} ${CFITSIO_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR}  include)
        include_directories(${PHOTOSPLINE_INCLUDE_DIR} ${Boost_INCLUDE_DIR} ${HEALPIX_INCLUDE_DIRS} ${CFITSIO_INCLUDE_DIRS} ${ROOT_INCLUDE_DIR} include)
        #set(LIBS ${LIBS} ${HEALPIX_LIBRARIES} ${CFITSIO_LIBRARIES} ${BOOST_LIBRARIES})
        #set(LIBS ${LIBS} ${HEALPIX_LIBRARIES} ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} ${PHOTOSPLINE_LIBRARIES} ${ASTRO_LIBRARIES})
        set(LIBS ${LIBS} ${HEALPIX_LIBRARIES} ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} ${PHOTOSPLINE_LIBRARIES}) 

        add_executable(make-local-maps
            src/make-local-maps.cc 
            src/esplines.cc
            src/pickle.cc
            )

        add_executable(get_azimuth
            src/get_azimuth.cc
            src/esplines.cc
            src/pickle.cc
            )



        add_executable(make-sidereal-maps
            src/make-sidereal-maps.cc 
            src/esplines.cc
            src/pickle.cc
            )
        add_executable(combine-local-maps
            src/combine-local-maps.cc 
            src/pickle.cc)

#        add_executable(time-scramble
#            src/TimeScramble.cpp
#            src/SimpleDST.cc)

        add_executable(make-random-local-maps
            src/make-random-local-maps.cc
            src/esplines.cc
            src/pickle.cc
            )

        add_executable(fitq
            src/esplines.cc
            src/pickle.cc
            src/fitq.cc
            )

        add_executable(fitd
            src/esplines.cc
            src/pickle.cc
            src/fitd.cc
            )

        add_executable(make-mc-hawc-local-maps
            src/make-hawc-mc-map.cc
            src/pickle.cc
            )
 

        target_link_libraries(make-local-maps ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )
        target_link_libraries(get_azimuth ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )
        target_link_libraries(make-mc-hawc-local-maps ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )
        target_link_libraries(make-sidereal-maps ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )
        target_link_libraries(combine-local-maps ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES})
        #target_link_libraries(time-scramble ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ASTRO_LIBRARIES})
        target_link_libraries(make-random-local-maps ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES})
        target_link_libraries(fitq ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} -lMinuit)
        target_link_libraries(fitd ${PHOTOSPLINE_LIBRARIES} ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} -lMinuit)

endif (SOFTWARE_NAME STREQUAL "AERIE")
