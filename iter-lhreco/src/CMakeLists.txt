#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
MESSAGE( STATUS "SOFTWARE_NAME: " ${SOFTWARE_NAME} )
if (SOFTWARE_NAME STREQUAL "AERIE")
   HAWC_PROJECT (iter-lhreco-proj CXX C)

   HAWC_ADD_LIBRARY (iter-lhreco-proj
     HEADERS include/iter-lhreco-proj/*.h
     SOURCES src/*.cc
     USE_PROJECTS hawcnest data-structures rng-service aerie-io
               astro-service detector-service rng-service track-fitter
     USE_PACKAGES CFITSIO HEALPix ROOT Boost)

   HAWC_ADD_EXECUTABLE (reconstruct
     SOURCES src/reconstruct.cc src/pickle.cc
     USE_PROJECTS hawcnest data-structures detector-service 
     USE_PACKAGES CFITSIO HEALPix Boost)

   HAWC_ADD_EXECUTABLE (multi-llh
     SOURCES src/multi-llh.cc src/pickle.cc
     USE_PROJECTS hawcnest data-structures detector-service 
     USE_PACKAGES CFITSIO HEALPix Boost)


   HAWC_ADD_EXECUTABLE (significance
     SOURCES src/significance.cc src/pickle.cc
     USE_PROJECTS hawcnest
     USE_PACKAGES CFITSIO HEALPix Boost)


   HAWC_ADD_EXECUTABLE (getevents
     SOURCES src/getevents.cc 
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

        #foreach(library ${HEALPIX_LIBRARIES})
        #  MESSAGE( STATUS "Boost libraries: " ${library} )
        #endforeach(library $HEALPIX_LIBRARIES})

        include_directories(${Boost_INCLUDE_DIR} ${HEALPIX_INCLUDE_DIRS} ${CFITSIO_INCLUDE_DIRS} include)
        set(LIBS ${LIBS} ${HEALPIX_LIBRARIES} ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES})
        

        add_executable(multi-llh
            src/multi-llh.cc 
            src/pickle.cc
            )
        add_executable(reconstruct
            src/reconstruct.cc 
            src/pickle.cc
            )
        add_executable(significance
            src/significance.cc 
            src/pickle.cc)



        target_link_libraries(multi-llh ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES}) 
        target_link_libraries(reconstruct ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES})
        target_link_libraries(significance ${HEALPIX_LIBRARIES} ${Boost_LIBRARIES})

endif (SOFTWARE_NAME STREQUAL "AERIE")
