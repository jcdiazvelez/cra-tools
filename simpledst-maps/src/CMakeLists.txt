cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -DBOOST_LOG_DYN_LINK")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

project(iter-lhreco-proj CXX C)

include(pretty)
include(config)
include(tools)
include(tooldef)
include(CheckIncludeFile)

FIND_PACKAGE(HEALPIX REQUIRED)
FIND_PACKAGE(CFITSIO REQUIRED)
FIND_PACKAGE(Boost 1.36.0 COMPONENTS system filesystem program_options iostreams REQUIRED)
FIND_PACKAGE(ROOT)
FIND_PACKAGE(photospline REQUIRED)
FIND_PACKAGE(PAL)
FIND_PACKAGE(GSL)

include_directories(
                ${PAL_INCLUDE_DIR} 
                ${Boost_INCLUDE_DIR} 
                ${CFITSIO_INCLUDE_DIRS} 
                ${HEALPIX-CXX_INCLUDE_DIRS} 
                ${GSL_INCLUDE_DIRS} 
                ${ROOT_INCLUDE_DIR} include)


set(LIBS ${LIBS} 
    ${PAL_LIBRARIES} ${HEALPIX-CXX_LIBRARIES} 
    ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} 
    ${PHOTOSPLINE_LIBRARIES} ) 

add_executable(dstcut
    private/dstcut.cc
    private/TDST.cxx
    TDSTDict.cxx
    private/Position.cxx
    private/Direction.cxx
    private/astro.cxx
    )

add_executable(splinechk
    private/splinechk.cc
    )



add_executable(make-local-maps
    private/make-local-maps.cc 
    private/cuts.cc
    private/SimpleDST.cc
    private/SimpleTrigger.cc
    private/pickle.cc
    private/iter-lhreco/illh-utils.cc
    private/solardipole.cc
    private/Direction.cxx
    private/astro.cxx
    )

add_executable(make-energy-histos
    private/make-energy-histos.cc 
    private/util/cnpy.cc
    private/histogram.cc
    private/cuts.cc
    private/SimpleDST.cc
    private/SimpleTrigger.cc
    private/pickle.cc
    private/iter-lhreco/illh-utils.cc
    private/solardipole.cc
    private/Direction.cxx
    private/astro.cxx
    )



add_executable(get_azimuth
    private/get_azimuth.cc
    private/cuts.cc
    private/SimpleDST.cc
    private/pickle.cc
    )



add_executable(make-sidereal-maps
    private/make-sidereal-maps.cc 
    private/cuts.cc
    private/SimpleDST.cc
    private/SimpleTrigger.cc
    private/pickle.cc
    private/Direction.cxx
    private/astro.cxx
    )
add_executable(combine-local-maps
    private/combine-local-maps.cc 
    private/pickle.cc)


add_executable(make-random-local-maps
    private/make-random-local-maps.cc
    private/cuts.cc
    private/pickle.cc
    private/SimpleDST.cc
    private/SimpleTrigger.cc
    private/pickle.cc
    private/solardipole.cc
    private/Direction.cxx
    private/astro.cxx
    )

add_executable(fitq
    private/pickle.cc
    private/fitq.cc
    )

add_executable(fitd
    private/pickle.cc
    private/fitd.cc
    )

add_executable(make-mc-hawc-local-maps
    private/make-hawc-mc-map.cc
    private/pickle.cc
    )

add_executable(TimeScramble
            private/TimeScramble.cc
            private/solardipole.cc
            private/SimpleTrigger.cc
            private/astro.cxx
            private/time.cxx
            private/Direction.cxx
            private/SimpleDST.cc)

add_executable(checksd
            private/checksd.cc
            private/solardipole.cc
            private/astro.cxx
            private/time.cxx)

add_executable(tsmc
            private/tsmc.cxx
            private/solardipole.cc)


target_link_libraries(dstcut
                ${PAL_LIBRARIES} 
                ${Boost_LIBRARIES} 
                ${GSL_LIBRARIES} 
                ${ROOT_LIBRARIES})

target_link_libraries(splinechk
                ${CFITSIO_LIBRARIES} 
                ${Boost_LIBRARIES} 
                ${PHOTOSPLINE_LIBRARIES} photospline
                )



target_link_libraries(TimeScramble
                ${PAL_LIBRARIES} 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES} 
                ${PHOTOSPLINE_LIBRARIES} photospline
                ${ROOT_LIBRARIES})

target_link_libraries(checksd 
                ${PAL_LIBRARIES} 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})

target_link_libraries(tsmc
                ${PAL_LIBRARIES} 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES} 
                ${PHOTOSPLINE_LIBRARIES} -lboost_log)


target_link_libraries(make-local-maps 
    ${PHOTOSPLINE_LIBRARIES} photospline 
    ${HEALPIX-CXX_LIBRARIES} 
    ${Boost_LIBRARIES} 
    ${ROOT_LIBRARIES} 
    ${PAL_LIBRARIES})

target_link_libraries(make-energy-histos
    ${PHOTOSPLINE_LIBRARIES} photospline 
    ${HEALPIX-CXX_LIBRARIES} 
    ${Boost_LIBRARIES} 
    ${ROOT_LIBRARIES} 
    ${PAL_LIBRARIES})



add_executable(multi-llh
            private/iter-lhreco/multi-llh.cc 
            private/iter-lhreco/pickle.cc
            )

add_executable(illh-reco
            private/iter-lhreco/illh-reco.cc 
            private/iter-lhreco/pickle.cc
            private/iter-lhreco/illh-utils.cc
            )



target_link_libraries(get_azimuth ${PHOTOSPLINE_LIBRARIES} photospline ${HEALPIX-CXX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

target_link_libraries(make-mc-hawc-local-maps ${PHOTOSPLINE_LIBRARIES} photospline ${HEALPIX-CXX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} )

target_link_libraries(make-sidereal-maps 
    ${PHOTOSPLINE_LIBRARIES} photospline 
    ${HEALPIX-CXX_LIBRARIES} 
    ${Boost_LIBRARIES} 
    ${ROOT_LIBRARIES} 
    ${PAL_LIBRARIES})

target_link_libraries(combine-local-maps ${HEALPIX-CXX_LIBRARIES} ${Boost_LIBRARIES})

target_link_libraries(make-random-local-maps 
    ${PHOTOSPLINE_LIBRARIES} photospline 
    ${HEALPIX-CXX_LIBRARIES} 
    ${Boost_LIBRARIES} 
    ${ROOT_LIBRARIES}
    ${PAL_LIBRARIES})

target_link_libraries(fitq ${PHOTOSPLINE_LIBRARIES} photospline ${HEALPIX-CXX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} -lMinuit)
target_link_libraries(fitd ${PHOTOSPLINE_LIBRARIES} photospline ${HEALPIX-CXX_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} -lMinuit)

target_link_libraries(multi-llh 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})


target_link_libraries(illh-reco
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/TDSTDict.cxx
  COMMAND rootcint -f ${CMAKE_CURRENT_SOURCE_DIR}/TDSTDict.cxx 
                ${CMAKE_CURRENT_SOURCE_DIR}/include/TDST.h 
                ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/TDST.h 
          ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h
  VERBATIM)
#add_library(myLib out.c)



configure_file(${CMAKE_SOURCE_DIR}/scripts/maker.py.in
           ${CMAKE_BINARY_DIR}/scripts/maker.py
           @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/scripts/missingFiles.py.in
           ${CMAKE_BINARY_DIR}/scripts/missingFiles.py
           @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/scripts/public_maker.py.in
           ${CMAKE_BINARY_DIR}/scripts/public_maker.py
           @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/scripts/purge.py.in
           ${CMAKE_BINARY_DIR}/scripts/purge.py
           @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/scripts/submit_all.py.in
           ${CMAKE_BINARY_DIR}/scripts/submit_all.py
           @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/scripts/make-local-maps.py.in
           ${CMAKE_BINARY_DIR}/scripts/make-local-maps.py
           @ONLY)

configure_file(${CMAKE_SOURCE_DIR}/scripts/combine-local-maps.py.in
           ${CMAKE_BINARY_DIR}/scripts/combine-local-maps.py
           @ONLY)

if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
       exec_program(test ARGS -h ${CMAKE_BINARY_DIR}/scripts/npx4 -a -e ${CMAKE_BINARY_DIR}/scripts/npxi4
            || ln -snf ${CMAKE_SOURCE_DIR}/scripts/npx4 ${CMAKE_BINARY_DIR}/scripts/npx4
            OUTPUT_VARIABLE DEV_NULL)
       install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
         DESTINATION ${CMAKE_BINARY_DIR}/scripts
         PATTERN ".svn" EXCLUDE
         PATTERN ".git" EXCLUDE
         PATTERN "*.py"
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)


#
#  Environment checking targets
#
configure_file(${CMAKE_SOURCE_DIR}/cmake/env-check.sh.in
  ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/env-check.sh
  @ONLY)

#
# env-shell.sh
#
foreach (file env-shell.sh)
  message(STATUS "Generating ${file}")
  configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/${file}.in
    ${CMAKE_CURRENT_BINARY_DIR}/${file}
    @ONLY
    )
  execute_process(COMMAND chmod 755 ${CMAKE_BINARY_DIR}/${file})
  execute_process(COMMAND cp ${CMAKE_BINARY_DIR}/${file} ${NOTES_DIR})
endforeach (file env-shell.sh)
