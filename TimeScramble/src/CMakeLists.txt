cmake_minimum_required(VERSION 3.0)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -DBOOST_LOG_DYN_LINK")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

project(time-scramble CXX C)

include(pretty)
include(config)
include(tools)
include(tooldef)
include(CheckIncludeFile)

FIND_PACKAGE(HEALPIX REQUIRED)
FIND_PACKAGE(CFITSIO REQUIRED)
#FIND_PACKAGE(Boost REQUIRED)
#FIND_PACKAGE(Boost COMPONENTS system filesystem program_options REQUIRED)
FIND_PACKAGE(Boost 1.36.0 COMPONENTS system filesystem program_options iostreams REQUIRED)
FIND_PACKAGE(ROOT)
#FIND_PACKAGE(PHOTOSPLINE)
FIND_PACKAGE(photospline REQUIRED)
FIND_PACKAGE(PAL)


include_directories(
                ${PAL_INCLUDE_DIR} 
                ${Boost_INCLUDE_DIR} 
                ${CFITSIO_INCLUDE_DIRS} 
                ${HEALPIX-CXX_INCLUDE_DIRS} 
                ${ROOT_INCLUDE_DIR} include)
set(LIBS ${LIBS} 
    ${PAL_LIBRARIES} ${HEALPIX-CXX_LIBRARIES} 
    ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES} ${ROOT_LIBRARIES} 
    ${PHOTOSPLINE_LIBRARIES} ) 



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
            private/time.cxx
            private/SimpleTrigger.cc
            private/SimpleDST.cc)

add_executable(tsmc
            private/tsmc.cxx
            private/solardipole.cc)


 

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
                ${Boost_LIBRARIES} 
                ${PHOTOSPLINE_LIBRARIES}
                ${ROOT_LIBRARIES})

target_link_libraries(tsmc
                ${PAL_LIBRARIES} 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES} 
                ${PHOTOSPLINE_LIBRARIES} -lboost_log)



if(IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    
       exec_program(test ARGS -h ${CMAKE_BINARY_DIR}/scripts -a -e ${CMAKE_BINARY_DIR}/scripts
            || ln -snf ${CMAKE_SOURCE_DIR}/scripts ${CMAKE_BINARY_DIR}/scripts
            OUTPUT_VARIABLE DEV_NULL)
       install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts
         DESTINATION ${CMAKE_BINARY_DIR}/scripts
         PATTERN ".svn" EXCLUDE
         PATTERN ".git" EXCLUDE
         PATTERN "*.py"
         PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ WORLD_EXECUTE WORLD_READ)
endif (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
