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
#FIND_PACKAGE(ROOT)
#FIND_PACKAGE(photospline REQUIRED)
#FIND_PACKAGE(PAL)

include_directories(
                ${Boost_INCLUDE_DIR} 
                ${CFITSIO_INCLUDE_DIRS} 
                ${HEALPIX-CXX_INCLUDE_DIRS} include)



set(LIBS ${LIBS} 
    ${PAL_LIBRARIES} ${HEALPIX-CXX_LIBRARIES} 
    ${CFITSIO_LIBRARIES} ${Boost_LIBRARIES} ) 

add_executable(multi-llh
            private/multi-llh.cc 
            private/pickle.cc
            )

add_executable(illh-reco
            private/illh-reco.cc 
            private/pickle.cc
            private/illh-utils.cc
            )

add_executable(illh-significance
            private/illh-significance.cc 
            private/pickle.cc)



target_link_libraries(multi-llh 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})


target_link_libraries(illh-reco
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})

target_link_libraries(illh-significance 
                ${CFITSIO_LIBRARIES} 
                ${HEALPIX-CXX_LIBRARIES} 
                ${Boost_LIBRARIES})

configure_file(${CMAKE_SOURCE_DIR}/scripts/illh-reco.py.in
           ${CMAKE_BINARY_DIR}/scripts/illh-reco.py
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


