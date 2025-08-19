cmake_minimum_required(VERSION 3.10)  # modern CMake features

project(time-scramble CXX C)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

include(pretty)
include(config)
include(tools)
include(tooldef)
include(CheckIncludeFile)

find_package(HEALPIX REQUIRED)
find_package(CFITSIO REQUIRED)
find_package(Boost 1.36.0 REQUIRED COMPONENTS
    system filesystem program_options iostreams log json)
find_package(photospline REQUIRED)
find_package(PAL REQUIRED)
find_package(OpenMP REQUIRED)

# Define the executable
add_executable(TimeScramble
    private/TimeScramble.cc
    private/solardipole.cc
    private/astro.cxx
    private/time.cxx
    private/Direction.cxx
)

# Attach include paths
target_include_directories(TimeScramble PRIVATE
    ${PAL_INCLUDE_DIR}
    ${Boost_INCLUDE_DIR}
    ${CFITSIO_INCLUDE_DIRS}
    ${HEALPIX-CXX_INCLUDE_DIRS}
    include
)

# Modern CMake: set language standard + definitions
target_compile_features(TimeScramble PRIVATE cxx_std_17)
target_compile_definitions(TimeScramble PRIVATE BOOST_LOG_DYN_LINK)

# Link libraries directly (no LIBS variable needed)
target_link_libraries(TimeScramble PRIVATE
    ${PAL_LIBRARIES}
    ${HEALPIX-CXX_LIBRARIES}
    ${CFITSIO_LIBRARIES}
    ${Boost_LIBRARIES}
    ${PHOTOSPLINE_LIBRARIES}
    photospline
    OpenMP::OpenMP_CXX
)

 
# Configure helper scripts
set(SCRIPTS
    maker.py
    missingFiles.py
    public_maker.py
    purge.py
)

foreach(script ${SCRIPTS})
    configure_file(
        ${CMAKE_SOURCE_DIR}/scripts/${script}.in
        ${CMAKE_BINARY_DIR}/scripts/${script}
        @ONLY
    )
endforeach()

# Install scripts directory, skipping VCS metadata
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/scripts)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/scripts/
        DESTINATION ${CMAKE_INSTALL_BINDIR}/scripts
        PATTERN ".svn" EXCLUDE
        PATTERN ".git" EXCLUDE
        PATTERN "*.py"
        PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ
                    GROUP_EXECUTE GROUP_READ
                    WORLD_EXECUTE WORLD_READ
    )
endif()



