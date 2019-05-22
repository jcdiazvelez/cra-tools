# - Try to find HEALPIX.
# Variables used by this module:
#  HEALPIX_ROOT_DIR     - HEALPIX root directory
# Variables defined by this module:
#  HEALPIX_FOUND        - system has HEALPIX
#  HEALPIX_INCLUDE_DIR  - the HEALPIX include directory (cached)
#  HEALPIX_INCLUDE_DIRS - the HEALPIX include directories
#                         (identical to HEALPIX_INCLUDE_DIR)
#  HEALPIX_LIBRARY      - the HEALPIX library (cached)
#  HEALPIX_LIBRARIES    - the HEALPIX libraries
#                         (identical to HEALPIX_LIBRARY)

# Copyright (C) 2009
# ASTRON (Netherlands Institute for Radio Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This file is part of the LOFAR software suite.
# The LOFAR software suite is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The LOFAR software suite is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
#
# $Id$

if(NOT HEALPIX_FOUND)

  message("healpix_root_dir" ${HEALPIX_ROOT_DIR})
  find_path( HEALPIX_INCLUDE_DIR
    NAMES healpix_cxx/healpix_base.h healpix_cxx/healpix_map.h
    PATHS ${HEALPIX_ROOT_DIR} PATH_SUFFIXES include)

  find_library(HEALPIX_LIBRARY healpix_cxx
    PATHS ${HEALPIX_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(HEALPIX_INCLUDE_DIR HEALPIX_LIBRARY )

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(HEALPIX DEFAULT_MSG HEALPIX_LIBRARY HEALPIX_INCLUDE_DIR)

  set(HEALPIX_INCLUDE_DIRS ${HEALPIX_INCLUDE_DIR})
  set(HEALPIX_LIBRARIES ${HEALPIX_LIBRARY})

endif(NOT HEALPIX_FOUND)
