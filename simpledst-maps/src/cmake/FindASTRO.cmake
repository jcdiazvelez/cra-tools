# - Try to find ASTRO.
# Variables used by this module:
#  ASTRO_ROOT_DIR     - ASTRO root directory
# Variables defined by this module:
#  ASTRO_FOUND        - system has ASTRO
#  ASTRO_INCLUDE_DIR  - the ASTRO include directory (cached)
#  ASTRO_INCLUDE_DIRS - the ASTRO include directories
#                         (identical to ASTRO_INCLUDE_DIR)
#  ASTRO_LIBRARY      - the ASTRO library (cached)
#  ASTRO_LIBRARIES    - the ASTRO libraries
#                         (identical to ASTRO_LIBRARY)

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

if(NOT ASTRO_FOUND)

  message("photospline_root_dir" ${ASTRO_ROOT_DIR})
  find_path(ASTRO_INCLUDE_DIR bspline.h
    PATHS ${ASTRO_ROOT_DIR} PATH_SUFFIXES include public/photospline)
  find_library(ASTRO_LIBRARY photospline
    PATHS ${ASTRO_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(ASTRO_INCLUDE_DIR ASTRO_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(ASTRO DEFAULT_MSG
    ASTRO_LIBRARY ASTRO_INCLUDE_DIR)

  set(ASTRO_INCLUDE_DIRS ${ASTRO_INCLUDE_DIR})
  set(ASTRO_LIBRARIES ${ASTRO_LIBRARY})

endif(NOT ASTRO_FOUND)
