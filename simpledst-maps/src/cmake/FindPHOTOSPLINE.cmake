# - Try to find PHOTOSPLINE.
# Variables used by this module:
#  PHOTOSPLINE_ROOT_DIR     - PHOTOSPLINE root directory
# Variables defined by this module:
#  PHOTOSPLINE_FOUND        - system has PHOTOSPLINE
#  PHOTOSPLINE_INCLUDE_DIR  - the PHOTOSPLINE include directory (cached)
#  PHOTOSPLINE_INCLUDE_DIRS - the PHOTOSPLINE include directories
#                         (identical to PHOTOSPLINE_INCLUDE_DIR)
#  PHOTOSPLINE_LIBRARY      - the PHOTOSPLINE library (cached)
#  PHOTOSPLINE_LIBRARIES    - the PHOTOSPLINE libraries
#                         (identical to PHOTOSPLINE_LIBRARY)

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

if(NOT PHOTOSPLINE_FOUND)

  message("photospline_root_dir" ${PHOTOSPLINE_ROOT_DIR})
  find_path(PHOTOSPLINE_INCLUDE_DIR bspline.h
    PATHS ${PHOTOSPLINE_ROOT_DIR} PATH_SUFFIXES include public/photospline)
  find_library(PHOTOSPLINE_LIBRARY photospline
    PATHS ${PHOTOSPLINE_ROOT_DIR} PATH_SUFFIXES lib)
  mark_as_advanced(PHOTOSPLINE_INCLUDE_DIR PHOTOSPLINE_LIBRARY)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PHOTOSPLINE DEFAULT_MSG
    PHOTOSPLINE_LIBRARY PHOTOSPLINE_INCLUDE_DIR)

  set(PHOTOSPLINE_INCLUDE_DIRS ${PHOTOSPLINE_INCLUDE_DIR})
  set(PHOTOSPLINE_LIBRARIES ${PHOTOSPLINE_LIBRARY})

endif(NOT PHOTOSPLINE_FOUND)
