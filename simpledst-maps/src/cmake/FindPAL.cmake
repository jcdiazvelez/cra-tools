# Debian packages this as starlink-pal
if (NOT PAL_FOUND)

  message("pal_root_dir" ${PAL_ROOT_DIR})
  find_path(PAL_INCLUDE_DIR star/fitsio.h
    PATHS ${PAL_ROOT_DIR} PATH_SUFFIXES include include/star/pal.h)
  find_library(PAL_LIBRARY starlink_pal
    PATHS ${PAL_ROOT_DIR} PATH_SUFFIXES lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(PAL DEFAULT_MSG
    PAL_LIBRARY PAL_INCLUDE_DIR)

  set(PAL_INCLUDE_DIRS ${CFITSIO_INCLUDE_DIR})
  set(PAL_LIBRARIES ${PAL_LIBRARY})

endif(NOT PAL_FOUND)
