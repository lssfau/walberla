if (PFFT_INCLUDE_DIR)
   # Already in cache, be silent
   set (PFFT_FIND_QUIETLY TRUE)
endif (PFFT_INCLUDE_DIR)

find_path (PFFT_INCLUDE_DIR pfft.h)
find_library (PFFT_LIBRARIES NAMES pfft)

# handle the QUIETLY and REQUIRED arguments and set PFFT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PFFT DEFAULT_MSG PFFT_LIBRARIES PFFT_INCLUDE_DIR)

mark_as_advanced (PFFT_LIBRARIES PFFT_INCLUDE_DIR)
