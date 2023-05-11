find_path (PARMETIS_INCLUDE_DIR parmetis.h)
find_library (PARMETIS_LIBRARY NAMES parmetis)

# handle the QUIETLY and REQUIRED arguments and set PFFT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Parmetis DEFAULT_MSG PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)

mark_as_advanced (PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)