find_path (PARMETIS_INCLUDE_DIR parmetis.h)
find_library (PARMETIS_LIBRARY NAMES parmetis)

# ParMetis requires MPI - find MPI as a dependency
if(NOT Parmetis_FOUND)
    find_package (MPI REQUIRED)
    if(MPI_FOUND)
        if(DEFINED MPI_CXX_INCLUDE_PATH)
            set (PARMETIS_INCLUDE_DIR
                 "${PARMETIS_INCLUDE_DIR};${MPI_CXX_INCLUDE_PATH}"
                 CACHE PATH "ParMetis include directories" FORCE
            )
        endif()
    endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set PFFT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (Parmetis DEFAULT_MSG PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)

mark_as_advanced (PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)
