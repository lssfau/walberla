if (FFTW3_INCLUDE_DIR)
   # Already in cache, be silent
   set (FFTW3_FIND_QUIETLY TRUE)
endif (FFTW3_INCLUDE_DIR)

find_path (FFTW3_INCLUDE_DIR fftw3.h)
find_library (FFTW3_LIBRARIES NAMES fftw3)
find_path (FFTW3_MPI_INCLUDE_DIR fftw3-mpi.h)
find_library (FFTW3_MPI_LIBRARIES NAMES fftw3_mpi)

# handle the QUIETLY and REQUIRED arguments and set FFTW3_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW3 DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)
set(FPHSA_NAME_MISMATCHED 1)
find_package_handle_standard_args (FFTW3_MPI DEFAULT_MSG FFTW3_MPI_LIBRARIES FFTW3_MPI_INCLUDE_DIR)
unset(FPHSA_NAME_MISMATCHED)

mark_as_advanced (FFTW3_LIBRARIES FFTW3_INCLUDE_DIR FFTW3_MPI_LIBRARIES FFTW3_MPI_INCLUDE_DIR)
