//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file MPIWrapper.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Feichtinger
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"



/// \cond internal

// CMake generated header, needed for WALBERLA_BUILD_WITH_MPI
#include "waLBerlaDefinitions.h"



// MPI SECTION //

#ifdef WALBERLA_BUILD_WITH_MPI

#define     WALBERLA_MPI_SECTION() if(true)
#define WALBERLA_NON_MPI_SECTION() if(false)

#else

#define     WALBERLA_MPI_SECTION() if(false)
#define WALBERLA_NON_MPI_SECTION() if(true)

#endif

namespace walberla {
namespace mpistubs {
    //empty namespace which can be used
} // namespace mpistubs
} // namespace walberla

#ifdef WALBERLA_BUILD_WITH_MPI



#ifdef _MSC_VER
#pragma warning ( push, 1 )
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wsign-conversion"
#endif
#include <mpi.h>
#if defined(OPEN_MPI) && OPEN_MPI
#include <mpi-ext.h>
#endif
#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning ( pop )
#endif



#else // WALBERLA_BUILD_WITH_MPI



// dummy definitions for MPI data types and MPI functions in order to guarantee successful compilation without MPI enabled

#include <stdexcept>

namespace walberla {

namespace mpistubs {


using MPI_Comm = int;
using MPI_Datatype = int;
using MPI_Group = int;
using MPI_Op = int;
using MPI_Request = int;
using MPI_File = int;
using MPI_Offset = int;
using MPI_Info = int;
using MPI_Aint = int;
using MPI_User_function = void (void *, void *, int *, MPI_Datatype *);

struct MPI_Status
{
   explicit inline MPI_Status()
      : count     ( 0 )
      , MPI_SOURCE( 0 )
      , MPI_TAG   ( 0 )
      , MPI_ERROR ( 0 )
   {}

   int count;
   int MPI_SOURCE;
   int MPI_TAG;
   int MPI_ERROR;
};



const int MPI_COMM_NULL  = 0;
const int MPI_COMM_WORLD = 1;

const int MPI_COMM_TYPE_SHARED = 0;

const int MPI_SUCCESS = 1;



const MPI_Datatype MPI_BYTE               =  0;
const MPI_Datatype MPI_CHAR               =  1;
const MPI_Datatype MPI_SHORT              =  2;
const MPI_Datatype MPI_INT                =  3;
const MPI_Datatype MPI_LONG               =  4;
const MPI_Datatype MPI_LONG_LONG          =  5;
const MPI_Datatype MPI_UNSIGNED_CHAR      =  6;
const MPI_Datatype MPI_UNSIGNED_SHORT     =  7;
const MPI_Datatype MPI_UNSIGNED           =  8;
const MPI_Datatype MPI_UNSIGNED_LONG      =  9;
const MPI_Datatype MPI_UNSIGNED_LONG_LONG = 10;
const MPI_Datatype MPI_FLOAT              = 11;
const MPI_Datatype MPI_DOUBLE             = 12;
const MPI_Datatype MPI_LONG_DOUBLE        = 13;

const int MPI_ORDER_C       = 0;
const int MPI_ORDER_FORTRAN = 0;

const MPI_Datatype MPI_ANY_SOURCE = -2;
const MPI_Datatype MPI_ANY_TAG    = -1;
const MPI_Datatype MPI_DATATYPE_NULL = 0;

const MPI_Op MPI_MIN  = 100;
const MPI_Op MPI_MAX  = 101;
const MPI_Op MPI_SUM  = 102;
const MPI_Op MPI_PROD = 103;
const MPI_Op MPI_LAND = 104;
const MPI_Op MPI_BAND = 105;
const MPI_Op MPI_LOR  = 106;
const MPI_Op MPI_BOR  = 107;
const MPI_Op MPI_LXOR = 108;
const MPI_Op MPI_BXOR = 109;

const int MPI_PACKED = 1;
const int MPI_UNDEFINED = -1;

#ifndef MPI_IN_PLACE
#define MPI_IN_PLACE ((void*)1)
#endif

const int MPI_MAX_PROCESSOR_NAME = 1;

const int MPI_MODE_WRONLY = 0;
const int MPI_MODE_CREATE = 1;
const int MPI_MODE_APPEND = 2;
const int MPI_MODE_RDONLY = 3;

const int MPI_INFO_NULL = 0;

const int MPI_REQUEST_NULL = 0;

const int MPI_FILE_NULL = 0;

static MPI_Status* MPI_STATUS_IGNORE   = nullptr;
static MPI_Status* MPI_STATUSES_IGNORE = nullptr;

namespace non_mpi_internal { inline void disableDefinedButNotUsedWarning() { MPI_STATUS_IGNORE = nullptr; MPI_STATUSES_IGNORE = nullptr; } }

const int MPI_MAX_OBJECT_NAME  = 255;
const int MPI_MAX_ERROR_STRING = 255;

#define WALBERLA_MPI_FUNCTION_ERROR WALBERLA_ABORT( "Invalid MPI function call! In case of compiling without MPI, MPI functions are not available and shouldn't be called!" );

inline int MPI_Init( int*, char*** )  { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Initialized( int *)    { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Finalize()             { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Abort( MPI_Comm, int ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Group_incl( MPI_Group, int, int*, MPI_Group * ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Group_excl( MPI_Group, int, int*, MPI_Group * ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Group_translate_ranks(MPI_Group, int, int*, MPI_Group, int*) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Group_free( MPI_Group *) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Comm_size( MPI_Comm, int* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_rank( MPI_Comm, int* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_get_name( MPI_Comm, char*, int* ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Info_create ( MPI_Info * ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Comm_group ( MPI_Comm, MPI_Group* )                         { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_create( MPI_Comm, MPI_Group, MPI_Comm* )               { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_free  ( MPI_Comm* )                                    { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_dup   ( MPI_Comm, MPI_Comm *)                          { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_split ( MPI_Comm, int, int, MPI_Comm *)                { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Comm_split_type ( MPI_Comm, int, int, MPI_Info, MPI_Comm *) { WALBERLA_MPI_FUNCTION_ERROR }


inline int MPI_Cart_create( MPI_Comm, int, int*, int*, int, MPI_Comm* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Dims_create( int, int, int* )                            { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Cart_coords( MPI_Comm, int, int, int* )                  { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Cart_rank  ( MPI_Comm, int*, int* )                      { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Cart_shift ( MPI_Comm, int, int, int*, int*)             { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Irecv( void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Isend( void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request* ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Recv( void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Send( void*, int, MPI_Datatype, int, int, MPI_Comm )              { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Sendrecv( void*, int, MPI_Datatype, int, int, void*, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Probe    ( int, int, MPI_Comm, MPI_Status* )       { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Iprobe   ( int, int, MPI_Comm, int*, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Get_count( MPI_Status*, MPI_Datatype, int* )       { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Start   ( MPI_Request* )      { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Startall( int, MPI_Request* ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Wait   ( MPI_Request*, MPI_Status* )            { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Waitall( int, MPI_Request*, MPI_Status* )       { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Waitany( int, MPI_Request*, int*, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Reduce   ( void*, void*, int, MPI_Datatype, MPI_Op, int, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Allreduce( void*, void*, int, MPI_Datatype, MPI_Op, MPI_Comm )      { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Scan  ( void*, void*, int, MPI_Datatype, MPI_Op, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Exscan( void*, void*, int, MPI_Datatype, MPI_Op, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Bcast     ( void*, int, MPI_Datatype, int, MPI_Comm )                      { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Allgather ( void*, int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Gather    ( void*, int, MPI_Datatype, void*, int, MPI_Datatype, int, MPI_Comm )  { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Allgatherv( void*, int, MPI_Datatype, void*, int*, int*, MPI_Datatype, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Gatherv   ( void*, int, MPI_Datatype, void*, int*, int*, MPI_Datatype, int, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Alltoall  ( void*, int, MPI_Datatype, void*, int, MPI_Datatype, MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Type_contiguous( int, MPI_Datatype, MPI_Datatype* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_create_subarray( int, const int*, const int*, const int*, int, MPI_Datatype, MPI_Datatype* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_create_indexed_block( int, int, const int*, MPI_Datatype, MPI_Datatype* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_commit( MPI_Datatype* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_free( MPI_Datatype* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_create_resized( MPI_Datatype, MPI_Aint, MPI_Aint, MPI_Datatype* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_size( MPI_Datatype, int * ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_get_extent(MPI_Datatype, MPI_Aint*, MPI_Aint*) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_Type_create_struct(int, const int[], const MPI_Aint[], const MPI_Datatype[], MPI_Datatype*) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Op_create(MPI_User_function*, int, MPI_Op*) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Get_processor_name( char*, int* ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Barrier( MPI_Comm ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_File_open        ( MPI_Comm, char*, int, int, MPI_File* )            { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_write_shared( MPI_File, void*, int, MPI_Datatype, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_write       ( MPI_File, void*, int, MPI_Datatype, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_write_all   ( MPI_File, void*, int, MPI_Datatype, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_write_at    ( MPI_File, MPI_Offset, void*, int, MPI_Datatype, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_read_all    ( MPI_File, void *, int, MPI_Datatype, MPI_Status * ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_read_at     ( MPI_File, MPI_Offset, void*, int, MPI_Datatype, MPI_Status* ) { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_close       ( MPI_File* )                                       { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_set_size    ( MPI_File, MPI_Offset )                            { WALBERLA_MPI_FUNCTION_ERROR }
inline int MPI_File_set_view    ( MPI_File, MPI_Offset, MPI_Datatype , MPI_Datatype, char*, MPI_Info ) { WALBERLA_MPI_FUNCTION_ERROR }

inline int MPI_Error_string ( int, char*, int* ) { WALBERLA_MPI_FUNCTION_ERROR }

inline double MPI_Wtime() { WALBERLA_MPI_FUNCTION_ERROR }

#undef WALBERLA_MPI_FUNCTION_ERROR

} // namespace mpistubs
using namespace mpistubs;

} // namespace walberla

#endif // WALBERLA_BUILD_WITH_MPI



namespace walberla {

namespace mpi {
   using MPIRank = int;
   using MPISize = int;
   const MPIRank INVALID_RANK = -1;
   const MPISize INVALID_SIZE = -1;
} // namespace MPI



/*!\class MPITrait
// \brief Base template for the MPITrait class
//
// The MPITrait class template offers a translation between the C++ built-in data types and
// the corresponding MPI data types required for send and receive operations. For a particular
// MPITrait instantiation, the corresponding MPI data type can be obtained by calling type()
// of the MPITrait. The following example demonstrates the application of the MPITrait class:

   \code
   // Initialization of the MPI communication
   int* pi;  // Integer buffer for the MPI send operation
   ...       // Initialization of the send buffer

   // Sending 50 integers to process 0
   MPI_Send( pi, 50, MPITrait< int >::type(), 0, 0, MPI_COMM_WORLD );
   \endcode
*/
template< typename T >
struct MPITrait;



/// Macro for the creation of MPITrait specializations for the supported data types.

#define WALBERLA_CREATE_MPITRAIT_SPECIALIZATION(CPP_TYPE,MPI_TYPE) \
   template<> \
   struct MPITrait< CPP_TYPE > \
   { \
      static inline MPI_Datatype type() { return (MPI_TYPE); } \
   }



// MPITRAIT SPECIALIZATIONS

WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( char               , MPI_CHAR               );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( signed char        , MPI_CHAR               );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( signed short int   , MPI_SHORT              );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( signed int         , MPI_INT                );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( signed long int    , MPI_LONG               );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( signed long long   , MPI_LONG_LONG          );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( unsigned char      , MPI_UNSIGNED_CHAR      );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( unsigned short int , MPI_UNSIGNED_SHORT     );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( unsigned int       , MPI_UNSIGNED           );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( unsigned long int  , MPI_UNSIGNED_LONG      );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( unsigned long long , MPI_UNSIGNED_LONG_LONG );
#ifdef WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT
   WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( walberla::float16  , MPI_WCHAR              );
#endif
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( float              , MPI_FLOAT              );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( double             , MPI_DOUBLE             );
WALBERLA_CREATE_MPITRAIT_SPECIALIZATION( long double        , MPI_LONG_DOUBLE        );

} // namespace walberla
/// \endcond

