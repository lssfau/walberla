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
//! \file ParMetisWrapper.cpp
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"

#ifdef WALBERLA_BUILD_WITH_PARMETIS

#ifdef _MSC_VER
#  pragma push_macro( "INT32_MIN" )
#  pragma push_macro( "INT32_MAX" )
#  pragma push_macro( "INT64_MIN" )
#  pragma push_macro( "INT64_MAX" )

#  ifdef INT32_MIN
#     undef INT32_MIN
#  endif

#  ifdef INT32_MAX
#     undef INT32_MAX
#  endif

#  ifdef INT64_MIN
#     undef INT64_MIN
#  endif

#  ifdef INT64_MAX
#     undef INT64_MAX
#  endif
#endif

#ifdef WALBERLA_BUILD_WITH_METIS
#  include "metis.h"
#endif

#ifdef _MSC_VER
#  pragma pop_macro( "INT64_MAX" )
#  pragma pop_macro( "INT64_MIN" )
#  pragma pop_macro( "INT32_MAX" )
#  pragma pop_macro( "INT32_MIN" )
#endif

#include <parmetis.h>
#endif

#include "ParMetisWrapper.h"

#include "core/Abort.h"

#include <type_traits>


namespace walberla {
namespace core {

#ifdef WALBERLA_BUILD_WITH_PARMETIS

int ParMETIS_V3_AdaptiveRepart(
   ::walberla::int64_t *vtxdist, ::walberla::int64_t *xadj, ::walberla::int64_t *adjncy, ::walberla::int64_t *vwgt,
   ::walberla::int64_t *vsize, ::walberla::int64_t *adjwgt, ::walberla::int64_t *wgtflag, ::walberla::int64_t *numflag, int64_t *ncon,
   ::walberla::int64_t *nparts, double *tpwgts, double *ubvec, double *ipc2redist,
   ::walberla::int64_t *options, ::walberla::int64_t *edgecut, ::walberla::int64_t *part, MPI_Comm *comm )
{

   static_assert( std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!"        );
   static_assert( std::is_same< double, ::real_t >::value,             "You have to compile the metis library with 64-bit wide floating-point type support!" );

   return ::ParMETIS_V3_AdaptiveRepart( vtxdist, xadj, adjncy, vwgt,
                                        vsize, adjwgt, wgtflag, numflag, ncon,
                                        nparts, tpwgts, ubvec, ipc2redist,
                                        options, edgecut, part, comm );
}

int ParMETIS_V3_PartKway(
   ::walberla::int64_t *vtxdist, ::walberla::int64_t *xadj, ::walberla::int64_t *adjncy, ::walberla::int64_t *vwgt,
   ::walberla::int64_t *adjwgt, ::walberla::int64_t *wgtflag, ::walberla::int64_t *numflag, ::walberla::int64_t *ncon, ::walberla::int64_t *nparts,
   double *tpwgts, double *ubvec, ::walberla::int64_t *options, ::walberla::int64_t *edgecut, ::walberla::int64_t *part,
   MPI_Comm *comm )
{
   static_assert( std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!" );
   static_assert( std::is_same< double, ::real_t >::value, "You have to compile the metis library with 64-bit wide floating-point type support!" );

   return ::ParMETIS_V3_PartKway( vtxdist, xadj, adjncy, vwgt,
                                  adjwgt, wgtflag, numflag, ncon, nparts,
                                  tpwgts, ubvec, options, edgecut, part, comm );

}

int ParMETIS_V3_PartGeomKway(
   ::walberla::int64_t *vtxdist, ::walberla::int64_t *xadj, ::walberla::int64_t *adjncy, ::walberla::int64_t *vwgt,
   ::walberla::int64_t *adjwgt, ::walberla::int64_t *wgtflag, ::walberla::int64_t *numflag, ::walberla::int64_t *ndims, double *xyz,
   ::walberla::int64_t *ncon, ::walberla::int64_t *nparts, double *tpwgts, double *ubvec, ::walberla::int64_t *options,
   ::walberla::int64_t *edgecut, ::walberla::int64_t *part, MPI_Comm *comm )
{
   static_assert( std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!" );
   static_assert( std::is_same< double, ::real_t >::value, "You have to compile the metis library with 64-bit wide floating-point type support!" );

   return ::ParMETIS_V3_PartGeomKway( vtxdist, xadj, adjncy, vwgt,
                                      adjwgt, wgtflag, numflag, ndims, xyz,
                                      ncon, nparts, tpwgts, ubvec, options,
                                      edgecut, part, comm );
}

int ParMETIS_V3_PartGeom(
   ::walberla::int64_t *vtxdist, ::walberla::int64_t *ndims, double *xyz, ::walberla::int64_t *part, MPI_Comm *comm )
{
   static_assert( std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!" );
   static_assert( std::is_same< double, ::real_t >::value, "You have to compile the metis library with 64-bit wide floating-point type support!" );

   return ::ParMETIS_V3_PartGeom( vtxdist, ndims, xyz, part, comm );
}

int ParMETIS_V3_RefineKway(
   ::walberla::int64_t *vtxdist, ::walberla::int64_t *xadj, ::walberla::int64_t *adjncy, ::walberla::int64_t *vwgt,
   ::walberla::int64_t *adjwgt, ::walberla::int64_t *wgtflag, ::walberla::int64_t *numflag, ::walberla::int64_t *ncon, ::walberla::int64_t *nparts,
   double *tpwgts, double *ubvec, ::walberla::int64_t *options, ::walberla::int64_t *edgecut,
   ::walberla::int64_t *part, MPI_Comm *comm )
{
   static_assert( std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!" );
   static_assert( std::is_same< double, ::real_t >::value, "You have to compile the metis library with 64-bit wide floating-point type support!" );

   return ::ParMETIS_V3_RefineKway( vtxdist, xadj, adjncy, vwgt,
                                    adjwgt, wgtflag, numflag, ncon, nparts,
                                    tpwgts, ubvec, options, edgecut,
                                    part, comm );
}

#else // build without ParMetis

int ParMETIS_V3_AdaptiveRepart(
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, int64_t *,
   ::walberla::int64_t *, double *, double *, double *,
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, MPI_Comm * )
{

   WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );
}

int ParMETIS_V3_PartKway(
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   double *, double *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   MPI_Comm * )
{
   WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );

}

int ParMETIS_V3_PartGeomKway(
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, double *,
   ::walberla::int64_t *, ::walberla::int64_t *, double *, double *, ::walberla::int64_t *,
   ::walberla::int64_t *, ::walberla::int64_t *, MPI_Comm * )
{
   WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );
}

int ParMETIS_V3_PartGeom(
   ::walberla::int64_t *, ::walberla::int64_t *, double *, ::walberla::int64_t *, MPI_Comm * )
{
   WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );
}

int ParMETIS_V3_RefineKway(
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *, ::walberla::int64_t *,
   double *, double *, ::walberla::int64_t *, ::walberla::int64_t *,
   ::walberla::int64_t *, MPI_Comm * )
{
   WALBERLA_ABORT( "You are trying to use ParMetis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_PARMETIS' to 'ON' in your CMake cache to build against an installed version of ParMetis!" );
}

#endif


} // namespace core
} // namespace walberla
