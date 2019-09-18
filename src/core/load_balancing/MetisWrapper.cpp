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
//! \file MetisWrapper.cpp
//! \ingroup blockforest
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "MetisWrapper.h"
#define WALBERLA_METIS_NOPTIONS METIS_NOPTIONS
#undef METIS_NOPTIONS

#include "waLBerlaDefinitions.h"

#include "core/Abort.h"

#include <type_traits>

#ifdef WALBERLA_BUILD_WITH_METIS

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

#include "metis.h"

#ifdef _MSC_VER
#  pragma pop_macro( "INT64_MAX" )
#  pragma pop_macro( "INT64_MIN" )
#  pragma pop_macro( "INT32_MAX" )
#  pragma pop_macro( "INT32_MIN" )
#endif

static_assert(WALBERLA_METIS_NOPTIONS == METIS_NOPTIONS, "The macro METIS_NOPTIONS defined in blockforest/loadbalancing/MetisWRapper.h does not match the value in metis.h!");

#endif // WALBERLA_BUILD_WITH_METIS

namespace walberla {
namespace core {

#ifdef WALBERLA_BUILD_WITH_METIS

int METIS_PartGraphKway( ::walberla::int64_t *nvtxs, ::walberla::int64_t *ncon, ::walberla::int64_t *xadj, ::walberla::int64_t *adjncy,
                         ::walberla::int64_t *vwgt, ::walberla::int64_t *vsize, ::walberla::int64_t *adjwgt, ::walberla::int64_t *nparts,
                         double *tpwgts, double *ubvec, ::walberla::int64_t *options, ::walberla::int64_t *edgecut, ::walberla::int64_t *part)
{
   static_assert(std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!");
   static_assert(std::is_same< double, ::real_t >::value, "You have to compile the metis library with 64-bit wide floating-point type support!");

   return ::METIS_PartGraphKway( nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, edgecut, part );
}


int METIS_PartGraphRecursive( ::walberla::int64_t *nvtxs, ::walberla::int64_t *ncon, ::walberla::int64_t *xadj, ::walberla::int64_t *adjncy,
                              ::walberla::int64_t *vwgt, ::walberla::int64_t *vsize, ::walberla::int64_t *adjwgt, ::walberla::int64_t *nparts,
                              double *tpwgts, double *ubvec, ::walberla::int64_t *options, ::walberla::int64_t *edgecut, ::walberla::int64_t *part)
{
   static_assert(std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!");
   static_assert(std::is_same< double, ::real_t >::value, "You have to compile the metis library with 64-bit wide floating-point type support!");

   return ::METIS_PartGraphRecursive( nvtxs, ncon, xadj, adjncy, vwgt, vsize, adjwgt, nparts, tpwgts, ubvec, options, edgecut, part );
}

int METIS_SetDefaultOptions( ::walberla::int64_t *options )
{
   static_assert(std::is_same< ::walberla::int64_t, ::idx_t >::value, "You have to compile the metis library with 64-bit wide integer type support!");
   
   return ::METIS_SetDefaultOptions( options );
}



const int METIS_OK           = ::METIS_OK;
const int METIS_ERROR        = ::METIS_ERROR;
const int METIS_ERROR_INPUT  = ::METIS_ERROR_INPUT;
const int METIS_ERROR_MEMORY = ::METIS_ERROR_MEMORY;

const int METIS_OPTION_PTYPE      = ::METIS_OPTION_PTYPE;   
const int METIS_OPTION_OBJTYPE    = ::METIS_OPTION_OBJTYPE;
const int METIS_OPTION_CTYPE      = ::METIS_OPTION_CTYPE; 
const int METIS_OPTION_IPTYPE     = ::METIS_OPTION_IPTYPE;   
const int METIS_OPTION_RTYPE      = ::METIS_OPTION_RTYPE;  
const int METIS_OPTION_DBGLVL     = ::METIS_OPTION_DBGLVL;   
const int METIS_OPTION_NITER      = ::METIS_OPTION_NITER;  
const int METIS_OPTION_NCUTS      = ::METIS_OPTION_NCUTS;   
const int METIS_OPTION_SEED       = ::METIS_OPTION_SEED;   
const int METIS_OPTION_NO2HOP     = ::METIS_OPTION_NO2HOP;   
const int METIS_OPTION_MINCONN    = ::METIS_OPTION_MINCONN;  
const int METIS_OPTION_CONTIG     = ::METIS_OPTION_CONTIG; 
const int METIS_OPTION_COMPRESS   = ::METIS_OPTION_COMPRESS; 
const int METIS_OPTION_CCORDER    = ::METIS_OPTION_CCORDER;
const int METIS_OPTION_PFACTOR    = ::METIS_OPTION_PFACTOR; 
const int METIS_OPTION_NSEPS      = ::METIS_OPTION_NSEPS; 
const int METIS_OPTION_UFACTOR    = ::METIS_OPTION_UFACTOR;  
const int METIS_OPTION_NUMBERING  = ::METIS_OPTION_NUMBERING;
const int METIS_OPTION_HELP       = ::METIS_OPTION_HELP;
const int METIS_OPTION_TPWGTS     = ::METIS_OPTION_TPWGTS;   
const int METIS_OPTION_NCOMMON    = ::METIS_OPTION_NCOMMON;  
const int METIS_OPTION_NOOUTPUT   = ::METIS_OPTION_NOOUTPUT; 
const int METIS_OPTION_BALANCE    = ::METIS_OPTION_BALANCE;
const int METIS_OPTION_GTYPE      = ::METIS_OPTION_GTYPE;    
const int METIS_OPTION_UBVEC      = ::METIS_OPTION_UBVEC;   


#else // build without Metis

int METIS_PartGraphKway( ::walberla::int64_t * /*nvtxs*/, ::walberla::int64_t * /*ncon*/, ::walberla::int64_t * /*xadj*/, ::walberla::int64_t * /*adjncy*/,
                         ::walberla::int64_t * /*vwgt*/, ::walberla::int64_t * /*vsize*/, ::walberla::int64_t * /*adjwgt*/, ::walberla::int64_t * /*nparts*/,
                         double * /*tpwgts*/, double * /*ubvec*/, ::walberla::int64_t * /*options*/, ::walberla::int64_t * /*edgecut*/, ::walberla::int64_t * /*part*/)
{
    WALBERLA_ABORT( "You are trying to use Metis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_METIS' to 'ON' in your CMake cache to build against an installed version of Metis!" );
}

int METIS_SetDefaultOptions( ::walberla::int64_t * /*options*/ )
{
   WALBERLA_ABORT("You are trying to use Metis functionality but waLBerla is not configured to use it. Set 'WALBERLA_BUILD_WITH_METIS' to 'ON' in your CMake cache to build against an installed version of Metis!");
}

const int METIS_OK           = 0;
const int METIS_ERROR        = 0;
const int METIS_ERROR_INPUT  = 0;
const int METIS_ERROR_MEMORY = 0;

const int METIS_OPTION_PTYPE      = 0;
const int METIS_OPTION_OBJTYPE    = 0;
const int METIS_OPTION_CTYPE      = 0;
const int METIS_OPTION_IPTYPE     = 0; 
const int METIS_OPTION_RTYPE      = 0;
const int METIS_OPTION_DBGLVL     = 0; 
const int METIS_OPTION_NITER      = 0;
const int METIS_OPTION_NCUTS      = 0;
const int METIS_OPTION_SEED       = 0;
const int METIS_OPTION_NO2HOP     = 0; 
const int METIS_OPTION_MINCONN    = 0; 
const int METIS_OPTION_CONTIG     = 0;
const int METIS_OPTION_COMPRESS   = 0; 
const int METIS_OPTION_CCORDER    = 0;
const int METIS_OPTION_PFACTOR    = 0;
const int METIS_OPTION_NSEPS      = 0;
const int METIS_OPTION_UFACTOR    = 0; 
const int METIS_OPTION_NUMBERING  = 0;
const int METIS_OPTION_HELP       = 0;
const int METIS_OPTION_TPWGTS     = 0; 
const int METIS_OPTION_NCOMMON    = 0; 
const int METIS_OPTION_NOOUTPUT   = 0; 
const int METIS_OPTION_BALANCE    = 0;
const int METIS_OPTION_GTYPE      = 0; 
const int METIS_OPTION_UBVEC      = 0;

#endif


} // namespace core
} // namespace walberla
