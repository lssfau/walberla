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
//! \file IteratorMacros.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//!
//! This file provides macros which make iterating over all cells of one or more fields easier.
//!
//! \section docFieldIteratorMacros Field Iterator Macros
//!
//! Iterator macros iterate cells, not all values in all cells, i.e., iteration includes x, y, and z, but **not** f values!
//!
//! ATTENTION: If OpenMP is enabled, iterating all cells is done in parallel!
//!            Meaning: Your code that is supposed to be executed for each cell must be free of any race condition!
//!
//! Existing macros are:
//!
//! - `WALBERLA_FOR_ALL_CELLS_XYZ`
//! - `WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ`
//! - `WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ`
//! - `WALBERLA_FOR_ALL_CELLS_YZ`
//! - `WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ`
//! - `WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ`
//! - `WALBERLA_FOR_ALL_CELLS`
//!
//! All macros also exist with `_OMP` at the end (for these macros, the OpenMP pragma can be specified by hand):
//!
//! - `WALBERLA_FOR_ALL_CELLS_XYZ_OMP`
//! - `WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ_OMP`
//! - `WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP`
//! - `WALBERLA_FOR_ALL_CELLS_YZ_OMP`
//! - `WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ_OMP`
//! - `WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP`
//! - `WALBERLA_FOR_ALL_CELLS_OMP`
//!
//! Some examples:
//!
//! Increase all values by '42':
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_XYZ( ptrToField,
//!
//!    ptrToField->get(x,y,z) += 42;
//! )
//! \endcode
//!
//! Increase all values by '42', including all cells in all ghost layers:
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( ptrToField,
//!
//!    ptrToField->get(x,y,z) += 42;
//! )
//! \endcode
//!
//! Increase all values by '42', including all cells in the first and second ghost layer:
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( ptrToField, 2,
//!
//!    ptrToField->get(x,y,z) += 42;
//! )
//! \endcode
//!
//! Increase all values by '42', using a field iterator:
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS( fieldIterator, ptrToField,
//!
//!    *fieldIterator += 42;
//! )
//! \endcode
//!
//! Add values from one field to another field:
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_XYZ( ptrToField,
//!
//!    ptrToAnotherField->get(x,y,z) += ptrToField->get(x,y,z);
//! )
//! \endcode
//!
//! Add values from one field to another field using field iterators:
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS( fieldIterator, ptrToField, anotherFieldIterator, ptrToAnotherField,
//!
//!    *anotherFieldIterator += *fieldIterator;
//! )
//! \endcode
//!
//! Only use every second x-coordinate for computation (aka: write the innermost loop yourself!):
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_YZ( ptrToField,
//!
//!    for( ::walberla::cell_idx_t x = ::walberla::cell_idx_t(0); x < ::walberla::cell_idx_c( ptrToField->xSize() ); x += ::walberla::cell_idx_t(2) )
//!       ptrToField->get(x,y,z) += 42;
//! )
//! \endcode
//!
//! If you want to customize to pragma used for OpenMP, you can do so by using the `_OMP` macros.
//! The OpenMP pragma is always the last parameter of the macro prior the the code itself.
//! Example (calculate scalar product):
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_XYZ_OMP( ptrToField, omp parallel for schedule(static) reduction(+:sum)
//!
//!    const auto value = ptrToField->get(x,y,z);
//!    sum += value * value;
//! )
//! \endcode
//!
//! The same code using a field iterator:
//!
//! \code
//! WALBERLA_FOR_ALL_CELLS_OMP( fieldIterator, ptrToField, omp parallel for schedule(static) reduction(+:sum)
//!
//!    const auto value = *fieldIterator;
//!    sum += value * value;
//! )
//! \endcode
//!
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "core/debug/Debug.h"

/////////////////////////////
//  WALBERLA_FOR_ALL_CELLS //
/////////////////////////////

#ifdef _OPENMP

#if (defined(_MSC_VER) && _MSC_VER < 1926)

#define WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (field)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = ::walberla::cell_idx_t(0); y < ySize__; ++y ) { \
            for( ::walberla::cell_idx_t x = ::walberla::cell_idx_t(0); x < xSize__; ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = ::walberla::cell_idx_t(0); z < zSize__; ++z ) { \
            for( ::walberla::cell_idx_t x = ::walberla::cell_idx_t(0); x < xSize__; ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ_OMP( interval, omp, CODE ) \
   { if( interval.zSize() >= interval.ySize() ) \
   { \
      const int izMax = ::walberla::int_c( interval.zMax() ); \
      __pragma(omp) \
      for( int iz = ::walberla::int_c( interval.zMin() ); iz <= izMax; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) { \
            for( ::walberla::cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iyMax = ::walberla::int_c( interval.yMax() ); \
      __pragma(omp) \
      for( int iy = ::walberla::int_c( interval.yMin() ); iy <= iyMax; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) { \
            for( ::walberla::cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } }
   
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_4( field, gl, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   WALBERLA_ASSERT_GREATER_EQUAL_2( (field)->nrOfGhostLayers(), gl ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (field)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   const ::walberla::cell_idx_t gl__ = ::walberla::cell_idx_c( gl ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         if( iz == 0 ) \
         { \
            for( ::walberla::cell_idx_t z = -gl__; z < ::walberla::cell_idx_t(0); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
            for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
               for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
               { \
                  CODE \
               } \
            } \
         } \
         if( iz == (izSize - 1) ) \
         { \
            for( ::walberla::cell_idx_t z = zSize__; z < (zSize__ + gl__); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         if( iy == 0 ) \
         { \
            for( ::walberla::cell_idx_t y = -gl__; y < ::walberla::cell_idx_t(0); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
            for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
               for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
               { \
                  CODE \
               } \
            } \
         } \
         if( iy == (iySize - 1) ) \
         { \
            for( ::walberla::cell_idx_t y = ySize__; y < (ySize__ + gl__); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
      } \
   } }
   
#define WALBERLA_FOR_ALL_CELLS_YZ_OMP( field, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = ::walberla::cell_idx_t(0); y < ySize__; ++y ) \
         { \
            CODE \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = ::walberla::cell_idx_t(0); z < zSize__; ++z ) \
         { \
            CODE \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ_OMP( interval, omp, CODE ) \
   { if( interval.zSize() >= interval.ySize() ) \
   { \
      const int izMax = ::walberla::int_c( interval.zMax() ); \
      __pragma(omp) \
      for( int iz = ::walberla::int_c( interval.zMin() ); iz <= izMax; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) \
         { \
            CODE \
         } \
      } \
   } \
   else \
   { \
      const int iyMax = ::walberla::int_c( interval.yMax() ); \
      __pragma(omp) \
      for( int iy = ::walberla::int_c( interval.yMin() ); iy <= iyMax; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) \
         { \
            CODE \
         } \
      } \
   } }
   
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_4( field, gl, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   WALBERLA_ASSERT_GREATER_EQUAL_2( (field)->nrOfGhostLayers(), gl ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   const ::walberla::cell_idx_t gl__ = ::walberla::cell_idx_c( gl ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         if( iz == 0 ) \
         { \
            for( ::walberla::cell_idx_t z = -gl__; z < ::walberla::cell_idx_t(0); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
               { \
                  CODE \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
            for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
            { \
               CODE \
            } \
         } \
         if( iz == (izSize - 1) ) \
         { \
            for( ::walberla::cell_idx_t z = zSize__; z < (zSize__ + gl__); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
               { \
                     CODE \
               } \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         if( iy == 0 ) \
         { \
            for( ::walberla::cell_idx_t y = -gl__; y < ::walberla::cell_idx_t(0); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) \
               { \
                  CODE \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
            for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) \
            { \
               CODE \
            } \
         } \
         if( iy == (iySize - 1) ) \
         { \
            for( ::walberla::cell_idx_t y = ySize__; y < (ySize__ + gl__); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) \
               { \
                  CODE \
               } \
            } \
         } \
      } \
   } }
   
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_4( it0, f0, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0) \
         { \
            CODE \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0) \
         { \
            CODE \
         } \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_6( it0, f0, it1, f1, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
      } \
   } }
   
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_8( it0, f0, it1, f1, it2, f2, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_10( it0, f0, it1, f1, it2, f2, it3, f3, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_12( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_14( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_16( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f6) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f6)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_18( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, it7, f7, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f6) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f7) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f6)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f7)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      __pragma(omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         auto it7 = (f7)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6, ++it7) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it7, (f7)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      __pragma(omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         auto it7 = (f7)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6, ++it7) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it7, (f7)->end() ); \
      } \
   } }

#else // == MSVC >= 2019 16.6 or not MSVC

#define WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (field)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = ::walberla::cell_idx_t(0); y < ySize__; ++y ) { \
            for( ::walberla::cell_idx_t x = ::walberla::cell_idx_t(0); x < xSize__; ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = ::walberla::cell_idx_t(0); z < zSize__; ++z ) { \
            for( ::walberla::cell_idx_t x = ::walberla::cell_idx_t(0); x < xSize__; ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ_OMP( interval, omp, CODE ) \
   { if( interval.zSize() >= interval.ySize() ) \
   { \
      const int izMax = ::walberla::int_c( interval.zMax() ); \
      _Pragma(#omp) \
      for( int iz = ::walberla::int_c( interval.zMin() ); iz <= izMax; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) { \
            for( ::walberla::cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iyMax = ::walberla::int_c( interval.yMax() ); \
      _Pragma(#omp) \
      for( int iy = ::walberla::int_c( interval.yMin() ); iy <= iyMax; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) { \
            for( ::walberla::cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x ) \
            { \
               CODE \
            } \
         } \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_4( field, gl, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   WALBERLA_ASSERT_GREATER_EQUAL_2( (field)->nrOfGhostLayers(), gl ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (field)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   const ::walberla::cell_idx_t gl__ = ::walberla::cell_idx_c( gl ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         if( iz == 0 ) \
         { \
            for( ::walberla::cell_idx_t z = -gl__; z < ::walberla::cell_idx_t(0); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
            for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
               for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
               { \
                  CODE \
               } \
            } \
         } \
         if( iz == (izSize - 1) ) \
         { \
            for( ::walberla::cell_idx_t z = zSize__; z < (zSize__ + gl__); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         if( iy == 0 ) \
         { \
            for( ::walberla::cell_idx_t y = -gl__; y < ::walberla::cell_idx_t(0); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
            for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
               for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
               { \
                  CODE \
               } \
            } \
         } \
         if( iy == (iySize - 1) ) \
         { \
            for( ::walberla::cell_idx_t y = ySize__; y < (ySize__ + gl__); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
                  for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
                  { \
                     CODE \
                  } \
               } \
            } \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_YZ_OMP( field, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = ::walberla::cell_idx_t(0); y < ySize__; ++y ) \
         { \
            CODE \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = ::walberla::cell_idx_t(0); z < zSize__; ++z ) \
         { \
            CODE \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ_OMP( interval, omp, CODE ) \
   { if( interval.zSize() >= interval.ySize() ) \
   { \
      const int izMax = ::walberla::int_c( interval.zMax() ); \
      _Pragma(#omp) \
      for( int iz = ::walberla::int_c( interval.zMin() ); iz <= izMax; ++iz ) { \
         ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         for( ::walberla::cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) \
         { \
            CODE \
         } \
      } \
   } \
   else \
   { \
      const int iyMax = ::walberla::int_c( interval.yMax() ); \
      _Pragma(#omp) \
      for( int iy = ::walberla::int_c( interval.yMin() ); iy <= iyMax; ++iy ) { \
         ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         for( ::walberla::cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) \
         { \
            CODE \
         } \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_4( field, gl, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   WALBERLA_ASSERT_GREATER_EQUAL_2( (field)->nrOfGhostLayers(), gl ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   const ::walberla::cell_idx_t gl__ = ::walberla::cell_idx_c( gl ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         if( iz == 0 ) \
         { \
            for( ::walberla::cell_idx_t z = -gl__; z < ::walberla::cell_idx_t(0); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
               { \
                  CODE \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
            for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
            { \
               CODE \
            } \
         } \
         if( iz == (izSize - 1) ) \
         { \
            for( ::walberla::cell_idx_t z = zSize__; z < (zSize__ + gl__); ++z ) { \
               for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
               { \
                     CODE \
               } \
            } \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         if( iy == 0 ) \
         { \
            for( ::walberla::cell_idx_t y = -gl__; y < ::walberla::cell_idx_t(0); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) \
               { \
                  CODE \
               } \
            } \
         } \
         { \
            ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
            for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) \
            { \
               CODE \
            } \
         } \
         if( iy == (iySize - 1) ) \
         { \
            for( ::walberla::cell_idx_t y = ySize__; y < (ySize__ + gl__); ++y ) { \
               for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) \
               { \
                  CODE \
               } \
            } \
         } \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_4( it0, f0, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0) \
         { \
            CODE \
         } \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0) \
         { \
            CODE \
         } \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_6( it0, f0, it1, f1, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_8( it0, f0, it1, f1, it2, f2, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_10( it0, f0, it1, f1, it2, f2, it3, f3, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_12( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_14( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_16( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f6) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f6)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_18( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, it7, f7, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f6) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f7) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f6)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f7)->xyzSize() ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (f0)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (f0)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (f0)->zSize() ); \
   if( zSize__ >= ySize__ ) \
   { \
      const int izSize = ::walberla::int_c( zSize__ ); \
      _Pragma(#omp) \
      for( int iz = 0; iz < izSize; ++iz ) { \
         const ::walberla::cell_idx_t z = ::walberla::cell_idx_c( iz ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), ::walberla::cell_idx_t(0), z, xSize__ - ::walberla::cell_idx_t(1), ySize__ - ::walberla::cell_idx_t(1), z ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         auto it7 = (f7)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6, ++it7) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it7, (f7)->end() ); \
      } \
   } \
   else \
   { \
      const int iySize = ::walberla::int_c( ySize__ ); \
      _Pragma(#omp) \
      for( int iy = 0; iy < iySize; ++iy ) { \
         const ::walberla::cell_idx_t y = ::walberla::cell_idx_c( iy ); \
         const CellInterval interval( ::walberla::cell_idx_t(0), y, ::walberla::cell_idx_t(0), xSize__ - ::walberla::cell_idx_t(1), y, zSize__ - ::walberla::cell_idx_t(1) ); \
         auto it0 = (f0)->beginSliceXYZ( interval ); \
         auto it1 = (f1)->beginSliceXYZ( interval ); \
         auto it2 = (f2)->beginSliceXYZ( interval ); \
         auto it3 = (f3)->beginSliceXYZ( interval ); \
         auto it4 = (f4)->beginSliceXYZ( interval ); \
         auto it5 = (f5)->beginSliceXYZ( interval ); \
         auto it6 = (f6)->beginSliceXYZ( interval ); \
         auto it7 = (f7)->beginSliceXYZ( interval ); \
         for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6, ++it7) \
         { \
            CODE \
         } \
         WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
         WALBERLA_ASSERT_EQUAL_2( it7, (f7)->end() ); \
      } \
   } }

#endif

#else // == no OpenMP

#define WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (field)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   for( ::walberla::cell_idx_t z = ::walberla::cell_idx_t(0); z < zSize__; ++z ) { \
      for( ::walberla::cell_idx_t y = ::walberla::cell_idx_t(0); y < ySize__; ++y ) { \
         for( ::walberla::cell_idx_t x = ::walberla::cell_idx_t(0); x < xSize__; ++x ) \
         { \
            CODE \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ_OMP( interval, omp, CODE ) \
   { for( ::walberla::cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) { \
      for( ::walberla::cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) { \
         for( ::walberla::cell_idx_t x = interval.xMin(); x <= interval.xMax(); ++x ) \
         { \
            CODE \
         } \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_4( field, gl, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   WALBERLA_ASSERT_GREATER_EQUAL_2( (field)->nrOfGhostLayers(), gl ); \
   const ::walberla::cell_idx_t xSize__ = ::walberla::cell_idx_c( (field)->xSize() ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   const ::walberla::cell_idx_t gl__ = ::walberla::cell_idx_c( gl ); \
   for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
      for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) { \
         for( ::walberla::cell_idx_t x = -gl__; x < (xSize__ + gl__); ++x ) \
         { \
            CODE \
         } \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_YZ_OMP( field, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   for( ::walberla::cell_idx_t z = ::walberla::cell_idx_t(0); z < zSize__; ++z ) { \
      for( ::walberla::cell_idx_t y = ::walberla::cell_idx_t(0); y < ySize__; ++y ) \
      { \
         CODE \
      } \
   } }

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ_OMP( interval, omp, CODE ) \
   { for( ::walberla::cell_idx_t z = interval.zMin(); z <= interval.zMax(); ++z ) { \
      for( ::walberla::cell_idx_t y = interval.yMin(); y <= interval.yMax(); ++y ) \
      { \
         CODE \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_4( field, gl, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (field) ); \
   WALBERLA_ASSERT_GREATER_EQUAL_2( (field)->nrOfGhostLayers(), gl ); \
   const ::walberla::cell_idx_t ySize__ = ::walberla::cell_idx_c( (field)->ySize() ); \
   const ::walberla::cell_idx_t zSize__ = ::walberla::cell_idx_c( (field)->zSize() ); \
   const ::walberla::cell_idx_t gl__ = ::walberla::cell_idx_c( gl ); \
   for( ::walberla::cell_idx_t z = -gl__; z < (zSize__ + gl__); ++z ) { \
      for( ::walberla::cell_idx_t y = -gl__; y < (ySize__ + gl__); ++y ) \
      { \
         CODE \
      } \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_4( it0, f0, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   auto it0 = (f0)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0) \
   { \
      CODE \
   } }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_6( it0, f0, it1, f1, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_8( it0, f0, it1, f1, it2, f2, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   auto it2 = (f2)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_10( it0, f0, it1, f1, it2, f2, it3, f3, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   auto it2 = (f2)->beginXYZ(); \
   auto it3 = (f3)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_12( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   auto it2 = (f2)->beginXYZ(); \
   auto it3 = (f3)->beginXYZ(); \
   auto it4 = (f4)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_14( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   auto it2 = (f2)->beginXYZ(); \
   auto it3 = (f3)->beginXYZ(); \
   auto it4 = (f4)->beginXYZ(); \
   auto it5 = (f5)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_16( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f6) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f6)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   auto it2 = (f2)->beginXYZ(); \
   auto it3 = (f3)->beginXYZ(); \
   auto it4 = (f4)->beginXYZ(); \
   auto it5 = (f5)->beginXYZ(); \
   auto it6 = (f6)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); }

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_OMP_18( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, it7, f7, omp, CODE ) \
   { WALBERLA_ASSERT_NOT_NULLPTR_1( (f0) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f1) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f2) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f3) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f4) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f5) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f6) ); \
   WALBERLA_ASSERT_NOT_NULLPTR_1( (f7) ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f1)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f2)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f3)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f4)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f5)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f6)->xyzSize() ); \
   WALBERLA_ASSERT_EQUAL_2( (f0)->xyzSize(), (f7)->xyzSize() ); \
   auto it0 = (f0)->beginXYZ(); \
   auto it1 = (f1)->beginXYZ(); \
   auto it2 = (f2)->beginXYZ(); \
   auto it3 = (f3)->beginXYZ(); \
   auto it4 = (f4)->beginXYZ(); \
   auto it5 = (f5)->beginXYZ(); \
   auto it6 = (f6)->beginXYZ(); \
   auto it7 = (f7)->beginXYZ(); \
   for(/* see above */; it0 != (f0)->end(); ++it0, ++it1, ++it2, ++it3, ++it4, ++it5, ++it6, ++it7) \
   { \
      CODE \
   } \
   WALBERLA_ASSERT_EQUAL_2( it1, (f1)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it2, (f2)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it3, (f3)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it4, (f4)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it5, (f5)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it6, (f6)->end() ); \
   WALBERLA_ASSERT_EQUAL_2( it7, (f7)->end() ); }

#endif // OpenMP



// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_3( field, omp, CODE ) \
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_4( field, (field)->nrOfGhostLayers(), omp, CODE ) \

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_3( field, gl, CODE ) \
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_4( field, gl, omp parallel for schedule(static), CODE )        

// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ' (using the same signature) instead        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_2( field, CODE ) \
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_4( field, (field)->nrOfGhostLayers(), omp parallel for schedule(static), CODE )          
    
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_3( field, omp, CODE ) \
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_4( field, (field)->nrOfGhostLayers(), omp, CODE ) \
        
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_3( field, gl, CODE ) \
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_4( field, gl, omp parallel for schedule(static), CODE )   
        
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ' (using the same signature) instead        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_2( field, CODE ) \
        WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_4( field, (field)->nrOfGhostLayers(), omp parallel for schedule(static), CODE )    
        
        
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_2(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_7(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_8(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_9(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_10(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO        
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_7(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_8(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_9(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_10(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_2(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_7(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_8(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_9(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_10(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO  
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_7(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_8(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_9(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_10(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ITERATOR_MACRO
        
        
        
#define WALBERLA_FOR_ALL_CELLS_XYZ( field, CODE ) \
        WALBERLA_FOR_ALL_CELLS_XYZ_OMP( field, omp parallel for schedule(static), CODE )

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ( interval, CODE ) \
        WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ_OMP( interval, omp parallel for schedule(static), CODE )
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(...) \
        WALBERLA_MACRO_OVERLOAD( WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP_, __VA_ARGS__ )    
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(...) \
        WALBERLA_MACRO_OVERLOAD( WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_, __VA_ARGS__ )
        
#define WALBERLA_FOR_ALL_CELLS_YZ( field, CODE ) \
        WALBERLA_FOR_ALL_CELLS_YZ_OMP( field, omp parallel for schedule(static), CODE )        

#define WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ( interval, CODE ) \
        WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_YZ_OMP( interval, omp parallel for schedule(static), CODE )
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP(...) \
        WALBERLA_MACRO_OVERLOAD( WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_OMP_, __VA_ARGS__ )    
        
#define WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ(...) \
        WALBERLA_MACRO_OVERLOAD( WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_YZ_, __VA_ARGS__ )
        


// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_3( it0, f0, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_4( it0, f0, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_5( it0, f0, it1, f1, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_6( it0, f0, it1, f1, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_7( it0, f0, it1, f1, it2, f2, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_8( it0, f0, it1, f1, it2, f2, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_9( it0, f0, it1, f1, it2, f2, it3, f3, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_10( it0, f0, it1, f1, it2, f2, it3, f3, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_11( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_12( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_13( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_14( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_15( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_16( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, omp parallel for schedule(static), CODE )
// Do not call this macro, call 'WALBERLA_FOR_ALL_CELLS' (using the same signature) instead
#define WALBERLA_FOR_ALL_CELLS_17( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, it7, f7, CODE ) \
        WALBERLA_FOR_ALL_CELLS_OMP_18( it0, f0, it1, f1, it2, f2, it3, f3, it4, f4, it5, f5, it6, f6, it7, f7, omp parallel for schedule(static), CODE )

        
#define WALBERLA_FOR_ALL_CELLS_OMP_1(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_1
#define WALBERLA_FOR_ALL_CELLS_OMP_2(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_2
#define WALBERLA_FOR_ALL_CELLS_OMP_3(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_3
#define WALBERLA_FOR_ALL_CELLS_OMP_5(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_5
#define WALBERLA_FOR_ALL_CELLS_OMP_7(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_7
#define WALBERLA_FOR_ALL_CELLS_OMP_9(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_9
#define WALBERLA_FOR_ALL_CELLS_OMP_11(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_11
#define WALBERLA_FOR_ALL_CELLS_OMP_13(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_13
#define WALBERLA_FOR_ALL_CELLS_OMP_15(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_15
#define WALBERLA_FOR_ALL_CELLS_OMP_17(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_17
#define WALBERLA_FOR_ALL_CELLS_OMP_19(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_19
#define WALBERLA_FOR_ALL_CELLS_OMP_20(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_OMP_20

#define WALBERLA_FOR_ALL_CELLS_1(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_1
#define WALBERLA_FOR_ALL_CELLS_2(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_2
#define WALBERLA_FOR_ALL_CELLS_4(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_4
#define WALBERLA_FOR_ALL_CELLS_6(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_6
#define WALBERLA_FOR_ALL_CELLS_8(...)  THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_8
#define WALBERLA_FOR_ALL_CELLS_10(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_10
#define WALBERLA_FOR_ALL_CELLS_12(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_12
#define WALBERLA_FOR_ALL_CELLS_14(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_14
#define WALBERLA_FOR_ALL_CELLS_16(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_16
#define WALBERLA_FOR_ALL_CELLS_18(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_18
#define WALBERLA_FOR_ALL_CELLS_19(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_19
#define WALBERLA_FOR_ALL_CELLS_20(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_ITERATOR_MACRO___WALBERLA_FOR_ALL_CELLS_20


#define WALBERLA_FOR_ALL_CELLS_OMP(...) \
        WALBERLA_MACRO_OVERLOAD( WALBERLA_FOR_ALL_CELLS_OMP_, __VA_ARGS__ )

#define WALBERLA_FOR_ALL_CELLS(...) \
        WALBERLA_MACRO_OVERLOAD( WALBERLA_FOR_ALL_CELLS_, __VA_ARGS__ )
