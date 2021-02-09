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
//! \file Boundary.h
//! \ingroup boundary
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "BoundaryUID.h"

#include "core/debug/Debug.h"
#include "core/DataTypes.h"

#include <type_traits>

namespace walberla {
namespace boundary {


/// Base class for all boundary configuration parameters

class BoundaryConfiguration {
public:
   virtual ~BoundaryConfiguration() = default;
   static const BoundaryConfiguration& null()                { return *boundaryNullPtr; }
   static const shared_ptr<BoundaryConfiguration> nullPtr()  { return  boundaryNullPtr; }
private:
   static const shared_ptr<BoundaryConfiguration> boundaryNullPtr;
};



//**********************************************************************************************************************
/*!
*   \brief Base class for all boundary classes
*
*   Every boundary class must be derived from this class and implement the following concept/member functions (for some
*   exemplary implementations look at classes NoSlip, FreeSlip, UBB, or SimpleUBB in module lbm/boundary):
*
*   1. "static shared_ptr<BoundaryConfiguration> createConfiguration( const Config::BlockHandle & config )"
*      This function is used to extract parameters from a configuration file and pass them to the boundary condition.
*      Meaning: This function defines how the parameter specification for this boundary condition must look like in the
*      configuration file.
*      If your boundary condition is stateless/does not need parameters, you can implement this function as follows:
*         "{ return make_shared<BoundaryConfiguration>(); }"
*      If you do need to pass parameters, these parameters must be implemented in terms of a derived class of
*      "BoundaryConfiguration" (as an example, see class "UBB").
*   2. "void pushFlags( std::vector< FlagUID >& uids )"
*      This function receives (by reference) a vector and must insert (using push_back!) all FlagUIDs into this vector
*      that mark cells as boundary cells that must be treated by this boundary condition.
*   3. "void beforeBoundaryTreatment()"
*      This function is called once before the boundary handler starts the boundary treatment. Of course, this function
*      is called every time the boundary treatment of the boundary handler is triggered (normally once per time step).
*   4. "void afterBoundaryTreatment()"
*      Just like "beforeBoundaryTreatment", this function is called once after the boundary handler has finished the
*      boundary treatment.
*   5. "template< typename Buffer_T >
*       void packCell( Buffer_T & buffer, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )"
*      This function is called every time a boundary cell that is handled by this boundary class is serialized. If
*      the boundary condition stores additional parameters for the boundary cell (x,y,z) these parameters must be
*      serialized and stored in the buffer "buffer". [5) serializes and 6) deserializes, both functions must match]
*   6. "template< typename Buffer_T >
*       void registerCell( Buffer_T & buffer,
*                          const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )"
*      This function is called every time a boundary cell that is handled by this boundary class is deserialized. The
*      flag which was set (must be part of this->mask_!), the cell (x,y,z), and a buffer that potentially stores
*      parameters are passed. [5) serializes and 6) deserializes, both functions must match]
*   7. "void registerCell( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
*                          const BoundaryConfiguration & parameter )"
*      This function is called every time a boundary that is handled by this boundary class is set at the boundary
*      handler. The flag which was set (must be part of this->mask_!), the cell (x,y,z), and a corresponding parameter
*      are passed.
*   8. "void registerCells( const flag_t flag, const CellInterval & cells, const BoundaryConfiguration & parameter )"
*      Just like "registerCell", only this function is called if a boundary that is handled by this boundary class is
*      set for each cell in a specific cell interval (-> "cells").
*   9. "template< typename CellIterator >
*       void registerCells( const flag_t flag, const CellIterator & begin, const CellIterator & end,
*       const BoundaryConfiguration & parameter )"
*      Just like the previous two functions, only this function is called if multiple boundary cells that are handled
*      by this boundary class are set at the boundary handler using CellIterators (= iterators the dereference to type
*      "Cell" - examples: iterators of classes CellVector, CellSet, ...).
*  10. "void unregisterCell( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )"
*      This function is called every time a boundary that is handled by this boundary class is removed at the boundary
*      handler. The flag which was removed (must be part of this->mask_!) and the corresponding cell (x,y,z) are passed.
*  11. "void treatDirection( const cell_idx_t  x, const cell_idx_t  y, const cell_idx_t  z, const stencil::Direction dir,
*                            const cell_idx_t nx, const cell_idx_t ny, const cell_idx_t nz, const flag_t mask )"
*      This function is called during boundary treatment, the parameters that are passed are as follows:
*      Coming from cell (x,y,z) and going in direction "dir" the cell (nx,ny,nz) was reached. Cell (nx,ny,nz) is a
*      boundary cell that is handled by this boundary class. All the flags that are currently set at cell (nx,ny,nz) in
*      the flag field are passed via "mask".
*  12. If the treatDirection member function is thread-safe, i.e., if this function can be called in parallel for
*      different input parameters, you have to add the following constant static boolean member variable to your class:
*      "static const bool threadsafe = true;"
*
*   Every boundary class must call the constructor of "Boundary" which expects a BoundaryUID (all boundary conditions
*   registered at the same boundary handler must have unique boundary UIDs).
*
*/
//**********************************************************************************************************************

template< typename flag_t >
class Boundary {
public:

   static_assert( std::is_unsigned<flag_t>::value, "You are trying to instantiate walberla::boundary::Boundary with "
                                                     "a flag_t which is not an unsigned integer!" );

#ifndef NDEBUG
   Boundary( const BoundaryUID & boundaryUID ) : boundaryUID_( boundaryUID ), mask_(0), maskSet_(false) {}

   void   setMask( const flag_t mask ) { WALBERLA_ASSERT( !maskSet_ ); mask_ = mask; maskSet_ = true; } // called by BoundaryHandler
   flag_t getMask() const { WALBERLA_ASSERT( maskSet_ ); return mask_; }
#else
   Boundary( const BoundaryUID & boundaryUID ) : boundaryUID_( boundaryUID ), mask_(0) {}

   void   setMask( const flag_t flag ) { mask_ = flag; } // called by BoundaryHandler
   flag_t getMask() const { return mask_; }
#endif

   const BoundaryUID & getUID() const { return boundaryUID_; }

protected:

   const BoundaryUID boundaryUID_;

   flag_t mask_; // can be just one flag, but can also be a mask that contains multiple flags
                 // If part of this mask is set for a specific cell, this boundary class/condition is responsible for the corresponding boundary treatment.

#ifndef NDEBUG
   bool maskSet_; // only used in debug mode!
#endif

}; // class Boundary



template< typename Boundary_T, class Enable = void >
struct isThreadSafe
{
   static const bool value = false;
};

template< typename Boundary_T >
struct isThreadSafe< Boundary_T, typename std::enable_if< Boundary_T::threadsafe >::type >
{
   static const bool value = Boundary_T::threadsafe;
};



} // namespace boundary

using boundary::BoundaryConfiguration;
using boundary::Boundary;

} // namespace walberla
