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
//! \file LatticeModelBase.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "CollisionModel.h"
#include "ForceModel.h"
#include "core/DataTypes.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace lbm {



//**********************************************************************************************************************
/*!
*   \brief Base class that is supposed (but not required) to be used for defining specific lattice models
*
*   The concept each lattice model has to follow (for an actual implementation see for example lbm::D3Q19) is
*   as follows:
*   1. It must define a type "Stencil" that refers to a stencil located in module "stencil"
*   2. It must define a type "CommunicationStencil" that is used to determine which neighbors are involved during
*      communication
*   3. There must be a static array named 'w' of size 'Stencil::Size' that contains the weights for each direction.
*      If you want to retrieve the weight that corresponds to a specific direction 'dir' you first have to map this
*      direction to the appropriate array index by using the underlying stencil: the weight of direction 'dir' is equal
*      to w[ Stencil::idx[dir] ]!
*      In Addition to 'w', there must also be a static array 'wInv' containing the inverse values of 'w'.
*   4. It must define a type "CollisionModel". An object of type "CollisionModel" must be returned when calling the
*      member function 'collisionModel()' - this function also needs to be implemented. Each "CollisionModel" must
*      define a type "tag" which may evaluate to SRT_tag, TRT_tag, or MRT_tag (see "CollisionModel.h"). For an exemplary
*      implementation of such a "CollisionModel" see classes lbm::collision_model::SRT or lbm::collision_model::TRT.
*   5. It must define a type "ForceModel". An object of type "ForceModel" must be returned when calling the
*      member function 'forceModel()' - this function also needs to be implemented. For details on force models
*      see \ref docForceModel in 'ForceModel.h'.
*   6. There must be a static boolean value named 'compressible' that must be set to false for incompressible lattice
*      models, otherwise it must be set to true.
*   7. There must be a member function "configure( IBlock & block, StructuredBlockStorage & sbs )" that returns nothing
*      and takes a block and a structured block storage as arguments. Everytime a PDF field is assigned to a specific
*      block, the "configure" function is called for the lattice model that is stored within this PDF field
*      (see lbm/field/AddToStorage.h).
*/
//**********************************************************************************************************************

template< typename CollisionModel_T, bool Compressible, typename ForceModel_T, int EquilibriumAccuracyOrder = 2 >
class LatticeModelBase
{
public:

   using CollisionModel = CollisionModel_T;
   using ForceModel = ForceModel_T;

   static const bool compressible = Compressible;
   static const int  equilibriumAccuracyOrder = EquilibriumAccuracyOrder;



   LatticeModelBase( const CollisionModel_T & cm, const ForceModel_T & fm ) :
      collisionModel_( cm ), forceModel_( fm ) {}

   virtual ~LatticeModelBase() = default;

   virtual void pack( mpi::SendBuffer & buffer ) const
   {
      collisionModel_.pack( buffer );
      forceModel_.pack( buffer );
   }

   virtual void unpack( mpi::RecvBuffer & buffer )
   {
      collisionModel_.unpack( buffer );
      forceModel_.unpack( buffer );
   }

   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      collisionModel_.configure( block, sbs );
      forceModel_.configure( block, sbs );

      config( block, sbs );
   }

   const CollisionModel_T & collisionModel() const { return collisionModel_; }
         CollisionModel_T & collisionModel()       { return collisionModel_; }

   const ForceModel_T & forceModel() const { return forceModel_; }
         ForceModel_T & forceModel()       { return forceModel_; }

protected:

   virtual void config( IBlock & block, StructuredBlockStorage & sbs ) = 0;
   
   

   CollisionModel_T collisionModel_;
   ForceModel_T     forceModel_;
};



} // namespace lbm
} // namespace walberla



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T,    // Element type of SendBuffer
          typename G,    // Growth policy of SendBuffer
          typename CM, bool CO, typename FM, int EQU >
mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const walberla::lbm::LatticeModelBase<CM,CO,FM,EQU> & lm )
{
   buffer.addDebugMarker( "lm" );
   lm.pack( buffer );
   return buffer;
}

template< typename T,    // Element type  of RecvBuffer
          typename CM, bool CO, typename FM, int EQU >
mpi::GenericRecvBuffer<T> & operator>>( mpi::GenericRecvBuffer<T> & buffer, walberla::lbm::LatticeModelBase<CM,CO,FM,EQU> & lm )
{
   buffer.readDebugMarker( "lm" );
   lm.unpack( buffer );
   return buffer;
}

template< typename CM, bool CO, typename FM, int EQU >
struct BufferSizeTrait< walberla::lbm::LatticeModelBase<CM,CO,FM,EQU> > { static const bool constantSize = false; };

}
}
