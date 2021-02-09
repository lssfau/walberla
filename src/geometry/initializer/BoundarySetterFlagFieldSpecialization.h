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
//! \file Helper.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "field/FlagFunctions.h"
#include "field/FlagField.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "BoundarySetter.h"

namespace walberla {
namespace geometry {
namespace initializer {

   /*
    * Specialization for FlagFields to make all initializers also work with pure FlagFields instead of
    * boundary handlings
    */
   template< typename Flag_T >
   class BoundarySetter< FlagField<Flag_T> >
   {
   public:

      //** Setup        ************************************************************************************************
      /*! \name Setup  */
      //@{
      void setConfigBlock( const Config::BlockHandle & blockHandle,
                           const std::set<std::string> & subblocksToIgnore = std::set<std::string> () );
      void setBoundaryConfigBlock( const BoundaryUID & boundaryUID, const Config::BlockHandle & blockHandle );
      void setBoundaryConfig( const BoundaryUID & boundaryUID, const shared_ptr<BoundaryConfiguration> & conf );
      void setFlagUID( const FlagUID & flag );
      //@}
      //****************************************************************************************************************


      //** Setting Boundaries        ***********************************************************************************
      /*! \name Setting Boundaries  */
      //@{
      void configure( IBlock & block,  BlockDataID boundaryHandlingID );
      void set( cell_idx_t x, cell_idx_t y, cell_idx_t z );
      void set( const CellInterval & ci );

      template<typename CellIterator>
      void set(  const CellIterator & begin, const CellIterator & end );

      FlagField<Flag_T> * getFlagField() const { return flagField_; }
      //@}
      //****************************************************************************************************************

      static const FlagField<Flag_T> * getFlagField(const IBlock & block, ConstBlockDataID  bdId) {
         return block.getData<FlagField< Flag_T > >(bdId);
      }
   private:
       FlagUID flagUID_;
       FlagField<Flag_T> * flagField_;
       Flag_T flag_;
   };

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T> >::setConfigBlock( const Config::BlockHandle & blockHandle,
                                                            const std::set<std::string> & subblocksToIgnore )
   {
      // Check that either there is one sub-block with boundary information  or a flag parameter
      Config::Blocks boundaryConfigBlocks;
      blockHandle.getBlocks ( boundaryConfigBlocks );

      for( auto blockHandleIt = boundaryConfigBlocks.begin(); blockHandleIt != boundaryConfigBlocks.end(); )
      {
         std::string key = blockHandleIt->getKey();
         if ( subblocksToIgnore.find( key ) != subblocksToIgnore.end() )
            boundaryConfigBlocks.erase( blockHandleIt );
         else
            ++blockHandleIt;
      }

      if (!boundaryConfigBlocks.empty() ) {
         WALBERLA_ABORT_NO_DEBUG_INFO( "No boundary setup blocks are allowed when configuring a flag field"
                                               << blockHandle.getKey() );
      }
      const bool hasFlagParameter = blockHandle.isDefined( "flag" );
      if( !hasFlagParameter ) {
         WALBERLA_ABORT_NO_DEBUG_INFO( "Geometry Block \"" << blockHandle.getKey() << ": missing flag parameter " );
      }
      flagUID_ = FlagUID( blockHandle.getParameter<std::string>("flag") );
   }

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T> >::setBoundaryConfigBlock( const BoundaryUID & boundaryUID, const Config::BlockHandle & blockHandle )
   {
      WALBERLA_ABORT("Passed boundary information to an initializer that sets up a pure flag field only");
   }

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T>>::setBoundaryConfig( const BoundaryUID & boundaryUID, const shared_ptr<BoundaryConfiguration> & conf )
   {
      WALBERLA_ABORT("Passed boundary information to an initializer that sets up a pure flag field only");
   }

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T>>::setFlagUID( const FlagUID & flag )
   {
      flagUID_ = flag;
   }

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T>>::configure( IBlock & block, BlockDataID flagFieldID )
   {
      flagField_ = block.getData< FlagField<Flag_T> >( flagFieldID );
      flag_ = flagField_->getOrRegisterFlag( flagUID_ );
   }

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T>>::set( cell_idx_t x, cell_idx_t y, cell_idx_t z )
   {
      flagField_->addFlag( x, y, z, flag_ );
   }

   template<typename Flag_T>
   void BoundarySetter<FlagField<Flag_T>>::set( const CellInterval & ci )
   {
      for( auto it = flagField_->beginSliceXYZ(ci); it != flagField_->end(); ++it )
         field::addFlag(it, flag_);
   }

   template<typename Flag_T>
   template< typename CellIterator >
   void BoundarySetter<FlagField<Flag_T> >::set( const CellIterator & begin, const CellIterator & end )
   {
      for(auto it = begin; it != end; ++it)
         flagField_->addFlag(it->x(), it->y(), it->z(), flag_);
   }

} // namespace initializer
} // namespace geometry
} // namespace walberla


