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
//! \file BoundarySetter.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "boundary/Boundary.h"
#include "boundary/BoundaryHandling.h"


#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/config/Config.h"

#include "field/FlagFunctions.h"

#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"



namespace walberla {
namespace geometry {
namespace initializer {
   /*
    * Usage Steps:
    *
    * 1) Setup
    *    - setConfigBlock: complete configuration block
    *       Example: a) setting boundary  { NoSlip {} } or { UBB{ vel 0.1; } }
    *                b) setting flag      { flag fluid; } also allowed for boundaries { flag NoSlip; }
    *
    *    - setBoundaryConfigBlock( BoundaryUID, config::block)
    *          setting boundaryUID with configuration as config::block
    *    - setBoundaryConfig ( BoundaryUID, BoundaryConfiguration )
    *          setting boundaryUID with already parsed configuration
    *    - setFlagUID
    *          setting FlagUID ( could also be fluid, or other non-boundary flags )
    *          no configuration is allowed here
    * 2) Configure:
    *    - pass the current  block
    *
    * 3) set() set the boundary with configuration , or a flags
    *
    */
   template< typename BoundaryHandling_T >
   class BoundarySetter
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

      typename BoundaryHandling_T::FlagField * getFlagField() const { return  boundaryHandling_->getFlagField(); }
      //@}
      //****************************************************************************************************************


       static const typename BoundaryHandling_T::FlagField * getFlagField(const IBlock & block, ConstBlockDataID bdId) {
          return block.getData<BoundaryHandling_T >(bdId)->getFlagField();
       }
   private:
       bool flagOrBoundaryCondition_;

       // When setup with boundary ( only valid when flagOrBoundaryCondition==false )
       BoundaryUID boundaryUID_;
       Config::BlockHandle boundaryConfigBlock_;

       // When setup with flag  ( only valid when flagOrBoundaryCondition==true )
       FlagUID flagUID_;

       typename BoundaryHandling_T::flag_t flag_;
       shared_ptr<BoundaryConfiguration> bcConfig_;
       BoundaryHandling_T * boundaryHandling_;
   };

   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::setConfigBlock( const Config::BlockHandle & blockHandle,
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

      if (boundaryConfigBlocks.size () > 1 ) {
         WALBERLA_ABORT_NO_DEBUG_INFO( "More than one boundary configuration sub-block in configuration block " << blockHandle.getKey() );
      }
      const bool hasBoundaryConfigSubBlock = ( boundaryConfigBlocks.size() == 1 );
      const bool hasFlagParameter = blockHandle.isDefined( "flag" );

      if ( hasBoundaryConfigSubBlock && hasFlagParameter ) {
         WALBERLA_ABORT_NO_DEBUG_INFO( "Boundary Block \"" << blockHandle.getKey() << "\" has a flag parameter and a boundary sub-block. Only one of them is allowed." );
      }
      if ( !hasBoundaryConfigSubBlock && !hasFlagParameter ) {
         WALBERLA_ABORT_NO_DEBUG_INFO( "Boundary Block \"" << blockHandle.getKey() << "\" has neither a boundary sub-block nor a flag parameter. "
                                        << "You have to specify which boundary/flag should be set." );
      }

      flagOrBoundaryCondition_ = hasFlagParameter;

      if ( !flagOrBoundaryCondition_ ) {
         boundaryConfigBlock_ = boundaryConfigBlocks.front();
         boundaryUID_ = boundaryConfigBlock_.getKey();
      }
      else
         flagUID_ = FlagUID( blockHandle.getParameter<std::string>("flag") );
   }

   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::setBoundaryConfigBlock( const BoundaryUID & boundaryUID, const Config::BlockHandle & blockHandle )
   {
      flagOrBoundaryCondition_ = false;
      boundaryUID_ = boundaryUID;
      boundaryConfigBlock_ = blockHandle;
      bcConfig_ = BoundaryConfiguration::nullPtr();
   }

   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::setBoundaryConfig( const BoundaryUID & boundaryUID, const shared_ptr<BoundaryConfiguration> & conf )
   {
      flagOrBoundaryCondition_ = false;
      boundaryUID_ = boundaryUID;
      boundaryConfigBlock_ = Config::BlockHandle();
      if ( !conf)
         bcConfig_ = BoundaryConfiguration::nullPtr();
      else
         bcConfig_ = conf;
   }

   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::setFlagUID( const FlagUID & flag )
   {
      flagOrBoundaryCondition_ = true;
      flagUID_ = flag;
   }


   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::configure( IBlock & block, BlockDataID boundaryHandlingID )
   {
      boundaryHandling_ = block.getData<BoundaryHandling_T>( boundaryHandlingID );

      if ( ! flagOrBoundaryCondition_ )
      {
         if ( !boundaryHandling_->containsBoundaryCondition( boundaryUID_ ) )
         {
            WALBERLA_ABORT_NO_DEBUG_INFO( "You tried to set a boundary with UID \"" << boundaryUID_ << "\". " <<
                                          "No boundary with this name was found is registered in this application." );
         }

         if ( boundaryConfigBlock_ )
         {
            bcConfig_ = boundaryHandling_->createBoundaryConfiguration( boundaryUID_, boundaryConfigBlock_ );
            boundaryConfigBlock_ = Config::BlockHandle(); // discard the config block so we don't unnecessarily run createBoundaryConfiguration more than once with the same arguments
         }

         flag_ = boundaryHandling_->getBoundaryMask( boundaryUID_ );
         if ( ! field::isFlag( flag_ ) )
         {
            WALBERLA_ABORT_NO_DEBUG_INFO (
                    "Setting boundary " << boundaryUID_ << " using a configuration file ( or GeometryInitializer)  ) "
                << " is not possible. \n This boundary is specified by a mask instead of a single flag "
                << " i.e. the boundary is handled whenever on of the flags of that mask is set. It is not clear "
                << " which of these flags that are part of the mask should be set by the initializer. You have to "
                << " write custom initialization code in this case. " );
         }
      }
      else
      {
         if ( ! boundaryHandling_->getFlagField()->flagExists( flagUID_ )  ) {
            WALBERLA_ABORT_NO_DEBUG_INFO( "You tried to set a flag named \"" << flagUID_ << "\". There is no such flag registered"
                     "at the BoundaryHandling." );
         }
         else
         {
            bcConfig_ = BoundaryConfiguration::nullPtr();
            flag_ = boundaryHandling_->getFlagField()->getFlag( flagUID_ );
         }
      }
   }

   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::set( cell_idx_t x, cell_idx_t y, cell_idx_t z )
   {
      boundaryHandling_->forceFlag( flag_, x,y,z, *bcConfig_ );
   }

   template<typename BoundaryHandling_T>
   void BoundarySetter<BoundaryHandling_T>::set( const CellInterval & ci )
   {
      boundaryHandling_->forceFlag( flag_, ci, *bcConfig_ );
   }

   template<typename BoundaryHandling_T>
   template< typename CellIterator >
   void BoundarySetter<BoundaryHandling_T>::set( const CellIterator & begin, const CellIterator & end )
   {
      boundaryHandling_->forceFlag( flag_, begin,end, *bcConfig_ );
   }

} // namespace initializer
} // namespace geometry
} // namespace walberla


#include "BoundarySetterFlagFieldSpecialization.h"
