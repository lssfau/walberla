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
//! \file BoundaryHandlingCollection.h
//! \ingroup boundary
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Boundary.h"
#include "BoundaryHandling.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/uid/UID.h"
#include "core/uid/UIDGenerators.h"

#include "domain_decomposition/IBlock.h"

#include "field/FlagField.h"

#include <tuple>

#include <ostream>
#include <string>
#include <type_traits>
#include <vector>


namespace walberla {
namespace boundary {



class BHCUIDGenerator : public uid::IndexGenerator< BHCUIDGenerator, uint_t >{};
typedef UID< BHCUIDGenerator > BoundaryHandlingCollectionUID;



template< typename FlagField_T, typename... Handlers > // Handlers: all the boundary handlers that are considered by this boundary handling collection
class BoundaryHandlingCollection
{
public:

   typedef FlagField_T                               FlagField;
   typedef typename FlagField_T::flag_t              flag_t;
   typedef typename FlagField_T::const_base_iterator ConstFlagFieldBaseIterator;



   class BlockSweep {
   public:
      BlockSweep( const BlockDataID & collection, const uint_t numberOfGhostLayersToInclude = 0 ) :
         collection_( collection ), numberOfGhostLayersToInclude_( numberOfGhostLayersToInclude ) {}
      void operator()( IBlock * block ) const {
         BoundaryHandlingCollection * collection = block->getData< BoundaryHandlingCollection >( collection_ );
         (*collection)( numberOfGhostLayersToInclude_ );
      }
   protected:
      const BlockDataID collection_;
      const uint_t numberOfGhostLayersToInclude_;
   };



   BoundaryHandlingCollection( const std::string & identifier, FlagField_T * const flagField, const Handlers & ... boundaryHandlers );

   bool operator==( const BoundaryHandlingCollection &     ) const { WALBERLA_ASSERT( false ); return false; } // For testing purposes, block data items must be comparable with operator "==".
                                                                                                      // Since instances of type "BoundaryHandlingCollection" are registered as block data items,
                                                                                                      // "BoundaryHandlingCollection" must implement operator "==". As of right now, comparing
                                                                                                      // two boundary handling collection instances will always fail... :-) TODO: fixit?
   bool operator!=( const BoundaryHandlingCollection & rhs ) const { return !operator==( rhs ); }

   const BoundaryHandlingCollectionUID & getUID() const { return uid_; }

   const FlagField_T * getFlagField() const { return flagField_; }
         FlagField_T * getFlagField()       { return flagField_; }

   inline bool isEmpty( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline bool isEmpty( const ConstFlagFieldBaseIterator & it ) const;

   inline bool consideredByAllHandlers( const uint_t numberOfGhostLayersToInclude = 0 ) const;               // These functions check if each
   inline bool consideredByAllHandlers( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;  // selected cell is either marked
   inline bool consideredByAllHandlers( const ConstFlagFieldBaseIterator & it ) const;                       // as domain or boundary
   inline bool consideredByAllHandlers( const CellInterval & cells ) const;                                  // in every boundary handling
   template< typename CellIterator >                                                                         // that belongs to this collection.
   inline bool consideredByAllHandlers( const CellIterator & begin, const CellIterator & end ) const;        // <--

   template< typename BoundaryHandling_T >
   inline const BoundaryHandling_T & getBoundaryHandling( const BoundaryHandlingUID & uid ) const; ///< You most likely have to call this function via "collection.template getBoundaryHandling< BoundaryHandling_T >(uid)".
   template< typename BoundaryHandling_T >
   inline       BoundaryHandling_T & getBoundaryHandling( const BoundaryHandlingUID & uid );       ///< You most likely have to call this function via "collection.template getBoundaryHandling< BoundaryHandling_T >(uid)".

   inline uint_t numberOfMatchingHandlers           ( const flag_t flag ) const;
   inline uint_t numberOfMatchingHandlersForDomain  ( const flag_t flag ) const;
   inline uint_t numberOfMatchingHandlersForBoundary( const flag_t flag ) const;

   inline bool containsBoundaryCondition( const BoundaryUID & uid ) const;
   inline bool containsBoundaryCondition( const FlagUID & flag ) const;
   inline bool containsBoundaryCondition( const flag_t flag ) const;

   inline flag_t getBoundaryMask( const BoundaryUID & uid ) const { return getBoundaryMask( boundaryHandlers_, uid ); }

   inline BoundaryUID getBoundaryUID( const FlagUID & flag ) const;
   inline BoundaryUID getBoundaryUID( const flag_t    flag ) const;

   inline shared_ptr<BoundaryConfiguration> createBoundaryConfiguration( const BoundaryUID & uid, const Config::BlockHandle & config ) const
            { return createBoundaryConfiguration( boundaryHandlers_, uid, config ); }

   inline bool checkConsistency( const uint_t numberOfGhostLayersToInclude = 0 ) const;
   inline bool checkConsistency( const CellInterval & cells ) const;

   inline void refresh( const uint_t numberOfGhostLayersToInclude = 0 ); // reset near boundary ...
   inline void refresh( const CellInterval & cells );                    // ... flags for all handlers

   inline void refreshOutermostLayer( cell_idx_t thickness = 1 ); // reset near boundary flags in the outermost "inner" layers

   //** General Flag Handling ******************************************************************************************
   /*! \name General Flag Handling */
   //@{
   inline void setFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void setFlag( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void setFlag( const FlagUID & flag, const CellInterval & cells,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void setFlag( const flag_t    flag, const CellInterval & cells,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   template< typename CellIterator >
   inline void setFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   template< typename CellIterator >
   inline void setFlag( const flag_t    flag, const CellIterator & begin, const CellIterator & end,
                        const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void forceFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void forceFlag( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void forceFlag( const FlagUID & flag, const CellInterval & cells,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   inline void forceFlag( const flag_t    flag, const CellInterval & cells,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   template< typename CellIterator >
   inline void forceFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );
   template< typename CellIterator >
   inline void forceFlag( const flag_t    flag, const CellIterator & begin, const CellIterator & end,
                          const BoundaryConfiguration & parameter = BoundaryConfiguration::null() );

   inline void removeFlag( const FlagUID & flag, const uint_t numberOfGhostLayersToInclude = 0 );
   inline void removeFlag( const flag_t    flag, const uint_t numberOfGhostLayersToInclude = 0 );

   inline void removeFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void removeFlag( const flag_t    flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );

   inline void removeFlag( const FlagUID & flag, const CellInterval & cells );
   inline void removeFlag( const flag_t    flag, const CellInterval & cells );

   template< typename CellIterator >
   inline void removeFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator >
   inline void removeFlag( const flag_t    flag, const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** Clear Cells ****************************************************************************************************
   /*! \name Clear Cells */
   //@{
   inline void clear( const uint_t numberOfGhostLayersToInclude = 0 );
   inline void clear( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void clear( const CellInterval & cells );
   template< typename CellIterator >
   inline void clear( const CellIterator & begin, const CellIterator & end );
   //@}
   //*******************************************************************************************************************

   //** Boundary Treatment *********************************************************************************************
   /*! \name Boundary Treatment */
   //@{
   static BlockSweep getBlockSweep( const BlockDataID handling, const uint_t numberOfGhostLayersToInclude = 0 )
      { return BlockSweep( handling, numberOfGhostLayersToInclude ); }

   inline void operator()( const uint_t numberOfGhostLayersToInclude = 0 );
   inline void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   inline void operator()( const CellInterval & cells );
   template< typename CellIterator >
   inline void operator()( const CellIterator & begin, const CellIterator & end );

   inline void beforeBoundaryTreatment() { beforeBoundaryTreatment( boundaryHandlers_ ); }
   inline void  afterBoundaryTreatment() {  afterBoundaryTreatment( boundaryHandlers_ ); }
   //@}
   //*******************************************************************************************************************

   //** Pack / Unpack boundary handling collection *********************************************************************
   /*! \name Pack / Unpack boundary handling collection */
   //@{
   template< typename Buffer_T >
   void   pack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers = 1, const bool assumeIdenticalFlagMapping = true ) const;
   template< typename Buffer_T >
   void unpack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers = 1, const bool assumeIdenticalFlagMapping = true );
   //@}
   //*******************************************************************************************************************

   inline void        toStream( std::ostream & os ) const;
   inline std::string toString() const;

private:

   CellInterval getGhostLayerCellInterval( const uint_t numberOfGhostLayersToInclude ) const;

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type isEmpty( const HandlersTuple & boundaryHandlers, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type isEmpty( const HandlersTuple &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { return true; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type isEmpty( const HandlersTuple & boundaryHandlers, const ConstFlagFieldBaseIterator & it ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type isEmpty( const HandlersTuple &, const ConstFlagFieldBaseIterator & ) const { return true; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type consideredByAllHandlers( const HandlersTuple & boundaryHandlers, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type consideredByAllHandlers( const HandlersTuple &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { return true; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type consideredByAllHandlers( const HandlersTuple & boundaryHandlers, const ConstFlagFieldBaseIterator & it ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type consideredByAllHandlers( const HandlersTuple &, const ConstFlagFieldBaseIterator & ) const { return true; }

   //** Get Boundary Handling (private helper functions) ***************************************************************
   /*! \name Get Boundary Handling (private helper functions) */
   //@{

   // matching type (-> BoundaryHandling_T) not yet found ...

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), BoundaryHandling_T>::type & getBoundaryHandling( const BoundaryHandlingUID & uid, const HandlersTuple & boundaryHandlers,
                                                          typename std::enable_if< std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::value >::type* /*dummy*/ = 0 ) const
   {
      if( uid == std::get<N>( boundaryHandlers ).getUID() )
         return std::get<N>( boundaryHandlers );
      else
         return getBoundaryHandling_TypeExists< BoundaryHandling_T, HandlersTuple, N-1 >( uid, boundaryHandlers );

   }

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), BoundaryHandling_T>::type & getBoundaryHandling( const BoundaryHandlingUID & uid, const HandlersTuple & boundaryHandlers,
                                                          typename std::enable_if< std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::value >::type* /*dummy*/ = 0 ) const
   {
      if( uid == std::get<N>( boundaryHandlers ).getUID() )
         return std::get<N>( boundaryHandlers );
      else
         WALBERLA_ABORT( "The requested boundary handler " << uid.getIdentifier() << " is not part of this boundary handling collection." );
   }

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), BoundaryHandling_T>::type & getBoundaryHandling( const BoundaryHandlingUID & uid, const HandlersTuple & boundaryHandlers,
                                                          typename std::enable_if< std::is_same< typename std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::type,
                                                                                                     std::false_type >::value >::type* /*dummy*/ = 0 ) const
   {
      return getBoundaryHandling< BoundaryHandling_T, HandlersTuple, N-1 >( uid, boundaryHandlers );
   }

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), BoundaryHandling_T>::type & getBoundaryHandling( const BoundaryHandlingUID & /*uid*/, const HandlersTuple & /*boundaryHandlers*/,
                                                          typename std::enable_if< std::is_same< typename std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::type,
                                                                                                     std::false_type >::value >::type* /*dummy*/ = 0 ) const
   {
      static_assert( sizeof(BoundaryHandling_T) == 0, "The requested boundary handling is not part of this boundary handling collection." );
   }

   // matching type (-> BoundaryHandling_T) exists!

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), BoundaryHandling_T>::type & getBoundaryHandling_TypeExists( const BoundaryHandlingUID & uid, const HandlersTuple & boundaryHandlers,
                                                                     typename std::enable_if< std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::value >::type* /*dummy*/ = 0 ) const
   {
      if( uid == std::get<N>( boundaryHandlers ).getUID() )
         return std::get<N>( boundaryHandlers );
      else
         return getBoundaryHandling_TypeExists< BoundaryHandling_T, HandlersTuple, N-1 >( uid, boundaryHandlers );

   }

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), BoundaryHandling_T>::type & getBoundaryHandling_TypeExists( const BoundaryHandlingUID & uid, const HandlersTuple & boundaryHandlers,
                                                                     typename std::enable_if< std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::value >::type* /*dummy*/ = 0 ) const
   {
      if( uid == std::get<N>( boundaryHandlers ).getUID() )
         return std::get<N>( boundaryHandlers );
      else
         WALBERLA_ABORT( "The requested boundary handler " << uid.getIdentifier() << " is not part of this boundary handling collection." );
   }

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N!=0), BoundaryHandling_T>::type & getBoundaryHandling_TypeExists( const BoundaryHandlingUID & uid, const HandlersTuple & boundaryHandlers,
                                                                     typename std::enable_if< std::is_same< typename std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::type,
                                                                                                                std::false_type >::value >::type* /*dummy*/ = 0 ) const
   {
      return getBoundaryHandling_TypeExists< BoundaryHandling_T, HandlersTuple, N-1 >( uid, boundaryHandlers );
   }

   template< typename BoundaryHandling_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline const typename std::enable_if<(N==0), BoundaryHandling_T>::type & getBoundaryHandling_TypeExists( const BoundaryHandlingUID & uid, const HandlersTuple & /*boundaryHandlers*/,
                                                                     typename std::enable_if< std::is_same< typename std::is_same< BoundaryHandling_T, typename std::tuple_element<N, HandlersTuple>::type >::type,
                                                                                                                std::false_type >::value >::type* /*dummy*/ = 0 ) const
   {
      WALBERLA_ABORT( "The requested boundary handler " << uid.getIdentifier() << " is not part of this boundary handling collection." );
   }

   //*******************************************************************************************************************

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type checkForUniqueBoundaryHandlingUIDs( const HandlersTuple & boundaryHandlers ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type checkForUniqueBoundaryHandlingUIDs( const HandlersTuple & ) const {}

   inline std::vector< BoundaryUID > getBoundaryUIDs() const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type getBoundaryUIDs( const HandlersTuple & boundaryHandlers, std::vector< BoundaryUID > & uids ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type getBoundaryUIDs( const HandlersTuple &, std::vector< BoundaryUID > & ) const {}

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type checkForIdenticalFlagFields( const HandlersTuple & boundaryHandlers ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type checkForIdenticalFlagFields( const HandlersTuple & ) const { return true; }

   inline uint_t numberOfMatchingBoundaryHandlers( const BoundaryHandlingUID & uid ) const;

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), uint_t>::type numberOfMatchingBoundaryHandlers( const HandlersTuple & boundaryHandlers, const BoundaryHandlingUID & uid ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), uint_t>::type numberOfMatchingBoundaryHandlers( const HandlersTuple &, const BoundaryHandlingUID & ) const { return uint_c(0); }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), uint_t>::type numberOfMatchingHandlers( const HandlersTuple & boundaryHandlers, const flag_t flag ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), uint_t>::type numberOfMatchingHandlers( const HandlersTuple &, const flag_t ) const { return uint_c(0); }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), uint_t>::type numberOfMatchingHandlersForDomain( const HandlersTuple & boundaryHandlers, const flag_t flag ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), uint_t>::type numberOfMatchingHandlersForDomain( const HandlersTuple &, const flag_t ) const { return uint_c(0); }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), uint_t>::type numberOfMatchingHandlersForBoundary( const HandlersTuple & boundaryHandlers, const flag_t flag ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), uint_t>::type numberOfMatchingHandlersForBoundary( const HandlersTuple &, const flag_t ) const { return uint_c(0); }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type containsBoundaryCondition( const HandlersTuple & boundaryHandlers, const BoundaryUID & uid ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type containsBoundaryCondition( const HandlersTuple &, const BoundaryUID & ) const { return false; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type containsBoundaryCondition( const HandlersTuple & boundaryHandlers, const flag_t flag ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type containsBoundaryCondition( const HandlersTuple &, const flag_t ) const { return false; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), flag_t>::type getBoundaryMask( const HandlersTuple & boundaryHandlers, const BoundaryUID & uid ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), flag_t>::type getBoundaryMask( const HandlersTuple &, const BoundaryUID & ) const { return numeric_cast<flag_t>(0); }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), BoundaryUID>::type getBoundaryUID( const HandlersTuple & boundaryHandlers, const flag_t flag ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), BoundaryUID>::type getBoundaryUID( const HandlersTuple &, const flag_t ) const;

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), shared_ptr<BoundaryConfiguration>>::type createBoundaryConfiguration( const HandlersTuple & boundaryHandlers,
                                                                         const BoundaryUID & uid, const Config::BlockHandle & config ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), shared_ptr<BoundaryConfiguration>>::type createBoundaryConfiguration( const HandlersTuple &, const BoundaryUID &,
                                                                         const Config::BlockHandle & ) const
                                                                                  { WALBERLA_ASSERT( false ); return make_shared<BoundaryConfiguration>(); }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), bool>::type checkConsistency( const HandlersTuple & boundaryHandlers, const CellInterval & cells ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), bool>::type checkConsistency( const HandlersTuple &, const CellInterval & ) const { return true; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type refresh(       HandlersTuple & boundaryHandlers, const CellInterval & cells );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type refresh( const HandlersTuple &, const CellInterval & ) const {}

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type refreshOutermostLayer(       HandlersTuple & boundaryHandlers, cell_idx_t thickness );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type refreshOutermostLayer( const HandlersTuple &, cell_idx_t ) const {}

   //** General Flag Handling (private helper functions) ***************************************************************
   /*! \name General Flag Handling (private helper functions) */
   //@{
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
          typename std::enable_if<(N!=-1), void>::type setFlag( HandlersTuple & boundaryHandlers, const flag_t flag,
                        const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & parameter );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type setFlag( const HandlersTuple &, const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                        const BoundaryConfiguration & );

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
          typename std::enable_if<(N!=-1), void>::type setFlag( HandlersTuple & boundaryHandlers, const flag_t flag,
                        const CellInterval & cells, const BoundaryConfiguration & parameter );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type setFlag( const HandlersTuple &, const flag_t flag, const CellInterval & cells, const BoundaryConfiguration & );

   void forceFlagHelper( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const BoundaryConfiguration & parameter );

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
          typename std::enable_if<(N!=-1), flag_t>::type flagsToRemove( HandlersTuple & boundaryHandlers, const flag_t flag,
                                const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), flag_t>::type flagsToRemove( const HandlersTuple &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { return 0; }

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
          typename std::enable_if<(N!=-1), void>::type removeFlag( HandlersTuple & boundaryHandlers, const flag_t flag,
                           const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type removeFlag( const HandlersTuple &, const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   //@}
   //*******************************************************************************************************************

   //** Clear Cells (private helper functions) *************************************************************************
   /*! \name Clear Cells (private helper functions) */
   //@{
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
          typename std::enable_if<(N!=-1), flag_t>::type clear( HandlersTuple & boundaryHandlers, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), flag_t>::type clear( const HandlersTuple &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { return 0; }
   //@}
   //*******************************************************************************************************************

   //** Boundary Treatment (private helper functions) ******************************************************************
   /*! \name Boundary Treatment (private helper functions) */
   //@{
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type execute( HandlersTuple & boundaryHandlers, const uint_t numberOfGhostLayersToInclude );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type execute( const HandlersTuple &, const uint_t ) const {}

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type execute( HandlersTuple & boundaryHandlers, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type execute( const HandlersTuple &, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type execute( HandlersTuple & boundaryHandlers, const CellInterval & cells );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type execute( const HandlersTuple &, const CellInterval & ) const {}

   template< typename CellIterator, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type execute( HandlersTuple & boundaryHandlers, const CellIterator & begin, const CellIterator & end );
   template< typename CellIterator, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type execute( const HandlersTuple &, const CellIterator &, const CellIterator & ) const {}

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type beforeBoundaryTreatment( HandlersTuple & boundaryHandlers );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type beforeBoundaryTreatment( const HandlersTuple & ) const {}

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type  afterBoundaryTreatment( HandlersTuple & boundaryHandlers );
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type  afterBoundaryTreatment( const HandlersTuple & ) const {}
   //@}
   //*******************************************************************************************************************

   //** Pack / Unpack boundary handling (private helper functions) *****************************************************
   /*! \name Pack / Unpack boundary handling (private helper functions) */
   //@{
   std::map< std::string, flag_t > getFlagMapping() const { return std::get<0>( boundaryHandlers_ ).getFlagMapping(); }

   template< typename Buffer_T >
   std::vector< flag_t > getNeighborFlagMapping( Buffer_T & buffer, const bool assumeIdenticalFlagMapping, bool & identicalFlagMapping ) const
   {
      return std::get<0>( boundaryHandlers_ ).getNeighborFlagMapping( buffer, assumeIdenticalFlagMapping, identicalFlagMapping );
   }

   inline void translateMask( flag_t & mask, const std::vector< flag_t > & flagMapping ) const;

   inline CellInterval   getPackingInterval( stencil::Direction direction, const uint_t numberOfLayers ) const;
   inline CellInterval getUnpackingInterval( stencil::Direction direction, const uint_t numberOfLayers ) const;

   template< typename Buffer_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type pack( const HandlersTuple & boundaryHandlers, Buffer_T & buffer, const flag_t mask,
                     const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   template< typename Buffer_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type pack( const HandlersTuple &, Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}

   template< typename Buffer_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type unpack( HandlersTuple & boundaryHandlers, Buffer_T & buffer, const flag_t mask,
                       const cell_idx_t x, const cell_idx_t y, const cell_idx_t z );
   template< typename Buffer_T, typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type unpack( const HandlersTuple &, Buffer_T &, const flag_t, const cell_idx_t, const cell_idx_t, const cell_idx_t ) const {}
   //@}
   //*******************************************************************************************************************

   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N!=-1), void>::type toStream( const HandlersTuple & boundaryHandlers, std::ostream & os ) const;
   template< typename HandlersTuple, int N = std::tuple_size<HandlersTuple>::value - 1 >
   inline typename std::enable_if<(N==-1), void>::type toStream( const HandlersTuple &, std::ostream & ) const {}



   const BoundaryHandlingCollectionUID uid_;

   FlagField_T * const flagField_;

   const CellInterval outerBB_;

   typedef std::tuple<Handlers...> Tuple;
   Tuple boundaryHandlers_;

}; // class BoundaryHandlingCollection



template< typename FlagField_T, typename... Handlers >
BoundaryHandlingCollection< FlagField_T, Handlers... >::BoundaryHandlingCollection( const std::string & identifier, FlagField_T * const flagField,
                                                                              const Handlers & ... boundaryHandlers ) :

   uid_( identifier ),
   flagField_( flagField ),
   outerBB_( -cell_idx_c( flagField_->nrOfGhostLayers() ), -cell_idx_c( flagField_->nrOfGhostLayers() ), -cell_idx_c( flagField_->nrOfGhostLayers() ),
             cell_idx_c( flagField_->xSize() + flagField_->nrOfGhostLayers() ) - 1, cell_idx_c( flagField_->ySize() + flagField_->nrOfGhostLayers() ) - 1,
             cell_idx_c( flagField_->zSize() + flagField_->nrOfGhostLayers() ) - 1 ),
   boundaryHandlers_( std::tuple< Handlers... >( boundaryHandlers... ) )
{
   if( flagField_->nrOfGhostLayers() < 1 )
      WALBERLA_ABORT( "The flag field passed to the boundary handling collection\"" << identifier << "\" must contain at least one ghost layer!" );

   if( !checkForIdenticalFlagFields( boundaryHandlers_ ) )
      WALBERLA_ABORT( "The flag field passed to the boundary handling collection\"" << identifier <<
                      "\" must be the same flag field that is registered at all boundary handlers!" );

   checkForUniqueBoundaryHandlingUIDs( boundaryHandlers_ );

   // check for unique boundaries
   std::vector< BoundaryUID > uids = getBoundaryUIDs();
   for( auto uid = uids.begin(); uid != uids.end(); ++uid )
      if( std::count( uids.begin(), uids.end(), *uid ) != 1 )
         WALBERLA_ABORT( "Every boundary condition registered at a boundary handler at the same boundary handling collection must have a unique boundary UID!\n"
                         "The boundary UID \"" << *uid << "\" is not unique for boundary handling collection \"" << uid_.getIdentifier() << "\"." );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::isEmpty( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   return isEmpty( boundaryHandlers_, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::isEmpty( const ConstFlagFieldBaseIterator & it ) const
{
   WALBERLA_ASSERT_EQUAL( it.getField(), flagField_ );
   WALBERLA_ASSERT( outerBB_.contains(it.x(),it.y(),it.z()) );

   return isEmpty( boundaryHandlers_, it );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const uint_t numberOfGhostLayersToInclude ) const
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   return consideredByAllHandlers( cells );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   return consideredByAllHandlers( boundaryHandlers_, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const ConstFlagFieldBaseIterator & it ) const
{
   WALBERLA_ASSERT_EQUAL( it.getField(), flagField_ );
   WALBERLA_ASSERT( outerBB_.contains(it.x(),it.y(),it.z()) );

   return consideredByAllHandlers( boundaryHandlers_, it );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const CellInterval & cells ) const
{
   WALBERLA_ASSERT( outerBB_.contains( cells ) );

   for( auto z = cells.zMin(); z <= cells.zMax(); ++z )
      for( auto y = cells.yMin(); y <= cells.yMax(); ++y )
         for( auto x = cells.xMin(); x <= cells.xMax(); ++x )
            if( !consideredByAllHandlers(x,y,z) ) return false;
   return true;
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const CellIterator & begin, const CellIterator & end ) const
{
   for( auto cell = begin; cell != end; ++cell )
      if( !consideredByAllHandlers( cell->x(), cell->y(), cell->z() ) )
         return false;
   return true;
}



template< typename FlagField_T, typename... Handlers >
template< typename BoundaryHandling_T >
inline const BoundaryHandling_T & BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryHandling( const BoundaryHandlingUID & uid ) const
{
   return getBoundaryHandling< BoundaryHandling_T, std::tuple< Handlers... >>( uid, boundaryHandlers_ );
}



template< typename FlagField_T, typename... Handlers >
template< typename BoundaryHandling_T >
inline BoundaryHandling_T & BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryHandling( const BoundaryHandlingUID & uid )
{
   return const_cast< BoundaryHandling_T & >( static_cast< const BoundaryHandlingCollection * >( this )->template getBoundaryHandling< BoundaryHandling_T >( uid ) );
}



template< typename FlagField_T, typename... Handlers >
inline uint_t BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingHandlers( const flag_t flag ) const
{
   WALBERLA_ASSERT( field::isFlag(flag) );
   return numberOfMatchingHandlers( boundaryHandlers_, flag );
}



template< typename FlagField_T, typename... Handlers >
inline uint_t BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingHandlersForDomain( const flag_t flag ) const
{
   WALBERLA_ASSERT( field::isFlag(flag) );
   return numberOfMatchingHandlersForDomain( boundaryHandlers_, flag );
}



template< typename FlagField_T, typename... Handlers >
inline uint_t BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingHandlersForBoundary( const flag_t flag ) const
{
   WALBERLA_ASSERT( field::isFlag(flag) );
   return numberOfMatchingHandlersForBoundary( boundaryHandlers_, flag );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::containsBoundaryCondition( const BoundaryUID & uid ) const
{
   return containsBoundaryCondition( boundaryHandlers_, uid );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::containsBoundaryCondition( const FlagUID & flag ) const
{
   if( flagField_->flagExists( flag ) )
      return containsBoundaryCondition( flagField_->getFlag( flag ) );
   return false;
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::containsBoundaryCondition( const flag_t flag ) const
{
   return containsBoundaryCondition( boundaryHandlers_, flag );
}



template< typename FlagField_T, typename... Handlers >
inline BoundaryUID BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryUID( const FlagUID & flag ) const
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   return getBoundaryUID( flagField_->getFlag( flag ) );
}



template< typename FlagField_T, typename... Handlers >
inline BoundaryUID BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryUID( const flag_t flag ) const
{
   WALBERLA_ASSERT( field::isFlag( flag ) );
   WALBERLA_ASSERT( flagField_->isRegistered( flag ) );

   return getBoundaryUID( boundaryHandlers_, flag );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::checkConsistency( const uint_t numberOfGhostLayersToInclude ) const
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   return checkConsistency( cells );
}



template< typename FlagField_T, typename... Handlers >
inline bool BoundaryHandlingCollection< FlagField_T, Handlers... >::checkConsistency( const CellInterval & cells ) const
{
   return checkConsistency( boundaryHandlers_, cells );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::refresh( const uint_t numberOfGhostLayersToInclude )
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   refresh( cells );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::refresh( const CellInterval & cells )
{
   refresh( boundaryHandlers_, cells );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::refreshOutermostLayer( cell_idx_t thickness  )
{
   refreshOutermostLayer( boundaryHandlers_, thickness );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                       const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setFlag( flagField_->getFlag( flag ), x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                       const BoundaryConfiguration & parameter )
{
   if( !outerBB_.contains(x,y,z) )
      return;

   WALBERLA_ASSERT( !flagField_->isFlagSet( x, y, z, flag ) );

   setFlag( boundaryHandlers_, flag, x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const FlagUID & flag, const CellInterval & cells,
                                                                       const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setFlag( flagField_->getFlag( flag ), cells, parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const flag_t flag, const CellInterval & cells,
                                                                       const BoundaryConfiguration & parameter )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

#ifndef NDEBUG
   for( auto cell = flagField_->beginSliceXYZ( localCells ); cell != flagField_->end(); ++cell )
      WALBERLA_ASSERT( !field::isFlagSet( cell, flag ) );
#endif

   setFlag( boundaryHandlers_, flag, localCells, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                                                                       const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   setFlag( flagField_->getFlag( flag ), begin, end, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                       const BoundaryConfiguration & parameter )
{
   for( auto cell = begin; cell != end; ++cell )
      setFlag( flag, cell->x(), cell->y(), cell->z(), parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                         const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceFlag( flagField_->getFlag( flag ), x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlag( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                         const BoundaryConfiguration & parameter )
{
   if( !outerBB_.contains(x,y,z) )
      return;

   forceFlagHelper( flag, x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlag( const FlagUID & flag, const CellInterval & cells,
                                                                         const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceFlag( flagField_->getFlag( flag ), cells, parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlag( const flag_t flag, const CellInterval & cells,
                                                                         const BoundaryConfiguration & parameter )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   for( auto z = localCells.zMin(); z <= localCells.zMax(); ++z )
      for( auto y = localCells.yMin(); y <= localCells.yMax(); ++y )
         for( auto x = localCells.xMin(); x <= localCells.xMax(); ++x )
            forceFlagHelper( flag, x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end,
                                                                         const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   forceFlag( flagField_->getFlag( flag ), begin, end, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlag( const flag_t flag, const CellIterator & begin, const CellIterator & end,
                                                                         const BoundaryConfiguration & parameter )
{
   for( auto cell = begin; cell != end; ++cell )
      forceFlag( flag, cell->x(), cell->y(), cell->z(), parameter );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const FlagUID & flag, const uint_t numberOfGhostLayersToInclude )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), numberOfGhostLayersToInclude );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const flag_t flag, const uint_t numberOfGhostLayersToInclude )
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   removeFlag( flag, cells );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const FlagUID & flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), x, y, z );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const flag_t flag, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( !outerBB_.contains(x,y,z) || !flagField_->isFlagSet( x, y, z, flag ) )
      return;

   removeFlag( boundaryHandlers_, flag, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const FlagUID & flag, const CellInterval & cells )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), cells );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const flag_t flag, const CellInterval & cells )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   for( auto z = localCells.zMin(); z <= localCells.zMax(); ++z )
      for( auto y = localCells.yMin(); y <= localCells.yMax(); ++y )
         for( auto x = localCells.xMin(); x <= localCells.xMax(); ++x )
            if( flagField_->isFlagSet( x, y, z, flag ) )
               removeFlag( boundaryHandlers_, flag, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const FlagUID & flag, const CellIterator & begin, const CellIterator & end )
{
   WALBERLA_ASSERT( flagField_->flagExists( flag ) );

   removeFlag( flagField_->getFlag( flag ), begin, end );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const flag_t flag, const CellIterator & begin, const CellIterator & end )
{
   for( auto cell = begin; cell != end; ++cell )
      removeFlag( flag, cell->x(), cell->y(), cell->z() );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::clear( const uint_t numberOfGhostLayersToInclude )
{
   CellInterval cells = getGhostLayerCellInterval( numberOfGhostLayersToInclude );
   clear( cells );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::clear( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   if( !outerBB_.contains(x,y,z) )
      return;

   flagField_->removeMask( x, y, z, clear( boundaryHandlers_, x, y, z ) );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::clear( const CellInterval & cells )
{
   CellInterval localCells( outerBB_ );
   localCells.intersect( cells );

   if( localCells.empty() )
      return;

   for( auto z = localCells.zMin(); z <= localCells.zMax(); ++z )
      for( auto y = localCells.yMin(); y <= localCells.yMax(); ++y )
         for( auto x = localCells.xMin(); x <= localCells.xMax(); ++x )
            flagField_->removeMask( x, y, z, clear( boundaryHandlers_, x, y, z ) );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::clear( const CellIterator & begin, const CellIterator & end )
{
   for( auto cell = begin; cell != end; ++cell )
      clear( cell->x(), cell->y(), cell->z() );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::operator()( const uint_t numberOfGhostLayersToInclude )
{
   execute( boundaryHandlers_, numberOfGhostLayersToInclude );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   execute( boundaryHandlers_, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::operator()( const CellInterval & cells )
{
   execute( boundaryHandlers_, cells );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::operator()( const CellIterator & begin, const CellIterator & end )
{
   execute( boundaryHandlers_, begin, end );
}



template< typename FlagField_T, typename... Handlers >
template< typename Buffer_T >
void BoundaryHandlingCollection< FlagField_T, Handlers... >::pack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers,
                                                             const bool assumeIdenticalFlagMapping ) const
{
#ifdef NDEBUG
   if( !assumeIdenticalFlagMapping )
#endif
      buffer << getFlagMapping();

   CellInterval interval = getPackingInterval( direction, numberOfLayers );

   for( auto z = interval.min()[2]; z <= interval.max()[2]; ++z ) {
      for( auto y = interval.min()[1]; y <= interval.max()[1]; ++y ) {
         for( auto x = interval.min()[0]; x <= interval.max()[0]; ++x )
         {
            const flag_t mask = flagField_->get(x,y,z);
            buffer << mask;
            pack( boundaryHandlers_, buffer, mask, x, y, z );
         }
      }
   }
}



template< typename FlagField_T, typename... Handlers >
template< typename Buffer_T >
void BoundaryHandlingCollection< FlagField_T, Handlers... >::unpack( Buffer_T & buffer, stencil::Direction direction, const uint_t numberOfLayers,
                                                               const bool assumeIdenticalFlagMapping )
{
   bool identicalFlagMapping = false;
   std::vector< flag_t > flagMapping = getNeighborFlagMapping( buffer, assumeIdenticalFlagMapping, identicalFlagMapping ); // neighbor-flag_t -> flag_t

   CellInterval interval = getUnpackingInterval( direction, numberOfLayers );
   clear( interval );

   for( auto z = interval.min()[2]; z <= interval.max()[2]; ++z ) {
      for( auto y = interval.min()[1]; y <= interval.max()[1]; ++y ) {
         for( auto x = interval.min()[0]; x <= interval.max()[0]; ++x )
         {
            flag_t mask;
            buffer >> mask;

            if( !identicalFlagMapping )
               translateMask( mask, flagMapping );

            (*flagField_)(x,y,z) = mask;
            unpack( boundaryHandlers_, buffer, mask, x, y, z );
         }
      }
   }
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::toStream( std::ostream & os ) const {

   os << "========================= BoundaryHandlingCollection =========================\n\n"
      << "Identifier: " << uid_.getIdentifier() << "\n\n"
      << "Included Boundary Handlers:\n\n";

   toStream( boundaryHandlers_, os );

   os << "\n========================= BoundaryHandlingCollection =========================\n";
}



template< typename FlagField_T, typename... Handlers >
inline std::string BoundaryHandlingCollection< FlagField_T, Handlers... >::toString() const {

   std::ostringstream oss;
   toStream( oss );
   return oss.str();
}



template< typename FlagField_T, typename... Handlers >
CellInterval BoundaryHandlingCollection< FlagField_T, Handlers... >::getGhostLayerCellInterval( const uint_t numberOfGhostLayersToInclude ) const
{
   CellInterval cells( -cell_idx_c( numberOfGhostLayersToInclude ),
                       -cell_idx_c( numberOfGhostLayersToInclude ),
                       -cell_idx_c( numberOfGhostLayersToInclude ),
                        cell_idx_c( flagField_->xSize() + numberOfGhostLayersToInclude ) - 1,
                        cell_idx_c( flagField_->ySize() + numberOfGhostLayersToInclude ) - 1,
                        cell_idx_c( flagField_->zSize() + numberOfGhostLayersToInclude ) - 1 );
   return cells;
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::isEmpty( const HandlersTuple & boundaryHandlers,
                                                                       const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   if( !(std::get<N>( boundaryHandlers ).isEmpty(x,y,z)) )
      return false;
   return isEmpty< HandlersTuple, N-1 >( boundaryHandlers, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::isEmpty( const HandlersTuple & boundaryHandlers,
                                                                       const ConstFlagFieldBaseIterator & it ) const
{
   if( !(std::get<N>( boundaryHandlers ).isEmpty(it)) )
      return false;
   return isEmpty< HandlersTuple, N-1 >( boundaryHandlers, it );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const HandlersTuple & boundaryHandlers,
                                                                                         const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   if( std::get<N>( boundaryHandlers ).isEmpty(x,y,z) )
      return false;
   return consideredByAllHandlers< HandlersTuple, N-1 >( boundaryHandlers, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::consideredByAllHandlers( const HandlersTuple & boundaryHandlers,
                                                                                         const ConstFlagFieldBaseIterator & it ) const
{
   if( std::get<N>( boundaryHandlers ).isEmpty(it) )
      return false;
   return consideredByAllHandlers< HandlersTuple, N-1 >( boundaryHandlers, it );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::checkForUniqueBoundaryHandlingUIDs( const HandlersTuple & boundaryHandlers ) const
{
   if( numberOfMatchingBoundaryHandlers( std::get<N>( boundaryHandlers ).getUID() ) != uint_c(1) )
      WALBERLA_ABORT( "Every boundary handler registered at the same boundary handling collection must have a unique boundary handling UID!\n"
                      "The boundary handling UID \"" << std::get<N>( boundaryHandlers ).getUID() << "\" is not unique for boundary handling collection \"" << uid_.getIdentifier() << "\"." );

   checkForUniqueBoundaryHandlingUIDs< HandlersTuple, N-1 >( boundaryHandlers );
}



template< typename FlagField_T, typename... Handlers >
inline std::vector< BoundaryUID > BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryUIDs() const
{
   std::vector< BoundaryUID > uids;
   getBoundaryUIDs( boundaryHandlers_, uids );
   return uids;
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryUIDs( const HandlersTuple & boundaryHandlers,
                                                                               std::vector< BoundaryUID > & uids ) const
{
   std::vector< BoundaryUID > handlerUIDs = std::get<N>( boundaryHandlers ).getBoundaryUIDs();
   uids.insert( uids.end(), handlerUIDs.begin(), handlerUIDs.end() );
   getBoundaryUIDs< HandlersTuple, N-1 >( boundaryHandlers, uids );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::checkForIdenticalFlagFields( const HandlersTuple & boundaryHandlers ) const
{
   return checkForIdenticalFlagFields< HandlersTuple, N-1 >( boundaryHandlers ) &&
          std::get<N>( boundaryHandlers ).getFlagField() == flagField_ && std::get<N>( boundaryHandlers ).outerBB_ == outerBB_;
}



template< typename FlagField_T, typename... Handlers >
inline uint_t BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingBoundaryHandlers( const BoundaryHandlingUID & uid ) const
{
   return numberOfMatchingBoundaryHandlers( boundaryHandlers_, uid );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), uint_t>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingBoundaryHandlers( const HandlersTuple & boundaryHandlers,
                                                                                                  const BoundaryHandlingUID & uid ) const
{
   return ( ( std::get<N>( boundaryHandlers ).getUID() == uid ) ? uint_c(1) : uint_c(0) ) +
          numberOfMatchingBoundaryHandlers< HandlersTuple, N-1 >( boundaryHandlers, uid );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), uint_t>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingHandlers( const HandlersTuple & boundaryHandlers,
                                                                                          const flag_t flag ) const
{
   return ( ( (std::get<N>( boundaryHandlers ).getBoundaryMask() | std::get<N>( boundaryHandlers ).getDomainMask()) & flag ) == flag ? 1 : 0 ) +
          numberOfMatchingHandlers< HandlersTuple, N-1 >( boundaryHandlers, flag );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), uint_t>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingHandlersForDomain( const HandlersTuple & boundaryHandlers,
                                                                                                   const flag_t flag ) const
{
   return ( ( std::get<N>( boundaryHandlers ).getDomainMask() & flag ) == flag ? 1 : 0 ) +
          numberOfMatchingHandlersForDomain< HandlersTuple, N-1 >( boundaryHandlers, flag );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), uint_t>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::numberOfMatchingHandlersForBoundary( const HandlersTuple & boundaryHandlers,
                                                                                                     const flag_t flag ) const
{
   return ( ( std::get<N>( boundaryHandlers ).getBoundaryMask() & flag ) == flag ? 1 : 0 ) +
          numberOfMatchingHandlersForBoundary< HandlersTuple, N-1 >( boundaryHandlers, flag );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::containsBoundaryCondition( const HandlersTuple & boundaryHandlers,
                                                                                         const BoundaryUID & uid ) const
{
   if( std::get<N>( boundaryHandlers ).containsBoundaryCondition( uid ) )
      return true;
   return containsBoundaryCondition< HandlersTuple, N-1 >( boundaryHandlers, uid );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::containsBoundaryCondition( const HandlersTuple & boundaryHandlers,
                                                                                         const flag_t flag ) const
{
   if( std::get<N>( boundaryHandlers ).containsBoundaryCondition( flag ) )
      return true;
   return containsBoundaryCondition< HandlersTuple, N-1 >( boundaryHandlers, flag );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), typename BoundaryHandlingCollection< FlagField_T, Handlers... >::flag_t>::type
   BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryMask( const HandlersTuple & boundaryHandlers,
                                                                      const BoundaryUID & uid ) const
{
   if( std::get<N>( boundaryHandlers ).containsBoundaryCondition(uid) )
      return std::get<N>( boundaryHandlers ).getBoundaryMask(uid);
   return getBoundaryMask< HandlersTuple, N-1 >( boundaryHandlers, uid );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), BoundaryUID>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryUID( const HandlersTuple & boundaryHandlers,
                                                                                     const flag_t flag ) const
{
   const auto & boundaryHandler = std::get<N>( boundaryHandlers );

   if( boundaryHandler.containsBoundaryCondition( flag ) )
   {
      return boundaryHandler.getBoundaryUID( flag );
   }
   else
   {
      return getBoundaryUID< HandlersTuple, N-1 >( boundaryHandlers, flag );
   }
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N==-1), BoundaryUID>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::getBoundaryUID( const HandlersTuple &, const flag_t flag ) const
{
   if( !flagField_->isRegistered( flag ) )
      WALBERLA_ABORT( "The requested flag with value " << flag << " is not registered at the flag field and is not handled "\
                      "by any boundary condition of boundary handling collection " << uid_.getIdentifier() << "!" );

   const FlagUID & flagUID = flagField_->getFlagUID( flag );
   WALBERLA_ABORT( "The requested flag " << flagUID.getIdentifier() << " is not handled by any boundary condition of "\
                   "boundary handling collection" << uid_.getIdentifier() << "!" );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), shared_ptr<BoundaryConfiguration>>::type
   BoundaryHandlingCollection< FlagField_T, Handlers... >::createBoundaryConfiguration( const HandlersTuple & boundaryHandlers,
                                                                                  const BoundaryUID & uid, const Config::BlockHandle & config ) const
{
   if( std::get<N>( boundaryHandlers ).containsBoundaryCondition(uid) )
      return std::get<N>( boundaryHandlers ).createBoundaryConfiguration( uid, config );
   return createBoundaryConfiguration< HandlersTuple, N-1 >( boundaryHandlers, uid, config );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), bool>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::checkConsistency( const HandlersTuple & boundaryHandlers,
                                                                                const CellInterval & cells ) const
{
   return checkConsistency< HandlersTuple, N-1 >( boundaryHandlers, cells ) && std::get<N>( boundaryHandlers ).checkConsistency(cells);
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::refresh( HandlersTuple & boundaryHandlers, const CellInterval & cells )
{
   std::get<N>( boundaryHandlers ).refresh( cells );
   refresh< HandlersTuple, N-1 >( boundaryHandlers, cells );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::refreshOutermostLayer( HandlersTuple & boundaryHandlers,
                                                                                     cell_idx_t thickness )
{
   std::get<N>( boundaryHandlers ).refreshOutermostLayer( thickness );
   refreshOutermostLayer< HandlersTuple, N-1 >( boundaryHandlers, thickness );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( HandlersTuple & boundaryHandlers, const flag_t flag,
                                                                const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   auto & handler = std::get<N>( boundaryHandlers );

   if( ( (handler.getBoundaryMask() | handler.getDomainMask()) & flag ) == flag )
   {
      flagField_->removeFlag( x, y, z, flag );
      handler.setFlag( flag, x, y, z, parameter );
   }

   setFlag< HandlersTuple, N-1 >( boundaryHandlers, flag, x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N==-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const HandlersTuple & /*boundaryHandlers*/, const flag_t flag,
                                                                       const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                       const BoundaryConfiguration & /*parameter*/ )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   flagField_->addFlag( x, y, z, flag );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( HandlersTuple & boundaryHandlers, const flag_t flag,
                                                                const CellInterval & cells, const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( outerBB_.contains(cells) );

   auto & handler = std::get<N>( boundaryHandlers );

   if( ( (handler.getBoundaryMask() | handler.getDomainMask()) & flag ) == flag )
   {
      if( flagField_->isFlagSet( cells.xMin(), cells.yMin(), cells.zMin(), flag ) )
      {
         for( auto cell = flagField_->beginSliceXYZ( cells ); cell != flagField_->end(); ++cell )
            field::removeFlag( cell, flag );
      }

      handler.setFlag( flag, cells, parameter );
   }

   setFlag< HandlersTuple, N-1 >( boundaryHandlers, flag, cells, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N==-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::setFlag( const HandlersTuple & /*boundaryHandlers*/, const flag_t flag,
                                                                       const CellInterval & cells, const BoundaryConfiguration & /*parameter*/ )
{
   WALBERLA_ASSERT( outerBB_.contains(cells) );

   for( auto cell = flagField_->beginSliceXYZ( cells ); cell != flagField_->end(); ++cell )
      field::addFlag( cell, flag );
}



template< typename FlagField_T, typename... Handlers >
void BoundaryHandlingCollection< FlagField_T, Handlers... >::forceFlagHelper( const flag_t flag, const cell_idx_t x, const cell_idx_t y,
                                                                        const cell_idx_t z, const BoundaryConfiguration & parameter )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   flag_t mask = flagsToRemove( boundaryHandlers_, flag, x, y, z );

   WALBERLA_ASSERT( flagField_->isMaskSet( x, y, z, mask ) );

   static const flag_t digits = numeric_cast< flag_t >( std::numeric_limits< flag_t >::digits );
   for( flag_t bit = 0; bit < digits; ++bit )
   {
      flag_t flagToRemove = numeric_cast< flag_t >( static_cast<flag_t>(1) << bit );
      if( ( flagToRemove & mask ) == flagToRemove )
         removeFlag( boundaryHandlers_, flagToRemove, x, y, z );
   }

   setFlag( boundaryHandlers_, flag, x, y, z, parameter );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), typename BoundaryHandlingCollection< FlagField_T, Handlers... >::flag_t>::type
BoundaryHandlingCollection< FlagField_T, Handlers... >::flagsToRemove( HandlersTuple & boundaryHandlers, const flag_t flag,
                                                                 const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   const auto & handler = std::get<N>( boundaryHandlers );
   flag_t mask = numeric_cast<flag_t>( handler.getBoundaryMask() | handler.getDomainMask() );

   if( ( mask & flag ) == flag )
      mask = numeric_cast<flag_t>( mask & flagField_->get(x,y,z) );
   else
      mask = numeric_cast<flag_t>(0);

   return numeric_cast<flag_t>( mask | flagsToRemove< HandlersTuple, N-1 >( boundaryHandlers, flag, x, y, z ) );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( HandlersTuple & boundaryHandlers, const flag_t flag,
                                                                   const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   auto & handler = std::get<N>( boundaryHandlers );

   if( ( (handler.getBoundaryMask() | handler.getDomainMask()) & flag ) == flag )
   {
      flagField_->addFlag( x, y, z, flag );
      handler.removeFlag( flag, x, y, z );
   }

   removeFlag< HandlersTuple, N-1 >( boundaryHandlers, flag, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N==-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::removeFlag( const HandlersTuple & /*boundaryHandlers*/, const flag_t flag,
                                                                          const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   flagField_->removeFlag( x, y, z, flag );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), typename BoundaryHandlingCollection< FlagField_T, Handlers... >::flag_t>::type
BoundaryHandlingCollection< FlagField_T, Handlers... >::clear( HandlersTuple & boundaryHandlers,
                                                         const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   const flag_t before = flagField_->get(x,y,z);
   std::get<N>( boundaryHandlers ).clear(x,y,z);
   const flag_t removedFlags = numeric_cast< flag_t >( before ^ flagField_->get(x,y,z) );
   flagField_->addMask( x, y, z, before );

   return numeric_cast< flag_t >( removedFlags | clear< HandlersTuple, N-1 >( boundaryHandlers, x, y, z ) );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::execute( HandlersTuple & boundaryHandlers,
                                                                       const uint_t numberOfGhostLayersToInclude )
{
   std::get<N>( boundaryHandlers )( numberOfGhostLayersToInclude );
   execute< HandlersTuple, N-1 >( boundaryHandlers, numberOfGhostLayersToInclude );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::execute( HandlersTuple & boundaryHandlers,
                                                                       const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   std::get<N>( boundaryHandlers )(x,y,z);
   execute< HandlersTuple, N-1 >( boundaryHandlers, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::execute( HandlersTuple & boundaryHandlers, const CellInterval & cells )
{
   std::get<N>( boundaryHandlers )( cells );
   execute< HandlersTuple, N-1 >( boundaryHandlers, cells );
}



template< typename FlagField_T, typename... Handlers >
template< typename CellIterator, typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::execute( HandlersTuple & boundaryHandlers,
                                                                       const CellIterator & begin, const CellIterator & end )
{
   std::get<N>( boundaryHandlers )( begin, end );
   execute< CellIterator, HandlersTuple, N-1 >( boundaryHandlers, begin, end );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::beforeBoundaryTreatment( HandlersTuple & boundaryHandlers )
{
   std::get<N>( boundaryHandlers ).beforeBoundaryTreatment();
   beforeBoundaryTreatment< HandlersTuple, N-1 >( boundaryHandlers );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::afterBoundaryTreatment( HandlersTuple & boundaryHandlers )
{
   std::get<N>( boundaryHandlers ).afterBoundaryTreatment();
   afterBoundaryTreatment< HandlersTuple, N-1 >( boundaryHandlers );
}



template< typename FlagField_T, typename... Handlers >
inline void BoundaryHandlingCollection< FlagField_T, Handlers... >::translateMask( flag_t & mask, const std::vector< flag_t > & flagMapping ) const
{
   std::get<0>( boundaryHandlers_ ).translateMask( mask, flagMapping );
}



template< typename FlagField_T, typename... Handlers >
inline CellInterval BoundaryHandlingCollection< FlagField_T, Handlers... >::getPackingInterval( stencil::Direction direction,
                                                                                          const uint_t numberOfLayers ) const
{
   return std::get<0>( boundaryHandlers_ ).getPackingInterval( direction, numberOfLayers );
}



template< typename FlagField_T, typename... Handlers >
inline CellInterval BoundaryHandlingCollection< FlagField_T, Handlers... >::getUnpackingInterval( stencil::Direction direction,
                                                                                            const uint_t numberOfLayers ) const
{
   return std::get<0>( boundaryHandlers_ ).getUnpackingInterval( direction, numberOfLayers );
}



template< typename FlagField_T, typename... Handlers >
template< typename Buffer_T, typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::pack( const HandlersTuple & boundaryHandlers, Buffer_T & buffer,
                                                                    const flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   std::get<0>( boundaryHandlers_ ).pack( buffer, mask, x, y, z );
   pack< Buffer_T, HandlersTuple, N-1 >( boundaryHandlers, buffer, mask, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename Buffer_T, typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::unpack( HandlersTuple & boundaryHandlers, Buffer_T & buffer,
                                                                      const flag_t mask, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   WALBERLA_ASSERT( outerBB_.contains(x,y,z) );

   auto & handler = std::get<N>( boundaryHandlers );

   flag_t flag = (handler.getBoundaryMask() | handler.getDomainMask()) & mask;
   if( flag )
   {
      flagField_->removeFlag( x, y, z, flag );
      handler.unpack( buffer, flag, x, y, z );
   }

   unpack< Buffer_T, HandlersTuple, N-1 >( boundaryHandlers, buffer, mask, x, y, z );
}



template< typename FlagField_T, typename... Handlers >
template< typename HandlersTuple, int N >
inline typename std::enable_if<(N!=-1), void>::type BoundaryHandlingCollection< FlagField_T, Handlers... >::toStream( const HandlersTuple & boundaryHandlers,
                                                                        std::ostream & os ) const
{
   os << std::get<N>( boundaryHandlers );
   toStream< HandlersTuple, N-1 >( boundaryHandlers, os );
}



} // namespace boundary



using boundary::BoundaryHandlingCollection;

template< typename FlagField_T, typename... Handlers >
inline std::ostream & operator<<( std::ostream & os, const BoundaryHandlingCollection< FlagField_T, Handlers... > & bhc ) {

   bhc.toStream( os );
   return os;
}



} // namespace walberla
