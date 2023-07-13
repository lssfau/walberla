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
//! \file NonUniformBufferedScheme.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LocalCommunicationMode.h"

#include "NonUniformPackInfo.h"
#include "blockforest/BlockNeighborhoodSection.h"
#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "core/selectable/IsSetSelected.h"
#include "core/uid/SUID.h"


#include <map>
#include <functional>
#include <set>
#include <vector>


namespace walberla::blockforest::communication {


template< typename Stencil >
class NonUniformBufferedScheme
{
public:

   enum INDEX { EQUAL_LEVEL = 0, COARSE_TO_FINE = 1, FINE_TO_COARSE = 2 };

   using SendBuffer = mpi::SendBuffer;
   using RecvBuffer = mpi::RecvBuffer;

   using PackInfo = shared_ptr<blockforest::communication::NonUniformPackInfo>;
   using VoidFunction = std::function<void ()>;
   using SendBufferFunction = std::function<void (SendBuffer &)>;

   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{
   explicit NonUniformBufferedScheme( const weak_ptr<StructuredBlockForest>& bf,
                                      int baseTag = 778 ); // waLBerla = 119+97+76+66+101+114+108+97

   NonUniformBufferedScheme( const weak_ptr<StructuredBlockForest>& bf,
                             const Set<SUID> & requiredBlockSelectors, 
                             const Set<SUID> & incompatibleBlockSelectors,
                             int baseTag = 778 ); // waLBerla = 119+97+76+66+101+114+108+97

   ~NonUniformBufferedScheme();
   //@}
   //*******************************************************************************************************************

   //** Pack Info Registration *****************************************************************************************
   /*! \name Pack Info Registration */
   //@{
   inline void addPackInfo( const PackInfo & packInfo );
   //@}
   //*******************************************************************************************************************

   //** Synchronous Communication **************************************************************************************
   /*! \name Synchronous Communication */
   //@{
   inline void operator() () { communicate(); }

   inline void communicate() {communicateEqualLevel(); communicateCoarseToFine(); communicateFineToCoarse();}

   inline void communicateEqualLevel();
   inline void communicateCoarseToFine();
   inline void communicateFineToCoarse();

   inline void communicateEqualLevel  ( uint_t level );
   inline void communicateCoarseToFine( uint_t fineLevel );
   inline void communicateFineToCoarse( uint_t fineLevel );

   std::function<void()>  communicateEqualLevelFunctor(const uint_t level) {
      return [level, this](){ NonUniformBufferedScheme::communicateEqualLevel(level);};
   }
   std::function<void()>  communicateCoarseToFineFunctor(const uint_t fineLevel) {
      return [fineLevel, this](){ NonUniformBufferedScheme::communicateCoarseToFine(fineLevel);};
   }
   std::function<void()>  communicateFineToCoarseFunctor(const uint_t fineLevel) {
      return [fineLevel, this](){ NonUniformBufferedScheme::communicateFineToCoarse(fineLevel);};
   }
   //@}
   //*******************************************************************************************************************
   

   [[nodiscard]] LocalCommunicationMode localMode() const { return localMode_; }
   inline void setLocalMode( const LocalCommunicationMode & mode );


   //** Asynchronous Communication *************************************************************************************
   /*! \name Asynchronous Communication */
   //@{

   void startCommunication() { startCommunicateEqualLevel(); startCommunicateCoarseToFine(); startCommunicateFineToCoarse(); }
   std::function<void()>  getStartCommunicateFunctor() { return std::bind( &NonUniformBufferedScheme::startCommunication, this ); }

   inline void startCommunicateEqualLevel();
   inline void startCommunicateCoarseToFine();
   inline void startCommunicateFineToCoarse();

   inline void startCommunicateEqualLevel  ( uint_t level );
   inline void startCommunicateCoarseToFine( uint_t fineLevel );
   inline void startCommunicateFineToCoarse( uint_t fineLevel );
   
   void wait() { waitCommunicateEqualLevel(); waitCommunicateCoarseToFine(); waitCommunicateFineToCoarse(); }
   std::function<void() >  getWaitFunctor() { return std::bind( &NonUniformBufferedScheme::wait, this ); }

   inline void waitCommunicateEqualLevel();
   inline void waitCommunicateCoarseToFine();
   inline void waitCommunicateFineToCoarse();

   inline void waitCommunicateEqualLevel  ( uint_t level );
   inline void waitCommunicateCoarseToFine( uint_t fineLevel );
   inline void waitCommunicateFineToCoarse( uint_t fineLevel );
   //@}
   //*******************************************************************************************************************

protected:

   void init();
   void refresh();

   void startCommunicationEqualLevel  ( uint_t index, std::set< uint_t > & participatingLevels );
   void startCommunicationCoarseToFine( uint_t index, uint_t coarsestLevel, uint_t finestLevel );
   void startCommunicationFineToCoarse( uint_t index, uint_t coarsestLevel, uint_t finestLevel );

   void resetBufferSystem( shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem );

   void start( INDEX i, uint_t j );
   void  wait( INDEX i, uint_t j );

   static void writeHeader( SendBuffer & buffer, const BlockID & sender, const BlockID & receiver, const stencil::Direction & dir );
   static void  readHeader( RecvBuffer & buffer,       BlockID & sender,       BlockID & receiver,       stencil::Direction & dir );

   static void send( SendBuffer & buffer, std::vector< SendBufferFunction > & functions );
          void receive( RecvBuffer & buffer );

   void localBufferPacking( INDEX i, uint_t j, uint_t bufferIndex, const PackInfo & packInfo,
                            const Block * sender, const Block * receiver, const stencil::Direction & dir );
   void localBufferUnpacking( INDEX i, uint_t j, uint_t bufferIndex, const PackInfo & packInfo,
                              Block * receiver, const Block * sender, const stencil::Direction & dir );

   [[nodiscard]] bool isAnyCommunicationInProgress() const;



   weak_ptr<StructuredBlockForest> blockForest_;
   uint_t forestModificationStamp_{uint_c(0)};

   std::vector< PackInfo > packInfos_;

   LocalCommunicationMode localMode_;

   int baseTag_;
   std::vector< std::vector< shared_ptr< mpi::OpenMPBufferSystem > > > bufferSystem_;

   std::vector< std::vector< std::vector< VoidFunction > > >           localCommunication_;
   std::vector< std::vector< std::vector< VoidFunction > > > threadsafeLocalCommunication_;

   std::vector< std::vector< std::vector< VoidFunction > > >           localCommunicationUnpack_;
   std::vector< std::vector< std::vector< VoidFunction > > > threadsafeLocalCommunicationUnpack_;

   std::vector< std::vector< std::vector< SendBuffer > > > localBuffers_;

   std::vector< std::vector< char > > setupBeforeNextCommunication_; // cannot be a vector of 'bool' since access by non-const reference is required
   std::vector< std::vector< bool > > communicationInProgress_;

   Set<SUID> requiredBlockSelectors_;
   Set<SUID> incompatibleBlockSelectors_;

}; // class NonUniformBufferedScheme



template< typename Stencil >
NonUniformBufferedScheme<Stencil>::NonUniformBufferedScheme( const weak_ptr<StructuredBlockForest>& bf, const int baseTag )
   : blockForest_( bf ), localMode_( START ), baseTag_( baseTag ),
     requiredBlockSelectors_( Set<SUID>::emptySet() ), incompatibleBlockSelectors_( Set<SUID>::emptySet() )
{
   init();
}



template< typename Stencil >
NonUniformBufferedScheme<Stencil>::NonUniformBufferedScheme( const weak_ptr<StructuredBlockForest>& bf,
                                                             const Set<SUID> & requiredBlockSelectors, 
                                                             const Set<SUID> & incompatibleBlockSelectors,
                                                             const int baseTag /*= 778*/ ) // waLBerla = 119+97+76+66+101+114+108+97
   : blockForest_( bf ), localMode_( START ), baseTag_( baseTag ),
     requiredBlockSelectors_( requiredBlockSelectors ), incompatibleBlockSelectors_( incompatibleBlockSelectors )
{
   init();
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::init()
{
   bufferSystem_.resize( 3 );

             localCommunication_.resize( 3 );
   threadsafeLocalCommunication_.resize( 3 );

             localCommunicationUnpack_.resize( 3 );
   threadsafeLocalCommunicationUnpack_.resize( 3 );

   localBuffers_.resize( 3 );

   setupBeforeNextCommunication_.resize( 3 );
   communicationInProgress_.resize( 3 );

   refresh();
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::refresh()
{
   WALBERLA_ASSERT( !isAnyCommunicationInProgress() )

   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levels = forest->getNumberOfLevels();

   for( uint_t i = 0; i != 3; ++i )
   {
      bufferSystem_[i].clear();
      localCommunication_[i].clear();
      threadsafeLocalCommunication_[i].clear();
      localCommunicationUnpack_[i].clear();
      threadsafeLocalCommunicationUnpack_[i].clear();
      localBuffers_[i].clear();
      setupBeforeNextCommunication_[i].clear();
      communicationInProgress_[i].clear();

      for( uint_t j = 0; j <= levels; ++j )
         bufferSystem_[i].push_back( make_shared< mpi::OpenMPBufferSystem >( mpi::MPIManager::instance()->comm(), baseTag_ + int_c( i * levels + j ) ) );

                localCommunication_[i].resize( levels + uint_t(1) );
      threadsafeLocalCommunication_[i].resize( levels + uint_t(1) );

                localCommunicationUnpack_[i].resize( levels + uint_t(1) );
      threadsafeLocalCommunicationUnpack_[i].resize( levels + uint_t(1) );

      localBuffers_[i].resize( levels + uint_t(1) );

      setupBeforeNextCommunication_[i].resize( levels + uint_t(1), char(1) );
      communicationInProgress_[i].resize( levels + uint_t(1), false );
   }
   
#ifndef NDEBUG
   for(auto & packInfo : packInfos_)
      packInfo->clearBufferSizeCheckMap();
#endif
   
   forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
}



template< typename Stencil >
NonUniformBufferedScheme<Stencil>::~NonUniformBufferedScheme()
{
   for( uint_t i = 0; i != bufferSystem_[EQUAL_LEVEL].size(); ++i )
   {
      wait( EQUAL_LEVEL, i );
      wait( FINE_TO_COARSE, i );
      wait( COARSE_TO_FINE, i );
   }
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::addPackInfo( const PackInfo & packInfo )
{
   if( isAnyCommunicationInProgress() )
   {
      WALBERLA_ABORT( "You may not add a PackInfo to a NonUniformBufferedScheme if any communication is in progress!" )
   }

   packInfos_.push_back( packInfo );

   for( uint_t i = 0; i != 3; ++i )
      for(char & level : setupBeforeNextCommunication_[i])
          level = char(1);
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::communicateEqualLevel()
{
   startCommunicateEqualLevel();
   waitCommunicateEqualLevel();
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::communicateCoarseToFine()
{
   startCommunicateCoarseToFine();
   waitCommunicateCoarseToFine();
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::communicateFineToCoarse()
{
   startCommunicateFineToCoarse();
   waitCommunicateFineToCoarse();
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::communicateEqualLevel( const uint_t level )
{
   startCommunicateEqualLevel( level );
   waitCommunicateEqualLevel( level );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::communicateCoarseToFine( const uint_t fineLevel )
{
   startCommunicateCoarseToFine( fineLevel );
   waitCommunicateCoarseToFine( fineLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::communicateFineToCoarse( const uint_t fineLevel )
{
   startCommunicateFineToCoarse( fineLevel );
   waitCommunicateFineToCoarse( fineLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::setLocalMode( const LocalCommunicationMode & mode )
{
   if( mode != localMode_ )
   {
      localMode_ = mode;

      for( uint_t i = 0; i != 3; ++i )
         for(char & level : setupBeforeNextCommunication_[i])
             level = char(1);
   }
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::startCommunicateEqualLevel()
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levelIndex = forest->getNumberOfLevels();

   if( forestModificationStamp_ != forest->getBlockForest().getModificationStamp() )
      refresh();

   std::set< uint_t > participatingLevels;
   for( uint_t i = 0; i != forest->getNumberOfLevels(); ++i )
      participatingLevels.insert(i);

   startCommunicationEqualLevel( levelIndex, participatingLevels );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::startCommunicateCoarseToFine()
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levelIndex = forest->getNumberOfLevels();

   if( levelIndex == 1 )
      return;

   if( forestModificationStamp_ != forest->getBlockForest().getModificationStamp() )
      refresh();

   const auto coarsestLevel = uint_t(0);
   const uint_t   finestLevel = levelIndex - uint_t(1);

   startCommunicationCoarseToFine( levelIndex, coarsestLevel, finestLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::startCommunicateFineToCoarse()
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levelIndex = forest->getNumberOfLevels();
   
   if( levelIndex == 1 )
      return;

   if( forestModificationStamp_ != forest->getBlockForest().getModificationStamp() )
      refresh();

   const auto coarsestLevel = uint_t(0);
   const uint_t   finestLevel = levelIndex - uint_t(1);

   startCommunicationFineToCoarse( levelIndex, coarsestLevel, finestLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::startCommunicateEqualLevel( const uint_t level )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   WALBERLA_ASSERT_LESS( level, forest->getNumberOfLevels() )

   if( forestModificationStamp_ != forest->getBlockForest().getModificationStamp() )
      refresh();

   std::set< uint_t > participatingLevels;
   participatingLevels.insert( level );

   startCommunicationEqualLevel( level, participatingLevels );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::startCommunicateCoarseToFine( const uint_t fineLevel )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   WALBERLA_ASSERT_GREATER( fineLevel, uint_t(0) )
   WALBERLA_ASSERT_LESS( fineLevel, forest->getNumberOfLevels() )

   if( forestModificationStamp_ != forest->getBlockForest().getModificationStamp() )
      refresh();

   const uint_t coarsestLevel = fineLevel - uint_t(1);
   const uint_t   finestLevel = fineLevel;

   startCommunicationCoarseToFine( fineLevel, coarsestLevel, finestLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::startCommunicateFineToCoarse( const uint_t fineLevel )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   WALBERLA_ASSERT_GREATER( fineLevel, uint_t(0) )
   WALBERLA_ASSERT_LESS( fineLevel, forest->getNumberOfLevels() )

   if( forestModificationStamp_ != forest->getBlockForest().getModificationStamp() )
      refresh();

   const uint_t coarsestLevel = fineLevel - uint_t(1);
   const uint_t   finestLevel = fineLevel;

   startCommunicationFineToCoarse( fineLevel, coarsestLevel, finestLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::waitCommunicateEqualLevel()
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levelIndex = forest->getNumberOfLevels();
   WALBERLA_ASSERT_EQUAL( levelIndex, bufferSystem_[EQUAL_LEVEL].size() - uint_t(1) )
   WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() )

   wait( EQUAL_LEVEL, levelIndex );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::waitCommunicateCoarseToFine()
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levelIndex = forest->getNumberOfLevels();
   WALBERLA_ASSERT_EQUAL( levelIndex, bufferSystem_[COARSE_TO_FINE].size() - uint_t(1) )
   WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() )

   if( levelIndex == 1 )
      return;

   wait( COARSE_TO_FINE, levelIndex );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::waitCommunicateFineToCoarse()
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   const uint_t levelIndex = forest->getNumberOfLevels();
   WALBERLA_ASSERT_EQUAL( levelIndex, bufferSystem_[FINE_TO_COARSE].size() - uint_t(1) )
   WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() )

   if( levelIndex == 1 )
      return;

   wait( FINE_TO_COARSE, levelIndex );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::waitCommunicateEqualLevel  ( const uint_t level )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   WALBERLA_ASSERT_LESS( level, forest->getNumberOfLevels() )
   WALBERLA_ASSERT_EQUAL( forest->getNumberOfLevels(), bufferSystem_[EQUAL_LEVEL].size() - uint_t(1) )
   WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() )

   wait( EQUAL_LEVEL, level );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::waitCommunicateCoarseToFine( const uint_t fineLevel )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   WALBERLA_ASSERT_GREATER( fineLevel, uint_t(0) )
   WALBERLA_ASSERT_LESS( fineLevel, forest->getNumberOfLevels() )
   WALBERLA_ASSERT_EQUAL( forest->getNumberOfLevels(), bufferSystem_[COARSE_TO_FINE].size() - uint_t(1) )
   WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() )

   wait( COARSE_TO_FINE, fineLevel );
}



template< typename Stencil >
inline void NonUniformBufferedScheme<Stencil>::waitCommunicateFineToCoarse( const uint_t fineLevel )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
   WALBERLA_ASSERT_GREATER( fineLevel, uint_t(0) )
   WALBERLA_ASSERT_LESS( fineLevel, forest->getNumberOfLevels() )
   WALBERLA_ASSERT_EQUAL( forest->getNumberOfLevels(), bufferSystem_[FINE_TO_COARSE].size() - uint_t(1) )
   WALBERLA_ASSERT_EQUAL( forestModificationStamp_, forest->getBlockForest().getModificationStamp() )

   wait( FINE_TO_COARSE, fineLevel );
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::startCommunicationEqualLevel( const uint_t index, std::set< uint_t > & participatingLevels )
{
   if( packInfos_.empty() )
      return;

   communicationInProgress_[ EQUAL_LEVEL ][ index ] = true;

   shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem = bufferSystem_[ EQUAL_LEVEL ][ index ];

   std::vector< VoidFunction > &           localCommunication =           localCommunication_[ EQUAL_LEVEL ][ index ];
   std::vector< VoidFunction > & threadsafeLocalCommunication = threadsafeLocalCommunication_[ EQUAL_LEVEL ][ index ];

   std::vector< VoidFunction > &           localCommunicationUnpack =           localCommunicationUnpack_[ EQUAL_LEVEL ][ index ];
   std::vector< VoidFunction > & threadsafeLocalCommunicationUnpack = threadsafeLocalCommunicationUnpack_[ EQUAL_LEVEL ][ index ];

   std::vector< SendBuffer > & localBuffers = localBuffers_[ EQUAL_LEVEL ][ index ];

   char & setupBeforeNextCommunication = setupBeforeNextCommunication_[ EQUAL_LEVEL ][ index ];

   if( setupBeforeNextCommunication == char(1) )
   {
      localCommunication.clear();
      threadsafeLocalCommunication.clear();

      localCommunicationUnpack.clear();
      threadsafeLocalCommunicationUnpack.clear();

      localBuffers.clear();

      std::map< uint_t, std::vector< SendBufferFunction > > sendFunctions;

      auto forest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )

      for( auto it = forest->begin(); it != forest->end(); ++it )
      {
         auto * block = dynamic_cast< Block * >( it.get() );

         if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
            continue;

         if( participatingLevels.find( block->getLevel() ) == participatingLevels.end() )
            continue;

         for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
         {
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

            if( !( block->neighborhoodSectionHasEquallySizedBlock(neighborIdx) ) )
               continue;

            WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize(neighborIdx), uint_t(1) )

            const BlockID & receiverId = block->getNeighborId( neighborIdx, uint_t(0) );

            if( !selectable::isSetSelected( block->getNeighborState( neighborIdx, uint_t(0) ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
               continue;

            if( block->neighborExistsLocally( neighborIdx, uint_t(0) ) )
            {
               auto neighbor = dynamic_cast< Block * >( forest->getBlock(receiverId) );
               WALBERLA_ASSERT_EQUAL( neighbor->getProcess(), block->getProcess() )

               for(auto & packInfo : packInfos_)
               {
                  if( localMode_ == BUFFER )
                  {
                     SendBuffer const buffer;
                     localBuffers.push_back( buffer );
                     const uint_t bufferIndex = uint_c( localBuffers.size() ) - uint_t(1);

                     VoidFunction pack = std::bind( &NonUniformBufferedScheme<Stencil>::localBufferPacking, this,
                                                      EQUAL_LEVEL, index, bufferIndex, std::cref( packInfo ), block, neighbor, *dir );

                     threadsafeLocalCommunication.push_back( pack );

                     VoidFunction unpack = std::bind( &NonUniformBufferedScheme<Stencil>::localBufferUnpacking, this,
                                                        EQUAL_LEVEL, index, bufferIndex, std::cref( packInfo ), neighbor, block, *dir  );

                     if( packInfo->threadsafeReceiving() )
                        threadsafeLocalCommunicationUnpack.push_back( unpack );
                     else
                        localCommunicationUnpack.push_back( unpack );
                  }
                  else
                  {
                     VoidFunction localCommunicationFunction = std::bind( &blockforest::communication::NonUniformPackInfo::communicateLocalEqualLevel,
                                                                            packInfo, block, neighbor, *dir );
                     if( packInfo->threadsafeReceiving() )
                        threadsafeLocalCommunication.push_back( localCommunicationFunction );
                     else
                        localCommunication.push_back( localCommunicationFunction );
                  }
               }
            }
            else
            {
               auto nProcess = block->getNeighborProcess( neighborIdx, uint_t(0) );

               if( !packInfos_.empty() )
                  sendFunctions[ nProcess ].push_back( std::bind( NonUniformBufferedScheme<Stencil>::writeHeader, std::placeholders::_1, block->getId(), receiverId, *dir ) );

               for(auto & packInfo : packInfos_)
                  sendFunctions[ nProcess ].push_back( std::bind( &blockforest::communication::NonUniformPackInfo::packDataEqualLevel, packInfo, block, *dir, std::placeholders::_1 ) );
            }
         }
      }

      resetBufferSystem( bufferSystem );

      for( auto & sender : sendFunctions)
      {
         bufferSystem->addSendingFunction  ( int_c(sender.first), std::bind(  NonUniformBufferedScheme<Stencil>::send, std::placeholders::_1, sender.second ) );
         bufferSystem->addReceivingFunction( int_c(sender.first), std::bind( &NonUniformBufferedScheme<Stencil>::receive, this, std::placeholders::_1 ) );
      }

      setupBeforeNextCommunication = char(0);
   }

   start( EQUAL_LEVEL, index );
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::startCommunicationCoarseToFine( const uint_t index, const uint_t coarsestLevel, const uint_t finestLevel )
{
   if( packInfos_.empty() )
      return;

   communicationInProgress_[ COARSE_TO_FINE ][ index ] = true;

   shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem = bufferSystem_[ COARSE_TO_FINE ][ index ];

   std::vector< VoidFunction > &           localCommunication =           localCommunication_[ COARSE_TO_FINE ][ index ];
   std::vector< VoidFunction > & threadsafeLocalCommunication = threadsafeLocalCommunication_[ COARSE_TO_FINE ][ index ];

   std::vector< VoidFunction > &           localCommunicationUnpack =           localCommunicationUnpack_[ COARSE_TO_FINE ][ index ];
   std::vector< VoidFunction > & threadsafeLocalCommunicationUnpack = threadsafeLocalCommunicationUnpack_[ COARSE_TO_FINE ][ index ];

   std::vector< SendBuffer > & localBuffers = localBuffers_[ COARSE_TO_FINE ][ index ];

   char & setupBeforeNextCommunication = setupBeforeNextCommunication_[ COARSE_TO_FINE ][ index ];

   if( setupBeforeNextCommunication == char(1) )
   {
      localCommunication.clear();
      threadsafeLocalCommunication.clear();

      localCommunicationUnpack.clear();
      threadsafeLocalCommunicationUnpack.clear();

      localBuffers.clear();

      std::map< uint_t, std::vector< SendBufferFunction > > sendFunctions;
      std::set< uint_t > ranksToReceiveFrom;

      auto forest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )
      for( auto it = forest->begin(); it != forest->end(); ++it )
      {
         auto * block = dynamic_cast< Block * >( it.get() );

         if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
            continue;

         const uint_t level = block->getLevel();

         if( level >= coarsestLevel && level < finestLevel ) // potential send
         {
            for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
            {
               const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

               if( !( block->neighborhoodSectionHasSmallerBlocks(neighborIdx) ) )
                  continue;

               for( uint_t n = 0; n != block->getNeighborhoodSectionSize(neighborIdx); ++n )
               {
                  const BlockID & receiverId = block->getNeighborId( neighborIdx, n );

                  if( !selectable::isSetSelected( block->getNeighborState( neighborIdx, n ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
                     continue;

                  if( block->neighborExistsLocally( neighborIdx, n ) )
                  {
                     auto neighbor = dynamic_cast< Block * >( forest->getBlock(receiverId) );
                     WALBERLA_ASSERT_EQUAL( neighbor->getProcess(), block->getProcess() )

                     for(auto & packInfo : packInfos_)
                     {
                        if( localMode_ == BUFFER )
                        {
                           SendBuffer const buffer;
                           localBuffers.push_back( buffer );
                           const uint_t bufferIndex = uint_c( localBuffers.size() ) - uint_t(1);

                           VoidFunction pack = std::bind( &NonUniformBufferedScheme<Stencil>::localBufferPacking, this,
                                                            COARSE_TO_FINE, index, bufferIndex, std::cref( packInfo ), block, neighbor, *dir );

                           threadsafeLocalCommunication.push_back( pack );

                           VoidFunction unpack = std::bind( &NonUniformBufferedScheme<Stencil>::localBufferUnpacking, this,
                                                              COARSE_TO_FINE, index, bufferIndex, std::cref( packInfo ), neighbor, block, *dir  );

                           if( packInfo->threadsafeReceiving() )
                              threadsafeLocalCommunicationUnpack.push_back( unpack );
                           else
                              localCommunicationUnpack.push_back( unpack );
                        }
                        else
                        {
                           VoidFunction localCommunicationFunction = std::bind( &blockforest::communication::NonUniformPackInfo::communicateLocalCoarseToFine,
                                                                                  packInfo, block, neighbor, *dir );
                           if( packInfo->threadsafeReceiving() )
                              threadsafeLocalCommunication.push_back( localCommunicationFunction );
                           else
                              localCommunication.push_back( localCommunicationFunction );
                        }
                     }
                  }
                  else
                  {
                     auto nProcess = block->getNeighborProcess( neighborIdx, n );
                     if( !packInfos_.empty() )
                        sendFunctions[ nProcess ].push_back( std::bind( NonUniformBufferedScheme<Stencil>::writeHeader, std::placeholders::_1, block->getId(), receiverId, *dir ) );

                     for(auto & packInfo : packInfos_)
                        sendFunctions[ nProcess ].push_back( std::bind( &blockforest::communication::NonUniformPackInfo::packDataCoarseToFine, packInfo, block, receiverId, *dir, std::placeholders::_1 ) );
                  }
               }
            }
         }

         if( level > coarsestLevel && level <= finestLevel ) // potential receive
         {
            for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
            {
               const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );
               if( block->neighborhoodSectionHasLargerBlock(neighborIdx) )
               {
                  WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize(neighborIdx), uint_t(1) )
                  if( block->neighborExistsRemotely( neighborIdx, uint_t(0) ) &&
                      selectable::isSetSelected( block->getNeighborState( neighborIdx, 0 ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
                  {
                     ranksToReceiveFrom.insert( block->getNeighborProcess( neighborIdx, uint_t(0) ) );
                  }
               }
            }
         }
      }

      resetBufferSystem( bufferSystem );

      for( auto sender : sendFunctions )
         bufferSystem->addSendingFunction( int_c(sender.first), std::bind(  NonUniformBufferedScheme<Stencil>::send, std::placeholders::_1, sender.second ) );

      for(auto receiver : ranksToReceiveFrom)
         bufferSystem->addReceivingFunction( int_c(receiver), std::bind( &NonUniformBufferedScheme<Stencil>::receive, this, std::placeholders::_1 ) );

      setupBeforeNextCommunication = char(0);
   }

   start( COARSE_TO_FINE, index );
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::startCommunicationFineToCoarse( const uint_t index, const uint_t coarsestLevel, const uint_t finestLevel )
{
   if( packInfos_.empty() )
      return;

   communicationInProgress_[ FINE_TO_COARSE ][ index ] = true;

   shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem = bufferSystem_[ FINE_TO_COARSE ][ index ];

   std::vector< VoidFunction > &           localCommunication =           localCommunication_[ FINE_TO_COARSE ][ index ];
   std::vector< VoidFunction > & threadsafeLocalCommunication = threadsafeLocalCommunication_[ FINE_TO_COARSE ][ index ];

   std::vector< VoidFunction > &           localCommunicationUnpack =           localCommunicationUnpack_[ FINE_TO_COARSE ][ index ];
   std::vector< VoidFunction > & threadsafeLocalCommunicationUnpack = threadsafeLocalCommunicationUnpack_[ FINE_TO_COARSE ][ index ];

   std::vector< SendBuffer > & localBuffers = localBuffers_[ FINE_TO_COARSE ][ index ];

   char & setupBeforeNextCommunication = setupBeforeNextCommunication_[ FINE_TO_COARSE ][ index ];

   if( setupBeforeNextCommunication == char(1) )
   {
      localCommunication.clear();
      threadsafeLocalCommunication.clear();

      localCommunicationUnpack.clear();
      threadsafeLocalCommunicationUnpack.clear();

      localBuffers.clear();

      std::map< uint_t, std::vector< SendBufferFunction > > sendFunctions;
      std::set< uint_t > ranksToReceiveFrom;

      auto forest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )

      for( auto it = forest->begin(); it != forest->end(); ++it )
      {
         auto * block = dynamic_cast< Block * >( it.get() );

         if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
            continue;

         const uint_t level = block->getLevel();

         if( level > coarsestLevel && level <= finestLevel ) // potential send
         {
            for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
            {
               const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

               if( !( block->neighborhoodSectionHasLargerBlock(neighborIdx) ) )
                  continue;

               WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize(neighborIdx), uint_t(1) )

               const BlockID & receiverId = block->getNeighborId( neighborIdx, uint_t(0) );

               if( !selectable::isSetSelected( block->getNeighborState( neighborIdx, uint_t(0) ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
                  continue;

               if( block->neighborExistsLocally( neighborIdx, uint_t(0) ) )
               {
                  auto neighbor = dynamic_cast< Block * >( forest->getBlock(receiverId) );
                  WALBERLA_ASSERT_EQUAL( neighbor->getProcess(), block->getProcess() )

                  for(auto & packInfo : packInfos_)
                  {
                     if( localMode_ == BUFFER )
                     {
                        SendBuffer const buffer;
                        localBuffers.push_back( buffer );
                        const uint_t bufferIndex = uint_c( localBuffers.size() ) - uint_t(1);

                        VoidFunction pack = std::bind( &NonUniformBufferedScheme<Stencil>::localBufferPacking, this,
                                                         FINE_TO_COARSE, index, bufferIndex, std::cref( packInfo ), block, neighbor, *dir );

                        threadsafeLocalCommunication.push_back( pack );

                        VoidFunction unpack = std::bind( &NonUniformBufferedScheme<Stencil>::localBufferUnpacking, this,
                                                           FINE_TO_COARSE, index, bufferIndex, std::cref( packInfo ), neighbor, block, *dir  );

                        if( packInfo->threadsafeReceiving() )
                           threadsafeLocalCommunicationUnpack.push_back( unpack );
                        else
                           localCommunicationUnpack.push_back( unpack );
                     }
                     else
                     {
                        VoidFunction localCommunicationFunction = std::bind( &blockforest::communication::NonUniformPackInfo::communicateLocalFineToCoarse,
                                                                               packInfo, block, neighbor, *dir );
                        if( packInfo->threadsafeReceiving() )
                           threadsafeLocalCommunication.push_back( localCommunicationFunction );
                        else
                           localCommunication.push_back( localCommunicationFunction );
                     }
                  }
               }
               else
               {
                  auto nProcess = block->getNeighborProcess( neighborIdx, uint_t(0) );

                  if( !packInfos_.empty() )
                     sendFunctions[ nProcess ].push_back( std::bind( NonUniformBufferedScheme<Stencil>::writeHeader, std::placeholders::_1, block->getId(), receiverId, *dir ) );

                  for(auto & packInfo : packInfos_)
                     sendFunctions[ nProcess ].push_back( std::bind( &blockforest::communication::NonUniformPackInfo::packDataFineToCoarse, packInfo, block, receiverId, *dir, std::placeholders::_1 ) );
               }
            }
         }

         if( level >= coarsestLevel && level < finestLevel ) // potential receive
         {
            for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
            {
               const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );
               if( block->neighborhoodSectionHasSmallerBlocks(neighborIdx) )
               {
                  for( uint_t n = 0; n != block->getNeighborhoodSectionSize(neighborIdx); ++n )
                     if( block->neighborExistsRemotely( neighborIdx, n ) &&
                         selectable::isSetSelected( block->getNeighborState( neighborIdx, n ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
                     {
                        ranksToReceiveFrom.insert( block->getNeighborProcess( neighborIdx, n ) );
                     }
               }
            }
         }
      }

      resetBufferSystem( bufferSystem );

      for( auto sender : sendFunctions )
         bufferSystem->addSendingFunction( int_c(sender.first), std::bind(  NonUniformBufferedScheme<Stencil>::send, std::placeholders::_1, sender.second ) );

      for(auto receiver : ranksToReceiveFrom)
         bufferSystem->addReceivingFunction( int_c(receiver), std::bind( &NonUniformBufferedScheme<Stencil>::receive, this, std::placeholders::_1 ) );

      setupBeforeNextCommunication = char(0);
   }

   start( FINE_TO_COARSE, index );
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::resetBufferSystem( shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem )
{
   bufferSystem->clearSendingFunctions();
   bufferSystem->clearReceivingFunctions();

   bool constantSizes     = true;
   bool threadsafeReceive = true;
   for(auto & packInfo : packInfos_)
   {
      if( !packInfo->constantDataExchange() ) constantSizes = false;
      if( !packInfo->threadsafeReceiving()  ) threadsafeReceive = false;
   }

   bufferSystem->setReceiverInfo( !constantSizes );
   bufferSystem->enforceSerialSends( false );
   bufferSystem->enforceSerialRecvs( !threadsafeReceive );
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::start( const INDEX i, const uint_t j )
{
   shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem = bufferSystem_[i][j];

   std::vector< VoidFunction > &           localCommunication =           localCommunication_[i][j];
   std::vector< VoidFunction > & threadsafeLocalCommunication = threadsafeLocalCommunication_[i][j];

   // MPI

   bufferSystem->startCommunication();

   // LOCAL

   if( localMode_ == START )
   {
      for(auto & function : localCommunication)
         function();

      const int threadsafeLocalCommunicationSize = int_c( threadsafeLocalCommunication.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int c = 0; c < threadsafeLocalCommunicationSize; ++c )
         threadsafeLocalCommunication[uint_c(c)]();
   }
   else if( localMode_ == BUFFER )
   {
      const int threadsafeLocalCommunicationSize = int_c( threadsafeLocalCommunication.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int c = 0; c < threadsafeLocalCommunicationSize; ++c )
         threadsafeLocalCommunication[uint_c(c)]();
   }
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::wait( const INDEX i, const uint_t j )
{
   if( !communicationInProgress_[i][j] || packInfos_.empty() )
      return;

   shared_ptr< mpi::OpenMPBufferSystem > & bufferSystem = bufferSystem_[i][j];

   std::vector< VoidFunction > &           localCommunication =           localCommunication_[i][j];
   std::vector< VoidFunction > & threadsafeLocalCommunication = threadsafeLocalCommunication_[i][j];

   std::vector< VoidFunction > &           localCommunicationUnpack =           localCommunicationUnpack_[i][j];
   std::vector< VoidFunction > & threadsafeLocalCommunicationUnpack = threadsafeLocalCommunicationUnpack_[i][j];

   // LOCAL

   if( localMode_ == WAIT )
   {
      for(auto & function : localCommunication)
         function();

      const int threadsafeLocalCommunicationSize = int_c( threadsafeLocalCommunication.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int c = 0; c < threadsafeLocalCommunicationSize; ++c )
         threadsafeLocalCommunication[uint_c(c)]();
   }
   else if( localMode_ == BUFFER )
   {
      for(auto & function : localCommunicationUnpack)
         function();

      const int threadsafeLocalCommunicationUnpackSize = int_c( threadsafeLocalCommunicationUnpack.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int c = 0; c < threadsafeLocalCommunicationUnpackSize; ++c )
         threadsafeLocalCommunicationUnpack[uint_c(c)]();
   }

   // MPI

   bufferSystem->wait();

   communicationInProgress_[i][j] = false;
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::writeHeader( SendBuffer & buffer, const BlockID & sender, const BlockID & receiver, const stencil::Direction & dir )
{
   sender.toBuffer( buffer );
   receiver.toBuffer( buffer );
   buffer << dir;
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::readHeader( RecvBuffer & buffer, BlockID & sender, BlockID & receiver, stencil::Direction & dir )
{
   sender.fromBuffer( buffer );
   receiver.fromBuffer( buffer );
   buffer >> dir;
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::send( SendBuffer & buffer, std::vector< SendBufferFunction > & functions )
{
   for(auto & function : functions)
      function( buffer );
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::receive( RecvBuffer & buffer )
{
   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" )

   while( !buffer.isEmpty() )
   {
      BlockID sender;
      BlockID receiver;
      stencil::Direction dir;

      readHeader( buffer, sender, receiver, dir );

      auto block = dynamic_cast< Block * >( forest->getBlock(receiver) );

      const uint_t   senderLevel = forest->getLevelFromBlockId( sender );
      const uint_t receiverLevel = block->getLevel();

      if( senderLevel == receiverLevel )
      {
         for(auto & packInfo : packInfos_)
            packInfo->unpackDataEqualLevel( block, stencil::inverseDir[dir], buffer );
      }
      else if( senderLevel < receiverLevel ) // coarse to fine
      {
         for(auto & packInfo : packInfos_)
            packInfo->unpackDataCoarseToFine( block, sender, stencil::inverseDir[dir], buffer );
      }
      else if( senderLevel > receiverLevel ) // fine to coarse
      {
         for(auto & packInfo : packInfos_)
            packInfo->unpackDataFineToCoarse( block, sender, stencil::inverseDir[dir], buffer );
      }
   }
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::localBufferPacking( const INDEX i, const uint_t j, const uint_t bufferIndex, const PackInfo & packInfo,
                                                            const Block * sender, const Block * receiver, const stencil::Direction & dir )
{
   WALBERLA_ASSERT_LESS( bufferIndex, localBuffers_[i][j].size() )

   SendBuffer & buffer = localBuffers_[i][j][ bufferIndex ];
   buffer.clear();

   if( i == EQUAL_LEVEL )
   {
      packInfo->packDataEqualLevel( sender, dir, buffer );
   }
   else if( i == COARSE_TO_FINE )
   {
      packInfo->packDataCoarseToFine( sender, receiver->getId(), dir, buffer );
   }
   else
   {
      WALBERLA_ASSERT( i == FINE_TO_COARSE )
      packInfo->packDataFineToCoarse( sender, receiver->getId(), dir, buffer );
   }
}



template< typename Stencil >
void NonUniformBufferedScheme<Stencil>::localBufferUnpacking( const INDEX i, const uint_t j, const uint_t bufferIndex, const PackInfo & packInfo,
                                                              Block * receiver, const Block * sender, const stencil::Direction & dir )
{
   WALBERLA_ASSERT_LESS( bufferIndex, localBuffers_[i][j].size() )

   SendBuffer & sendBuffer = localBuffers_[i][j][ bufferIndex ];
   RecvBuffer recvBuffer( sendBuffer );

   if( i == EQUAL_LEVEL )
   {
      packInfo->unpackDataEqualLevel( receiver, stencil::inverseDir[dir], recvBuffer );
   }
   else if( i == COARSE_TO_FINE )
   {
      packInfo->unpackDataCoarseToFine( receiver, sender->getId(), stencil::inverseDir[dir], recvBuffer );
   }
   else
   {
      WALBERLA_ASSERT( i == FINE_TO_COARSE )
      packInfo->unpackDataFineToCoarse( receiver, sender->getId(), stencil::inverseDir[dir], recvBuffer );
   }
}



template< typename Stencil >
bool NonUniformBufferedScheme<Stencil>::isAnyCommunicationInProgress() const
{
   for(const auto & communicationInProgres : communicationInProgress_)
      for(bool inProgress : communicationInProgres)
         if( inProgress )
            return true;

   return false;
}


} // namespace walberla
