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
//! \file UniformBufferedScheme.h
//! \ingroup blockforest
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LocalCommunicationMode.h"

#include "blockforest/BlockNeighborhoodSection.h"
#include "blockforest/StructuredBlockForest.h"
#include "communication/UniformPackInfo.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/OpenMPBufferSystem.h"
#include "core/selectable/IsSetSelected.h"
#include "core/uid/SUID.h"

#include <functional>
#include <map>
#include <vector>


namespace walberla {
namespace blockforest {
namespace communication {

//*******************************************************************************************************************
/*! Communication scheme for buffered communication in uniform block grids.
*
* Most common use case: Synchronize a set of GhostLayerFields with the neighboring processes
*     - when multiple fields have been changed they can be synchronized at once, using one MPI message per
*       communication partner
*   \code
*      UniformBufferedScheme<stencil::D3Q19> scheme;  // the stencil defines the communication neighbors
*      scheme.addPackInfo( make_shared<field::PackInfo<FieldType> >( idOfFirstField ) );
*      scheme.addPackInfo( make_shared<field::PackInfo<FieldType> >( idOfSecondField ) );
*
*      // either synchronous communication...
*      scheme();
*
*      // .. or asynchronous:
*      scheme.startCommunication();
*      functionWhichDoesNotNeedCommunicatedValues();
*      scheme.wait();
*   \endcode
*
* This scheme sends one message per communication step and neighbor process.
* Therefore all contents that have to be sent, are packed into a buffer before.
* Multiple PackInfos can be registered to send their contents in a single step.
* Another option is to omit the buffering step and send multiple messages.
* This strategy is implemented in blockforest::communication::UniformDirectScheme
*
* When running multiple Schemes concurrently different MPI tags have to be used
* for the schemes: the tag can be passed in the constructor.
*/
//*******************************************************************************************************************
template< typename Stencil_T >
class UniformBufferedScheme
{
public:
   using Stencil = Stencil_T;
   using SendBuffer = mpi::SendBuffer;
   using RecvBuffer = mpi::RecvBuffer;

   using PackInfo = shared_ptr<walberla::communication::UniformPackInfo>;

   using VoidFunction = std::function<void ()>;
   using SendBufferFunction = std::function<void (SendBuffer &)>;

   using CommunicationItemInfo = walberla::communication::UniformPackInfo;

   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{

   explicit UniformBufferedScheme( weak_ptr<StructuredBlockForest> bf,
                                   const int tag = 778 ) // waLBerla = 119+97+76+66+101+114+108+97
      : blockForest_( bf ),
        localMode_( START ),
        bufferSystem_( mpi::MPIManager::instance()->comm(), tag ),
        setupBeforeNextCommunication_( true ),
        communicationInProgress_( false ),
        requiredBlockSelectors_( Set<SUID>::emptySet() ),
        incompatibleBlockSelectors_( Set<SUID>::emptySet() )
   {
      auto forest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" );
      forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
   }

   UniformBufferedScheme( weak_ptr<StructuredBlockForest> bf,
                          const Set<SUID> & requiredBlockSelectors,
                          const Set<SUID> & incompatibleBlockSelectors,
                          const int tag = 778 ) // waLBerla = 119+97+76+66+101+114+108+97
      : blockForest_( bf ),
        localMode_( START ),
        bufferSystem_( mpi::MPIManager::instance()->comm(), tag ),
        setupBeforeNextCommunication_( true ),
        communicationInProgress_( false ),
        requiredBlockSelectors_( requiredBlockSelectors ),
        incompatibleBlockSelectors_( incompatibleBlockSelectors )
   {
      auto forest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to access communication for a block storage object that doesn't exist anymore" );
      forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
   }

   ~UniformBufferedScheme() { wait(); }
   //@}
   //*******************************************************************************************************************

   //** Pack Info Registration *****************************************************************************************
   /*! \name Pack Info Registration */
   //@{
   inline void addPackInfo( const PackInfo & packInfo );
   inline void addDataToCommunicate(  const PackInfo & packInfo ) { addPackInfo( packInfo); }
   //@}
   //*******************************************************************************************************************

   //** Synchronous Communication **************************************************************************************
   /*! \name Synchronous Communication */
   //@{
   inline void operator() () { communicate(); }
   inline void communicate();
   //@}
   //*******************************************************************************************************************


   LocalCommunicationMode localMode() const { return localMode_; }
   inline void setLocalMode( const LocalCommunicationMode & mode );


   //** Asynchronous Communication *************************************************************************************
   /*! \name Asynchronous Communication */
   //@{
   void startCommunication();
   void wait();

   std::function<void()> getCommunicateFunctor();
   std::function<void()> getStartCommunicateFunctor();
   std::function<void()> getWaitFunctor();
   //@}
   //*******************************************************************************************************************

protected:

   static void writeHeader( SendBuffer & buffer, const BlockID & id, const stencil::Direction & dir );
   static void  readHeader( RecvBuffer & buffer,       BlockID & id,       stencil::Direction & dir );

   static void send   ( SendBuffer & buffer, const std::vector< SendBufferFunction > & functions );
          void receive( RecvBuffer & buffer );

   void localBufferPacking  ( const uint_t index, const PackInfo & packInfo, const Block * sender,   const stencil::Direction & dir );
   void localBufferUnpacking( const uint_t index, const PackInfo & packInfo,       Block * receiver, const stencil::Direction & dir );



   weak_ptr<StructuredBlockForest> blockForest_;
   uint_t forestModificationStamp_;

   std::vector< PackInfo > packInfos_;
   
   LocalCommunicationMode localMode_;

   mpi::OpenMPBufferSystem bufferSystem_;

   std::vector< VoidFunction >           localCommunication_;
   std::vector< VoidFunction > threadsafeLocalCommunication_;

   std::vector< VoidFunction >           localCommunicationUnpack_;
   std::vector< VoidFunction > threadsafeLocalCommunicationUnpack_;

   std::vector< SendBuffer > localBuffers_;

   bool setupBeforeNextCommunication_;
   bool communicationInProgress_;

   Set<SUID> requiredBlockSelectors_;
   Set<SUID> incompatibleBlockSelectors_;

}; // class UniformBufferedScheme






template< typename Stencil >
inline void UniformBufferedScheme<Stencil>::addPackInfo( const PackInfo & packInfo )
{
   WALBERLA_ASSERT( !communicationInProgress_, "Cannot add pack info while communication is in progress");

   packInfos_.push_back( packInfo );
   setupBeforeNextCommunication_ = true;
}



template< typename Stencil >
inline void UniformBufferedScheme<Stencil>::communicate()
{
   startCommunication();
   wait();
}



template< typename Stencil >
inline void UniformBufferedScheme<Stencil>::setLocalMode( const LocalCommunicationMode & mode )
{
   if( mode != localMode_ )
   {
      localMode_ = mode;
      setupBeforeNextCommunication_ = true;
   }
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::startCommunication()
{
   if( packInfos_.empty() )
      return;

   auto forest = blockForest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to execute communication for a block storage object that doesn't exist anymore" );

   if ( forest->getBlockForest().getModificationStamp() != forestModificationStamp_ )
      setupBeforeNextCommunication_ = true;

   communicationInProgress_ = true;

   for(auto & packInfo : packInfos_)
      packInfo->beforeStartCommunication();

   bool constantSizes = true;
   bool threadsafeReceive = true;
   for(auto & packInfo : packInfos_)
   {
      if( !packInfo->constantDataExchange() ) constantSizes = false;
      if( !packInfo->threadsafeReceiving()  ) threadsafeReceive = false;
   }

   // Redo setup if a PackInfo has changed its requirements
   if( constantSizes == bufferSystem_.sizeChangesEverytime() )
      setupBeforeNextCommunication_ = true;

   if( threadsafeReceive == bufferSystem_.serialRecvs() )
      setupBeforeNextCommunication_ = true;

   if( setupBeforeNextCommunication_ )
   {
      localCommunication_.clear();
      threadsafeLocalCommunication_.clear();

      localCommunicationUnpack_.clear();
      threadsafeLocalCommunicationUnpack_.clear();

      localBuffers_.clear();

      std::map< uint_t, std::vector< SendBufferFunction > > sendFunctions;

      for( auto it = forest->begin(); it != forest->end(); ++it )
      {
         auto * block = dynamic_cast< Block * >( it.get() );

         if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
            continue;


         for( auto dir = Stencil::beginNoCenter(); dir != Stencil::end(); ++dir )
         {
            const auto neighborIdx = blockforest::getBlockNeighborhoodSectionIndex( *dir );

            if( block->getNeighborhoodSectionSize(neighborIdx) == uint_t(0) )
               continue;

            WALBERLA_ASSERT( block->neighborhoodSectionHasEquallySizedBlock(neighborIdx) );
            WALBERLA_ASSERT_EQUAL( block->getNeighborhoodSectionSize(neighborIdx), uint_t(1) );

            const BlockID nBlockId = block->getNeighborId( neighborIdx, uint_t(0) );

            if( !selectable::isSetSelected( block->getNeighborState( neighborIdx, uint_t(0) ), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
               continue;

            if( block->neighborExistsLocally( neighborIdx, uint_t(0) ) && localMode_ != NO_OPTIMIZATION )
            {
               auto neighbor = dynamic_cast< Block * >( forest->getBlock(nBlockId) );
               WALBERLA_ASSERT_EQUAL( neighbor->getProcess(), block->getProcess() );

               for(auto & packInfo : packInfos_)
               {
                  if( localMode_ == BUFFER )
                  {
                     SendBuffer const buffer;
                     localBuffers_.push_back( buffer );
                     const uint_t index = uint_c( localBuffers_.size() ) - uint_t(1);

                     VoidFunction pack = [this, index, pi=packInfo, block, direction = *dir]() {
                        this->localBufferPacking(index, pi, block, direction);
                     };

                     threadsafeLocalCommunication_.push_back( pack );

                     VoidFunction unpack = [this, index, pi=packInfo, neighbor, direction=*dir]() {
                        this->localBufferUnpacking(index, pi, neighbor, direction);
                     };

                     if( packInfo->threadsafeReceiving() )
                        threadsafeLocalCommunicationUnpack_.push_back( unpack );
                     else
                        localCommunicationUnpack_.push_back( unpack );
                  }
                  else
                  {
                     VoidFunction localCommunicationFunction = [pi=packInfo, block, neighbor, direction = *dir](){
                        pi->communicateLocal(block, neighbor, direction);
                     };

                     if( packInfo->threadsafeReceiving() )
                        threadsafeLocalCommunication_.push_back( localCommunicationFunction );
                     else
                        localCommunication_.push_back( localCommunicationFunction );
                  }
               }
            }
            else
            {
               auto nProcess = block->getNeighborProcess( neighborIdx, uint_t(0) );

               if( !packInfos_.empty() ){
                  auto writeHeader = [bId=nBlockId, direction = *dir](SendBuffer & buf){
                     UniformBufferedScheme<Stencil>::writeHeader(buf, bId, direction);
                  };
                  sendFunctions[ nProcess ].push_back( writeHeader );
               }
                  

               for(auto & packInfo : packInfos_){
                  auto packData = [pi = packInfo, block, direction = *dir](SendBuffer & buf){
                     pi->packData(block, direction, buf);
                  };
                  sendFunctions[ nProcess ].push_back( packData );
               }
                  
            }
         }
      }

      bufferSystem_.clearSendingFunctions();
      bufferSystem_.clearReceivingFunctions();

      bufferSystem_.setReceiverInfo( !constantSizes );
      bufferSystem_.enforceSerialSends( false );
      bufferSystem_.enforceSerialRecvs( !threadsafeReceive );

      for( const auto& sIt : sendFunctions )
      {
         auto sendingFunc = [sfunc = sIt.second](auto & sbuf) { UniformBufferedScheme< Stencil >::send(sbuf, sfunc); };
         bufferSystem_.addSendingFunction  (int_c(sIt.first), sendingFunc );

         auto receivingFunc = [this](auto & rbuf) { this->receive(rbuf); };
         bufferSystem_.addReceivingFunction( int_c(sIt.first), receivingFunc );
      }

      setupBeforeNextCommunication_ = false;
      forestModificationStamp_ = forest->getBlockForest().getModificationStamp();
   }
   
   // MPI

   bufferSystem_.startCommunication();

   // LOCAL
   
   if( localMode_ == START )
   {
      for(auto & function : localCommunication_)
         function();

      const int threadsafeLocalCommunicationSize = int_c( threadsafeLocalCommunication_.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int i = 0; i < threadsafeLocalCommunicationSize; ++i )
         threadsafeLocalCommunication_[uint_c(i)]();
   }
   else if( localMode_ == BUFFER )
   {
      const int threadsafeLocalCommunicationSize = int_c( threadsafeLocalCommunication_.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int i = 0; i < threadsafeLocalCommunicationSize; ++i )
         threadsafeLocalCommunication_[uint_c(i)]();
   }

   for(auto & packInfo : packInfos_)
      packInfo->afterStartCommunication();
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::wait()
{
   if( packInfos_.empty() || !communicationInProgress_ )
      return;

   for(auto & packInfo : packInfos_)
      packInfo->beforeWait();

   // LOCAL

   if( localMode_ == WAIT )
   {
      for(auto & function : localCommunication_)
         function();

      const int threadsafeLocalCommunicationSize = int_c( threadsafeLocalCommunication_.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int i = 0; i < threadsafeLocalCommunicationSize; ++i )
         threadsafeLocalCommunication_[uint_c(i)]();
   }
   else if( localMode_ == BUFFER )
   {
      for(auto & function : localCommunicationUnpack_)
         function();

      const int threadsafeLocalCommunicationUnpackSize = int_c( threadsafeLocalCommunicationUnpack_.size() );
#ifdef _OPENMP
      #pragma omp parallel for schedule(dynamic)
#endif
      for( int i = 0; i < threadsafeLocalCommunicationUnpackSize; ++i )
         threadsafeLocalCommunicationUnpack_[uint_c(i)]();
   }
   
   // MPI

   bufferSystem_.wait();

   for(auto & packInfo : packInfos_)
      packInfo->afterWait();

   communicationInProgress_ = false;
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::writeHeader( SendBuffer & buffer, const BlockID & id, const stencil::Direction & dir )
{
   id.toBuffer( buffer );
   buffer << dir;
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::readHeader( RecvBuffer & buffer, BlockID & id, stencil::Direction & dir )
{
   id.fromBuffer( buffer );
   buffer >> dir;
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::send( SendBuffer & buffer, const std::vector< SendBufferFunction > & functions )
{
   for(auto & function : functions)
      function( buffer );
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::receive( RecvBuffer & buffer )
{
   if( !buffer.isEmpty() )
   {
      auto forest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( forest, "Trying to execute communication for a block storage object that doesn't exist anymore" );

      while( !buffer.isEmpty() )
      {
         BlockID blockID;
         stencil::Direction dir;

         readHeader( buffer, blockID, dir );

         auto block = dynamic_cast< Block * >( forest->getBlock(blockID) );

         for(auto & packInfo : packInfos_)
            packInfo->unpackData( block, stencil::inverseDir[dir], buffer );
      }
   }
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::localBufferPacking( const uint_t index, const PackInfo & packInfo,
                                                         const Block * sender, const stencil::Direction & dir )
{
   WALBERLA_ASSERT_LESS( index, localBuffers_.size() );

   SendBuffer & buffer = localBuffers_[ index ];
   buffer.clear();

   packInfo->packData( sender, dir, buffer );
}



template< typename Stencil >
void UniformBufferedScheme<Stencil>::localBufferUnpacking( const uint_t index, const PackInfo & packInfo,
                                                           Block * receiver, const stencil::Direction & dir )
{
   WALBERLA_ASSERT_LESS( index, localBuffers_.size() );

   SendBuffer & sendBuffer = localBuffers_[ index ];
   RecvBuffer recvBuffer( sendBuffer );

   packInfo->unpackData( receiver, stencil::inverseDir[dir], recvBuffer );
}

template< typename Stencil >
std::function<void()> UniformBufferedScheme<Stencil>::getCommunicateFunctor()
{
   return [this] { communicate(); };
}

template< typename Stencil >
std::function<void()> UniformBufferedScheme<Stencil>::getStartCommunicateFunctor()
{
   return [this] { startCommunication(); };
}

template< typename Stencil >
std::function<void()> UniformBufferedScheme<Stencil>::getWaitFunctor()
{
   return [this] { wait(); };
}

} // namespace communication
} // namespace blockforest
} // namespace walberla
