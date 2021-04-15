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
//! \file DirectionBasedReduceScheme.h
//! \ingroup blockforest
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Scheme for direction based ghost layer communication
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "communication/ReducePackInfo.h"
#include "core/Abort.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"

#include <vector>


namespace walberla {
namespace blockforest {

//**********************************************************************************************************************
/*! Communication Scheme class for ghost layer based communication with the direction neighbors
*
* Most common use case: Reduce a GhostLayerField with all processes in one direction, e.g. summing up all values in one direction
*     - when multiple fields have been changed they can be synchronized at once, using one MPI message per communication partner
*   \code
*      DirectionBasedReduceScheme<stencil::N> scheme;  // the direction defines the communication direction, here from S to N
*      scheme.addPackInfo( make_shared<ReduceFieldPackInfo<std::plus,real_t,1> >( idOfFirstField  ) ); // field is of type real_t with F_SIZE = 1 and reduce operation = plus
*      scheme.addPackInfo( make_shared<ReduceFieldPackInfo<std::plus,real_t,3> >( idOfSecondField ) ); // field is of type real_t with F_SIZE = 3 and reduce operation = plus
*
*      // synchronous communication...
*      scheme();
*   \endcode
*/
//**********************************************************************************************************************
template< stencil::Direction dir_ >
class DirectionBasedReduceScheme
{
static_assert(
   dir_ == stencil::B || dir_ == stencil::T || dir_ == stencil::N ||
   dir_ == stencil::S || dir_ == stencil::E || dir_ == stencil::W,
   "DirectionBasedReduceScheme only implemented for directions of D3Q6 stencil: N,S,E,W,B,T" );

public:
   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{
   explicit DirectionBasedReduceScheme( const shared_ptr<StructuredBlockForest> & bf )
      : blockForest_(bf),
        bufferSystem_( make_shared<mpi::BufferSystem>( MPIManager::instance()->comm() ) )
   {
      init();
   }

   explicit DirectionBasedReduceScheme( const shared_ptr<StructuredBlockForest> & bf,
                                        const shared_ptr<mpi::BufferSystem> & bs )
      : blockForest_(bf), bufferSystem_(bs){ init(); }
   //@}
   //*******************************************************************************************************************

   //** Pack Info Registration *****************************************************************************************
   /*! \name Pack Info Registration */
   //@{
   size_t addPackInfo( const shared_ptr<walberla::communication::ReducePackInfo> & pi );
   //@}
   //*******************************************************************************************************************

   //** Synchronous Communication **************************************************************************************
   /*! \name Synchronous Communication */
   //@{
   inline void operator() (){ communicate(); }
   inline void communicate();
   //@}
   //*******************************************************************************************************************

protected:
   shared_ptr<StructuredBlockForest>                        blockForest_;
   shared_ptr<mpi::BufferSystem>                            bufferSystem_;
   std::vector< shared_ptr<walberla::communication::ReducePackInfo> > packInfos_;

   std::vector< Block* >                 localBlocks_;
   std::vector< std::vector< Block*  > > localCopy_;
   std::vector< std::vector< BlockID > > remoteSend_;

   /// Uses PackInfo's of the given step to copy data into the MPI Buffers
   void copyIntoMPIBuffers();

   /// Communicates data from MPI buffers using a BufferSystem
   void mpiCommunication(){ bufferSystem_->sendAll(); }

   /// Copies received data back into own ghost layers
   void copyFromMPIBuffers();

   void readHeader ( mpi::RecvBuffer & buf,       BlockID & id ) const { id.fromBuffer( buf ); }
   void writeHeader( mpi::SendBuffer & buf, const BlockID & id ) const { id.toBuffer  ( buf ); }

private:
   void init();
};




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////     PackInfo Registration        ///////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< stencil::Direction dir_ >
size_t DirectionBasedReduceScheme<dir_>::addPackInfo( const shared_ptr<walberla::communication::ReducePackInfo>& pi )
{
   packInfos_.push_back( pi );
   return packInfos_.size() - 1u;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Synchronous Communication
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< stencil::Direction dir_ >
void DirectionBasedReduceScheme<dir_>::communicate()
{
   copyIntoMPIBuffers();
   mpiCommunication();
   copyFromMPIBuffers();
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Private member functions
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
template< stencil::Direction dir_ >
void DirectionBasedReduceScheme<dir_>::init()
{
   if( !blockForest_->containsGlobalBlockInformation() )
      WALBERLA_ABORT( "For this communication scheme global block information is needed" );

   //  -------------------  Block Loop ------------------------------------------------
   uint_t x,y,z;
   int nx,ny,nz;

   blockForest_->getBlocks(localBlocks_,0);

   localCopy_.resize( localBlocks_.size() );
   remoteSend_.resize( localBlocks_.size() );

   std::set < int > remoteRecv;

   for(size_t b = 0u; b<localBlocks_.size(); ++b)
   {
      Block * block = localBlocks_[b];
      blockForest_->getBlockForest().getRootBlockCoordinates( x, y, z, block->getId() );

      // Loop over all neighbors in direction dir_
      for( int step = 1; true; ++step )
      {
         nx = int_c(x) + step * stencil::cx[dir_];
         ny = int_c(y) + step * stencil::cy[dir_];
         nz = int_c(z) + step * stencil::cz[dir_];

         if( nx < 0 || ny < 0 || nz < 0 ||
               !blockForest_->getBlockForest().rootBlockExists( uint_c(nx), uint_c(ny), uint_c(nz) ) ) break;

         BlockID nBlockId;
         blockForest_->getBlockForest().getRootBlockID( nBlockId, uint_c(nx), uint_c(ny), uint_c(nz) );

         if ( blockForest_->getBlockForest().blockExistsLocally(nBlockId) )
         {
            localCopy_[b].push_back( blockForest_->getBlockForest().getBlock(nBlockId) );
         }
         else
         {
            remoteSend_[b].push_back( nBlockId );
         }
      }

      // Loop over all neighbors in inverse direction dir_

      for( int step = 1; true; ++step )
      {
         nx = int_c(x) - step * stencil::cx[dir_];
         ny = int_c(y) - step * stencil::cy[dir_];
         nz = int_c(z) - step * stencil::cz[dir_];

         if( nx < 0 || ny < 0 || nz < 0 ||
               !blockForest_->getBlockForest().rootBlockExists( uint_c(nx), uint_c(ny), uint_c(nz) ) ) break;

         BlockID nBlockId;
         blockForest_->getBlockForest().getRootBlockID( nBlockId, uint_c(nx), uint_c(ny), uint_c(nz) );

         if( !blockForest_->getBlockForest().blockExistsLocally(nBlockId) )
         {
            uint_t rank;
            blockForest_->getBlockForest().getRootBlockProcessRank( rank, uint_c(nx), uint_c(ny), uint_c(nz) );
            remoteRecv.insert( int_c(rank) );
         }
      }
   }

   //TODO Martin to Matthias: is size changing from step to step?
   // Assumed yes to be safe
   bufferSystem_->setReceiverInfo( remoteRecv, true );
}


template< stencil::Direction dir_ >
void DirectionBasedReduceScheme<dir_>::copyIntoMPIBuffers()
{
   //  -------------------  Block Loop ------------------------------------------------
   for(size_t b = 0u; b < localBlocks_.size(); ++b)
   {
      for( size_t nb = 0u; nb < localCopy_[b].size(); ++nb )
      {
         for(size_t s = 0u; s < packInfos_.size(); ++s)
         {
            packInfos_[s]->communicateLocal( localBlocks_[b], localCopy_[b][nb], dir_ );
         }
      }

      for( size_t nb = 0u; nb < remoteSend_[b].size(); ++nb )
      {
         bool headerWritten = false;
         for(size_t s = 0u; s < packInfos_.size(); ++s)
         {
            //Packing of MPI Buffer
            uint_t rank,nx,ny,nz;
            blockForest_->getBlockForest().getRootBlockCoordinates( nx, ny, nz, remoteSend_[b][nb] );
            blockForest_->getBlockForest().getRootBlockProcessRank( rank, nx, ny, nz );
            mpi::SendBuffer & buffer = bufferSystem_->sendBuffer( rank );

            if( !headerWritten )
            {
               writeHeader( buffer, remoteSend_[b][nb] );
               headerWritten = true;
            }

            // Call Packer
            packInfos_[s]->packData( localBlocks_[b], dir_, buffer );
         }
      }
   }
}



template< stencil::Direction dir_ >
void DirectionBasedReduceScheme<dir_>::copyFromMPIBuffers()
{
   for(size_t s = 0u; s < packInfos_.size(); ++s){
      packInfos_[s]->reset();
   }

   //  --------  Loop over all incoming messages  -------------------------------------
   for( auto i = bufferSystem_->begin() ; i != bufferSystem_->end(); ++i )
   {
      mpi::RecvBuffer & buffer = i.buffer();
      while( ! buffer.isEmpty() )
      {
         // Parse Header
         BlockID rBlockID;
         readHeader(buffer, rBlockID);

         // Get receiving block
         Block * rBlock = blockForest_->getBlockForest().getBlock(rBlockID);

         for(size_t s = 0u; s < packInfos_.size(); ++s)
         {
            // Call Unpacker
            packInfos_[s]->unpackData( rBlock, stencil::inverseDir[dir_], buffer );
         }
      }
   }
}


} // namespace blockforest
} // namespace walberla
