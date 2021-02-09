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
//! \file HandlingPackInfo.h
//! \ingroup boundary
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "communication/UniformPackInfo.h"

#include "core/debug/Debug.h"


namespace walberla {
namespace boundary {



template< typename Handling_T >
class HandlingPackInfo : public communication::UniformPackInfo
{
public:

   HandlingPackInfo( const BlockDataID & bdId, const bool assumeIdenticalFlagMapping = true, const uint_t numberOfLayers = 0 ) :
      bdId_( bdId ), numberOfLayers_( numberOfLayers ), assumeIdenticalFlagMapping_( assumeIdenticalFlagMapping ) {}

   ~HandlingPackInfo() override = default;

   bool constantDataExchange() const override { return false; }
   bool threadsafeReceiving()  const override { return false; }

   void unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) override;

   void communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir ) override;

protected:

   void packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const override;



   const BlockDataID bdId_;

   const uint_t numberOfLayers_;
   const bool assumeIdenticalFlagMapping_;

}; // class HandlingPackInfo



template< typename Handling_T >
void HandlingPackInfo< Handling_T >::unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer )
{
   Handling_T * handling = receiver->getData< Handling_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( handling );

   const uint_t layers = ( numberOfLayers_ == 0 ) ? handling->getFlagField()->nrOfGhostLayers() : numberOfLayers_;

   handling->unpack( buffer, dir, layers, assumeIdenticalFlagMapping_ );
}



template< typename Handling_T >
void HandlingPackInfo< Handling_T >::communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir )
{
   mpi::SendBuffer sBuffer;
   packData( sender, dir, sBuffer );
   mpi::RecvBuffer rBuffer( sBuffer );
   unpackData( receiver, stencil::inverseDir[dir], rBuffer );
}



template< typename Handling_T >
void HandlingPackInfo< Handling_T >::packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
   const Handling_T * handling = sender->getData< Handling_T >( bdId_ );
   WALBERLA_ASSERT_NOT_NULLPTR( handling );

   const uint_t layers = ( numberOfLayers_ == 0 ) ? handling->getFlagField()->nrOfGhostLayers() : numberOfLayers_;

   handling->pack( buffer, dir, layers, assumeIdenticalFlagMapping_ );
}



} // namespace boundary
} // namespace walberla
