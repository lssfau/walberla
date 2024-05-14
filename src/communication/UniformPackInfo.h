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
//! \file UniformPackInfo.h
//! \ingroup communication
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/IBlock.h"
#include "stencil/Directions.h"


namespace walberla {
namespace communication {



/**
 * \brief Data packing/unpacking for ghost layer based communication of a field.
 *
 * Encapsulate information on how to extract data from blocks that should be
 * communicated to neighboring blocks (see \ref packData())
 * and how to inject this data in a receiving block (see \ref unpackData()).
 * This involves a memory buffer and two memory copy operations.
 *
 * A special method exists for communication between two blocks which are
 * allocated on the same process (see \ref communicateLocal()).
 * In this case the data does not have be communicated via a buffer,
 * but can be copied directly.
 *
 * Data that is packed in direction "dir" at one block is unpacked in
 * direction "stencil::inverseDir[dir]" at the neighboring block. This
 * behavior must be implemented in \ref communicateLocal()!
 *
 * \ingroup communication
 */
class UniformPackInfo
{
public:

   //**Construction & Destruction***************************************************************************************
   /*! \name Construction & Destruction */
   //@{
            UniformPackInfo() = default;
   virtual ~UniformPackInfo() = default;
   //@}
   //*******************************************************************************************************************

   /**
    * Should return true if the amount of data that is packed for a given block in direction
    * "dir" is guaranteed to remain constant over time. False otherwise.
    * If you are not sure what to return, return false! Returning false is always safe.
    * Falsely returning true will lead to errors! However, if the data can be guaranteed to remain
    * constant over time, returning true enables performance optimizations during the communication.
    */
   virtual bool constantDataExchange() const = 0;

   /**
    * Must return false if calling \ref unpackData and/or \ref communicateLocal is not thread-safe.
    * True otherwise.
    * If you are not sure what to return, return false! Returning false is always safe.
    * Falsely returning true will most likely lead to errors! However, if both \ref unpackData AND
    * \ref communicateLocal are thread-safe, returning true can lead to performance improvements.
    */
   virtual bool threadsafeReceiving() const = 0;

   /**
    * \brief Pack data from a block into a send buffer.
    *
    * Must be thread-safe! Calls \ref packDataImpl.
    *
    * @param sender     the block whose data should be packed into a buffer
    * @param dir        pack data for neighbor in this direction
    * @param buffer     buffer for writing the data into
    *
    */
   inline void packData( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const;

   /**
    * \brief Unpack received Data.
    *
    * If NOT thread-safe, \ref threadsafeReceiving must return false!
    *
    * @param receiver the block where the unpacked data should be stored into
    * @param dir      receive data from neighbor in this direction
    * @param buffer   buffer for reading the data from
    */
   virtual void unpackData( IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer ) = 0;

   /**
    * \brief Copy data from one local block to another local block.
    *
    * Both blocks are allocated on the current process.
    * If NOT thread-safe, \ref threadsafeReceiving must return false!
    *
    * @param sender    id of block where the data should be copied from
    * @param receiver  id of block where the data should be copied to
    * @param dir       the direction of the communication ( from  sender to receiver )
    */
   virtual void communicateLocal( const IBlock * sender, IBlock * receiver, stencil::Direction dir ) = 0;

   /**
   * This function is called once before the communication is started by the UniformBufferedScheme
   */
   virtual void beforeStartCommunication() { };

   /**
   * This function is called once after the communication has been started by the UniformBufferedScheme
   */
   virtual void  afterStartCommunication() { };

   /**
   * This function is called once before the UniformBufferedScheme waits for the communication to finish
   */
   virtual void beforeWait() { };

   /**
   * This function is called once after the communication has been finished by the UniformBufferedScheme
   */
   virtual void afterWait() { };

protected:

   /**
    * \brief Pack data from a block into a send buffer.
    *
    * Must be thread-safe!
    *
    * @param sender     the block whose data should be packed into a buffer
    * @param dir        pack data for neighbor in this direction
    * @param buffer     buffer for writing the data into
    *
    */
   virtual void packDataImpl( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const = 0;

#ifndef NDEBUG
   mutable std::map< const IBlock *, std::map< stencil::Direction, size_t > > bufferSize_;
#endif
};



inline void UniformPackInfo::packData( const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & buffer ) const
{
#ifndef NDEBUG
   size_t const sizeBefore = buffer.size();
#endif

   packDataImpl( sender, dir, buffer );

#ifndef NDEBUG
   size_t const sizeAfter = buffer.size();
   if( constantDataExchange() )
   {
#ifdef _OPENMP
      #pragma omp critical (UniformPackInfo_packData) // -> packData must be thread-safe!
      {
#endif
      auto & blockMap = bufferSize_[ sender ];
      auto dirEntry = blockMap.find( dir );
      if( dirEntry == blockMap.end() )
         blockMap[ dir ] = sizeAfter - sizeBefore;
      else
         WALBERLA_ASSERT_EQUAL( blockMap[ dir ], (sizeAfter - sizeBefore) )
#ifdef _OPENMP
      }
#endif
   }
#endif
}



} // namespace communication
} // namespace walberla
