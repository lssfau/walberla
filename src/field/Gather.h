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
//! \file Gather.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Collect field from all blocks
//
//======================================================================================================================

#pragma once


#include "domain_decomposition/StructuredBlockStorage.h"
#include "core/logging/Logging.h"
#include "core/mpi/Gatherv.h"

namespace walberla {
namespace field {



   //*******************************************************************************************************************
   /*! Gathers (part of) a field on a single process
   *
   *  \param gatheredField  output parameter, on targetRank this field is resized and overwritten with the
   *                        gathered data. On other processes this field is left unchanged
   *  \param blocks         the block storage where the field is stored
   *  \param fieldID        the block data id of the field
   *  \param boundingBox    cell bounding box in global coordinates of the interval which is gathered
   *  \param targetRank     rank of the process where field is gathered, if negative rank is passed, field is gathered
    *                       on all processes
   *  \param comm           MPI communicator
   *
   *  \tparam Field_T       the type of field which is stored in the given fieldID
   *  \tparam ResultField_T the type of the given result field, can be different to Field_T ( for example
   *                        Field_T may be an adaptor class, ResultField_T a normal Field
   */
   //*******************************************************************************************************************
   template<typename Field_T, typename ResultField_T >
   void gather ( ResultField_T & gatheredField,
                 const shared_ptr<StructuredBlockStorage> & blocks, ConstBlockDataID fieldID,
                 CellInterval boundingBox = CellInterval(),  int targetRank = 0, MPI_Comm comm = MPI_COMM_WORLD )
   {
      // If no boundingBox was specified use the whole domain
      if ( boundingBox == CellInterval() ) {
         boundingBox = blocks->getDomainCellBB();
      }

      // -- Packing --
      mpi::SendBuffer sendBuffer;
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         IBlock & block = *blockIt;
         auto field = blockIt->template getData<Field_T>( fieldID );

         // intersection between block cell bounding box and the interval to collect - in global coordinates
         CellInterval intersection = blocks->getBlockCellBB( block );
         intersection.intersect(  boundingBox );


         // Target interval describes the coordinates where this part is located in the collected field
         CellInterval targetInterval = intersection;
         const Cell offset = boundingBox.min();
         targetInterval.min() -= offset;
         targetInterval.max() -= offset;
         sendBuffer << targetInterval;


         // Pack unpack in layout zyxf - to be on the safe side in case the field is stored in different
         // layouts on different blocks
         blocks->transformGlobalToBlockLocalCellInterval( intersection, block );
         for( auto it = field->beginSliceXYZ(intersection); it != field->end(); ++it )
            for( uint_t f = 0; f < Field_T::F_SIZE; ++f )
               sendBuffer << it.getF(f);
      }

      // -- Gather message sizes --
      mpi::RecvBuffer recvBuffer;
      if ( targetRank >= 0 )
         mpi::gathervBuffer( sendBuffer, recvBuffer, targetRank, comm );
      else
         mpi::allGathervBuffer( sendBuffer, recvBuffer, comm );

      // -- Unpacking --
      if ( recvBuffer.size() > 0 )
      {
         gatheredField.resize( boundingBox.size(0), boundingBox.size(1), boundingBox.size(2) );

         while ( ! recvBuffer.isEmpty() )
         {
            CellInterval targetInterval;
            recvBuffer >> targetInterval;

            for( auto it = gatheredField.beginSliceXYZ(targetInterval); it != gatheredField.end(); ++it )
               for( uint_t f = 0; f < Field_T::F_SIZE; ++f )
                  recvBuffer >> it.getF(f);
         }
      }
   }


   template<typename Field_T, typename ResultField_T>
   void gatherSlice( ResultField_T & gatheredField, const shared_ptr<StructuredBlockStorage> & blocks,
                     ConstBlockDataID fieldID, int sliceDim, cell_idx_t sliceValue, int targetRank = 0,
                     MPI_Comm comm = MPI_COMM_WORLD )
   {
      CellInterval boundingBox = blocks->getDomainCellBB();
      boundingBox.min()[ uint_c( sliceDim ) ] = sliceValue;
      boundingBox.max()[ uint_c( sliceDim ) ] = sliceValue;

      gather<Field_T,ResultField_T>( gatheredField, blocks, fieldID, boundingBox, targetRank, comm );
   }

} // namespace field
} // namespace walberla





