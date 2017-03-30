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
//! \file CellGatherPackInfo.impl.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "CellGatherPackInfo.h"


namespace walberla {
namespace gather {


template< typename Field_T, typename CC>
CellGatherPackInfo<Field_T,CC>::CellGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                                    ConstBlockDataID fieldID,
                                                    const CC & containerOfGlobalCells,
                                                    const shared_ptr<DataProcessor> & dp )
     : blocks_(bs), fieldID_( fieldID ), dataProcessor_(dp)
{
   for ( auto blockIt = bs->begin(); blockIt != bs->end(); ++blockIt )
   {
      uint_t cellNr = 0;
      Samples & s = localSamplePoints_[ & (*blockIt) ];
      for( auto cellIt = containerOfGlobalCells.begin(); cellIt != containerOfGlobalCells.end(); ++ cellIt )
      {
         Cell localCell;
         bs->transformGlobalToBlockLocalCell( localCell, *blockIt, *cellIt );
         if ( localCell[0] >= 0 &&
              localCell[1] >= 0 &&
              localCell[2] >= 0 &&
              localCell[0] < cell_idx_c( bs->getNumberOfXCells( *blockIt  ) ) &&
              localCell[1] < cell_idx_c( bs->getNumberOfYCells( *blockIt  ) ) &&
              localCell[2] < cell_idx_c( bs->getNumberOfZCells( *blockIt  ) )    )
         {
            s.globalCellNr.push_back( cellNr );
            s.cells.push_back( localCell );
         }
         cellNr++;
      }
   }
}

template< typename Field_T, typename CC>
void CellGatherPackInfo<Field_T,CC>::packData  ( const IBlock * sender,
                                                 mpi::SendBuffer & outBuffer )
{
   auto i = localSamplePoints_.find( sender );

   if(i == localSamplePoints_.end() ) //nothing to send
      return;

   const std::vector<uint_t> & cellNumbers = i->second.globalCellNr;
   const CellVector & cells                = i->second.cells;
   WALBERLA_ASSERT_EQUAL( cells.size(), cellNumbers.size() );

   auto field = sender->getData<Field_T>( fieldID_ );

   const size_t fieldSize = Field_T::F_SIZE;

   outBuffer << cells.size();
   outBuffer << fieldSize;

   // Write data for each component of the field
   for(unsigned int j=0; j< cellNumbers.size(); ++j)
   {
      outBuffer << cellNumbers[j];
      for(uint_t f=0; f< fieldSize; ++f )
         outBuffer << real_c( field->get( cells[j][0], cells[j][1], cells[j][2], f ) );
   }

   // Message format ( example with fieldSize=#f=4 and  points.size()=#Points=2  )
   // #Points |  #f  | t0 (f0 f1 f2 f3) | t1 (f0 f1 f2 f3)
}


template< typename Field_T, typename CC>
void CellGatherPackInfo<Field_T,CC>::unpackData( mpi::RecvBuffer & buffer )
{
   size_t nrPoints;
   buffer >> nrPoints;
   size_t fieldSize;
   buffer >> fieldSize;

   for( size_t i=0; i< nrPoints; ++i )
   {
      receivedData.push_back(std::vector<real_t>(fieldSize+1)); //+1 because we also store t value as first entry
      std::vector<real_t> & pointVec = receivedData[receivedData.size()-1];

      uint_t t;
      real_t val;

      buffer >> t;
      pointVec[0] = real_c( t );
      for( size_t f=0; f<fieldSize; ++f )
      {
         buffer >> val;
         pointVec[f+1] = val;
      }
   }
}

template< typename Field_T, typename CC>
void CellGatherPackInfo<Field_T,CC>::gatherFinished()
{
   //sort according to "t" value, which is the first entry of the inner vector
   std::sort(receivedData.begin(), receivedData.end(), Compare() );
   dataProcessor_->process( receivedData );
   receivedData.clear();
}




} // namespace gather
} // namespace walberla




