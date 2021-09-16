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
//! \file CellGatherPackInfo.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "DataProcessor.h"
#include "GatherPackInfo.h"
#include "core/cell/CellVector.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <type_traits>


namespace walberla {
namespace gather {




//**********************************************************************************************************************
/*! Packs data in cells given by a CellContainter ( CellInterval, CellVector..)
*
* Template Parameters:
*  - Field_T        Field or Field-Adaptor that should be packed
*  - CellContainer  for example CellInterval, CellVector
*
* \ingroup gather
*/
//**********************************************************************************************************************
template< typename Field_T, typename CellContainer>
class CellGatherPackInfo : public GatherPackInfo
{
public:

   //*******************************************************************************************************************
   /*! Constructor
   *
   * This PackInfo packs data of a field in all cells that are provided by the CellContainer.
   * The output to the DataProcessor contains in first column the cellNr in the container, then the field value(s)
   * in the following columns.
   *
   * Example: Collect Data along a line through the domain
   *     - CellContainer would be a CellInterval, that has size 1 in two dimensions, and the length of the line
   *       is given by the third size
   *
   * \param bs          StructuredBlockStorage containing the field
   * \param fieldID     BlockDataID of the field that should be collected
   * \param containerOfGlobalCells   Container of cells in global coordinates, order of cells is important for output
   * \param dp          data processor where gathered data is delivered to
   */
   //*******************************************************************************************************************
   CellGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                       ConstBlockDataID fieldID,
                       const CellContainer & containerOfGlobalCells,
                       const shared_ptr<DataProcessor> & dp );



   //** Packing Interface  *******************************************************************************
   /*! \name Packing Interface  */
   //@{

   void packData  ( const IBlock * sender,
                            mpi::SendBuffer & outBuffer ) override;

   void unpackData( mpi::RecvBuffer & buffer ) override;

   void gatherFinished() override;
   //@}
   //*******************************************************************************************************************


protected:
   static_assert( std::is_fundamental<typename Field_T::value_type >::value,
                  "CellGatherPackInfo supports fields of build in datatypes"  );


   shared_ptr<StructuredBlockStorage> blocks_;

   /// DataInterpolator acting as source, for the data that has to be packed
   ConstBlockDataID fieldID_;

   struct Samples {
         std::vector<uint_t>    globalCellNr;
         CellVector             cells;
   };

   /// For every LOCAL block, where data has to be packed, a vector of
   /// local cell coordinates is stored
   std::map<const IBlock*, Samples > localSamplePoints_;


   //** Members for Receiving  ***************************************************************************
   /*! \name Members for Receiving  */
   //@{

   /// Two dimensional array of received data:
   /// the outer vector has one entry for every sample point
   /// the inner vector represents one sample point
   ///  - the first entry (receivedData[i][0] ) is the t-value
   ///  - subsequent entries are for the f values
   std::vector<std::vector<real_t> > receivedData;

   /// Channel for output of gathered data
   shared_ptr<DataProcessor> dataProcessor_;

   /// Helper class for sorting the receivedData array according to t-value
   struct Compare : public std::function<bool(std::vector<real_t> , std::vector<real_t>)>
   {
      inline bool operator()(const std::vector<real_t> & v1, const std::vector<real_t> & v2) const {
         return v1[0] < v2[0];
      }
   };
   //@}
   //*******************************************************************************************************************
};




} // namespace gather
} // namespace walberla

#include "CellGatherPackInfo.impl.h"

