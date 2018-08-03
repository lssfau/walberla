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
//! \file StructuredBlockStorage.cpp
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockDataHandling.h"
#include "StructuredBlockStorage.h"
#include "core/debug/Debug.h"


namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   This function can be used to transform any cell on level "level" into the periodic cell simulation space. For
*   example, if the simulation is periodic in x direction and the simulation domain spans from x = 0 to x = 9, then a
*   cell located at x = 38 is mapped to x = 8, and a cell located at x = -13 is mapped to x = 7.
*/
//**********************************************************************************************************************
void StructuredBlockStorage::mapToPeriodicDomain( Cell& cell, const uint_t level ) const {

   for( uint_t i = 0; i != 3; ++i ) {
      if( isPeriodic(i) ) {

         const cell_idx_t width = cell_idx_c( getNumberOfCells( i, level ) );

         if( cell[i] >= cell_idx_t(0) )
            cell[i] = cell[i] % width;
         else
            cell[i] = ( width - cell_idx_t(1) ) - ( ( -( cell[i] + cell_idx_t(1) ) ) % width );

         WALBERLA_ASSERT_GREATER_EQUAL( cell[i], domainCellBB_[level].min()[i] );
         WALBERLA_ASSERT_LESS_EQUAL(    cell[i], domainCellBB_[level].max()[i] );
      }
   }
}



//**********************************************************************************************************************
/*!
*   \brief Transforms any axis-aligned bounding box within the 3D simulation space into a cell bounding box with respect
*          to the grid refinement level "level".
*
*   Attention: If you know that your axis-aligned bounding box is perfectly aligned with the cell grid you should use
*              the member function "getCellBBFromCellAlignedAABB" in order to be save from floating point inaccuracies.
*/
//**********************************************************************************************************************
void StructuredBlockStorage::getCellBBFromAABB( CellInterval& cellBB, const AABB& aabb, const uint_t level ) const {

   const AABB& domain = getDomain();

   cellBB.min().x() = cell_idx_c( std::floor( ( aabb.xMin() - domain.xMin() ) / dx( level ) ) );
   cellBB.min().y() = cell_idx_c( std::floor( ( aabb.yMin() - domain.yMin() ) / dy( level ) ) );
   cellBB.min().z() = cell_idx_c( std::floor( ( aabb.zMin() - domain.zMin() ) / dz( level ) ) );

   cellBB.max().x() = cell_idx_c( std::ceil( ( aabb.xMax() - domain.xMin() ) / dx( level ) - real_c(1) ) );
   cellBB.max().y() = cell_idx_c( std::ceil( ( aabb.yMax() - domain.yMin() ) / dy( level ) - real_c(1) ) );
   cellBB.max().z() = cell_idx_c( std::ceil( ( aabb.zMax() - domain.zMin() ) / dz( level ) - real_c(1) ) );
}



//**********************************************************************************************************************
/*!
*   \brief Transforms an axis-aligned bounding box within the 3D simulation space that is perfectly aligned with the
*          overlaying cell grid on grid refinement level "level" into a cell bounding box.
*
*   Attention: This function should only be used if the corners of the axis-aligned bounding box "aabb" are located
*              precisely on the border between eight cells (see member function "isCellAlignedAABB" for checking whether
*              a given axis-aligned bounding box is aligned with the overlaying cell grid or not). For transforming
*              arbitrary axis-aligned bounding boxes see member function "getCellBBFromAABB".
*/
//**********************************************************************************************************************
void StructuredBlockStorage::getCellBBFromCellAlignedAABB( CellInterval& cellBB, const AABB& aabb, const uint_t level ) const {

   const real_t dx2 = real_c(0.5) * dx( level );
   const real_t dy2 = real_c(0.5) * dy( level );
   const real_t dz2 = real_c(0.5) * dz( level );

   getCell( cellBB.min(), aabb.xMin() + dx2, aabb.yMin() + dy2, aabb.zMin() + dz2, level );
   getCell( cellBB.max(), aabb.xMax() - dx2, aabb.yMax() - dy2, aabb.zMax() - dz2, level );
}



//**********************************************************************************************************************
/*!
*   Checks whether a given axis-aligned bounding box is aligned with the overlaying cell grid on level "level".
*/
//**********************************************************************************************************************
bool StructuredBlockStorage::isCellAlignedAABB( const AABB& aabb, const uint_t level ) const
{
   const AABB& domain = getDomain();

   const real_t _dx = dx( level );
   const real_t _dy = dy( level );
   const real_t _dz = dz( level );

   const real_t xMinModDx = real_c( fmod( aabb.xMin() - domain.xMin(), _dx ) );
   const real_t yMinModDy = real_c( fmod( aabb.yMin() - domain.yMin(), _dy ) );
   const real_t zMinModDz = real_c( fmod( aabb.zMin() - domain.zMin(), _dz ) );
   const real_t xMaxModDx = real_c( fmod( aabb.xMax() - domain.xMin(), _dx ) );
   const real_t yMaxModDy = real_c( fmod( aabb.yMax() - domain.yMin(), _dy ) );
   const real_t zMaxModDz = real_c( fmod( aabb.zMax() - domain.zMin(), _dz ) );

   return ( realIsEqual( xMinModDx, real_c(0) ) || realIsEqual( xMinModDx, _dx) ) &&
          ( realIsEqual( yMinModDy, real_c(0) ) || realIsEqual( yMinModDy, _dy) ) &&
          ( realIsEqual( zMinModDz, real_c(0) ) || realIsEqual( zMinModDz, _dz) ) &&
          ( realIsEqual( xMaxModDx, real_c(0) ) || realIsEqual( xMaxModDx, _dx) ) &&
          ( realIsEqual( yMaxModDy, real_c(0) ) || realIsEqual( yMaxModDy, _dy) ) &&
          ( realIsEqual( zMaxModDz, real_c(0) ) || realIsEqual( zMaxModDz, _dz) ) ;
}



//**********************************************************************************************************************
/*!
*   Returns an axis-aligned bounding box within the 3D simulation space that covers the area defined by the cell
*   bounding box "cellBB" on level "level".
*/
//**********************************************************************************************************************
void StructuredBlockStorage::getAABBFromCellBB( AABB& aabb, const CellInterval& cellBB, const uint_t level ) const {

   if( cellBB.empty() )
      aabb.init(0,0,0,0,0,0);
   else {

      const AABB& domain = getDomain();

      aabb.init( domain.xMin() + real_c( cellBB.xMin()     ) * dx( level ),
                 domain.yMin() + real_c( cellBB.yMin()     ) * dy( level ),
                 domain.zMin() + real_c( cellBB.zMin()     ) * dz( level ),
                 domain.xMin() + real_c( cellBB.xMax() + 1 ) * dx( level ),
                 domain.yMin() + real_c( cellBB.yMax() + 1 ) * dy( level ),
                 domain.zMin() + real_c( cellBB.zMax() + 1 ) * dz( level ) );
   }
}



StructuredBlockStorage::StructuredBlockStorage( const shared_ptr<BlockStorage>& blockStorage, const std::vector< uint_t >& xCells,
                                                const std::vector< uint_t >& yCells, const std::vector< uint_t >& zCells ) :

   blockStorage_( blockStorage ), blockCellBBCreated_( false ), blockCellBBId_( 0 )
{
   resetCellDecomposition( xCells, yCells, zCells );
}



void StructuredBlockStorage::resetCellDecomposition( const std::vector< uint_t > & xCells, const std::vector< uint_t > & yCells,
                                                     const std::vector< uint_t > & zCells )
{
   WALBERLA_ASSERT_GREATER( xCells.size(), 0 );
   WALBERLA_ASSERT_EQUAL( xCells.size(), yCells.size() );
   WALBERLA_ASSERT_EQUAL( yCells.size(), zCells.size() );

   levels_ = xCells.size();

   const AABB& domain = getDomain();

   const real_t xWidth = domain.xMax() - domain.xMin();
   const real_t yWidth = domain.yMax() - domain.yMin();
   const real_t zWidth = domain.zMax() - domain.zMin();

   domainCellBB_.clear();
   dx_.clear();
   dy_.clear();
   dz_.clear();

   for( uint_t i = 0; i != levels_; ++i )
   {
      WALBERLA_ASSERT_GREATER( xCells[i], 0 );
      WALBERLA_ASSERT_GREATER( yCells[i], 0 );
      WALBERLA_ASSERT_GREATER( zCells[i], 0 );

      domainCellBB_.emplace_back( 0, 0, 0, cell_idx_c( xCells[i]-1 ), cell_idx_c( yCells[i]-1 ), cell_idx_c( zCells[i]-1 ) );

      dx_.push_back( xWidth / real_c( xCells[i] ) );
      dy_.push_back( yWidth / real_c( yCells[i] ) );
      dz_.push_back( zWidth / real_c( zCells[i] ) );
   }
}



//**********************************************************************************************************************
/*!
*   Adds a cell bounding box to every block as a block data "item". For accessing this cell bounding box see member
*   functions "getBlockCellBBId()" and "getBlockCellBB( [const] IBlock& block )".
*/
//**********************************************************************************************************************
void StructuredBlockStorage::createCellBoundingBoxes() {

   if( !blockCellBBCreated_ )
   {
      blockCellBBId_ = addCellBoundingBoxesAsBlockData( "cell bounding box" );
      blockCellBBCreated_ = true;
   }
}



CellInterval * StructuredBlockStorage::initializeCellBoundingBox( IBlock * const block ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_EQUAL( &(getBlockStorage()), &(block->getBlockStorage()) );

   const uint_t level = getLevel( *block );
   const AABB & aabb  = block->getAABB();

   WALBERLA_ASSERT( isCellAlignedAABB( aabb, level ), "AABB = " << aabb << ", level = " << level );

   CellInterval* cellBB = new CellInterval();

   getCellBBFromCellAlignedAABB( *cellBB, aabb, level );

   return cellBB;
}



} // namespace domain_decomposition
} // namespace walberla
