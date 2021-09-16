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
//! \file HashGrids.cpp
//! \author Klaus Iglberger
//! \author Florian Schornbaum
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "HashGrids.h"
#include <pe/rigidbody/BodyStorage.h>

#include "core/timing/TimingTree.h"

namespace walberla{
namespace pe{
namespace ccd {

//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for the HashGrid class.
 */
HashGrids::HashGrid::HashGrid( real_t cellSpan )
{
   // Initialization of all member variables and ...
   xCellCount_   = powerOfTwo( xCellCount ) ? xCellCount : 16;
   yCellCount_   = powerOfTwo( yCellCount ) ? yCellCount : 16;
   zCellCount_   = powerOfTwo( zCellCount ) ? zCellCount : 16;

   xHashMask_    = xCellCount_ - 1;
   yHashMask_    = yCellCount_ - 1;
   zHashMask_    = zCellCount_ - 1;

   xyCellCount_  = xCellCount_ * yCellCount_;
   xyzCellCount_ = xyCellCount_ * zCellCount_;

   enlargementThreshold_ = xyzCellCount_ / minimalGridDensity;

   // ... allocation of the linear array that is representing the hash grid.
   cell_ = new Cell[ xyzCellCount_ ];

   // Initialization of each cell - i.e., initially setting the pointer to the body container to
   // NULL (=> no bodies are assigned to this hash grid yet!) and ...
   for( Cell* c  = cell_; c < cell_ + xyzCellCount_; ++c ) {
      c->bodies_ = nullptr;
   }
   // ... setting up the neighborhood relationship (using the offset array).
   initializeNeighborOffsets();

   cellSpan_        = cellSpan;
   inverseCellSpan_ = real_c( 1 ) / cellSpan;

   occupiedCells_.reserve( occupiedCellsVectorSize );

   bodyCount_ = 0;
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the HashGrid class.
 */
HashGrids::HashGrid::~HashGrid()
{
   clear();

   for( Cell* c  = cell_; c < cell_ + xyzCellCount_; ++c ) {
      if( c->neighborOffset_ != stdNeighborOffset_ ) delete[] c->neighborOffset_;
   }
   delete[] cell_;
}
//*************************************************************************************************




//=================================================================================================
//
//  ADD/REMOVE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adds a body to this hash grid.
 *
 * \param body The body that is about to be added to this hash grid.
 * \return void
 *
 * This function is called every time a new rigid body is added to this grid. If adding the body
 * will cause the total number of bodies assigned to this grid to exceed the enlargement threshold,
 * the size of this hash grid is increased in order to maintain a fixed minimal grid density (=
 * ratio of cells to bodies). This function may also be called during the update phase.
 */
void HashGrids::HashGrid::add( BodyID body )
{
   // If adding the body will cause the total number of bodies assigned to this grid to exceed the
   // enlargement threshold, the size of this hash grid must be increased.
   if( bodyCount_ == enlargementThreshold_ ) enlarge();

   // Calculate (and store) the hash value (= the body's cell association) and ...
   size_t h = hash( body );
   body->setHash( h );

   // ... insert the body into the corresponding cell.
   Cell* cell = cell_ + h;
   add( body, cell );

   ++bodyCount_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removes a body from this hash grid.
 *
 * \param body The body that is about to be removed from this hash grid.
 * \return void
 *
 * This function is called every time a rigid body is removed from this grid. It may also be called
 * during the update phase.
 */
void HashGrids::HashGrid::remove( BodyID body )
{
   // The stored hash value (= the body's cell association) is used in order to directly access the
   // cell from which this body will be removed.
   Cell* cell = cell_ + body->getHash();
   remove( body, cell );

   --bodyCount_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Updates the cell association of a body that is assigned to this grid.
 *
 * \param body The body whose cell association is being updated.
 * \return void
 *
 * Checks if a given rigid body can remain stored in its currently assigned cell - and if not,
 * removes the body from this cell and reassigns it to another cell (which is identified by
 * re-evaluating the body's hash value).
 */
void HashGrids::HashGrid::update( BodyID body )
{
   // The hash value is recomputed based on the body's current spatial location.
   size_t newHash = hash( body );
   size_t oldHash = body->getHash();

   // If this new hash value is identical to the hash value of the previous time step, the body
   // remains assigned to its current grid cell.
   if( newHash == oldHash ) return;

   // Only if the hash value changes, the cell association has to be changed, too - meaning, the
   // body has to be removed from its currently assigned cell and ...
   Cell* cell = cell_ + oldHash;
   remove( body, cell );

   body->setHash( newHash );

   // ... stored in the cell that corresponds to the new hash value.
   cell = cell_ + newHash;
   add( body, cell );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Clears the hash grid.
 *
 * \return void
 *
 * This function removes all (the handles to) the bodies that are assigned to this hash grid from
 * the body containers of the grid's cells.
 */
void HashGrids::HashGrid::clear()
{
   for( auto cellIt = occupiedCells_.begin(); cellIt < occupiedCells_.end(); ++cellIt ) {
      delete (*cellIt)->bodies_;
      (*cellIt)->bodies_ = nullptr;
   }
   occupiedCells_.clear();
   bodyCount_ = 0;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Sets up the neighborhood relationships for all grid cells.
 *
 * \return void
 *
 * This function is used to initialize the offset arrays of all grid cells. The offsets are required
 * for ensuring fast direct access to all directly adjacent cells for each cell in the hash grid.
 */
void HashGrids::HashGrid::initializeNeighborOffsets()
{
   offset_t xc   = static_cast<offset_t>( xCellCount_ );
   offset_t yc   = static_cast<offset_t>( yCellCount_ );
   offset_t zc   = static_cast<offset_t>( zCellCount_ );
   offset_t xyc  = static_cast<offset_t>( xyCellCount_ );
   offset_t xyzc = static_cast<offset_t>( xyzCellCount_ );

   // Initialization of the grid-global offset array that is valid for all inner cells in the hash grid.
   unsigned int i = 0;
   for( offset_t zz = -xyc; zz <= xyc; zz += xyc ) {
      for( offset_t yy = -xc; yy <= xc; yy += xc ) {
         for( offset_t xx = -1; xx <= 1; ++xx, ++i ) {
            stdNeighborOffset_[i] = xx + yy + zz;
         }
      }
   }

   // Allocation and initialization of the offset arrays of all the border cells. All inner cells
   // are set to point to the grid-global offset array.
   Cell* c = cell_;
   for( offset_t z = 0; z < zc; ++z ) {
      for( offset_t y = 0; y < yc; ++y ) {
         for( offset_t x = 0; x < xc; ++x, ++c ) {

            /* border cell */
            if( x == 0 || x == (xc - 1) || y == 0 || y == (yc - 1) || z == 0 || z == (zc - 1) ) {

               c->neighborOffset_ = new offset_t[27];

               i = 0;
               for( offset_t zz = -xyc; zz <= xyc; zz += xyc )
               {
                  offset_t zo = zz;
                  if( z == 0 && zz == -xyc ) {
                     zo = xyzc - xyc;
                  }
                  else if( z == (zc - 1) && zz == xyc ) {
                     zo = xyc - xyzc;
                  }

                  for( offset_t yy = -xc; yy <= xc; yy += xc )
                  {
                     offset_t yo = yy;
                     if( y == 0 && yy == -xc ) {
                        yo = xyc - xc;
                     }
                     else if( y == (yc - 1) && yy == xc ) {
                        yo = xc - xyc;
                     }

                     for( offset_t xx = -1; xx <= 1; ++xx, ++i ) {

                        offset_t xo = xx;
                        if( x == 0 && xx == -1 ) {
                           xo = xc - 1;
                        }
                        else if( x == (xc - 1) && xx == 1 ) {
                           xo = 1 - xc;
                        }

                        c->neighborOffset_[i] = xo + yo + zo;
                     }
                  }
               }
            }
            /* inner cell */
            else {
               c->neighborOffset_ = stdNeighborOffset_;
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hash value (= cell association) of a given rigid body.
 *
 * \param body The body whose hash value is about to be calculated.
 * \return The hash value (=cell association) of the body.
 */
size_t HashGrids::HashGrid::hash( BodyID body ) const
{
   const AABB& bodyAABB = body->getAABB();
   return hashPoint(bodyAABB.xMin(), bodyAABB.yMin(), bodyAABB.zMin());
}
//*************************************************************************************************

   
/*!\brief Computes the hash for a given point.
 *
 * \param x X value of the point.
 * \param y Y value of the point.
 * \param z Z value of the point.
 * \return The hash value (=cell association) of the point.
 *
 * The hash calculation uses modulo operations in order to spatially map entire blocks of connected
 * cells to the origin of the coordinate system. This block of cells at the origin of the coordinate
 * system that is filled with all the bodies of the simulation is referred to as the hash grid. The
 * key feature, and ultimately the advantage, of hash grids is the fact that two adjacent cells that
 * are located anywhere in the simulation space are mapped to two cells that are still adjacent in
 * the hashed storage structure.
 *
 * Note that the modulo calculations are replaced with fast bitwise AND operations - hence, the
 * spatial dimensions of the hash grid must be restricted to powers of two!
 */
size_t HashGrids::HashGrid::hashPoint(real_t x, real_t y, real_t z) const {
   size_t xHash;
   size_t yHash;
   size_t zHash;
   
   if( x < 0 ) {
      real_t i = ( -x ) * inverseCellSpan_;
      xHash  = xCellCount_ - 1 - ( static_cast<size_t>( i ) & xHashMask_ );
   }
   else {
      real_t i = x * inverseCellSpan_;
      xHash  = static_cast<size_t>( i ) & xHashMask_;
   }
   
   if( y < 0 ) {
      real_t i = ( -y ) * inverseCellSpan_;
      yHash  = yCellCount_ - 1 - ( static_cast<size_t>( i ) & yHashMask_ );
   }
   else {
      real_t i = y * inverseCellSpan_;
      yHash  = static_cast<size_t>( i ) & yHashMask_;
   }
   
   if( z < 0 ) {
      real_t i = ( -z ) * inverseCellSpan_;
      zHash  = zCellCount_ - 1 - ( static_cast<size_t>( i ) & zHashMask_ );
   }
   else {
      real_t i = z * inverseCellSpan_;
      zHash  = static_cast<size_t>( i ) & zHashMask_;
   }
   
   return xHash + yHash * xCellCount_ + zHash * xyCellCount_;
}

//*************************************************************************************************
/*!\brief Adds a body to a specific cell in this hash grid.
 *
 * \param body The body that is about to be added to this hash grid.
 * \param cell The cell the body is assigned to.
 * \return void
 */
void HashGrids::HashGrid::add( BodyID body, Cell* cell )
{
   // If this cell is already occupied by other bodies, which means the pointer to the body
   // container holds a valid address and thus the container itself is properly initialized, then
   // the body is simply added to this already existing body container. Note that the index position
   // is memorized (=> "body->setCellId()") in order to ensure constant time removal.

   // If, however, the cell is still empty, then the object container, first of all, must be created
   // (i.e., allocated) and properly initialized (i.e., sufficient initial storage capacity must be
   // reserved). Furthermore, the cell must be inserted into the grid-global vector 'occupiedCells_'
   // in which all cells that are currently occupied by bodies are recorded.
   if( cell->bodies_ == nullptr )
   {
      cell->bodies_ = new BodyVector;
      cell->bodies_->reserve( cellVectorSize );

      cell->lastNonFixedBody_ = -1;

      cell->occupiedCellsId_ = occupiedCells_.size();
      occupiedCells_.push_back( cell );
   }

   // If body is not fixed add it at the end. If body is not fixed add it at the end and swap it with
   // lastNonFixedBody_

   body->setCellId( cell->bodies_->size() );
   cell->bodies_->push_back( body );

   if ( !body->isFixed() )
   {
      ++(cell->lastNonFixedBody_);
      BodyID bd2 = cell->bodies_->at(static_cast<size_t>(cell->lastNonFixedBody_));
      auto tmp = body->getCellId();
      body->setCellId(bd2->getCellId());
      bd2->setCellId(tmp);
      std::swap(cell->bodies_->at(static_cast<size_t>(cell->lastNonFixedBody_)), cell->bodies_->at( cell->bodies_->size()-1 ));
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removes a body from a specific cell in this hash grid.
 *
 * \param body The body that is about to be removed from this hash grid.
 * \param cell The cell the body is removed from.
 * \return void
 */
void HashGrids::HashGrid::remove( BodyID body, Cell* cell )
{
   // If the body is the last body that is stored in this cell ...
   if( cell->bodies_->size() == 1 ) {

      // ... the cell's body container is destroyed and ...
      delete cell->bodies_;
      cell->bodies_ = nullptr;
      cell->lastNonFixedBody_ = -1;

      // ... the cell is removed from the grid-global vector 'occupiedCells_' that records all
      // body-occupied cells. Since the cell memorized its index (=> 'occupiedCellsId_') in this
      // vector, it can be removed in constant time, O(1).
      if( cell->occupiedCellsId_ == occupiedCells_.size() - 1 ) {
         occupiedCells_.pop_back();
      }
      else {
         Cell* lastCell = occupiedCells_.back();
         occupiedCells_.pop_back();
         lastCell->occupiedCellsId_ = cell->occupiedCellsId_;
         occupiedCells_[ cell->occupiedCellsId_ ] = lastCell;
      }
   }
   // If the body is *not* the last body that is stored in this cell ...
   else
   {
      size_t cellId = body->getCellId();

      if (!body->isFixed())
      {
         //copy last nonfixed body to body that should be removed
         BodyID bd2 = cell->bodies_->at(static_cast<size_t>(cell->lastNonFixedBody_));
         bd2->setCellId( body->getCellId() );
         cell->bodies_->at( cellId ) = cell->bodies_->at(static_cast<size_t>(cell->lastNonFixedBody_));

         //mark lastNonFixedBody_ to be removed
         cellId = static_cast<size_t>(cell->lastNonFixedBody_);

         --(cell->lastNonFixedBody_);
      }

      // ... the body is removed from the cell's body container. Since the body memorized its
      // index (=> 'cellId') in this container, it can be removed in constant time, O(1).
      if( cellId == cell->bodies_->size() - 1 )
      {
         cell->bodies_->pop_back();
      }
      else
      {
         BodyID lastElement = cell->bodies_->back();
         cell->bodies_->pop_back();
         lastElement->setCellId( cellId );
         (*cell->bodies_)[ cellId ] = lastElement;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Increases the number of cells used by this hash grid.
 *
 * \return void
 *
 * In order to handle an initially unknown and ultimately arbitrary number of bodies, the hash grid,
 * starting with a rather small number of cells at the time of its creation, must have the ability
 * to grow as new bodies are inserted. Therefore, if by inserting a body into this hash grid the
 * associated grid density - that is the ratio of cells to bodies - drops below the threshold
 * specified by \a minimalGridDensity, or in other words, if the number of bodies grows larger than
 * specified by \a enlargementThreshold_, the number of cells in each coordinate direction is
 * doubled (thus the total number of grid cells is increased by a factor of 8).
 */
void HashGrids::HashGrid::enlarge()
{
   BodyID* bodies = new BodyID[ bodyCount_ ];
   BodyID* body   = bodies;

   // All bodies that are assigned to this grid are temporarily removed, ...
   for( auto cellIt = occupiedCells_.begin(); cellIt < occupiedCells_.end(); ++cellIt ) {
      BodyVector* cellBodies = (*cellIt)->bodies_;
      for( auto eIt = cellBodies->begin(); eIt < cellBodies->end(); ++eIt ) {
         *(body++) = *eIt;
      }
   }

   // ... the grid's current data structures are deleted, ...
   clear();

   for( Cell* c  = cell_; c < cell_ + xyzCellCount_; ++c ) {
      if( c->neighborOffset_ != stdNeighborOffset_ ) delete[] c->neighborOffset_;
   }
   delete[] cell_;

   // ... the number of cells is doubled in each coordinate direction, ...
   xCellCount_  *= 2;
   yCellCount_  *= 2;
   zCellCount_  *= 2;

   xHashMask_    = xCellCount_ - 1;
   yHashMask_    = yCellCount_ - 1;
   zHashMask_    = zCellCount_ - 1;

   xyCellCount_  = xCellCount_ * yCellCount_;
   xyzCellCount_ = xyCellCount_ * zCellCount_;

   // ... a new threshold for enlarging this hash grid is set, ...
   enlargementThreshold_ = xyzCellCount_ / minimalGridDensity;

   // ... a new linear array of cells representing this enlarged hash grid is allocated and ...
   cell_ = new Cell[ xyzCellCount_ ];

   // ... initialized, and finally ...
   for( Cell* c  = cell_; c < cell_ + xyzCellCount_; ++c ) {
      c->bodies_ = nullptr;
      c->lastNonFixedBody_ = -1;
   }
   initializeNeighborOffsets();

   // ... all previously removed bodies are reinserted.
   for( BodyID* p = bodies; p < body; ++p ) {
      add(*p);
   }
   delete[] bodies;
}
//*************************************************************************************************








//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
//\brief Constructor for the HashGrids class.
//*
//* \param bodystorage Reference to the general body storage.
//*
//* Note that all (local, global and remote) must be contained in the body storage.
//
//HashGrids::HashGrids( BodyStorage& bodystorage )
//   : bodystorage_( bodystorage )
//   , bodystorageShadowCopies_( bodystorage )
//{
//   nonGridBodies_.reserve( gridActivationThreshold );

//   // As long as the number of bodies is less-or-equal than specified by 'gridActivationThreshold',
//   // no hierarchy of hash grids is constructed. Instead, all bodies are simply stored in
//   // 'nonGridBodies_' and pairwise collision checks are performed.
//   gridActive_ = false;

//   bodystorage_.registerAddCallback( "HashGrids", std::bind(&HashGrids::add, this, std::placeholders::_1) );
//   bodystorage_.registerRemoveCallback( "HashGrids", std::bind(&HashGrids::remove, this, std::placeholders::_1) );
//}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the HashGrids class.
 *
 * \param bodystorage Reference to the central body storage.
 * \param bodystorageShadowCopies Reference to the body storage containing all body shadow copies.
 */
HashGrids::HashGrids( BodyStorage& globalStorage, BodyStorage& bodystorage, BodyStorage& bodystorageShadowCopies )
   : globalStorage_(globalStorage)
   , bodystorage_( bodystorage )
   , bodystorageShadowCopies_( bodystorageShadowCopies )
   , observedBodyCount_(0)
{
   nonGridBodies_.reserve( gridActivationThreshold );

   // As long as the number of bodies is less-or-equal than specified by 'gridActivationThreshold',
   // no hierarchy of hash grids is constructed. Instead, all bodies are simply stored in
   // 'nonGridBodies_' and pairwise collision checks are performed.
   gridActive_ = false;

   bodystorage_.registerAddCallback( "HashGrids", std::bind(&HashGrids::add, this, std::placeholders::_1) );
   bodystorage_.registerRemoveCallback( "HashGrids", std::bind(&HashGrids::remove, this, std::placeholders::_1) );

   bodystorageShadowCopies_.registerAddCallback( "HashGrids", std::bind(&HashGrids::add, this, std::placeholders::_1) );
   bodystorageShadowCopies_.registerRemoveCallback( "HashGrids", std::bind(&HashGrids::remove, this, std::placeholders::_1) );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destructor for the HashGrids class.
 */
HashGrids::~HashGrids()
{
   // Delete all grids that are stored in the grid hierarchy (=> gridList_).
   for( auto gridIt = gridList_.begin(); gridIt != gridList_.end(); ++gridIt ) {
      delete (*gridIt);
   }

   bodystorage_.deregisterAddCallback( "HashGrids" );
   bodystorage_.deregisterRemoveCallback( "HashGrids" );

   if( &bodystorage_ != &bodystorageShadowCopies_ )
   {
      bodystorageShadowCopies_.deregisterAddCallback( "HashGrids" );
      bodystorageShadowCopies_.deregisterRemoveCallback( "HashGrids" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  ADD/REMOVE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Registering a rigid body with the coarse collision detector.
 *
 * \param body The body to be registered with the coarse collision detector.
 * \return void
 *
 * This function is a requirement for all coarse collision detectors. It is automatically
 * called every time a new rigid body is added to the simulation world.
 */
void HashGrids::add( BodyID body )
{
   ++observedBodyCount_;
   // The body is marked as being added to 'bodiesToAdd_' by setting the grid pointer to NULL and
   // setting the cell-ID to '0'. Additionally, the hash value is used to memorize the body's
   // index position in the 'bodiesToAdd_' vector.
   body->setGrid  ( nullptr );
   body->setHash  ( bodiesToAdd_.size() );
   body->setCellId( 0 );

   // Temporarily add the body to 'bodiesToAdd_'. As soon as "findContacts()" is called, all
   // bodies stored in 'bodiesToAdd_' are finally inserted into the data structure.
   bodiesToAdd_.push_back( body );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deregistering a rigid body from the coarse collision detector.
 *
 * \param body The body to be deregistered from the coarse collision detector.
 * \return void
 *
 * This function is a requirement for all coarse collision detectors. It is automatically
 * called every time a rigid body is removed from the simulation world.
 */
void HashGrids::remove( BodyID body )
{
   --observedBodyCount_;

   HashGrid* grid = static_cast<HashGrid*>( body->getGrid() );

   // The body is stored in a hash grid from which it must be removed.
   if( grid != nullptr ) {
      grid->remove( body );
   }
   // The body's grid pointer is equal to NULL.
   // => The body is either stored in 'bodiesToAdd_' (-> cell-ID = 0) or 'nonGridBodies_' (-> cell-ID = 1).
   else {
      if( body->getCellId() == 0 ) {
         // the body's hash value => index of this body in 'bodiesToAdd_'
         if( body->getHash() == bodiesToAdd_.size() - 1 ) {
            bodiesToAdd_.pop_back();
         }
         else {
            BodyID lastElement = bodiesToAdd_.back();
            bodiesToAdd_.pop_back();
            lastElement->setHash( body->getHash() );
            bodiesToAdd_[ body->getHash() ] = lastElement;
         }
      }
      else {
         // the body's hash value => index of this body in 'nonGridBodies_'
         if( body->getHash() == nonGridBodies_.size() - 1 ) {
            nonGridBodies_.pop_back();
         }
         else {
            BodyID lastElement = nonGridBodies_.back();
            nonGridBodies_.pop_back();
            lastElement->setHash( body->getHash() );
            nonGridBodies_[ body->getHash() ] = lastElement;
         }
      }
   }
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Clears the coarse collision detector.
 *
 * \return void
 */
void HashGrids::clear()
{
   for( auto gridIt = gridList_.begin(); gridIt != gridList_.end(); ++gridIt ) {
      delete (*gridIt);
   }
   gridList_.clear();

   gridActive_ = false;

   nonGridBodies_.clear();

   bodiesToAdd_.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief clears all bodies from the hash grid and reloads bodies from bodystorage and shadowbodystorage
 *
 * \return void
 */
void HashGrids::reloadBodies()
{
   clear();

   for (auto& body : bodystorage_)
   {
      add( &body );
   }
   if( &bodystorage_ != &bodystorageShadowCopies_ )
   {
      for (auto& body : bodystorageShadowCopies_)
      {
         add( &body );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Updates all hash grids and reassigns bodies.
 */
void HashGrids::update(WcTimingTree* tt)
{
   // ----- UPDATE PHASE ----- //

   if (tt != nullptr) tt->start("AddNewBodies");
   // Finally add all bodies that were temporarily stored in 'bodiesToAdd_' to the data structure.
   if( !bodiesToAdd_.empty() )
   {
      for( auto bodyIt = bodiesToAdd_.begin(); bodyIt < bodiesToAdd_.end(); ++bodyIt ) {

         if( gridActive_ ) addGrid( *bodyIt );
         else              addList( *bodyIt );
      }
      bodiesToAdd_.clear();
   }
   if (tt != nullptr) tt->stop("AddNewBodies");

   if (tt != nullptr) tt->start("Update");
   // Update the data structure (=> adapt to the current body distribution) by taking care of
   // moved, rotated and/or deformed bodies.
   if( gridActive_ )
   {
      // 1) Update strategy: All bodies are simultaneously removed from the data structure, and
      //                     immediately afterwards they are reinserted based on their current sizes
      //                     and spatial locations - which may or may not have changed since the
      //                     last time step.
      //
      // for( typename GridList::iterator grid = gridList_.begin(); grid != gridList_.end(); ++grid ) {
      //    (*grid)->clear();
      // }
      //
      // {
      //    const Iterator end( bodystorage_.end() );
      //    for( Iterator b = bodystorage_.begin(); b != end; ++b ) {
      //
      //       BodyID    body = *b;
      //       HashGrid* grid = static_cast<HashGrid*>( body->getGrid() );
      //
      //       if( grid != NULL ) {
      //
      //          real_t size     = body->getAABBSize();
      //          real_t cellSpan = grid->getCellSpan();
      //
      //          if( size >= cellSpan || size < ( cellSpan / hierarchyFactor ) ) {
      //             addGrid( body );
      //          }
      //          else {
      //             grid->add( body );
      //          }
      //       }
      //    }
      // }
      //
      // if( &bodystorage_ != &bodystorageShadowCopies_ ){
      //    const Iterator end( bodystorageShadowCopies_.end() );
      //    for( Iterator b = bodystorageShadowCopies_.begin(); b != end; ++b ) {
      //
      //       BodyID    body = *b;
      //       HashGrid* grid = static_cast<HashGrid*>( body->getGrid() );
      //
      //       if( grid != NULL ) {
      //
      //          real_t size     = body->getAABBSize();
      //          real_t cellSpan = grid->getCellSpan();
      //
      //          if( size >= cellSpan || size < ( cellSpan / hierarchyFactor ) ) {
      //             addGrid( body );
      //          }
      //          else {
      //             grid->add( body );
      //          }
      //       }
      //    }
      // }

      // 2) Update strategy: For every body, the grid association is recomputed based on the body's
      //                     current size (=> size of its AABB -> length of the longest edge of the
      //                     AABB, which most likely changes every time the body is rotated or
      //                     deformed):
      //
      //                      - If the new grid association is identical to the current grid
      //                        association, the body may remain assigned to its current grid, yet
      //                        still the body's cell association within this grid has to be checked
      //                        and potentially changed if the body was moved (=> "grid->update()").
      //
      //                      - If the grid association has changed, the body must be removed from
      //                        its current grid (=> "grid->remove()") and reassigned to a grid with
      //                        suitably sized cells (=> "addGrid()").

      {
         for( auto& body : bodystorage_ )
         {
            HashGrid* grid = static_cast<HashGrid*>( body.getGrid() );

            if( grid != nullptr )
            {
               real_t size     = body.getAABBSize();
               real_t cellSpan = grid->getCellSpan();

               if( size >= cellSpan || size < ( cellSpan / hierarchyFactor ) ) {
                  grid->remove( &body );
                  addGrid( &body );
               }
               else {
                  grid->update( &body );
               }
            }
         }
      }

      if( &bodystorage_ != &bodystorageShadowCopies_ ) {
         for( auto& body : bodystorageShadowCopies_ )
         {
            HashGrid* grid = static_cast<HashGrid*>( body.getGrid() );

            if( grid != nullptr )
            {
               real_t size     = body.getAABBSize();
               real_t cellSpan = grid->getCellSpan();

               if( size >= cellSpan || size < ( cellSpan / hierarchyFactor ) ) {
                  grid->remove( &body );
                  addGrid( &body );
               }
               else {
                  grid->update( &body );
               }
            }
         }
      }
   }
   if (tt != nullptr) tt->stop("Update");
}

//**Implementation of ICCD interface ********************************************************
//*************************************************************************************************
/*!\brief Contact generation between colliding rigid bodies.
 *
 * \return Vector of possible contacts.
 *
 * This function generates all contacts between all registered rigid bodies. The contacts are
 * added to the contact container which can be retrieved via getPossibleContacts().
 */
PossibleContacts& HashGrids::generatePossibleContacts( WcTimingTree* tt )
{
   if (tt != nullptr) tt->start("CCD");

   contacts_.clear();

   WALBERLA_LOG_DETAIL( "   Finding the contacts via the hierarchical hash grids algorithm...");

   update(tt);

   if (tt != nullptr) tt->start("Detection");
   // ----- DETECTION STEP ----- //

   // Contact generation by traversing through all hash grids (which are sorted in ascending order
   // with respect to the size of their cells).
   for( auto gridIt = gridList_.begin(); gridIt != gridList_.end(); ++gridIt ) {

      // Contact generation for all bodies stored in the currently processed grid 'grid'.
      BodyID* bodies     = nullptr;
      size_t  bodyCount = (*gridIt)->process( &bodies, contacts_ );

      if( bodyCount > 0 ) {

         // Test all bodies stored in 'grid' against bodies stored in grids with larger sized cells.
         auto nextGridIt = gridIt;
         for( ++nextGridIt; nextGridIt != gridList_.end(); ++nextGridIt ) {
            (*nextGridIt)->processBodies( bodies, bodyCount, contacts_ );
         }

         BodyID* bodiesEnd = bodies + bodyCount;
         for( BodyID* a = bodies; a < bodiesEnd; ++a ) {
            // Test all bodies stored in 'grid' against all bodies stored in 'nonGridBodies_'.
            for( auto bIt = nonGridBodies_.begin(); bIt < nonGridBodies_.end(); ++bIt ) {
               collide( *a, *bIt, contacts_ );
            }
            // Test all bodies stored in 'grid' against all bodies stored in 'globalStorage_'.
            for( auto bIt = globalStorage_.begin(); bIt < globalStorage_.end(); ++bIt ) {
               collide( *a, &(*bIt), contacts_ );
            }
         }
      }

      delete[] bodies;
   }

   for( auto aIt = nonGridBodies_.begin(); aIt < nonGridBodies_.end(); ++aIt ) {
      // Pairwise test (=> contact generation) for all bodies that are stored in 'nonGridBodies_'.
      for( auto bIt = aIt + 1; bIt < nonGridBodies_.end(); ++bIt ) {
         collide( *aIt, *bIt, contacts_ );
      }

      // Pairwise test (=> contact generation) for all bodies that are stored in 'nonGridBodies_' with global bodies.
      for( auto bIt = globalStorage_.begin(); bIt < globalStorage_.end(); ++bIt ) {
         collide( *aIt, &(*bIt), contacts_ );
      }
   }
   if (tt != nullptr) tt->stop("Detection");

   WALBERLA_LOG_DETAIL_SECTION()
   {
      std::stringstream log;
      if( contacts_.empty() )
         log << "      no contacts found!\n";
      else {
         log << "      State of the contacts:\n";
         for( auto cIt=contacts_.begin(); cIt!=contacts_.end(); ++cIt ) {
            log << "possible Contact " << cIt->first->getSystemID() << " : " << cIt->second->getSystemID() << "\n";
         }
      }
      WALBERLA_LOG_DETAIL( log.str() );
   }

   if (tt != nullptr) tt->stop("CCD");

   return contacts_;
}


//=================================================================================================
//
//  ADD FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adds a body to the hierarchical hash grids data structure.
 *
 * \param body The body that is about to be added to the hierarchical hash grids data structure.
 * \return void
 *
 * If the usage of hierarchical hash grids is activated, this function is called every time a new
 * rigid body is added to the data structure. It may also be called during the update phase.
 */
void HashGrids::addGrid( BodyID body )
{
   real_t size = ( body->isFinite() ) ? body->getAABBSize() : real_c( -1 );

   // If the body is finite in size, it must be assigned to a grid with suitably sized cells.
   if( size > 0 )
   {
      HashGrid* grid = nullptr;

      if( gridList_.empty() )
      {
         // If no hash grid yet exists in the hierarchy, an initial hash grid is created
         // based on the body's size.

         grid = new HashGrid( size * std::sqrt( hierarchyFactor ) );
      }
      else
      {
         // Check the hierarchy for a hash grid with suitably sized cells - if such a grid does not
         // yet exist, it will be created.

         real_t cellSpan = 0;
         for( auto gIt = gridList_.begin(); gIt != gridList_.end(); ++gIt )
         {
            grid     = *gIt;
            cellSpan = grid->getCellSpan();

            if( size < cellSpan )
            {
               cellSpan /= hierarchyFactor;
               if( size < cellSpan ) {
                  while( size < cellSpan ) cellSpan /= hierarchyFactor;
                  grid = new HashGrid( cellSpan * hierarchyFactor );
                  gridList_.insert( gIt, grid );
               }

               grid->add( body );
               body->setGrid( static_cast<void*>( grid ) );

               return;
            }
         }

         while( size >= cellSpan) cellSpan *= hierarchyFactor;
         grid = new HashGrid( cellSpan );
      }

      grid->add( body );
      body->setGrid( static_cast<void*>( grid ) );

      gridList_.push_back( grid );

      return;
   }

   // The body - which is infinite in size - is marked as being added to 'nonGridBodies_' by setting
   // the grid pointer to NULL and setting the cell-ID to '1'. Additionally, the hash value is used
   // to memorize the body's index position in the 'nonGridBodies_' vector.

   body->setGrid  ( nullptr );
   body->setHash  ( nonGridBodies_.size() );
   body->setCellId( 1 );

   nonGridBodies_.push_back( body );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Adds a body to the data structure (the usage of hash grids is not yet activated!).
 *
 * \param body The body that is about to be added to the data structure.
 * \return void
 *
 * As long as the number of bodies is less-or-equal than specified by \a gridActivationThreshold,
 * this function is called every time a new rigid body is added to the data structure. The body then
 * is stored in \a nonGridBodies_ - so that later during the detection step pairwise collision
 * checks for all bodies stored in \a nonGridBodies_ can be performed.
 *
 * However, the moment the threshold is exceeded, all bodies are removed from the array and finally
 * added to the grid hierarchy - thereby creating the initial set of hash grids which are adapted to
 * the just inserted bodies. Once this has happened, the usage of the hierarchical hash grids is
 * activated irrevocably, which means the coarse collision detection phase will continue to use the
 * hierarchical hash grids, even if removing bodies from the simulation might cause the total number
 * of bodies to again drop below the threshold.
 */
void HashGrids::addList( BodyID body )
{
   // If the threshold is exceeded ...
   if( nonGridBodies_.size() == gridActivationThreshold )
   {
      if( gridActivationThreshold > 0 )
      {
         BodyID* bodies = new BodyID[ gridActivationThreshold ];

         // ... all bodies stored in 'nonGridBodies_' are temporarily copied ...
         for( size_t i = 0; i < gridActivationThreshold; ++i ) {
            bodies[i] = nonGridBodies_[i];
         }

         // ... the 'nonGridBodies_' vector is cleared ...
         nonGridBodies_.clear();

         // ... the bodies are reinserted into the data structure - yet this time they are
         // inserted into the grid hierarchy (=> "addGrid()") ...
         for( size_t i = 0; i < gridActivationThreshold; ++i ) {
            addGrid( bodies[i] );
         }

         delete[] bodies;
      }

      addGrid( body );

      // ... and the usage of the hierarchical hash grids is activated irrevocably.
      gridActive_ = true;

      return;
   }

   // The body is marked as being added to 'nonGridBodies_' by setting the grid pointer to NULL and
   // setting the cell-ID to '1'. Additionally, the hash value is used to memorize the body's index
   // position in the 'nonGridBodies_' vector.

   body->setGrid  ( nullptr );
   body->setHash  ( nonGridBodies_.size() );
   body->setCellId( 1 );

   nonGridBodies_.push_back( body );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks whether a certain number is equal to a power of two.
 *
 * \param number The number that is about to be checked.
 * \return \a true if the number is equal to a power of two, \a false if not.
 *
 * This function is used to ensure that the number of cells in each coordinate dimension of a hash
 * grid is equal to a power of two so that the modulo calculations required for the hash value
 * computation can be replaced by bitwise AND operations.
 */
bool HashGrids::powerOfTwo( size_t number )
{
   return ( ( number > 0 ) && ( ( number & ( number - 1 ) ) == 0 ) );
}
//*************************************************************************************************


uint64_t HashGrids::intersectionTestCount = 0; //ToDo remove again
   
//=================================================================================================
//
//  CONSTANTS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief The initial number of cells in x-direction of a newly created hash grid.
 *
 * This value represents the initial number of cells of a newly created hash grid in x-direction.
 * The larger the value (i.e. the greater the number of cells of every newly created hash grid),
 * the more memory is required for the storage of the hash grid. Since the size of a hash grid is
 * increased at runtime in order to adapt to the number of currently inserted bodies, 16x16x16
 * is a suitable choice for the initial size of a newly created hash grid - it already consists
 * of four thousand cells, yet only requires a few hundred kilobytes of memory. Note that the
 * initial number of cells must both be greater-or-equal to 4 and equal to a power of two. Also
 * note that the initial number of cells does not necessarily have to be equal for all three
 * coordinate directions.
 */
const size_t HashGrids::xCellCount = 16;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The initial number of cells in y-direction of a newly created hash grid.
 *
 * This value represents the initial number of cells of a newly created hash grid in y-direction.
 * The larger the value (i.e. the greater the number of cells of every newly created hash grid),
 * the more memory is required for the storage of the hash grid. Since the size of a hash grid is
 * increased at runtime in order to adapt to the number of currently inserted bodies, 16x16x16
 * is a suitable choice for the initial size of a newly created hash grid - it already consists
 * of four thousand cells, yet only requires a few hundred kilobytes of memory. Note that the
 * initial number of cells must both be greater-or-equal to 4 and equal to a power of two. Also
 * note that the initial number of cells does not necessarily have to be equal for all three
 * coordinate directions.
 */
const size_t HashGrids::yCellCount = 16;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The initial number of cells in z-direction of a newly created hash grid.
 *
 * This value represents the initial number of cells of a newly created hash grid in z-direction.
 * The larger the value (i.e. the greater the number of cells of every newly created hash grid),
 * the more memory is required for the storage of the hash grid. Since the size of a hash grid is
 * increased at runtime in order to adapt to the number of currently inserted bodies, 16x16x16
 * is a suitable choice for the initial size of a newly created hash grid - it already consists
 * of four thousand cells, yet only requires a few hundred kilobytes of memory. Note that the
 * initial number of cells must both be greater-or-equal to 4 and equal to a power of two. Also
 * note that the initial number of cells does not necessarily have to be equal for all three
 * coordinate directions.
 */
const size_t HashGrids::zCellCount = 16;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The initial storage capacity of a newly created grid cell body container.
 *
 * This value specifies the initial storage capacity reserved for every grid cell body container,
 * i.e., the number of bodies that can initially be assigned to a grid cell with the need to
 * increase the storage capacity. The smaller this number, the more likely the storage capacity
 * of a body container must be increased, leading to potentially costly reallocation operations,
 * which generally involve the entire storage space to be copied to a new location. The greater
 * this number, the more memory is required. Rule of thumb:
 *
 *                        \f$ cellVectorSize = 2 \cdot hierarchyFactor^3 \f$
 */
const size_t HashGrids::cellVectorSize = 16;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The initial storage capacity of the grid-global vector.
 *
 * This value specifies the initial storage capacity of the grid-global vector that keeps track
 * of all body-occupied cells. As long as at least one body is assigned to a certain cell, this
 * cell is recorded in a grid-global list that keeps track of all body-occupied cells in order to
 * avoid iterating through all grid cells whenever all bodies that are stored in the grid need
 * to be addressed.
 */
const size_t HashGrids::occupiedCellsVectorSize = 256;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The minimal ratio of cells to bodies that must be maintained at any time.
 *
 * This \a minimalGridDensity specifies the minimal ratio of cells to bodies that is allowed
 * before a grid grows.\n
 * In order to handle an initially unknown and ultimately arbitrary number of bodies, each hash
 * grid, starting with a rather small number of cells at the time of its creation, must have the
 * ability to grow as new bodies are inserted. Therefore, if by inserting a body into a hash grid
 * the associated grid density - that is the ratio of cells to bodies - drops below the threshold
 * specified by \a minimalGridDensity, the number of cells in each coordinate direction is doubled
 * (thus the total number of grid cells is increased by a factor of 8).
 *
 * Possible settings: any integral value greater than 0.
 */
const size_t HashGrids::minimalGridDensity = 8;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Activation threshold for the hierarchical hash grids coarse collision detection algorithm.
 *
 * If the simulation only consists of a very small number of bodies, simply checking each body
 * against each other body proves to be faster than involving the far more complex mechanisms
 * of the hierarchical hash grids. In other words, despite its quadratic complexity, as long as
 * there are only a couple of bodies present in the simulation, the naive approach of conducting
 * pairwise checks for all existing bodies will always result in the best runtime performance. As
 * a consequence, a threshold is introduced, and as long as the number of bodies is less-or-equal
 * than specified by this threshold value, no hierarchy of hash grids is constructed and thus no
 * detection algorithm based on grids is applied.
 *
 * Possible settings: any integral value greater-or-equal to 0.
 */
const size_t HashGrids::gridActivationThreshold = 32;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The constant factor by which the cell size of any two successive grids differs.
 *
 * This factor specifies the size difference of two successive grid levels of the hierarchical
 * hash grids. The grid hierarchy is constructed such that the cell size of any two successive
 * grids differs by a constant factor - the hierarchy factor \a hierarchyFactor. As a result,
 * the cell size \f$ c_k \f$ of grid \f$ k \f$ can be expressed as:
 *
 *                          \f$ c_k = c_0 \cdot hierarchyFactor^k \f$.
 *
 * Note that the hierarchy does not have to be dense, which means, if not every valid cell size
 * that can be generated is required, some in-between grids are not created. Consequently, the
 * cell size of two successive grids differs by a factor of \f$ hierarchyFactor^x \f$, with x
 * being an integral value that is not necessarily equal to 1.
 *
 * The larger the ratio between the cell size of two successive grids, the more bodies are
 * potentially assigned to one single cell, but overall fewer grids have to be used. On the other
 * hand, the smaller the ratio between the cell size of two successive grids, the fewer bodies
 * are assigned to one single cell, but overall more grids have to be created. Hence, the number
 * of bodies that are stored in one single cell is inversely proportional to the number of grids
 * which are in use. Unfortunately, minimizing the number of bodies that are potentially assigned
 * to the same cell and at the same time also minimizing the number of grids in the hierarchy are
 * two opposing goals. In general - based on the evaluation of a number of different scenarios -
 * the best choice seems to be a hierarchy factor that is equal to 2.0.
 *
 * Possible settings: any floating point value that is greater than 1.0.
 */
const real_t HashGrids::hierarchyFactor = real_c(2);
//*************************************************************************************************

}  // namespace ccd

}  // namespace pe

}  // namespace walberla
