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
//! \file HashGridsBodyTrait.h
//! \author Klaus Iglberger
//! \author Florian Schornbaum
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla {
namespace pe {
namespace ccd {

//=================================================================================================
//
//  SPECIALIZATION FOR THE HIERARCHICAL HASH GRID ALGORITHM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Specialization of the BodyTrait class template for the 'hierarchical hash-grid' algorithm.
 *
 * This specialization of the BodyTrait class template adapts rigid bodies to the hierarchical hash-grid
 * coarse collision detection algorithm.
 */
class HashGridsBodyTrait
{
protected:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   explicit inline HashGridsBodyTrait( );
   //@}
   //**********************************************************************************************

public:
   //**Get functions*******************************************************************************
   /*!\name Get functions */
   //@{
   inline void*  getGrid    () const;
   inline size_t getHash    () const;
   inline size_t getCellId  () const;
   //@}
   //**********************************************************************************************

   //**Set functions*******************************************************************************
   /*!\name Set functions */
   //@{
   inline void setGrid  ( void*  grid );
   inline void setHash  ( size_t hash );
   inline void setCellId( size_t cell );
   //@}
   //**********************************************************************************************

protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   void*  grid_;    //!< Pointer to the hash grid this rigid body is currently assigned to.
   size_t hash_;    //!< Current hash value of this rigid body.
   size_t cellId_;  //!< The body's index in the body container of the grid cell this rigid body is currently assigned to.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the BodyTrait<HashGrids> specialization.
 *
 * \param body The rigid body containing this bounding box.
 */
inline HashGridsBodyTrait::HashGridsBodyTrait( )
   : grid_  (nullptr)  // Pointer to the hash grid this rigid body is currently assigned to
   , hash_  (0)  // Current hash value of this rigid body
   , cellId_(0)  // The body's index in the body container of the grid cell this rigid body is currently assigned to
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the grid the rigid body is currently assigned to.
 *
 * \return The grid the rigid body is currently assigned to.
 */
inline void* HashGridsBodyTrait::getGrid() const
{
   return grid_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current hash value of the rigid body.
 *
 * \return The current hash value of the rigid body.
 */
inline size_t HashGridsBodyTrait::getHash() const
{
   return hash_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current body container index within the cell the body is currently assigned to.
 *
 * \return The current body container index.
 */
inline size_t HashGridsBodyTrait::getCellId() const
{
   return cellId_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the grid the rigid body is associated with.
 *
 * \param grid The grid the rigid body is assigned to.
 * \return void
 */
inline void HashGridsBodyTrait::setGrid( void* grid )
{
   grid_ = grid;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the hash value of the rigid body.
 *
 * \param hash The new hash value of the rigid body.
 * \return void
 */
inline void HashGridsBodyTrait::setHash( size_t hash )
{
   hash_ = hash;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the body container index within the cell the body is currently assigned to.
 *
 * \param cell The new body container index.
 * \return void
 */
inline void HashGridsBodyTrait::setCellId( size_t cell )
{
   cellId_ = cell;
}
//*************************************************************************************************


} // namespace ccd

} // namespace pe

} // namespace walberla
