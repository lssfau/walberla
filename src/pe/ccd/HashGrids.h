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
//! \file ICCD.h
//! \author Klaus Iglberger
//! \author Florian Schornbaum
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include "ICCD.h"
#include <pe/rigidbody/BodyStorage.h>
#include <pe/Types.h>

#include <core/logging/Logging.h>
#include <core/debug/Debug.h>
#include <core/NonCopyable.h>

#include <cmath>
#include <list>
#include <sstream>
#include <vector>


namespace walberla{
namespace pe{

class BodyStorage;

namespace ccd {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of the 'Hierarchical Hash Grids' coarse collision detection algorithm.
 *
 * The 'Hierarchical Hash Grids' coarse collision detection algorithm is based on a spatial
 * partitioning scheme that uses a hierarchy of uniform, hash storage based grids. Uniform grids
 * subdivide the simulation space into cubic cells. A hierarchy of spatially superimposed grids
 * provides a set of grids with each grid subdividing the very same space, however possessing
 * different and thus unique-sized cells. Spatial partitioning is achieved by assigning every
 * rigid body to exactly one cell in exactly one grid - that is the grid with the smallest cells
 * that are larger than the rigid body (more precisely cells that are larger than the longest
 * edge of the rigid body's axis-aligned bounding box). As a consequence, the collision detection
 * is reduced to only checking body's that are assigned to spatially directly adjacent cells.
 *
 * Key features of this implementation are not only an <b>average-case computational complexity
 * of order O(N)</b> as well as a space complexity of order O(N), but also a short actual runtime
 * combined with low memory consumption. Moreover, the data structure representing the hierarchy
 * of grids has the ability to, at runtime, automatically and perfectly adapt to the bodies of the
 * underlying simulation. Also, arbitrary bodies can be added and removed in constant time, O(1).
 * Consequently, the coarse collision detection algorithm based on these hierarchical hash grids
 * is <b>especially well-suited for large-scale simulations</b> involving very large numbers of
 * bodies.
 *
 * For further information and a much more detailed explanation of this algorithm see
 *     http://www10.informatik.uni-erlangen.de/Publications/Theses/2009/Schornbaum_SA09.pdf
 */
class HashGrids : public ICCD
{
public:
   //**Constants **********************************************************************************
   // for documentation see cpp file
   static const size_t xCellCount;
   static const size_t yCellCount;
   static const size_t zCellCount;
   static const size_t cellVectorSize;
   static const size_t occupiedCellsVectorSize;
   static const size_t minimalGridDensity;
   static const size_t gridActivationThreshold;
   static const real_t hierarchyFactor;
   //**********************************************************************************************

private:
   //**Type definitions****************************************************************************
   //! Vector for storing (handles to) rigid bodies.
   typedef std::vector<BodyID>  BodyVector;
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Implementation of the hash grid data structure.
   */
   class HashGrid
   {
    private:
      //**Type definitions*************************************************************************
      //! The signed integer type that is used for storing offsets to neighboring cells.
      typedef long offset_t;
      //*******************************************************************************************

      //*******************************************************************************************
      /*!\brief Data structure for representing a cell in the hash grid (used by the 'Hierarchical
      //        Hash Grids' coarse collision detection algorithm).
      */
      struct Cell
      {
         BodyVector* bodies_;           /*!< \brief The cell's body container that stores (handles to)
                                             all the bodies that are assigned to this cell. */
                                        /*!< Note that only a pointer to such a body container is
                                             stored: in order to save memory, this container object
                                             is only allocated as long as there are bodies assigned
                                             to this cell. */
         offset_t*   neighborOffset_;   /*!< \brief Pointer to an array that is storing offsets that
                                             can be used to directly access all the neighboring
                                             cells in the hash grid. */
         size_t      occupiedCellsId_;  //!< The cell's index in the \a occupiedCells_ vector.
         int         lastNonFixedBody_; //!< marks the last body in the array which is not fixed
      };
      //*******************************************************************************************

      //**Type definitions*************************************************************************
      //! Vector for storing pointers to (body-occupied) cells.
      typedef std::vector<Cell*>  CellVector;
      //*******************************************************************************************

    public:
      //**Constructor******************************************************************************
      /*!\name Constructor */
      //@{
      explicit HashGrid( real_t cellSpan );
      //@}
      //*******************************************************************************************

      //**Destructor*******************************************************************************
      /*!\name Destructor */
      //@{
      ~HashGrid();
      //@}
      //*******************************************************************************************

      //**Getter functions*************************************************************************
      /*!\name Getter functions */
      //@{
      real_t getCellSpan() const { return cellSpan_; }  //!< Getter for \a cellSpan_.
      //@}
      //*******************************************************************************************

      //**Add/remove functions*********************************************************************
      /*!\name Add/remove functions */
      //@{
      inline void add   ( BodyID body );
      inline void remove( BodyID body );
      //@}
      //*******************************************************************************************

      //**Utility functions************************************************************************
      /*!\name Utility functions */
      //@{
      void update( BodyID body );

      template< typename Contacts >
      size_t process      ( BodyID** gridBodies, Contacts& contacts ) const;

      template< typename Contacts >
      void   processBodies( BodyID* bodies, size_t bodyCount, Contacts& contacts ) const;

      void clear();
      //@}
      //*******************************************************************************************

    private:
      //**Utility functions************************************************************************
      /*!\name Utility functions */
      //@{
      void initializeNeighborOffsets();

      size_t hash( BodyID body ) const;

      void add   ( BodyID body, Cell* cell );
      void remove( BodyID body, Cell* cell );

      void enlarge();
      //@}
      //*******************************************************************************************

      //**Member variables*************************************************************************
      /*!\name Member variables */
      //@{
      Cell* cell_;  //!< Linear array of cells representing the three-dimensional hash grid.

      size_t xCellCount_;  //!< Number of cells allocated in x-axis direction.
      size_t yCellCount_;  //!< Number of cells allocated in y-axis direction.
      size_t zCellCount_;  //!< Number of cells allocated in z-axis direction.

      size_t xHashMask_;  //!< Bit-mask required for the hash calculation in x-axis direction (\a xCellCount_ - 1).
      size_t yHashMask_;  //!< Bit-mask required for the hash calculation in y-axis direction (\a yCellCount_ - 1).
      size_t zHashMask_;  //!< Bit-mask required for the hash calculation in z-axis direction (\a zCellCount_ - 1).

      size_t xyCellCount_;   /*!< \brief Number of cells allocated in x-axis direction multiplied by
                                         the number of cells allocated in y-axis direction. */
      size_t xyzCellCount_;  //!< Total number of allocated cells.

      size_t enlargementThreshold_;  /*!< \brief The enlargement threshold - the moment the number
                                                 of assigned bodies exceeds this threshold, the
                                                 size of this hash grid is increased. */

      real_t cellSpan_;         //!< The grid's cell size (edge length of the underlying cubic grid cells).
      real_t inverseCellSpan_;  //!< The inverse cell size.
                              /*!< Required for replacing floating point divisions with multiplications
                                   during the hash computation. */

      CellVector occupiedCells_;  //!< A grid-global list that keeps track of all body-occupied cells.
                                  /*!< The list is required in order to avoid iterating through all
                                       grid cells whenever all bodies that are stored in the grid
                                       need to be addressed.  */

      size_t bodyCount_;  //!< Number of bodies assigned to this hash grid.

      offset_t stdNeighborOffset_[27];  /*!< \brief Array of offsets to the neighboring cells (valid
                                                    for all the inner cells of the hash grid). */
      //@}
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! List for storing all the hash grids that are in use.
   /*! This data structure is used to represent the grid hierarchy. All hash grids are stored in
       ascending order by the size of their cells. */
   typedef std::list<HashGrid*>  GridList;
   //**********************************************************************************************

public:
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
//   HashGrids( BodyStorage& bodystorage );
   HashGrids( BodyStorage& globalStorage, BodyStorage& bodystorage, BodyStorage& bodystorageShadowCopies );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~HashGrids();
   //@}
   //**********************************************************************************************

   //**Add/remove functions************************************************************************
   /*!\name Add/remove functions */
   //@{
   inline void add   ( BodyID body );
          void remove( BodyID body );
   inline int getObservedBodyCount() const { return observedBodyCount_ + static_cast<int> (globalStorage_.size()); }
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
          void clear       ();
   //@}
   //**********************************************************************************************

   //**Implementation of ICCD interface ********************************************************
   virtual PossibleContacts& generatePossibleContacts( WcTimingTree* tt = NULL );

   bool active() const { return gridActive_; }

protected:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename Contacts >
   static inline void collide( BodyID a, BodyID b, Contacts& contacts );
   //@}
   //**********************************************************************************************

private:
   //**Add functions*******************************************************************************
   /*!\name Add functions */
   //@{
   void addGrid( BodyID body );
   void addList( BodyID body );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static inline bool powerOfTwo( size_t number );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   BodyStorage& globalStorage_;    //!< Reference to the global body storage.
   BodyStorage& bodystorage_;    //!< Reference to the central body storage.
   BodyStorage& bodystorageShadowCopies_;    //!< Reference to the body storage containing body shadow copies.
   BodyVector   bodiesToAdd_;    //!< Vector of all bodies to be added to the HHG data structure.
                                 /*!< This vector contains all bodies that were registered with
                                      the coarse collision detector and have not yet been added
                                      to the hierarchical hash grids data structure. */
   BodyVector   nonGridBodies_;  //!< Vector of all unassigned rigid bodies.
                                 /*!< This vector contains all bodies that are not assigned to any
                                      grid in the hierarchy. A body is not assigned to any grid if
                                      the body is infinite in size or if the detection algorithm
                                      based on hierarchical hash grids has not yet been activated
                                      (see \a gridActive_). */
   GridList     gridList_;       //!< List of all grids that form the hierarchy.
                                 /*!< The grids in this list are sorted in ascending order by the
                                      size of their cells. */
   bool         gridActive_;     //!< Grid activation flag.
                                 /*!< This flag indicates whether the detection algorithm based on
                                      hierarchical hash grids is activated. If the simulation only
                                      consists of a very small number of bodies, each body is simply
                                      checked against every other body. This strategy proves to be
                                      faster than involving the far more complex mechanisms of the
                                      hierarchical hash grids. */
   int          observedBodyCount_;  /// number of bodies currently tracked by this hashgrid
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================


//*************************************************************************************************
/*!\brief Processes this hash grid for colliding bodies.
 *
 * \param gridBodies Linear array of (handles to) all the bodies stored in this hash grid.
 * \param contacts Contact container for the generated contacts.
 * \return The number of bodies stored in this grid - which is the number of bodies stored in \a gridBodies.
 *
 * This function generates all contacts between all rigid bodies that are assigned to this grid. The
 * contacts are added to the contact container \a contacts. Moreover, a linear array that contains
 * (handles to) all bodies that are stored in this grid is returned in order to being able to check
 * these bodies against other bodies that are stored in grids with larger sized cells.
 */
template< typename Contacts >  // Contact container type
size_t HashGrids::HashGrid::process( BodyID** gridBodies, Contacts& contacts ) const
{
   BodyID* bodies = new BodyID[ bodyCount_ ];
   *gridBodies    = bodies;

   // Iterate through all cells that are occupied by bodies (=> 'occupiedCells_') and ...
   for( typename CellVector::const_iterator cell = occupiedCells_.begin(); cell < occupiedCells_.end(); ++cell )
   {
      BodyVector* cellBodies = (*cell)->bodies_;

      // ... perform pairwise collision checks within each of these cells.
      for( auto aIt = cellBodies->begin(); aIt < cellBodies->end(); ++aIt ) {
         auto end = cellBodies->begin();
         if ((*aIt)->isFixed())
         {
            end = cellBodies->begin() + ((*cell)->lastNonFixedBody_ + 1);
         } else
         {
            end = cellBodies->end();
         }
         for( auto bIt = aIt + 1; bIt < end; ++bIt ) {
            WALBERLA_ASSERT( !((*aIt)->isFixed() && (*bIt)->isFixed()), "collision between two fixed bodies" );
            HashGrids::collide( *aIt, *bIt, contacts );
         }
         *(bodies++) = *aIt;
      }

      // Moreover, check all the bodies that are stored in the currently processed cell against all
      // bodies that are stored in the first half of all directly adjacent cells.
      for( unsigned int i = 0; i < 13; ++i )
      {
         Cell*       nbCell   = (*cell) + (*cell)->neighborOffset_[i];
         BodyVector* nbBodies = nbCell->bodies_;

         if( nbBodies != NULL )
         {
            for( auto aIt = cellBodies->begin(); aIt < cellBodies->end(); ++aIt ) {
               auto endNeighbour = nbBodies->begin();
               if ((*aIt)->isFixed())
               {
                  endNeighbour = nbBodies->begin() + (nbCell->lastNonFixedBody_ + 1);
               } else
               {
                  endNeighbour = nbBodies->end();
               }
               for( auto bIt = nbBodies->begin(); bIt < endNeighbour; ++bIt ) {
                  WALBERLA_ASSERT( !((*aIt)->isFixed() && (*bIt)->isFixed()), "collision between two fixed bodies" );
                  HashGrids::collide( *aIt, *bIt, contacts );
               }
            }
         }
      }
   }

   return bodyCount_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks all the bodies that are stored in \a bodies for collisions with the bodies that are
 *        stored in this grid.
 *
 * \param bodies Linear array of (handles to) all the bodies that are about to be checked for
 *        collisions with the bodies that are stored in this grid.
 * \param bodyCount The number of bodies that are stored in \a bodies.
 * \param contacts Contact container for the generated contacts.
 * \return void
 *
 * This function generates all contacts between the rigid bodies that are stored in \a bodies and
 * all the rigid bodies that are assigned to this grid. The contacts are added to the contact
 * container \a contacts.
 */
template< typename Contacts >  // Contact container type
void HashGrids::HashGrid::processBodies( BodyID* bodies, size_t bodyCount, Contacts& contacts ) const
{
   // For each body 'a' that is stored in 'bodies' ...
   for( BodyID* aIt = bodies; aIt < bodies + bodyCount; ++aIt )
   {
      // ... calculate the body's cell association (=> "hash()") within this hash grid and ...
      Cell* cell = cell_ + hash( *aIt );

      // ... check 'a' against every body that is stored in this or in any of the directly adjacent
      // cells. Note: one entry in the offset array of a cell is always referring back to the cell
      // itself. As a consequence, a specific cell X and all of its neighbors can be addressed by
      // simply iterating through all entries of X's offset array!
      for( unsigned int i = 0; i < 27; ++i )
      {
         Cell*       nbCell   = cell + cell->neighborOffset_[i];
         BodyVector* nbBodies = nbCell->bodies_;

         if( nbBodies != NULL ) {
            auto endNeighbour = nbBodies->begin();
            if ((*aIt)->isFixed())
            {
               endNeighbour = nbBodies->begin() + (nbCell->lastNonFixedBody_ + 1);
            } else
            {
               endNeighbour = nbBodies->end();
            }
            for( auto bIt = nbBodies->begin(); bIt != endNeighbour; ++bIt ) {
               WALBERLA_ASSERT( !((*aIt)->isFixed() && (*bIt)->isFixed()), "collision between two fixed bodies" );
               HashGrids::collide( *aIt, *bIt, contacts );
            }
         }
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Checks whether two bodies are colliding or not, and if so generates the related points
 *        of contact.
 *
 * \param a The first body.
 * \param b The second body.
 * \param contacts Contact container for the generated contacts.
 * \return void
 */
template< typename Contacts >  // Contact container type
void HashGrids::collide( BodyID a, BodyID b, Contacts& contacts )
{
   if( ( !a->hasInfiniteMass() || !b->hasInfiniteMass() ) &&        // Ignoring contacts between two fixed bodies
       ( a->getAABB().intersects( b->getAABB() ) ) )  // Testing for overlapping bounding boxes
   {
      //FineDetector::collide( a, b, contacts );
      contacts.push_back( std::make_pair(a, b) );
   }
}
//*************************************************************************************************


}  // namespace ccd

}  // namespace pe

}  // namespace walberla
