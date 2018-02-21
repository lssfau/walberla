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

#include <unordered_set>
#include <pe/raytracing/Ray.h>
#include <pe/utility/BodyCast.h>
#include <pe/raytracing/Intersects.h>

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
   
   static size_t intersectionTestCount; // ToDo remove again

//private: // toDo uncomment to change to private again
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
      
      template<typename BodyTuple>
      BodyID getRayIntersectingBody(const raytracing::Ray& ray, const AABB& blockAABB, real_t& t, Vec3& n) const;
      
      void clear();
      //@}
      //*******************************************************************************************


    //private: // ToDo uncomment to change to public again
      //**Utility functions************************************************************************
      /*!\name Utility functions */
      //@{
      void initializeNeighborOffsets();

      size_t hash( BodyID body ) const;
      size_t hashPoint(real_t x, real_t y, real_t z) const;

      void add   ( BodyID body, Cell* cell );
      void remove( BodyID body, Cell* cell );

      void enlarge();
      
      template<typename BodyTuple>
      bool getBodyIntersectionForCellCenter(real_t x, real_t y, real_t z, const AABB& blockAABB,
                                            const raytracing::Ray& ray, size_t c,
                                            real_t& t, Vec3& n, BodyID& body) const;
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
          void reloadBodies();
   //@}
   //**********************************************************************************************

   //**Implementation of ICCD interface ********************************************************
   virtual PossibleContacts& generatePossibleContacts( WcTimingTree* tt = NULL );
   void update(WcTimingTree* tt = NULL);
   
   bool active() const { return gridActive_; }
   
   template<typename BodyTuple>
   BodyID getClosestBodyIntersectingWithRay(const raytracing::Ray& ray, const AABB& blockAABB,
                                            real_t& t, Vec3& n);
   
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

/*!\brief Computes closest ray-body intersection of cell with center at point x,y,z and neighboring ones.
 *
 * \param x x value of center of cell to be processed.
 * \param y y value of center of cell to be processed.
 * \param z z value of center of cell to be processed.
 * \param blockAABB AABB of the block this grid corresponds to.
 * \param ray Ray being casted trough grid.
 * \param t Reference for calculated distance from ray origin.
 * \param n Reference for normal of intersection point.
 * \param body Reference for closest body.
 *
 * Inserts bodies of specified cell into bodiesContainer and additionally considers bodies in neighboring cells
 * in all negative coordinate direction to take bodies in consideration, which protrude from their
 * cell into the intersected cell (and thus, possibly intersect with the ray as well).
 */
template<typename BodyTuple>
bool HashGrids::HashGrid::getBodyIntersectionForCellCenter(real_t x, real_t y, real_t z, const AABB& blockAABB,
                                                           const raytracing::Ray& ray, size_t c,
                                                           real_t& t, Vec3& n, BodyID& body) const {
   //const Vec3& minCorner = blockAABB.minCorner();
   //const Vec3& maxCorner = blockAABB.maxCorner();
   
   real_t t_local;
   Vec3 n_local;
   bool intersected = false;
   
   raytracing::IntersectsFunctor intersectsFunc(ray, t_local, n_local);

   //const Vec3& dir = ray.getDirection();
   
   int minX = -1, minY = -1, minZ = -1;
   int maxX = 0, maxY = 0, maxZ = 0;
   
   /*if (c == 0) {
      minX = -2;
      maxY = 1;
      maxZ = 1;
   } else if (c == 1) {
      minY = -2;
      maxX = 1;
      maxZ = 1;
   } else if (c == 2) {
      minZ = -2;
      maxY = 1;
      maxX = 1;
   }*/

   for (int i = minX; i <= maxX; ++i) {
      const real_t x_shifted = x + i*cellSpan_;
      for (int j = minY; j <= maxY; ++j) {
         const real_t y_shifted = y + j*cellSpan_;
         for (int k = minZ; k <= maxZ; ++k) {
            const real_t z_shifted = z + k*cellSpan_;
            // commenting upper line of condition fixes a bug where objects with their minCorner
            // out of the blocks bounds would not get considered for intersections.
            //if (/*x_shifted > minCorner[0] && y_shifted > minCorner[1] && z_shifted > minCorner[2] &&*/
            //    /*x_shifted < maxCorner[0] && y_shifted < maxCorner[1] && z_shifted < maxCorner[2]*/ true) {
               size_t hash = hashPoint(x_shifted, y_shifted, z_shifted);
               
               const Cell& cell = cell_[hash];
               if (cell.bodies_ != NULL) {
                  for (const BodyID& cellBody: *cell.bodies_) {
                     HashGrids::intersectionTestCount++;
                     bool intersects = SingleCast<BodyTuple, raytracing::IntersectsFunctor, bool>::execute(cellBody, intersectsFunc);
                     if (intersects && t_local < t) {
                        body = cellBody;
                        t = t_local;
                        n = n_local;
                        intersected = true;
                     }
                  }
               }
            //}
         }
      }
   }
   
   return intersected;
}

/*!\brief Calculates ray-cell intersections and determines the closest body from the ray origin.
 *
 * \param ray Ray getting shot through this grid.
 * \param blockAABB AABB of the block this grid corresponds to.
 * \param t Reference for the distance.
 * \param n Reference for the intersetion normal.
 * \return BodyID of closest body, NULL if none found.
 *
 * This function calculates ray-cell intersections and the closest body in those cells. Also, neighboring
 * cells in all negative coordinate directions are regarded to take bodies in consideration, which
 * protrude from their cell into the intersected cell (and thus, possibly intersect with the ray as well).
 */
template<typename BodyTuple>
BodyID HashGrids::HashGrid::getRayIntersectingBody(const raytracing::Ray& ray, const AABB& blockAABB,
                                                   real_t& t_closest, Vec3& n_closest) const {
   //real_t inf = std::numeric_limits<real_t>::max();
   
   const Vec3& rayDirection = ray.getDirection();
   const Vec3& minCorner = blockAABB.minCorner();
   const Vec3& maxCorner = blockAABB.maxCorner();
   const Vec3& rayOrigin = ray.getOrigin();
   const Vec3& rayInvDirection = ray.getInvDirection();
   
   real_t xCenterOffsetFactor = rayDirection[0] < 0 ? -0.5 : 0.5;
   real_t yCenterOffsetFactor = rayDirection[1] < 0 ? -0.5 : 0.5;
   real_t zCenterOffsetFactor = rayDirection[2] < 0 ? -0.5 : 0.5;
   
   //real_t t_closest = inf;
   //Vec3 n_closest;
   BodyID body_closest = NULL;
   
   int i_first_intersection = -1;
   for (size_t i = 0; i <= xCellCount_; i++) {
      size_t i_dir = i;
      if (rayDirection[0] < 0) {
         i_dir = xCellCount_ - i;
      }
      real_t xValue = minCorner[0] + i_dir*cellSpan_; // calculate x value
      // get lambda in ray equation for which ray intersects with yz plane with this x value
      real_t lambda = (xValue - rayOrigin[0]) * rayInvDirection[0];
      // lambda is the distance from the ray origin to the cell intersection point
      if (lambda < 0) {
         // intersection is in rays negative direction
         continue;
      }
      if (rayDirection[0] >= 0 && lambda > t_closest) {
         // cell is too far away to generate useful body intersections
         // this condition only works for checks in positive direction
         break;
      }
      if (rayDirection[0] < 0 && lambda > t_closest+cellSpan_) {
         break;
      }
      real_t yValue = rayOrigin[1] + lambda * rayDirection[1]; // get y and z values for this intersection
      real_t zValue = rayOrigin[2] + lambda * rayDirection[2];
      if (yValue > maxCorner[1]+cellSpan_ || yValue < minCorner[1]-cellSpan_ ||
          zValue > maxCorner[2]+cellSpan_ || zValue < minCorner[2]-cellSpan_ ||
          lambda != lambda) {
         // calculated intersection with yz plane is out of the blocks bounds
         continue;
      }
      real_t centeredXValue = xValue + cellSpan_*xCenterOffsetFactor;
      if (centeredXValue > maxCorner[0] || centeredXValue < minCorner[0]) {
         // intersection was with one of the outer bounds of the block and in this direction no more cells exist
         continue;
      }
      // calculate the y and z values of the center of the cell
      real_t centeredYValue = (int((yValue - minCorner[1])*inverseCellSpan_) + real_t(0.5)) * cellSpan_ + minCorner[1];
      real_t centeredZValue = (int((zValue - minCorner[2])*inverseCellSpan_) + real_t(0.5)) * cellSpan_ + minCorner[2];
      bool intersected = getBodyIntersectionForCellCenter<BodyTuple>(centeredXValue, centeredYValue, centeredZValue, blockAABB,
                                                                     ray, 0,
                                                                     t_closest, n_closest, body_closest);
      if (intersected && i_first_intersection == -1) {
         i_first_intersection = int(i);
      }
      if (i_first_intersection != -1 && i_first_intersection == i-2) {
         break;
      }
   }
   
   i_first_intersection = -1;
   for (size_t i = 0; i < yCellCount_; i++) {
      size_t i_dir = i;
      if (rayDirection[1] < 0) {
         i_dir = yCellCount_ - i;
      }
      real_t yValue = minCorner[1] + i_dir*cellSpan_;
      real_t lambda = (yValue - rayOrigin[1]) * rayInvDirection[1];
      if (lambda < 0) {
         continue;
      }
      if (rayDirection[1] >= 0 && lambda > t_closest) {
         break;
      }
      if (rayDirection[1] < 0 && lambda > t_closest+cellSpan_) {
         break;
      }
      real_t xValue = rayOrigin[0] + lambda * rayDirection[0];
      real_t zValue = rayOrigin[2] + lambda * rayDirection[2];
      if (xValue >= maxCorner[0]+cellSpan_ || xValue < minCorner[0]-cellSpan_ ||
          zValue >= maxCorner[2]+cellSpan_ || zValue < minCorner[2]-cellSpan_ ||
          lambda != lambda) {
         continue;
      }
      real_t centeredXValue = (int((xValue - minCorner[0])*inverseCellSpan_) + real_t(0.5)) * cellSpan_ + minCorner[0];
      real_t centeredYValue = yValue + cellSpan_*yCenterOffsetFactor;
      if (centeredYValue > maxCorner[1] || centeredYValue < minCorner[1]) {
         continue;
      }
      real_t centeredZValue = (int((zValue - minCorner[2])*inverseCellSpan_) + real_t(0.5)) * cellSpan_ + minCorner[2];
      bool intersected = getBodyIntersectionForCellCenter<BodyTuple>(centeredXValue, centeredYValue, centeredZValue, blockAABB,
                                                                     ray, 1,
                                                                     t_closest, n_closest, body_closest);
      if (intersected && i_first_intersection == -1) {
         i_first_intersection = int(i);
      }
      if (i_first_intersection != -1 && i_first_intersection == i-2) {
         break;
      }
   }
   
   i_first_intersection = -1;
   for (size_t i = 0; i < zCellCount_; i++) {
      size_t i_dir = i;
      if (rayDirection[2] < 0) {
         i_dir = zCellCount_ - i;
      }
      real_t zValue = minCorner[2] + i_dir*cellSpan_;
      real_t lambda = (zValue - rayOrigin[2]) * rayInvDirection[2];
      if (lambda < 0) {
         continue;
      }
      if (rayDirection[2] >= 0 && lambda > t_closest) {
         break;
      }
      if (rayDirection[2] < 0 && lambda > t_closest+cellSpan_) {
         break;
      }
      real_t xValue = rayOrigin[0] + lambda * rayDirection[0];
      real_t yValue = rayOrigin[1] + lambda * rayDirection[1];
      if (xValue > maxCorner[0]+cellSpan_ || xValue < minCorner[0]-cellSpan_ ||
          yValue > maxCorner[1]+cellSpan_ || yValue < minCorner[1]-cellSpan_ ||
          lambda != lambda) {
         continue;
      }
      real_t centeredXValue = (int((xValue - minCorner[0])*inverseCellSpan_) + real_t(0.5)) * cellSpan_ + minCorner[0];
      real_t centeredYValue = (int((yValue - minCorner[1])*inverseCellSpan_) + real_t(0.5)) * cellSpan_ + minCorner[1];
      real_t centeredZValue = zValue + cellSpan_*zCenterOffsetFactor;
      if (centeredZValue > maxCorner[2] || centeredZValue < minCorner[2]) {
         continue;
      }
      bool intersected = getBodyIntersectionForCellCenter<BodyTuple>(centeredXValue, centeredYValue, centeredZValue, blockAABB,
                                                                     ray, 2,
                                                                     t_closest, n_closest, body_closest);
      if (intersected && i_first_intersection == -1) {
         i_first_intersection = int(i);
      }
      if (i_first_intersection != -1 && i_first_intersection == i-2) {
         break;
      }
   }
   
   //n = n_closest;
   //t = t_closest;
   
   return body_closest;
}

/*!\brief Determines the closest intersecting body in the underlying hash grids.
 *
 * \param ray Ray getting shot through this grid.
 * \param blockAABB AABB of the block this grid corresponds to.
 * \param t Reference for the distance.
 * \param n Reference for the intersetion normal.
 * \return BodyID of closest body, NULL if none found.
 */
template<typename BodyTuple>
BodyID HashGrids::getClosestBodyIntersectingWithRay(const raytracing::Ray& ray, const AABB& blockAABB, 
                                                    real_t& t, Vec3& n) {
   real_t inf = std::numeric_limits<real_t>::max();

   real_t t_closest = inf;
   Vec3 n_closest;
   BodyID body_closest = NULL;
   
   real_t t_local;
   Vec3 n_local;
   
   for(auto grid: gridList_) {
      BodyID body = grid->getRayIntersectingBody<BodyTuple>(ray, blockAABB, t_closest, n_closest);
      if (body != NULL) {
         body_closest = body;
      }
   }
   
   raytracing::IntersectsFunctor intersectsFunc(ray, t_local, n_local);
   for(auto body: nonGridBodies_) {
      bool intersects = SingleCast<BodyTuple, raytracing::IntersectsFunctor, bool>::execute(body, intersectsFunc);
      if (intersects && t_local < t_closest) {
         body_closest = body;
         t_closest = t_local;
         n_closest = n_local;
      }
   }
   
   n = n_closest;
   t = t_closest;
   
   return body_closest;
}


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
   //make sure to always check in the correct order (a<b)
   if (a->getSystemID() > b->getSystemID())
      std::swap(a, b);

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
