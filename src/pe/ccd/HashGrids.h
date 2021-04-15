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
   
   static uint64_t intersectionTestCount; // ToDo remove again
   
private:
   //**Type definitions****************************************************************************
   //! Vector for storing (handles to) rigid bodies.
   using BodyVector = std::vector<BodyID>;
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Implementation of the hash grid data structure.
   */
   class HashGrid
   {
    private:
      //**Type definitions*************************************************************************
      //! The signed integer type that is used for storing offsets to neighboring cells.
      using offset_t = long;
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
      using CellVector = std::vector<Cell *>;
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
      BodyID getRayIntersectingBody(const raytracing::Ray& ray, const AABB& blockAABB, real_t& t, Vec3& n,
                                    std::function<bool (const BodyID body)> isBodyVisibleFunc) const;
      
      static const int BLOCKCELL_NORMAL_INDETERMINATE = 3;

      template<typename BodyTuple>
      BodyID getBodyIntersectionForBlockCell(const Vector3<int32_t>& blockCell,
                                             const int8_t cellNormalAxis, const int8_t cellNormalDir,
                                             const raytracing::Ray& ray,
                                             real_t& t_closest, Vec3& n_closest,
                                             std::function<bool (const BodyID body)> isBodyVisibleFunc) const;
      
      void clear();
      //@}
      //*******************************************************************************************


    private:
      //**Utility functions************************************************************************
      /*!\name Utility functions */
      //@{
      void initializeNeighborOffsets();

      size_t hash( BodyID body ) const;
      size_t hashPoint(real_t x, real_t y, real_t z) const;

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
   public: size_t xyzCellCount_;  //!< Total number of allocated cells.

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
   using GridList = std::list<HashGrid *>;
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
   ~HashGrids() override;
   //@}
   //**********************************************************************************************

   //**Add/remove functions************************************************************************
   /*!\name Add/remove functions */
   //@{
   inline void add   ( BodyID body );
          void remove( BodyID body );
   inline int getObservedBodyCount() const override { return observedBodyCount_ + static_cast<int> (globalStorage_.size()); }
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
          void clear       ();
          void reloadBodies() override;
   //@}
   //**********************************************************************************************

   //**Implementation of ICCD interface ********************************************************
   PossibleContacts& generatePossibleContacts( WcTimingTree* tt = nullptr ) override;
   void update(WcTimingTree* tt = nullptr);
   
   bool active() const { return gridActive_; }
   
   template<typename BodyTuple>
   BodyID getClosestBodyIntersectingWithRay(const raytracing::Ray& ray, const AABB& blockAABB,
                                            real_t& t, Vec3& n,
                                            std::function<bool (const BodyID body)> isBodyVisibleFunc) const;
   
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

         if( nbBodies != nullptr )
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

         if( nbBodies != nullptr ) {
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
 * \param blockCell index of cell within block grid.
 * \param blockAABB AABB of the block this grid corresponds to.
 * \param ray Ray being casted trough grid.
 * \param t_closest Distance of closest object from ray origin. Will be updated if closer body found.
 * \param n_closest Normal of intersection point.
 *
 * \return BodyID of intersected body, NULL if no intersection found.
 *
 * Inserts bodies of specified cell into bodiesContainer and additionally considers bodies in neighboring cells
 * in all negative coordinate direction to take bodies in consideration, which protrude from their
 * cell into the intersected cell (and thus, possibly intersect with the ray as well).
 */
template<typename BodyTuple>
BodyID HashGrids::HashGrid::getBodyIntersectionForBlockCell(const Vector3<int32_t>& blockCell,
                                                            const int8_t cellNormalAxis, const int8_t cellNormalDir,
                                                            const raytracing::Ray& ray,
                                                            real_t& t_closest, Vec3& n_closest,
                                                            std::function<bool (const BodyID body)> isBodyVisibleFunc) const {
   real_t t_local;
   Vec3 n_local;
   BodyID body = nullptr;
   
   raytracing::IntersectsFunctor intersectsFunc(ray, t_local, n_local);
   
   // calculate center coordinates of the cell in the block
   real_t x = (real_c(blockCell[0]) + real_t(0.5)) * cellSpan_;
   real_t y = (real_c(blockCell[1]) + real_t(0.5)) * cellSpan_;
   real_t z = (real_c(blockCell[2]) + real_t(0.5)) * cellSpan_;
   
   // hash of cell in the hashgrid
   size_t cellHash = hashPoint(x, y, z);
   
   const Cell& centerCell = cell_[cellHash];
   
   std::vector<offset_t> relevantNeighborIndices;
   
   if (cellNormalAxis == 0) {
      if (cellNormalDir == -1) {
         relevantNeighborIndices = {2,5,8, 11,14,17, 20,23,26};
      } else {
         relevantNeighborIndices = {0,3,6, 9,12,15, 18,21,24};
      }
   } else if (cellNormalAxis == 1) {
      if (cellNormalDir == -1) {
         relevantNeighborIndices = {6,7,8, 15,16,17, 24,25,26};
      } else {
         relevantNeighborIndices = {0,1,2, 9,10,11, 18,19,20};
      }
   } else if (cellNormalAxis == 2) {
      if (cellNormalDir == -1) {
         relevantNeighborIndices = {18,19,20,21,22,23,24,25,26};
      } else {
         relevantNeighborIndices = {0,1,2,3,4,5,6,7,8};
      }
   } else {
      // cellNormalAxis == BLOCKCELL_NORMAL_INDETERMINATE
      relevantNeighborIndices = {
         0, 1, 2, 3, 4, 5, 6, 7, 8,
         9, 10, 11, 12, 13, 14, 15, 16, 17,
         18, 19, 20, 21, 22, 23, 24, 25, 26
      };
   }
   
#ifdef HASHGRIDS_RAYTRACING_CHECK_ALL_NEIGHBORS
   relevantNeighborIndices = {
      0, 1, 2, 3, 4, 5, 6, 7, 8,
      9, 10, 11, 12, 13, 14, 15, 16, 17,
      18, 19, 20, 21, 22, 23, 24, 25, 26
   };
#endif
   
   for (uint_t i = 0; i < relevantNeighborIndices.size(); ++i) {
      const offset_t neighborIndex = relevantNeighborIndices[i];
      const Cell* nbCell = &centerCell + centerCell.neighborOffset_[neighborIndex];
      const BodyVector* nbBodies = nbCell->bodies_;
      
      if (nbBodies != nullptr) {
         for (const BodyID& cellBody: *nbBodies) {
            if (cellBody->isRemote()) {
               continue;
            }
            if (!isBodyVisibleFunc(cellBody)) {
               continue;
            }
               
            HashGrids::intersectionTestCount++;
            bool intersects = SingleCast<BodyTuple, raytracing::IntersectsFunctor, bool>::execute(cellBody, intersectsFunc);
            if (intersects && t_local < t_closest) {
               body = cellBody;
               t_closest = t_local;
               n_closest = n_local;
            }
         }
      }
   }
      
   return body;
}

/*!\brief Calculates ray-cell intersections and determines the closest body to the ray origin.
 *
 * \param ray Ray getting shot through this grid.
 * \param blockAABB AABB of the block this grid corresponds to.
 * \param t Reference for the distance.
 * \param n Reference for the intersection normal.
 * \return BodyID of closest body, NULL if none found.
 *
 * This function calculates ray-cell intersections and the closest body in those cells. Also, neighboring
 * cells in all negative coordinate directions are regarded to take bodies in consideration, which
 * protrude from their cell into the intersected cell (and thus, possibly intersect with the ray as well).
 */
template<typename BodyTuple>
BodyID HashGrids::HashGrid::getRayIntersectingBody(const raytracing::Ray& ray, const AABB& blockAABB,
                                                   real_t& t_closest, Vec3& n_closest,
                                                   std::function<bool (const BodyID body)> isBodyVisibleFunc) const {
   const real_t realMax = std::numeric_limits<real_t>::max();
   
   BodyID body_local = nullptr;
   BodyID body_closest = nullptr;
   
   int32_t blockXCellCountMin = int32_c(blockAABB.xMin() * inverseCellSpan_) - 1;
   int32_t blockXCellCountMax = int32_c(std::ceil(blockAABB.xMax() * inverseCellSpan_)) + 1;
   int32_t blockYCellCountMin = int32_c(blockAABB.yMin() * inverseCellSpan_) - 1;
   int32_t blockYCellCountMax = int32_c(std::ceil(blockAABB.yMax() * inverseCellSpan_)) + 1;
   int32_t blockZCellCountMin = int32_c(blockAABB.zMin() * inverseCellSpan_) - 1;
   int32_t blockZCellCountMax = int32_c(std::ceil(blockAABB.zMax() * inverseCellSpan_)) + 1;
   
   Vec3 firstPoint;
   Vec3 firstPointCenteredInCell;
   real_t tRayOriginToGrid = 0;
   if (blockAABB.contains(ray.getOrigin(), cellSpan_)) {
      firstPoint = ray.getOrigin();
      firstPointCenteredInCell = firstPoint;
   } else {
      real_t t_start;
      Vec3 firstPointNormal;
      if (intersects(blockAABB, ray, t_start, cellSpan_, &firstPointNormal)) {
         firstPoint = ray.getOrigin() + ray.getDirection()*t_start;
         firstPointCenteredInCell = firstPoint - firstPointNormal * (cellSpan_/real_t(2));
         tRayOriginToGrid = (ray.getOrigin() - firstPoint).length();
      } else {
         return nullptr;
      }
   }
   
   Vector3<int32_t> firstCell(int32_c(std::floor(firstPointCenteredInCell[0]*inverseCellSpan_)),
                              int32_c(std::floor(firstPointCenteredInCell[1]*inverseCellSpan_)),
                              int32_c(std::floor(firstPointCenteredInCell[2]*inverseCellSpan_)));
   
   const int8_t stepX = ray.xDir() >= 0 ? 1 : -1;
   const int8_t stepY = ray.yDir() >= 0 ? 1 : -1;
   const int8_t stepZ = ray.zDir() >= 0 ? 1 : -1;
   
   Vec3 nearPoint((stepX >= 0) ? real_c(firstCell[0]+1)*cellSpan_-firstPoint[0] : firstPoint[0]-real_c(firstCell[0])*cellSpan_,
                  (stepY >= 0) ? real_c(firstCell[1]+1)*cellSpan_-firstPoint[1] : firstPoint[1]-real_c(firstCell[1])*cellSpan_,
                  (stepZ >= 0) ? real_c(firstCell[2]+1)*cellSpan_-firstPoint[2] : firstPoint[2]-real_c(firstCell[2])*cellSpan_);
   
   // tMax: distance along the ray to the next cell change in the axis direction
   real_t tMaxX = (!realIsEqual(ray.xDir(), 0)) ? std::abs(nearPoint[0]*ray.xInvDir()) : realMax;
   real_t tMaxY = (!realIsEqual(ray.yDir(), 0)) ? std::abs(nearPoint[1]*ray.yInvDir()) : realMax;
   real_t tMaxZ = (!realIsEqual(ray.zDir(), 0)) ? std::abs(nearPoint[2]*ray.zInvDir()) : realMax;
   
   // tDelta: how far along the ray must be moved to encounter a new cell in the specified axis direction
   real_t tDeltaX = (!realIsEqual(ray.xDir(), 0)) ? std::abs(cellSpan_*ray.xInvDir()) : realMax;
   real_t tDeltaY = (!realIsEqual(ray.yDir(), 0)) ? std::abs(cellSpan_*ray.yInvDir()) : realMax;
   real_t tDeltaZ = (!realIsEqual(ray.zDir(), 0)) ? std::abs(cellSpan_*ray.zInvDir()) : realMax;
   
   Vector3<int32_t> currentCell = firstCell;
   
   // First cell needs extra treatment, as it might lay out of the blocks upper bounds
   // due to the nature of how it is calculated: If the first point lies on a upper border
   // it maps to the cell "above" the grid.
   if (currentCell[0] < blockXCellCountMax &&
       currentCell[1] < blockYCellCountMax &&
       currentCell[2] < blockZCellCountMax) {
      body_local = getBodyIntersectionForBlockCell<BodyTuple>(currentCell, BLOCKCELL_NORMAL_INDETERMINATE, 0,
                                                              ray, t_closest, n_closest,
                                                              isBodyVisibleFunc);
      if (body_local != nullptr) {
         body_closest = body_local;
      }
   }
   
   int8_t blockCellNormalAxis;
   int8_t blockCellNormalDir;
   
   while (true) {
      if (tMaxX < tMaxY) {
         if (tMaxX < tMaxZ) {
#if !defined(HASHGRIDS_DISABLE_EARLY_CUTOFF)
            if (tRayOriginToGrid+tMaxX-tDeltaX > t_closest) {
               break;
            }
#endif
            tMaxX += tDeltaX;
            currentCell[0] += stepX;
            blockCellNormalAxis = 0;
            blockCellNormalDir = int8_c(-stepX);
            if (currentCell[0] >= blockXCellCountMax || currentCell[0] < blockXCellCountMin) {
               break;
            }
         } else {
#if !defined(HASHGRIDS_DISABLE_EARLY_CUTOFF)
            if (tRayOriginToGrid+tMaxZ-tDeltaZ > t_closest) {
               break;
            }
#endif
            tMaxZ += tDeltaZ;
            currentCell[2] += stepZ;
            blockCellNormalAxis = 2;
            blockCellNormalDir = int8_c(-stepZ);
            if (currentCell[2] >= blockZCellCountMax || currentCell[2] < blockZCellCountMin) {
               break;
            }
         }
      } else {
         if (tMaxY < tMaxZ) {
#if !defined(HASHGRIDS_DISABLE_EARLY_CUTOFF)
            if (tRayOriginToGrid+tMaxY-tDeltaY > t_closest) {
               break;
            }
#endif
            tMaxY += tDeltaY;
            currentCell[1] += stepY;
            blockCellNormalAxis = 1;
            blockCellNormalDir = int8_c(-stepY);
            if (currentCell[1] >= blockYCellCountMax || currentCell[1] < blockYCellCountMin) {
               break;
            }
         } else {
#if !defined(HASHGRIDS_DISABLE_EARLY_CUTOFF)
            if (tRayOriginToGrid+tMaxZ-tDeltaZ > t_closest) {
               break;
            }
#endif
            tMaxZ += tDeltaZ;
            currentCell[2] += stepZ;
            blockCellNormalAxis = 2;
            blockCellNormalDir = int8_c(-stepZ);
            if (currentCell[2] >= blockZCellCountMax || currentCell[2] < blockZCellCountMin) {
               break;
            }
         }
      }
      
      body_local = getBodyIntersectionForBlockCell<BodyTuple>(currentCell, blockCellNormalAxis, blockCellNormalDir,
                                                              ray, t_closest, n_closest,
                                                              isBodyVisibleFunc);
      if (body_local != nullptr) {
         body_closest = body_local;
      }
   }
   
   return body_closest;
}

/*!\brief Determines the closest intersecting body in the underlying hash grids.
 *
 * \param ray Ray getting shot through those grids.
 * \param blockAABB AABB of the block the grids correspond to.
 * \param t Reference for the distance.
 * \param n Reference for the intersection normal.
 * \return BodyID of closest body, NULL if none found.
 */
template<typename BodyTuple>
BodyID HashGrids::getClosestBodyIntersectingWithRay(const raytracing::Ray& ray, const AABB& blockAABB,
                                                    real_t& t, Vec3& n,
                                                    std::function<bool (const BodyID body)> isBodyVisibleFunc) const {
   const real_t realMax = std::numeric_limits<real_t>::max();

   BodyID body_closest = nullptr;
   real_t t_closest = realMax;
   Vec3 n_closest;
   
   BodyID body_local;
   real_t t_local = realMax;
   Vec3 n_local;
   
   raytracing::IntersectsFunctor intersectsFunc(ray, t_local, n_local);
   for(auto body: nonGridBodies_) {
      bool intersects = SingleCast<BodyTuple, raytracing::IntersectsFunctor, bool>::execute(body, intersectsFunc);
      if (intersects && t_local < t_closest) {
         body_closest = body;
         t_closest = t_local;
         n_closest = n_local;
      }
   }
   
   for(auto grid: gridList_) {
      body_local = grid->getRayIntersectingBody<BodyTuple>(ray, blockAABB, t_closest, n_closest, isBodyVisibleFunc);
      if (body_local != nullptr){
         body_closest = body_local;
      }
   }
   
   t = t_closest;
   n = n_closest;

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
