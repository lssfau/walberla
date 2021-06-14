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
//! \file
//! \author Christoph Rettinger christoph.rettinger@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>
#include <mesa_pd/data/ParticleStorage.h>

#include <core/Abort.h>
#include <core/debug/Debug.h>

#include <atomic>
#include <cmath>
#include <vector>
#include <array>

namespace walberla {
namespace mesa_pd {
namespace data {

//*************************************************************************************************
/*!\brief Implementation of the 'Hierarchical Hash Grids' coarse collision detection algorithm.
 *
 * This is a port of the pe implementation (src/pe/ccd/HashGrids.h) for mesa_pd.
 *
 * The 'Hierarchical Hash Grids' coarse collision detection algorithm is based on a spatial
 * partitioning scheme that uses a hierarchy of uniform, hash storage based grids. Uniform grids
 * subdivide the simulation space into cubic cells. A hierarchy of spatially superimposed grids
 * provides a set of grids with each grid subdividing the very same space, however possessing
 * different and thus unique-sized cells. Spatial partitioning is achieved by assigning every
 * rigid particle to exactly one cell in exactly one grid - that is the grid with the smallest cells
 * that are larger than the rigid particle (more precisely cells that are larger than the longest
 * edge of the rigid particle's axis-aligned bounding box). As a consequence, the collision detection
 * is reduced to only checking particle's that are assigned to spatially directly adjacent cells.
 *
 * Key features of this implementation are not only an <b>average-case computational complexity
 * of order O(N)</b> as well as a space complexity of order O(N), but also a short actual runtime
 * combined with low memory consumption. Moreover, the data structure representing the hierarchy
 * of grids has the ability to, at runtime, automatically and perfectly adapt to the particles of the
 * underlying simulation. Also, arbitrary particles can be added and removed in constant time, O(1).
 * Consequently, the coarse collision detection algorithm based on these hierarchical hash grids
 * is <b>especially well-suited for large-scale simulations</b> involving very large numbers of
 * particles.
 *
 * For further information and a much more detailed explanation of this algorithm see
 *     https://www10.cs.fau.de/publications/theses/2009/Schornbaum_SA_2009.pdf
 */
class HashGrids
{

public:
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
    * increased at runtime in order to adapt to the number of currently inserted particles, 16x16x16
    * is a suitable choice for the initial size of a newly created hash grid - it already consists
    * of four thousand cells, yet only requires a few hundred kilobytes of memory. Note that the
    * initial number of cells must both be greater-or-equal to 4 and equal to a power of two. Also
    * note that the initial number of cells does not necessarily have to be equal for all three
    * coordinate directions.
    */
   static constexpr size_t xCellCount = 16;
   //*************************************************************************************************


   //*************************************************************************************************
   /*!\brief The initial number of cells in y-direction of a newly created hash grid.
    *
    * See HashGrids::xCellCount for more infos.
    */
   static constexpr size_t yCellCount = 16;
   //*************************************************************************************************


   //*************************************************************************************************
   /*!\brief The initial number of cells in z-direction of a newly created hash grid.
    *
    * See HashGrids::xCellCount for more infos.
    */
   static constexpr size_t zCellCount = 16;
   //*************************************************************************************************


   //*************************************************************************************************
   /*!\brief The initial storage capacity of a newly created grid cell particle container.
    *
    * This value specifies the initial storage capacity reserved for every grid cell particle container,
    * i.e., the number of particles that can initially be assigned to a grid cell with the need to
    * increase the storage capacity. The smaller this number, the more likely the storage capacity
    * of a particle container must be increased, leading to potentially costly reallocation operations,
    * which generally involve the entire storage space to be copied to a new location. The greater
    * this number, the more memory is required. Rule of thumb:
    *
    *                        \f$ cellVectorSize = 2 \cdot hierarchyFactor^3 \f$
    */
   static constexpr size_t cellVectorSize = 16;
   //*************************************************************************************************


   //*************************************************************************************************
   /*!\brief The initial storage capacity of the grid-global vector.
    *
    * This value specifies the initial storage capacity of the grid-global vector that keeps track
    * of all particle-occupied cells. As long as at least one particle is assigned to a certain cell, this
    * cell is recorded in a grid-global list that keeps track of all particle-occupied cells in order to
    * avoid iterating through all grid cells whenever all particles that are stored in the grid need
    * to be addressed.
    */
   static constexpr size_t occupiedCellsVectorSize = 256;
   //*************************************************************************************************


   //*************************************************************************************************
   /*!\brief The minimal ratio of cells to particles that must be maintained at any time.
    *
    * This \a minimalGridDensity specifies the minimal ratio of cells to particles that is allowed
    * before a grid grows.\n
    * In order to handle an initially unknown and ultimately arbitrary number of particles, each hash
    * grid, starting with a rather small number of cells at the time of its creation, must have the
    * ability to grow as new particles are inserted. Therefore, if by inserting a particle into a hash grid
    * the associated grid density - that is the ratio of cells to particles - drops below the threshold
    * specified by \a minimalGridDensity, the number of cells in each coordinate direction is doubled
    * (thus the total number of grid cells is increased by a factor of 8).
    *
    * Possible settings: any integral value greater than 0.
    */
   static constexpr size_t minimalGridDensity = 8;
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
    * The larger the ratio between the cell size of two successive grids, the more particles are
    * potentially assigned to one single cell, but overall fewer grids have to be used. On the other
    * hand, the smaller the ratio between the cell size of two successive grids, the fewer particles
    * are assigned to one single cell, but overall more grids have to be created. Hence, the number
    * of particles that are stored in one single cell is inversely proportional to the number of grids
    * which are in use. Unfortunately, minimizing the number of particles that are potentially assigned
    * to the same cell and at the same time also minimizing the number of grids in the hierarchy are
    * two opposing goals. In general - based on the evaluation of a number of different scenarios -
    * the best choice seems to be a hierarchy factor that is equal to 2.0.
    *
    * Possible settings: any floating point value that is greater than 1.0.
    */
   static constexpr real_t hierarchyFactor = real_t(2);
   //*************************************************************************************************

private:
   //**Type definitions****************************************************************************
   //! Vector for storing (handles to) rigid particles.
   using ParticleIdxVector = std::vector<size_t>;
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
         ParticleIdxVector* particles_;           /*!< \brief The cell's particle container that stores (handles to)
                                             all the particles that are assigned to this cell. */
         /*!< Note that only a pointer to such a particle container is
              stored: in order to save memory, this container object
              is only allocated as long as there are particles assigned
              to this cell. */
         offset_t*   neighborOffset_;   /*!< \brief Pointer to an array that is storing offsets that
                                             can be used to directly access all the neighboring
                                             cells in the hash grid. */
         size_t      occupiedCellsId_;  //!< The cell's index in the \a occupiedCells_ vector.
      };
      //*******************************************************************************************

      //**Type definitions*************************************************************************
      //! Vector for storing pointers to (particle-occupied) cells.
      using CellVector = std::vector<Cell *>;
      //*******************************************************************************************

   public:

      explicit HashGrid( real_t cellSpan );

      ~HashGrid();

      real_t getCellSpan() const { return cellSpan_; }

      template <typename Accessor>
      inline void addParticle( size_t p_idx, Accessor& ac );

      template <typename Selector, typename Accessor, typename Func, typename... Args>
      void checkEachParticlePairHalfAndStore( ParticleIdxVector& particlesOnGrid,
                                              const bool openmp, const Selector& selector, Accessor& ac, Func&& func, Args&&... args ) const;

      template <typename Selector, typename Accessor, typename Func, typename... Args>
      void checkAgainstVectorEachParticlePairHalf( const ParticleIdxVector& particlesOnGrid,
                                                   const bool openmp, const Selector& selector, Accessor& ac, Func&& func, Args&&... args ) const;

      void clear();


   private:

      void initializeNeighborOffsets();

      template <typename Accessor>
      size_t hashOfParticle( size_t p_idx, Accessor& ac ) const;
      size_t hashPoint(real_t x, real_t y, real_t z) const;

      void addParticleToCell( size_t p_idx, Cell* cell );

      template <typename Accessor>
      void enlarge(Accessor& ac);

      Cell* cell_;  //!< Linear array of cells representing the three-dimensional hash grid.

      size_t xCellCount_;  //!< Number of cells allocated in x-axis direction.
      size_t yCellCount_;  //!< Number of cells allocated in y-axis direction.
      size_t zCellCount_;  //!< Number of cells allocated in z-axis direction.

      size_t xHashMask_;  //!< Bit-mask required for the hash calculation in x-axis direction (\a xCellCount_ - 1).
      size_t yHashMask_;  //!< Bit-mask required for the hash calculation in y-axis direction (\a yCellCount_ - 1).
      size_t zHashMask_;  //!< Bit-mask required for the hash calculation in z-axis direction (\a zCellCount_ - 1).

      size_t xyCellCount_;   /*!< \brief Number of cells allocated in x-axis direction multiplied by
                                         the number of cells allocated in y-axis direction. */
   public:
      size_t xyzCellCount_;  //!< Total number of allocated cells.

      size_t enlargementThreshold_;  /*!< \brief The enlargement threshold - the moment the number
                                                 of assigned particles exceeds this threshold, the
                                                 size of this hash grid is increased. */

      real_t cellSpan_;         //!< The grid's cell size (edge length of the underlying cubic grid cells).
      real_t inverseCellSpan_;  //!< The inverse cell size.
      /*!< Required for replacing floating point divisions with multiplications
           during the hash computation. */

      CellVector occupiedCells_;  //!< A grid-global list that keeps track of all particle-occupied cells.
      /*!< The list is required in order to avoid iterating through all
           grid cells whenever all particles that are stored in the grid
           need to be addressed.  */

      size_t particleCount_;  //!< Number of particles assigned to this hash grid.

      offset_t stdNeighborOffset_[27];  /*!< \brief Array of offsets to the neighboring cells (valid
                                                    for all the inner cells of the hash grid). */

   };

   //**Type definitions****************************************************************************
   //! List for storing all the hash grids that are in use.
   /*! This data structure is used to represent the grid hierarchy. All hash grids are stored in
       ascending order by the size of their cells. */
   using GridList = std::list<shared_ptr<HashGrid>>;

public:
   explicit HashGrids() = default;
   ~HashGrids() = default;

   // initializes Hash Grid with given particle
   template <typename Accessor>
   inline void operator()( size_t p_idx, Accessor& ac );

   void clear();
   void clearAll();

   /**
    * Calls the provided functor \p func for all particle pairs.
    *
    * Additional arguments can be provided. No pairs with twice the same particle are generated.
    * No pair is called twice!
    * Call syntax for the provided functor
    * \code
    * func( *this, i, j, std::forward<Args>(args)... );
    * \endcode
    * \param openmp enables/disables OpenMP parallelization of the kernel call
    */
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   void forEachParticlePairHalf(const bool openmp,
                                const Selector& selector,
                                Accessor& acForLC,
                                Func&& func,
                                Args&&... args) const;

private:

   inline void addInfiniteParticle(size_t p_idx);

   ParticleIdxVector infiniteParticles_;
   GridList     gridList_;       //!< List of all grids that form the hierarchy.
   /*!< The grids in this list are sorted in ascending order by the
        size of their cells. */

};

// HashGrids::HashGrid member function implementations

HashGrids::HashGrid::HashGrid( real_t cellSpan )
{
   // Initialization of all member variables and ...
   xCellCount_   = math::uintIsPowerOfTwo( xCellCount ) ? xCellCount : 16;
   yCellCount_   = math::uintIsPowerOfTwo( yCellCount ) ? yCellCount : 16;
   zCellCount_   = math::uintIsPowerOfTwo( zCellCount ) ? zCellCount : 16;

   xHashMask_    = xCellCount_ - 1;
   yHashMask_    = yCellCount_ - 1;
   zHashMask_    = zCellCount_ - 1;

   xyCellCount_  = xCellCount_ * yCellCount_;
   xyzCellCount_ = xyCellCount_ * zCellCount_;

   enlargementThreshold_ = xyzCellCount_ / minimalGridDensity;

   // ... allocation of the linear array that is representing the hash grid.
   cell_ = new Cell[ xyzCellCount_ ];

   // Initialization of each cell - i.e., initially setting the pointer to the particle container to
   // NULL (=> no particles are assigned to this hash grid yet!) and ...
   for( Cell* c  = cell_; c < cell_ + xyzCellCount_; ++c ) {
      c->particles_ = nullptr;
   }
   // ... setting up the neighborhood relationship (using the offset array).
   initializeNeighborOffsets();

   cellSpan_        = cellSpan;
   inverseCellSpan_ = real_c( 1 ) / cellSpan;

   occupiedCells_.reserve( occupiedCellsVectorSize );

   particleCount_ = 0;
}

HashGrids::HashGrid::~HashGrid()
{
   clear();

   for( Cell* c  = cell_; c < cell_ + xyzCellCount_; ++c ) {
      if( c->neighborOffset_ != stdNeighborOffset_ ) delete[] c->neighborOffset_;
   }
   delete[] cell_;
}

void HashGrids::HashGrid::clear()
{
   for( auto cellIt = occupiedCells_.begin(); cellIt < occupiedCells_.end(); ++cellIt ) {
      delete (*cellIt)->particles_;
      (*cellIt)->particles_ = nullptr;
   }
   occupiedCells_.clear();
   particleCount_ = 0;
}

//*************************************************************************************************
/*!\brief Sets up the neighborhood relationships for all grid cells.
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
/*!\brief Adds a particle to this hash grid.
 *
 * This function is called every time a new rigid particle is added to this grid. If adding the particle
 * will cause the total number of particles assigned to this grid to exceed the enlargement threshold,
 * the size of this hash grid is increased in order to maintain a fixed minimal grid density (=
 * ratio of cells to particles). This function may also be called during the update phase.
 */
template <typename Accessor>
void HashGrids::HashGrid::addParticle( size_t p_idx, Accessor& ac )
{
   // If adding the particle will cause the total number of particles assigned to this grid to exceed the
   // enlargement threshold, the size of this hash grid must be increased.
   if( particleCount_ == enlargementThreshold_ ) enlarge(ac);

   // Calculate (and store) the hash value (= the particle's cell association) and ...
   size_t hash = hashOfParticle( p_idx, ac );

   // ... insert the particle into the corresponding cell.
   Cell* cell = cell_ + hash;
   addParticleToCell( p_idx, cell );

   ++particleCount_;
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Adds a particle to a specific cell in this hash grid.
 */
void HashGrids::HashGrid::addParticleToCell( size_t p_idx, Cell* cell )
{
   // If this cell is already occupied by other particles, which means the pointer to the particle
   // container holds a valid address and thus the container itself is properly initialized, then
   // the particle is simply added to this already existing particle container.

   // If, however, the cell is still empty, then the object container, first of all, must be created
   // (i.e., allocated) and properly initialized (i.e., sufficient initial storage capacity must be
   // reserved). Furthermore, the cell must be inserted into the grid-global vector 'occupiedCells_'
   // in which all cells that are currently occupied by particles are recorded.
   if( cell->particles_ == nullptr )
   {
      cell->particles_ = new ParticleIdxVector ;
      cell->particles_->reserve( cellVectorSize );

      cell->occupiedCellsId_ = occupiedCells_.size();
      occupiedCells_.push_back( cell );
   }
   cell->particles_->push_back( p_idx );
}
//*************************************************************************************************



//*************************************************************************************************
/*!\brief Computes the hash value (= cell association) of a given particle.
 */
template <typename Accessor>
size_t HashGrids::HashGrid::hashOfParticle( size_t p_idx, Accessor& ac ) const
{
   auto particlePosition = ac.getPosition(p_idx);
   return hashPoint(particlePosition[0], particlePosition[1], particlePosition[2]);
}
//*************************************************************************************************


/*!\brief Computes the hash for a given point.
 *
 * The hash calculation uses modulo operations in order to spatially map entire blocks of connected
 * cells to the origin of the coordinate system. This block of cells at the origin of the coordinate
 * system that is filled with all the particles of the simulation is referred to as the hash grid. The
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
/*!\brief Increases the number of cells used by this hash grid.
 *
 * In order to handle an initially unknown and ultimately arbitrary number of particles, the hash grid,
 * starting with a rather small number of cells at the time of its creation, must have the ability
 * to grow as new particles are inserted. Therefore, if by inserting a particle into this hash grid the
 * associated grid density - that is the ratio of cells to particles - drops below the threshold
 * specified by \a minimalGridDensity, or in other words, if the number of particles grows larger than
 * specified by \a enlargementThreshold_, the number of cells in each coordinate direction is
 * doubled (thus the total number of grid cells is increased by a factor of 8).
 */
template <typename Accessor>
void HashGrids::HashGrid::enlarge(Accessor& ac)
{
   ParticleIdxVector particles(particleCount_);
   size_t curIdx = 0;

   // All particles that are assigned to this grid are temporarily removed, ...
   for( auto cellIt = occupiedCells_.begin(); cellIt < occupiedCells_.end(); ++cellIt ) {
      auto cellParticles = (*cellIt)->particles_;
      for( auto eIt = cellParticles->begin(); eIt < cellParticles->end(); ++eIt ) {
         particles[curIdx++] = *eIt;
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
      c->particles_ = nullptr;
   }
   initializeNeighborOffsets();

   // ... all previously removed particles are reinserted.
   for( auto& p : particles ) {
      addParticle(p, ac);
   }
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Processes this hash grid for colliding particles.
 *
 * This function generates all contacts between all rigid particles that are assigned to this grid.
 * Moreover, a linear array that contains
 * (handles to) all particles that are stored in this grid is returned in order to being able to check
 * these particles against other particles that are stored in grids with larger sized cells.
 */
template <typename Selector, typename Accessor, typename Func, typename... Args>
void HashGrids::HashGrid::checkEachParticlePairHalfAndStore( ParticleIdxVector& particlesOnGrid,
                                                             const bool openmp, const Selector& selector, Accessor& ac, Func&& func, Args&&... args ) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");
   WALBERLA_UNUSED(openmp);

   particlesOnGrid = ParticleIdxVector(particleCount_);
   size_t curVecIdx = 0;

   // Iterate through all cells that are occupied by particles (=> 'occupiedCells_') and ...
   for( typename CellVector::const_iterator cell = occupiedCells_.begin(); cell < occupiedCells_.end(); ++cell )
   {
      ParticleIdxVector* cellParticles = (*cell)->particles_;

      // ... perform pairwise interactions within each of these cells.
      for( auto aIt = cellParticles->begin(); aIt < cellParticles->end(); ++aIt ) {
         for( auto bIt = aIt + 1; bIt < cellParticles->end(); ++bIt ) {
            if (selector(uint_c(*aIt), uint_c(*bIt), ac))
            {
               func(uint_c(*aIt), uint_c(*bIt), std::forward<Args>(args)...);
            }
         }
         particlesOnGrid[curVecIdx++] = *aIt;
      }

      // Moreover, check all the particles that are stored in the currently processed cell against all
      // particles that are stored in the first half of all directly adjacent cells.
      for( unsigned int i = 0; i < 13; ++i )
      {
         Cell*       nbCell   = (*cell) + (*cell)->neighborOffset_[i];
         ParticleIdxVector* nbParticles = nbCell->particles_;

         if( nbParticles != nullptr )
         {
            for( auto aIt = cellParticles->begin(); aIt < cellParticles->end(); ++aIt ) {
               for( auto bIt = nbParticles->begin(); bIt < nbParticles->end(); ++bIt ) {
                  if (selector(uint_c(*aIt), uint_c(*bIt), ac))
                  {
                     func(uint_c(*aIt), uint_c(*bIt), std::forward<Args>(args)...);
                  }
               }
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks all the particles that are stored in \a particles for collisions with the particles that are
 *        stored in this grid.
 *
 * This function generates all contacts between the rigid particles that are stored in \a particles and
 * all the rigid particles that are assigned to this grid. The contacts are added to the contact
 * container \a contacts.
 */
template <typename Selector, typename Accessor, typename Func, typename... Args>
void HashGrids::HashGrid::checkAgainstVectorEachParticlePairHalf( const ParticleIdxVector& particlesOnGrid,
                                                                  const bool openmp, const Selector& selector, Accessor& ac, Func&& func, Args&&... args ) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");
   WALBERLA_UNUSED(openmp);

   // For each particle 'a' that is stored in 'particles' ...
   for( auto& aIt : particlesOnGrid )
   {
      // ... calculate the particle's cell association (=> "hash()") within this hash grid and ...
      Cell* cell = cell_ + hashOfParticle( aIt, ac );

      // ... check 'a' against every particle that is stored in this or in any of the directly adjacent
      // cells. Note: one entry in the offset array of a cell is always referring back to the cell
      // itself. As a consequence, a specific cell X and all of its neighbors can be addressed by
      // simply iterating through all entries of X's offset array!
      for( unsigned int i = 0; i < 27; ++i )
      {
         Cell* nbCell = cell + cell->neighborOffset_[i];
         auto nbParticles = nbCell->particles_;

         if( nbParticles != nullptr ) {
            for( auto& bIt : *nbParticles ) {
               if (selector(uint_c(aIt), uint_c(bIt), ac))
               {
                  func(uint_c(aIt), uint_c(bIt), std::forward<Args>(args)...);
               }
            }
         }
      }
   }
}



//*************************************************************************************************
//
// implementation of HashGrids member functions
//
//*************************************************************************************************

// clear all particles and grids
void HashGrids::clearAll()
{
   gridList_.clear();
   infiniteParticles_.clear();
}

// clear only all particles from the hash grids, but maintain the overall grid hierarchy
// useful for "updating" the data structure in each time step, but also prevents clean-up of unnecessary grids
void HashGrids::clear()
{
   for( auto gridIt = gridList_.begin(); gridIt != gridList_.end(); ++gridIt ) {
      (*gridIt)->clear();
   }
   infiniteParticles_.clear();
}

template <typename Accessor>
void HashGrids::operator()(const size_t p_idx, Accessor& ac)
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (data::particle_flags::isSet(ac.getFlags(p_idx), data::particle_flags::INFINITE))
   {
      addInfiniteParticle(p_idx);
   } else
   {
      WALBERLA_ASSERT_GREATER(ac.getInteractionRadius(p_idx), 0_r, "Did you forget to set the interaction radius?");

      auto particleSize = 2_r * ac.getInteractionRadius(p_idx);

      shared_ptr<HashGrid> grid = nullptr;

      if( gridList_.empty() )
      {
         // If no hash grid yet exists in the hierarchy, an initial hash grid is created
         // based on the particle's size.

         grid = std::make_shared<HashGrid>( particleSize * std::sqrt( hierarchyFactor ) );
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

            if( particleSize < cellSpan )
            {
               cellSpan /= hierarchyFactor;
               if( particleSize < cellSpan ) {
                  while( particleSize < cellSpan ) cellSpan /= hierarchyFactor;
                  grid = std::make_shared<HashGrid>(cellSpan * hierarchyFactor );
                  gridList_.insert( gIt, grid );
               }

               grid->addParticle( p_idx, ac );
               return;
            }
         }

         while( particleSize >= cellSpan) cellSpan *= hierarchyFactor;
         grid = std::make_shared<HashGrid>( cellSpan );
      }

      grid->addParticle( p_idx, ac );
      gridList_.push_back( grid );

      return;

   }
}

void HashGrids::addInfiniteParticle(size_t p_idx)
{
   infiniteParticles_.push_back(p_idx);
}


template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void HashGrids::forEachParticlePairHalf(const bool openmp, const Selector& selector, Accessor& ac, Func&& func, Args&&... args) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");
   WALBERLA_UNUSED(openmp);

   // Pair generation by traversing through all hash grids (which are sorted in ascending order
   // with respect to the size of their cells).
   for( auto gridIt = gridList_.begin(); gridIt != gridList_.end(); ++gridIt ) {

      // Contact generation for all particles stored in the currently processed grid 'grid'.
      ParticleIdxVector particlesOnGrid;
      (*gridIt)->checkEachParticlePairHalfAndStore( particlesOnGrid, openmp, selector, ac, func, std::forward<Args>(args)... );

      if( !particlesOnGrid.empty() ) {

         // Test all particles stored in 'grid' against particles stored in grids with larger sized cells.
         auto nextGridIt = gridIt;
         for( ++nextGridIt; nextGridIt != gridList_.end(); ++nextGridIt ) {
            (*nextGridIt)->checkAgainstVectorEachParticlePairHalf( particlesOnGrid, openmp, selector, ac, func, std::forward<Args>(args)... );
         }

         for( auto& aIt : particlesOnGrid ) {
            // Test all particles stored in 'grid' against all particles stored in 'nonGridParticles_'.
            for( auto& infIt : infiniteParticles_ ) {
               if (selector(uint_c(aIt), uint_c(infIt), ac))
               {
                  func(uint_c(aIt), uint_c(infIt), std::forward<Args>(args)...);
               }
            }
         }
      }
   }

   // exclude infinite - infinite interactions
}


} //namespace data
} //namespace mesa_pd
} //namespace walberla