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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
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
#include <core/math/AABB.h>
#include <stencil/D3Q27.h>

#include <atomic>
#include <cmath>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace data {

struct LinkedCells
{
   LinkedCells(const math::AABB& domain, const real_t cellDiameter) : LinkedCells(domain, Vec3(cellDiameter,cellDiameter,cellDiameter)) {}
   LinkedCells(const math::AABB& domain, const Vec3& cellDiameter);

   void clear();

   /**
    * Calls the provided functor \p func for all particle pairs.
    *
    * Additional arguments can be provided. No pairs with twice the same particle.
    * Call syntax for the provided functor
    * \code
    * func( *this, i, j, std::forward<Args>(args)... );
    * \endcode
    * \param openmp enables/disables OpenMP parallelization of the kernel call
    */
   template <typename Selector, typename Accessor, typename Func, typename... Args>
   void forEachParticlePair(const bool openmp,
                            const Selector& selector,
                            Accessor& acForLC,
                            Func&& func,
                            Args&&... args) const;
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

   math::AABB   domain_ {};
   Vector3<int> numCellsPerDim_ {};
   Vec3         cellDiameter_ {};
   Vec3         invCellDiameter_ {};
   std::atomic<int> infiniteParticles_ {};
   std::vector< std::atomic<int> > cells_ {};
};

inline
math::AABB getCellAABB(const LinkedCells& ll,
                       const int64_t hash0,
                       const int64_t hash1,
                       const int64_t hash2)
{
   WALBERLA_ASSERT_GREATER_EQUAL(hash0, 0);
   WALBERLA_ASSERT_LESS(hash0, ll.numCellsPerDim_[0]);
   WALBERLA_ASSERT_GREATER_EQUAL(hash1, 0);
   WALBERLA_ASSERT_LESS(hash1, ll.numCellsPerDim_[1]);
   WALBERLA_ASSERT_GREATER_EQUAL(hash2, 0);
   WALBERLA_ASSERT_LESS(hash2, ll.numCellsPerDim_[2]);
   const auto& minCorner = ll.domain_.minCorner();
   const real_t xMin = ll.cellDiameter_[0] * real_c(hash0) + minCorner[0];
   const real_t yMin = ll.cellDiameter_[1] * real_c(hash1) + minCorner[1];
   const real_t zMin = ll.cellDiameter_[2] * real_c(hash2) + minCorner[2];
   const real_t xMax = ll.cellDiameter_[0] * real_c(hash0 + 1) + minCorner[0];
   const real_t yMax = ll.cellDiameter_[1] * real_c(hash1 + 1) + minCorner[1];
   const real_t zMax = ll.cellDiameter_[2] * real_c(hash2 + 1) + minCorner[2];
   return math::AABB(xMin, yMin, zMin, xMax, yMax, zMax);
}

inline
uint_t getCellIdx(const LinkedCells& ll,
                  const int64_t hash0,
                  const int64_t hash1,
                  const int64_t hash2)
{
   WALBERLA_ASSERT_GREATER_EQUAL(hash0, 0);
   WALBERLA_ASSERT_LESS(hash0, ll.numCellsPerDim_[0]);
   WALBERLA_ASSERT_GREATER_EQUAL(hash1, 0);
   WALBERLA_ASSERT_LESS(hash1, ll.numCellsPerDim_[1]);
   WALBERLA_ASSERT_GREATER_EQUAL(hash2, 0);
   WALBERLA_ASSERT_LESS(hash2, ll.numCellsPerDim_[2]);
   return uint_c(hash2 * ll.numCellsPerDim_[1] * ll.numCellsPerDim_[0] + hash1 * ll.numCellsPerDim_[0] + hash0);
}

inline
void getCellCoordinates(const LinkedCells& ll,
                        const uint64_t idx,
                        int64_t& hash0,
                        int64_t& hash1,
                        int64_t& hash2)
{
   hash2 = int64_c(idx) / (ll.numCellsPerDim_[1] * ll.numCellsPerDim_[0]);
   hash1 = (int64_c(idx) - (hash2 * ll.numCellsPerDim_[1] * ll.numCellsPerDim_[0])) / (ll.numCellsPerDim_[0]);
   hash0 = int64_c(idx) - hash2 * ll.numCellsPerDim_[1] * ll.numCellsPerDim_[0] - hash1 * ll.numCellsPerDim_[0];

   WALBERLA_ASSERT_GREATER_EQUAL(hash0, 0);
   WALBERLA_ASSERT_LESS(hash0, ll.numCellsPerDim_[0]);
   WALBERLA_ASSERT_GREATER_EQUAL(hash1, 0);
   WALBERLA_ASSERT_LESS(hash1, ll.numCellsPerDim_[1]);
   WALBERLA_ASSERT_GREATER_EQUAL(hash2, 0);
   WALBERLA_ASSERT_LESS(hash2, ll.numCellsPerDim_[2]);
}

inline
LinkedCells::LinkedCells(const math::AABB& domain, const Vec3& cellDiameter)
   : domain_(domain)
   , numCellsPerDim_( static_cast<int>(std::ceil( domain.sizes()[0] / cellDiameter[0])),
     static_cast<int>(std::ceil( domain.sizes()[1] / cellDiameter[1])),
     static_cast<int>(std::ceil( domain.sizes()[2] / cellDiameter[2])) )
   , cellDiameter_( cellDiameter)
   , invCellDiameter_( real_t(1) / cellDiameter_[0], real_t(1) / cellDiameter_[1], real_t(1) / cellDiameter_[2] )
   , cells_(uint_c(numCellsPerDim_[0]*numCellsPerDim_[1]*numCellsPerDim_[2]))
{
   //precondition
   WALBERLA_CHECK_GREATER_EQUAL(cellDiameter[0], real_t(0));
   WALBERLA_CHECK_GREATER_EQUAL(cellDiameter[1], real_t(0));
   WALBERLA_CHECK_GREATER_EQUAL(cellDiameter[2], real_t(0));

   //postcondition
   WALBERLA_CHECK_GREATER_EQUAL(real_c(numCellsPerDim_[0]) * cellDiameter_[0], domain.size(0));
   WALBERLA_CHECK_GREATER_EQUAL(real_c(numCellsPerDim_[1]) * cellDiameter_[1], domain.size(1));
   WALBERLA_CHECK_GREATER_EQUAL(real_c(numCellsPerDim_[2]) * cellDiameter_[2], domain.size(2));

   std::fill(cells_.begin(), cells_.end(), -1);
}

void LinkedCells::clear()
{
   const uint64_t cellsSize = cells_.size();
   //clear existing linked cells
   for (int64_t i = 0; i < int64_c(cellsSize); ++i)
      cells_[uint64_c(i)] = -1;
   infiniteParticles_ = -1;
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void LinkedCells::forEachParticlePair(const bool openmp, const Selector& selector, Accessor& acForLC, Func&& func, Args&&... args) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");
   WALBERLA_UNUSED(openmp);
   for (int z = 0; z < numCellsPerDim_[2]; ++z)
   {
      for (int y = 0; y < numCellsPerDim_[1]; ++y)
      {
         for (int x = 0; x < numCellsPerDim_[0]; ++x)
         {
            const uint_t cell_idx = getCellIdx(*this, x, y, z); ///< current cell index
            int p_idx = cells_[cell_idx]; ///< current particle index
            int np_idx = -1; ///< particle to be checked against

            while (p_idx != -1)
            {
               WALBERLA_ASSERT_GREATER_EQUAL(p_idx, 0);
               WALBERLA_ASSERT_LESS(p_idx, acForLC.size());

               // check particles in own cell
               np_idx = acForLC.getNextParticle(uint_c(p_idx)); ///< neighbor particle index
               while (np_idx != -1)
               {
                  WALBERLA_ASSERT_GREATER_EQUAL(np_idx, 0);
                  WALBERLA_ASSERT_LESS(np_idx, acForLC.size());

                  if (selector(uint_c(p_idx), uint_c(np_idx), acForLC))
                  {
                     func(uint_c(p_idx), uint_c(np_idx), std::forward<Args>(args)...);
                     func(uint_c(np_idx), uint_c(p_idx), std::forward<Args>(args)...);
                  }

                  np_idx = acForLC.getNextParticle(uint_c(np_idx));
               }

               // check particles in neighboring cells (only positive ones)
               for (auto dir : stencil::D3Q27::dir_pos)
               {
                  const int nx = x + stencil::cx[dir];
                  const int ny = y + stencil::cy[dir];
                  const int nz = z + stencil::cz[dir];
                  if (nx < 0) continue;
                  if (ny < 0) continue;
                  if (nz < 0) continue;
                  if (nx >= numCellsPerDim_[0]) continue;
                  if (ny >= numCellsPerDim_[1]) continue;
                  if (nz >= numCellsPerDim_[2]) continue;

                  const uint_t ncell_idx = getCellIdx(*this, nx, ny, nz); ///< neighbor cell index

                  WALBERLA_ASSERT_GREATER_EQUAL(p_idx, 0);
                  WALBERLA_ASSERT_LESS(p_idx, acForLC.size());
                  np_idx = cells_[ncell_idx]; ///< neighbor particle index
                  while (np_idx != -1)
                  {
                     WALBERLA_ASSERT_GREATER_EQUAL(np_idx, 0);
                     WALBERLA_ASSERT_LESS(np_idx, acForLC.size());

                     if (selector(uint_c(p_idx), uint_c(np_idx), acForLC))
                     {
                        func(uint_c(p_idx), uint_c(np_idx), std::forward<Args>(args)...);
                        func(uint_c(np_idx), uint_c(p_idx), std::forward<Args>(args)...);
                     }

                     np_idx = acForLC.getNextParticle(uint_c(np_idx));
                  }
               }

               // check particles in infiniteParticles list
               np_idx = infiniteParticles_; ///< neighbor particle index
               while (np_idx != -1)
               {
                  WALBERLA_ASSERT_GREATER_EQUAL(np_idx, 0);
                  WALBERLA_ASSERT_LESS(np_idx, acForLC.size());

                  if (selector(uint_c(p_idx), uint_c(np_idx), acForLC))
                  {
                     func(uint_c(p_idx), uint_c(np_idx), std::forward<Args>(args)...);
                     func(uint_c(np_idx), uint_c(p_idx), std::forward<Args>(args)...);
                  }

                  np_idx = acForLC.getNextParticle(uint_c(np_idx));
               }

               // go to next particle
               p_idx = acForLC.getNextParticle(uint_c(p_idx));
            }
         }
      }
   }
}
template <typename Selector, typename Accessor, typename Func, typename... Args>
inline void LinkedCells::forEachParticlePairHalf(const bool openmp, const Selector& selector, Accessor& acForLC, Func&& func, Args&&... args) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");
   WALBERLA_UNUSED(openmp);
   for (int z = 0; z < numCellsPerDim_[2]; ++z)
   {
      for (int y = 0; y < numCellsPerDim_[1]; ++y)
      {
         for (int x = 0; x < numCellsPerDim_[0]; ++x)
         {
            const uint_t cell_idx = getCellIdx(*this, x, y, z); ///< current cell index
            int p_idx = cells_[cell_idx]; ///< current particle index
            int np_idx = -1; ///< particle to be checked against

            while (p_idx != -1)
            {
               WALBERLA_ASSERT_GREATER_EQUAL(p_idx, 0);
               WALBERLA_ASSERT_LESS(p_idx, acForLC.size());

               // check particles in own cell
               np_idx = acForLC.getNextParticle(uint_c(p_idx)); ///< neighbor particle index
               while (np_idx != -1)
               {
                  WALBERLA_ASSERT_GREATER_EQUAL(np_idx, 0);
                  WALBERLA_ASSERT_LESS(np_idx, acForLC.size());

                  if (selector(uint_c(p_idx), uint_c(np_idx), acForLC))
                  {
                     func(uint_c(p_idx), uint_c(np_idx), std::forward<Args>(args)...);
                  }

                  np_idx = acForLC.getNextParticle(uint_c(np_idx));
               }

               // check particles in neighboring cells (only positive ones)
               for (auto dir : stencil::D3Q27::dir_pos)
               {
                  const int nx = x + stencil::cx[dir];
                  const int ny = y + stencil::cy[dir];
                  const int nz = z + stencil::cz[dir];
                  if (nx < 0) continue;
                  if (ny < 0) continue;
                  if (nz < 0) continue;
                  if (nx >= numCellsPerDim_[0]) continue;
                  if (ny >= numCellsPerDim_[1]) continue;
                  if (nz >= numCellsPerDim_[2]) continue;

                  const uint_t ncell_idx = getCellIdx(*this, nx, ny, nz); ///< neighbor cell index

                  WALBERLA_ASSERT_GREATER_EQUAL(p_idx, 0);
                  WALBERLA_ASSERT_LESS(p_idx, acForLC.size());
                  np_idx = cells_[ncell_idx]; ///< neighbor particle index
                  while (np_idx != -1)
                  {
                     WALBERLA_ASSERT_GREATER_EQUAL(np_idx, 0);
                     WALBERLA_ASSERT_LESS(np_idx, acForLC.size());

                     if (selector(uint_c(p_idx), uint_c(np_idx), acForLC))
                     {
                        func(uint_c(p_idx), uint_c(np_idx), std::forward<Args>(args)...);
                     }

                     np_idx = acForLC.getNextParticle(uint_c(np_idx));
                  }
               }

               // check particles in infiniteParticles list
               np_idx = infiniteParticles_; ///< neighbor particle index
               while (np_idx != -1)
               {
                  WALBERLA_ASSERT_GREATER_EQUAL(np_idx, 0);
                  WALBERLA_ASSERT_LESS(np_idx, acForLC.size());

                  if (selector(uint_c(p_idx), uint_c(np_idx), acForLC))
                  {
                     func(uint_c(p_idx), uint_c(np_idx), std::forward<Args>(args)...);
                  }

                  np_idx = acForLC.getNextParticle(uint_c(np_idx));
               }

               // go to next particle
               p_idx = acForLC.getNextParticle(uint_c(p_idx));
            }
         }
      }
   }
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla