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
#include <mesa_pd/data/SparseLinkedCells.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Inserts a particle into the data::SparseLinkedCells data structure
 *
 * \attention Make sure to data::SparseLinkedCells::clear() the data structure before
 * reinserting new particles.
 *
 * This kernel requires the following particle accessor interface
 * \code
 * const walberla::mesa_pd::Vec3& getPosition(const size_t p_idx) const;
 *
 * const walberla::mesa_pd::data::particle_flags::FlagT& getFlags(const size_t p_idx) const;
 *
 * const size_t& getNextParticle(const size_t p_idx) const;
 * void setNextParticle(const size_t p_idx, const size_t& v);
 *
 * \endcode
 * \ingroup mesa_pd_kernel
 */
class InsertParticleIntoSparseLinkedCells
{
public:
   template <typename Accessor>
   void operator()(const size_t p_idx, Accessor& ac, data::SparseLinkedCells& lc) const;
};

template <typename Accessor>
inline void InsertParticleIntoSparseLinkedCells::operator()(const size_t p_idx, Accessor& ac, data::SparseLinkedCells& lc) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   const auto& minCorner = lc.domain_.minCorner();
   if (data::particle_flags::isSet(ac.getFlags(p_idx), data::particle_flags::INFINITE))
   {
      ac.setNextParticle(p_idx, lc.infiniteParticles_.exchange(int_c(p_idx)));
   } else
   {
      WALBERLA_ASSERT_GREATER(ac.getInteractionRadius(p_idx), 0_r, "Did you forget to set the interaction radius?");
      WALBERLA_ASSERT_LESS(2_r * ac.getInteractionRadius(p_idx), lc.cellDiameter_[0], "Interaction radius is to large for this cell size. Contacts might get lost.");
      int hash0 = static_cast<int>(std::floor((ac.getPosition(p_idx)[0] - minCorner[0]) * lc.invCellDiameter_[0]));
      WALBERLA_ASSERT_LESS(2_r * ac.getInteractionRadius(p_idx), lc.cellDiameter_[0], "Interaction radius is to large for this cell size. Contacts might get lost.");
      int hash1 = static_cast<int>(std::floor((ac.getPosition(p_idx)[1] - minCorner[1]) * lc.invCellDiameter_[1]));
      WALBERLA_ASSERT_LESS(2_r * ac.getInteractionRadius(p_idx), lc.cellDiameter_[0], "Interaction radius is to large for this cell size. Contacts might get lost.");
      int hash2 = static_cast<int>(std::floor((ac.getPosition(p_idx)[2] - minCorner[2]) * lc.invCellDiameter_[2]));
      if (hash0 < 0) hash0 = 0;
      if (hash0 >= lc.numCellsPerDim_[0]) hash0 = lc.numCellsPerDim_[0] - 1;
      if (hash1 < 0) hash1 = 0;
      if (hash1 >= lc.numCellsPerDim_[1]) hash1 = lc.numCellsPerDim_[1] - 1;
      if (hash2 < 0) hash2 = 0;
      if (hash2 >= lc.numCellsPerDim_[2]) hash2 = lc.numCellsPerDim_[2] - 1;
      uint64_t cell_idx = getCellIdx(lc, hash0, hash1, hash2);
      ac.setNextParticle(p_idx, lc.cells_[cell_idx].exchange(int_c(p_idx)));
      if (ac.getNextParticle(p_idx) == -1)
      {
         lc.nonEmptyCells_.emplace_back(cell_idx);
      }
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla