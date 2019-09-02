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
//! \file InsertParticleIntoSparseLinkedCells.h
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
   {%- for prop in interface %}
   {%- if 'g' in prop.access %}
 * const {{prop.type}}& get{{prop.name | capFirst}}(const size_t p_idx) const;
   {%- endif %}
   {%- if 's' in prop.access %}
 * void set{{prop.name | capFirst}}(const size_t p_idx, const {{prop.type}}& v);
   {%- endif %}
   {%- if 'r' in prop.access %}
 * {{prop.type}}& get{{prop.name | capFirst}}Ref(const size_t p_idx);
   {%- endif %}
 *
   {%- endfor %}
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
      {%- for dim in range(3) %}
      int hash{{dim}} = static_cast<int>(std::floor((ac.getPosition(p_idx)[{{dim}}] - minCorner[{{dim}}]) * lc.invCellDiameter_[{{dim}}]));
      {%- endfor %}
      {%- for dim in range(3) %}
      if (hash{{dim}} < 0) hash{{dim}} = 0;
      if (hash{{dim}} >= lc.numCellsPerDim_[{{dim}}]) hash{{dim}} = lc.numCellsPerDim_[{{dim}}] - 1;
      {%- endfor %}
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
