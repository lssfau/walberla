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
//! \file Bubble.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Describes a bubble as gas volume via volume and density.
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <iostream>

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
// forward declarations of friend classes
template< typename Stencil_T >
class BubbleModel;

class MergeInformation;

/***********************************************************************************************************************
 * Describes a bubble as gas volume via volume and density. A bubble can be located on multiple blocks, i.e., processes.
 **********************************************************************************************************************/
class Bubble
{
 public:
   explicit Bubble(real_t initVolume)
      : initVolume_(initVolume), currentVolume_(initVolume), volumeDiff_(real_c(0)), rho_(real_c(1.0))
   {}

   Bubble(real_t initVolume, real_t density)
      : initVolume_(initVolume), currentVolume_(density * initVolume), volumeDiff_(real_c(0)), rho_(density)
   {}

   // dummy constructor with meaningless default values
   Bubble() : initVolume_(real_c(-1.0)), currentVolume_(real_c(-1.0)), volumeDiff_(real_c(0)), rho_(real_c(-1.0)) {}

   real_t getInitVolume() const { return initVolume_; }
   real_t getCurrentVolume() const { return currentVolume_; }
   real_t getDensity() const { return rho_; }

   bool hasConstantDensity() const { return hasConstantDensity_; }

   void setConstantDensity(real_t density = real_c(1.0))
   {
      hasConstantDensity_ = true;
      rho_                = density;
   }

   // update the bubble volume change
   void updateVolumeDiff(real_t diff) { volumeDiff_ += diff; }

   void setDensity(real_t rho)
   {
      initVolume_ = rho * currentVolume_;
      updateDensity();
   }

 private:
   template< typename Stencil_T >
   friend class BubbleModel;

   friend class MergeInformation;

   // merge bubbles by adding volumes, and update density (see dissertation of S. Bogner, 2017, section 4.3)
   void merge(const Bubble& other)
   {
      initVolume_ += other.initVolume_;
      currentVolume_ += other.currentVolume_;
      updateDensity();
   }

   // update bubble volume and density using the change in the bubble's volume (see dissertation of S. Bogner, 2017,
   // section 4.3)
   void applyVolumeDiff(real_t diff)
   {
      WALBERLA_ASSERT(volumeDiff_ <= real_c(0.0));
      currentVolume_ += diff;
      updateDensity();
   }

   // return and reset the bubble's volume change
   real_t getAndResetVolumeDiff()
   {
      real_t ret  = volumeDiff_;
      volumeDiff_ = real_c(0);
      return ret;
   }

   // update bubble density (see dissertation of S. Bogner, 2017, section 4.3)
   void updateDensity()
   {
      if (hasConstantDensity_) return;
      rho_ = initVolume_ / currentVolume_;
   }

   real_t initVolume_;
   real_t currentVolume_;
   real_t volumeDiff_; // bubble's volume change (caused by interface fill level or movement)
   real_t rho_;

   bool hasConstantDensity_{ false };

   // function for packing a bubble into SendBuffer
   friend mpi::SendBuffer& operator<<(mpi::SendBuffer& buf, const Bubble& b);

   // function for unpacking a bubble from RecvBuffer
   friend mpi::RecvBuffer& operator>>(mpi::RecvBuffer& buf, Bubble& b);
}; // class Bubble

inline mpi::SendBuffer& operator<<(mpi::SendBuffer& buf, const Bubble& b)
{
   return buf << b.initVolume_ << b.currentVolume_;
}

inline mpi::RecvBuffer& operator>>(mpi::RecvBuffer& buf, Bubble& b)
{
   buf >> b.initVolume_ >> b.currentVolume_;
   b.updateDensity();
   return buf;
}

inline std::ostream& operator<<(std::ostream& os, const Bubble& b)
{
   os << "Bubble (" << b.getInitVolume() << "," << b.getCurrentVolume() << "," << b.getDensity() << ")";

   return os;
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

namespace walberla
{
namespace mpi
{
template<> // value type
struct BufferSizeTrait< free_surface::bubble_model::Bubble >
{
   static const bool constantSize = true;
   // 2 * real_t since the buffers above are filled with initVolume_ and currentVolume_
   static const uint_t size = 2 * sizeof(real_t) + mpi::BUFFER_DEBUG_OVERHEAD;
};
} // namespace mpi
} // namespace walberla
