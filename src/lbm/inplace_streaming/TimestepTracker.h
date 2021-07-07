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
//! \\file TimestepTracker.h
//! \\ingroup lbm
//! \\author Frederik Hennig <frederik.hennig@fau.de>
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla
{
namespace lbm
{
class TimestepTracker
{
 private:
   uint8_t counter_{ 0 };

 public:
   TimestepTracker() = default;
   TimestepTracker(uint8_t zeroth_timestep) : counter_(zeroth_timestep & 1) {}

   void advance() { counter_ = (counter_ + 1) & 1; }

   std::function< void() > getAdvancementFunction()
   {
      return [this]() { this->advance(); };
   }

   uint8_t getCounter() const { return counter_; }
   uint8_t getCounterPlusOne() const { return (counter_ + 1) & 1; }

}; // class TimestepTracker

} // namespace lbm
} // namespace walberla
