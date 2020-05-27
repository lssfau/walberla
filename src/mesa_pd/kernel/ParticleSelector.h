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
//! \file ParticleSelector.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/ParticleStorage.h>

namespace walberla {
namespace mesa_pd {
namespace kernel {

class SelectAll
{
public:
   template <typename Accessor>
   bool operator()(const size_t /*idx*/, Accessor& /*ac*/) const { return true; }

   template <typename Accessor>
   bool operator()(const size_t /*idx*/, const size_t /*jdx*/, Accessor& /*ac*/) const { return true; }
};

class SelectMaster
{
public:
   template <typename Accessor>
   bool operator()(const size_t idx, Accessor& ac) const
   {
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(ac.getFlags(idx), GHOST)) return false;
      if (isSet(ac.getFlags(idx), GLOBAL)) return false;
      return true;
   }
};

/// prefer SelectMaster over SelectLocal
using SelectLocal = SelectMaster;

class SelectGhost
{
public:
   template <typename Accessor>
   bool operator()(const size_t idx, Accessor& ac) const
   {
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(ac.getFlags(idx), GHOST)) return true;
      return false;
   }
};

class ExcludeInfiniteInfinite
{
public:
   template <typename Accessor>
   bool operator()(const size_t idx, const size_t jdx, Accessor& ac) const
   {
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(ac.getFlags(idx), INFINITE) && isSet(ac.getFlags(jdx), INFINITE)) return false;
      return true;
   }
};

} //namespace data
} //namespace mesa_pd
} //namespace walberla
