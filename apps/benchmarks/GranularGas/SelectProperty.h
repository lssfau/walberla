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
//! \file   SelectProperty.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>

#include <core/mpi/MPIManager.h>

namespace walberla {
namespace mesa_pd {

class SelectRank
{
public:
   using return_type = int;
   int operator()(const data::Particle& /*p*/) const { return rank_; }
   int operator()(const data::Particle&& /*p*/) const { return rank_; }
private:
   int rank_ = walberla::mpi::MPIManager::instance()->rank();
};

class SelectIdx
{
public:
   using return_type = int;
   auto operator()(const data::Particle& p) const { return p.getIdx(); }
   auto operator()(const data::Particle&& p) const { return p.getIdx(); }
private:
};

} // namespace mesa_pd
} // namespace walberla
