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
//! \file   CreateParticles.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/ShapeStorage.h>

#include <core/mpi/MPIManager.h>

namespace walberla {
namespace mesa_pd {

data::ParticleStorage::iterator createPlane( data::ParticleStorage& ps,
                                             data::ShapeStorage& ss,
                                             const Vec3& pos,
                                             const Vec3& normal );

data::ParticleStorage::iterator createSphere( data::ParticleStorage& ps,
                                              const Vec3& pos,
                                              const real_t& radius,
                                              const uint64_t shapeID);

} // namespace mesa_pd
} // namespace walberla
