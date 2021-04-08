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
//! \file   CreateParticles.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "CreateParticles.h"

namespace walberla {
namespace mesa_pd {

data::ParticleStorage::iterator createPlane( data::ParticleStorage& ps,
                                             data::ShapeStorage& ss,
                                             const Vec3& pos,
                                             const Vec3& normal )
{
   auto p0              = ps.create(true);
   p0->getPositionRef() = pos;
   p0->getInteractionRadiusRef() = std::numeric_limits<real_t>::infinity();
   p0->getShapeIDRef()  = ss.create<data::HalfSpace>( normal );
   p0->getOwnerRef()    = walberla::mpi::MPIManager::instance()->rank();
   p0->getTypeRef()     = 0;
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::FIXED);
   data::particle_flags::set(p0->getFlagsRef(), data::particle_flags::NON_COMMUNICATING);
   return p0;
}

data::ParticleStorage::iterator createSphere( data::ParticleStorage& ps,
                                              const Vec3& pos,
                                              const real_t& radius,
                                              const uint64_t shapeID)
{
   auto p                       = ps.create();
   p->getPositionRef()          = pos;
   p->getInteractionRadiusRef() = radius;
   p->getShapeIDRef()           = shapeID;
   p->getOwnerRef()             = walberla::MPIManager::instance()->rank();
   p->getTypeRef()              = 0;
   return p;
}

} // namespace mesa_pd
} // namespace walberla
