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
//! \file   SpherePile.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/common/ParticleFunctions.h"

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"

#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/SemiImplicitEuler.h"
#include "mesa_pd/kernel/SpringDashpot.h"
#include "mesa_pd/kernel/SpringDashpotSpring.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "vtk/VTKOutput.h"

#include <iostream>

namespace walberla {
namespace mesa_pd {

auto createSphere(data::ParticleStorage &ps, const Vec3 &pos)
{
   auto p = ps.create();
   p->setPosition(pos);
   p->setType(0);
   return p;
}

/*
 * Simulates oblique sphere-wall collision and checks rebound angle, i.e. the tangential part of the collision model.
 *
 */
int main(int argc, char **argv)
{
   walberla::mpi::Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   auto dt = real_t(0.001);
   auto cor = real_t(0.2);
   auto ct = dt * real_t(20);
   auto radius = real_t(1);
   auto density = real_t(2700);
   auto simSteps = 100000;

   //init data structures
   auto ps = walberla::make_shared<data::ParticleStorage>(2);
   auto ss = walberla::make_shared<data::ShapeStorage>();
   auto ac = walberla::make_shared<data::ParticleAccessorWithShape>(ps, ss);

   auto sphereShape = ss->create<data::Sphere>(radius);
   ss->shapes[sphereShape]->updateMassAndInertia(density);

   // create sphere
   createSphere(*ps, Vec3(0, 0, 0));
   createSphere(*ps, Vec3(2, 0, 0));
   auto sp = createSphere(*ps, Vec3(1, 0, std::sqrt(real_t(4) - real_t(1))));

   // create plane
   data::Particle &&p0 = *ps->create(true);
   p0.setPosition(Vec3(0, 0, -1));
   p0.setShapeID(ss->create<data::HalfSpace>(Vector3<real_t>(0, 0, 1)));
   p0.setType(0);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::INFINITE);
   data::particle_flags::set(p0.getFlagsRef(), data::particle_flags::FIXED);

   WALBERLA_LOG_INFO_ON_ROOT("*** VTK ***");
   auto vtkOutput       = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps) ;
   auto vtkWriter       = walberla::vtk::createVTKOutput_PointData(vtkOutput,
         "Bodies",
         1,
         "vtk",
         "simulation_step",
         false,
         false);

   // explicit euler
   kernel::SemiImplicitEuler implEuler(dt);
   collision_detection::AnalyticContactDetection acd;
   kernel::DoubleCast double_cast;
   kernel::SpringDashpot sd(1);
   kernel::SpringDashpotSpring sds(1);
   mpi::ReduceContactHistory rch;

   WALBERLA_LOG_DEVEL_VAR(sp->getPosition());

   for (auto i = 0; i < simSteps; ++i)
   {
      ps->forEachParticle(false,
                          kernel::SelectLocal(),
                          *ac,
                          [&](const size_t idx1)
                          {
                             if (ac->getShape(idx1)->getShapeType() == data::Sphere::SHAPE_TYPE)
                             {
                                ac->setForce(idx1, Vec3(0,0,-1) * real_t(9.81) / ac->getInvMass(idx1));
                             }
                          });
      ps->forEachParticlePairHalf(false,
                                  kernel::SelectAll(),
                                  *ac,
                                  [&](const size_t idx1, const size_t idx2)
                                  {
                                     if ((ac->getShape(idx1)->getShapeType() != data::Sphere::SHAPE_TYPE) &&
                                         (ac->getShape(idx2)->getShapeType() != data::Sphere::SHAPE_TYPE))
                                     {
                                        //skip plane - plane collision
                                     } else
                                     {
                                        if (double_cast(idx1, idx2, *ac, acd, *ac))
                                        {
                                           //WALBERLA_LOG_DEVEL_VAR(acd.getPenetrationDepth());
                                           auto meff = real_t(1) / (ac->getInvMass(idx1) + ac->getInvMass(idx2));
                                           sds.setParametersFromCOR(0, 0, cor, ct, meff);
                                           sds.setCoefficientOfFriction(0,0, real_t(0.4));
                                           sds.setStiffnessT(0,0, sds.getStiffnessN(0,0));
                                           sds(acd.getIdx1(), acd.getIdx2(), *ac, acd.getContactPoint(),
                                                      acd.getContactNormal(), acd.getPenetrationDepth(), dt);
                                        }
                                     }
                                  });

      rch(*ps);
      ps->forEachParticle(false,
                    kernel::SelectLocal(),
                    *ac,
                    implEuler,
                    *ac);
   }
   WALBERLA_LOG_DEVEL_VAR(sp->getPosition());
   WALBERLA_CHECK_GREATER(sp->getPosition()[2], real_t(1));

   for (auto i = 0; i < simSteps; ++i)
   {
      ps->forEachParticle(false,
                          kernel::SelectLocal(),
                          *ac,
                          [&](const size_t idx1)
                          {
                             if (ac->getShape(idx1)->getShapeType() == data::Sphere::SHAPE_TYPE)
                             {
                                ac->setForce(idx1, Vec3(0,0,-1) * real_t(9.81) / ac->getInvMass(idx1));
                             }
                          });
      ps->forEachParticlePairHalf(false,
                                  kernel::SelectAll(),
                                  *ac,
                                  [&](const size_t idx1, const size_t idx2)
                                  {
                                     if ((ac->getShape(idx1)->getShapeType() != data::Sphere::SHAPE_TYPE) &&
                                         (ac->getShape(idx2)->getShapeType() != data::Sphere::SHAPE_TYPE))
                                     {
                                        //skip plane - plane collision
                                     } else
                                     {
                                        if (double_cast(idx1, idx2, *ac, acd, *ac))
                                        {
                                           //WALBERLA_LOG_DEVEL_VAR(acd.getPenetrationDepth());
                                           auto meff = real_t(1) / (ac->getInvMass(idx1) + ac->getInvMass(idx2));
                                           sd.setParametersFromCOR(0, 0, cor, ct, meff);
                                           sd.setFriction(0,0, real_t(0.4));
                                           sd.setDampingT(0,0, sd.getDampingN(0,0));
                                           sd(acd.getIdx1(), acd.getIdx2(), *ac, acd.getContactPoint(),
                                              acd.getContactNormal(), acd.getPenetrationDepth());
                                        }
                                     }
                                  });

      rch(*ps);
      ps->forEachParticle(false,
                    kernel::SelectLocal(),
                    *ac,
                    implEuler,
                    *ac);
   }
   WALBERLA_LOG_DEVEL_VAR(sp->getPosition());
   WALBERLA_CHECK_LESS(sp->getPosition()[2], real_t(1));

   return EXIT_SUCCESS;
}

} //namespace mesa_pd
} //namespace walberla

int main(int argc, char **argv)
{
   return walberla::mesa_pd::main(argc, argv);
}
