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
//! \file LubricationCorrection.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"

#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/LinkedCells.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/InsertParticleIntoLinkedCells.h"
#include "mesa_pd/kernel/ParticleSelector.h"

namespace lubrication_correction_test
{
using namespace walberla;


// from Brenner - "The slow motion of a sphere through a viscous fluid towards a plane surface" (1961), Gl. 3.13
// from Cooley & O'Neill - "On the slow motion generated in a viscous fluid by the approach of a sphere to a plane wall or stationary sphere" (1969) -> Brenners work for free surface is identical to two spheres
//NOTE: Brenner uses a sphere with velocity U and thus calculates the force as F = 6 PI visc R U beta
// thus, the normalization to get beta is F / (..U..). However, we have two approaching spheres here with U and -U as velocities, thus relative velocity of 2U. Therefore, we have to divide this value by 2.
// see also Tab.2 of Cooley & O'Neill and discussion above.
real_t analyticalNonDimForceSphSph( real_t radius, real_t gap )
{
   uint_t numberOfAdditions = 50;

   real_t h = radius + real_t(0.5)*gap; // half gap because the formula was derived for sph-free surface where the gap is half of the gap here
   real_t alpha = std::acosh(h/radius);

   real_t temp = real_t(0);

   for( uint_t n = 1; n < numberOfAdditions; ++n )
   {
      real_t nr = real_t(n);
      real_t temp1 = real_t(4) * std::cosh( ( nr + real_t(0.5) ) * alpha ) * std::cosh( ( nr + real_t(0.5) ) * alpha ) +
                     ( real_t(2) * nr + real_t(1) ) * ( real_t(2) * nr + real_t(1) ) * std::sinh( alpha ) * std::sinh( alpha );
      real_t temp2 = real_t(2) * std::sinh( ( real_t(2) * nr + real_t(1) ) * alpha) -
                     (real_t(2) * nr + real_t(1) ) * std::sinh( real_t(2) * alpha );
      temp += nr * ( nr + real_t(1) ) / ( (real_t(2) * nr - real_t(1) ) * ( real_t(2) * nr + real_t(3) ) )  * ( temp1 / temp2 - real_t(1));
   }

   real_t beta = real_t(4)/real_t(3) * std::sinh(alpha) * temp;

   return beta * real_t(0.5);
}

// from Brenner - "The slow motion of a sphere through a viscous fluid towards a plane surface" (1961), Gl. 2.19
real_t analyticalNonDimForceSphWall( real_t radius, real_t gap )
{
   uint_t numberOfAdditions = 50;

   real_t h = radius + gap;
   real_t alpha = std::acosh(h/radius);

   real_t temp = real_t(0);

   for( uint_t n = 1; n < numberOfAdditions; ++n )
   {
      real_t nr = real_t(n);
      real_t temp1 = real_t(2) * std::sinh( ( real_t(2) * nr + real_t(1) ) * alpha ) + ( real_t(2) * nr + real_t(1) ) * std::sinh(real_t(2) * alpha );
      real_t temp2 = real_t(4) * std::sinh( ( nr + real_t(0.5) ) * alpha) * std::sinh( ( nr + real_t(0.5) ) * alpha) - (real_t(2) * nr + real_t(1)) * (real_t(2) * nr + real_t(1)) * std::sinh(alpha) * std::sinh(alpha);
      temp += nr * ( nr + real_t(1) ) / ( (real_t(2) * nr - real_t(1) ) * ( real_t(2) * nr + real_t(3) ) )  * ( temp1 / temp2 - real_t(1));
   }

   real_t lambda = real_t(4)/real_t(3) * std::sinh(alpha) * temp;

   return lambda;

}

/*!\brief Checks the implementation of the lubrication correction
 *
 * the resulting forces of various scenarios are checked for plausibility (sign, fi = -fj, etc) and compared to the analytical function
 * note that the lubrication formula of the correction is typically only a simple approximation of the analytical value (an infinite series) which is why differences will be present
 *
 * Currently, only the normal force is checked
 */

//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto shapeStorage = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor = mesa_pd::data::ParticleAccessorWithShape;
   ParticleAccessor accessor(ps, shapeStorage);

   // note: these are just arbitrary test values and are not necessarily physically meaningful!
   real_t cutOffDistance = real_t(0.8);
   real_t minimalGapSize = real_t(0.1);
   real_t dynamicViscosity = real_t(1.7);

   mesa_pd::kernel::DoubleCast doubleCast;
   mesa_pd::collision_detection::AnalyticContactDetection acd;
   acd.getContactThreshold() = cutOffDistance;


   /////////////////////
   // SPHERE - SPHERE //
   /////////////////////
   {
      real_t sphereRadius = real_t(2);
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      real_t sphereGapSize = real_t(0.3);
      real_t relativeVelocity = real_t(0.7);

      real_t normalizationFactor = real_t(6) * math::pi * dynamicViscosity * sphereRadius * relativeVelocity;

      real_t analyticalLubForce = (analyticalNonDimForceSphSph(sphereRadius, sphereGapSize) - analyticalNonDimForceSphSph(sphereRadius, cutOffDistance)) * normalizationFactor;

      Vector3<real_t> position1(real_t(0), real_t(0), real_t(0));
      Vector3<real_t> position2(real_t(2) * sphereRadius + sphereGapSize, real_t(0), real_t(0));
      Vector3<real_t> position3(real_t(0), real_t(2) * sphereRadius + sphereGapSize, real_t(0));
      Vector3<real_t> position4(real_t(0), real_t(0), real_t(2) * sphereRadius + sphereGapSize);

      mesa_pd::data::Particle&& p1 = *ps->create();
      p1.setPosition(position1);
      p1.setShapeID(sphereShape);
      p1.setLinearVelocity(Vector3<real_t>(relativeVelocity,-relativeVelocity,real_t(0.5)*relativeVelocity));
      auto idxSph1 = p1.getIdx();

      mesa_pd::data::Particle&& p2 = *ps->create();
      p2.setPosition(position2);
      p2.setShapeID(sphereShape);
      p2.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(2),real_t(0)));
      auto idxSph2 = p2.getIdx();

      mesa_pd::data::Particle&& p3 = *ps->create();
      p3.setPosition(position3);
      p3.setShapeID(sphereShape);
      p3.setLinearVelocity(Vector3<real_t>(real_t(2),real_t(0),real_t(2)));
      auto idxSph3 = p3.getIdx();

      mesa_pd::data::Particle&& p4 = *ps->create();
      p4.setPosition(position4);
      p4.setShapeID(sphereShape);
      p4.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),real_t(-0.5)*relativeVelocity));
      auto idxSph4 = p4.getIdx();

      lbm_mesapd_coupling::LubricationCorrectionKernel lubCorrFctr(dynamicViscosity, [minimalGapSize](real_t /*r*/){return minimalGapSize;}, cutOffDistance, real_t(0), real_t(0));

      ps->forEachParticlePairHalf(false, mesa_pd::kernel::ExcludeInfiniteInfinite(), accessor,
                                 [&](const size_t idx1, const size_t idx2)
                                 {
                                    if (doubleCast(idx1, idx2, accessor, acd, accessor ))
                                    {
                                       //note: in parallel simulations, first use the contact filter to assign the pairs to a single process
                                       doubleCast(idx1, idx2, accessor, lubCorrFctr, accessor, acd.getContactNormal(), acd.getPenetrationDepth());
                                    }
                                 }
      );

      auto force1 = accessor.getForce(idxSph1);
      auto force2 = accessor.getForce(idxSph2);
      auto force3 = accessor.getForce(idxSph3);
      auto force4 = accessor.getForce(idxSph4);

      WALBERLA_CHECK_LESS(force1[0], real_t(0), "Sph1, dir 0, sign check");
      WALBERLA_CHECK_LESS((std::abs(force1[0]) - std::abs(analyticalLubForce))/ std::abs(analyticalLubForce), real_t(0.2), "Sph1, dir 0, rel error check");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[0], -force2[0], "Sph1 - Sph2, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force2[1], real_t(0), "Sph2, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(force2[2], real_t(0), "Sph2, dir 2");

      WALBERLA_CHECK_GREATER(force1[1], real_t(0), "Sph1, dir 1, sign check");
      WALBERLA_CHECK_LESS((std::abs(force1[1]) - std::abs(analyticalLubForce))/ std::abs(analyticalLubForce), real_t(0.2), "Sph1, dir 1, rel error check");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[1], -force3[1], "Sph1 - Sph3, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(force3[0], real_t(0), "Sph3, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force3[2], real_t(0), "Sph3, dir 2");

      WALBERLA_CHECK_LESS(force1[2], real_t(0), "Sph1, dir 2, sign check");
      WALBERLA_CHECK_LESS((std::abs(force1[2]) - std::abs(analyticalLubForce))/ std::abs(analyticalLubForce), real_t(0.2), "Sph1, dir 2, rel error check");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[2], -force4[2], "Sph1 - Sph4, dir 2");
      WALBERLA_CHECK_FLOAT_EQUAL(force4[0], real_t(0), "Sph4, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force4[1], real_t(0), "Sph4, dir 1");

      // clean up
      ps->clear();

   }

   ///////////////////////////////////
   // SPHERE - SPHERE Special Cases //
   ///////////////////////////////////
   {
      real_t sphereRadius = real_t(2);
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      real_t relativeVelocity = real_t(0.7);

      real_t normalizationFactor = real_t(6) * math::pi * dynamicViscosity * sphereRadius * relativeVelocity;

      Vector3<real_t> position1(real_t(0), real_t(0), real_t(0));
      Vector3<real_t> position2(real_t(2) * sphereRadius - real_t(0.1), real_t(0), real_t(0)); // with overlap
      Vector3<real_t> position3(real_t(0), real_t(2) * sphereRadius + real_t(0.3) * minimalGapSize, real_t(0)); // below minimal gap size
      Vector3<real_t> position4(real_t(0), real_t(0), real_t(2) * sphereRadius + real_t(1.2) * cutOffDistance); // larger than cutoff distance

      mesa_pd::data::Particle&& p1 = *ps->create();
      p1.setPosition(position1);
      p1.setShapeID(sphereShape);
      p1.setLinearVelocity(Vector3<real_t>(relativeVelocity,-relativeVelocity,real_t(0.5)*relativeVelocity));
      auto idxSph1 = p1.getIdx();

      mesa_pd::data::Particle&& p2 = *ps->create();
      p2.setPosition(position2);
      p2.setShapeID(sphereShape);
      p2.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(2),real_t(0)));
      auto idxSph2 = p2.getIdx();

      mesa_pd::data::Particle&& p3 = *ps->create();
      p3.setPosition(position3);
      p3.setShapeID(sphereShape);
      p3.setLinearVelocity(Vector3<real_t>(real_t(2),real_t(0),real_t(2)));
      auto idxSph3 = p3.getIdx();

      mesa_pd::data::Particle&& p4 = *ps->create();
      p4.setPosition(position4);
      p4.setShapeID(sphereShape);
      p4.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),real_t(-0.5)*relativeVelocity));
      auto idxSph4 = p4.getIdx();

      // variant with linked cells
      // they can be used to break the O(N^2) complexity

      math::AABB domain(real_t(0), real_t(0), real_t(0), real_t(10), real_t(10), real_t(10));
      real_t spacing = real_t(2.1) * ( sphereRadius + cutOffDistance );
      mesa_pd::data::LinkedCells lc(domain.getExtended(spacing), spacing );
      mesa_pd::kernel::InsertParticleIntoLinkedCells ipilc;

      lc.clear();
      ps->forEachParticle( false, mesa_pd::kernel::SelectAll(), accessor, ipilc, accessor, lc);

      lbm_mesapd_coupling::LubricationCorrectionKernel lubCorrFctr(dynamicViscosity, [minimalGapSize](real_t /*r*/){return minimalGapSize;}, cutOffDistance, real_t(0), real_t(0));

      lc.forEachParticlePairHalf( false, mesa_pd::kernel::ExcludeInfiniteInfinite(), accessor,
                                  [&](const size_t idx1, const size_t idx2, auto& ac)
                                  {
                                     if (doubleCast(idx1, idx2, ac, acd, ac ))
                                     {
                                        //note: in parallel simulations, first use the contact filter to assign the pairs to a single process
                                        doubleCast(idx1, idx2, ac, lubCorrFctr, ac, acd.getContactNormal(), acd.getPenetrationDepth());
                                     }
                                  },
      accessor
      );

      auto force1 = accessor.getForce(idxSph1);
      auto force2 = accessor.getForce(idxSph2);
      auto force3 = accessor.getForce(idxSph3);
      auto force4 = accessor.getForce(idxSph4);

      WALBERLA_CHECK_FLOAT_EQUAL(force1[0], real_t(0), "SpecialCase: Sph1, dir 0, sign check");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[0], -force2[0], "SpecialCase:  Sph1 - Sph2, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force2[1], real_t(0), "SpecialCase:  Sph2, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(force2[2], real_t(0), "SpecialCase:  Sph2, dir 2");

      real_t analyticalLubForce = (analyticalNonDimForceSphSph(sphereRadius, minimalGapSize) - analyticalNonDimForceSphSph(sphereRadius, cutOffDistance)) * normalizationFactor;
      WALBERLA_CHECK_GREATER(force1[1], real_t(0), "SpecialCase: Sph1, dir 1, sign check");
      WALBERLA_CHECK_LESS(std::abs((std::abs(force1[1]) - std::abs(analyticalLubForce))/ analyticalLubForce), real_t(0.2), "SpecialCase: Sph1, dir 1, rel error check");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[1], -force3[1], "SpecialCase: Sph1 - Sph3, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(force3[0], real_t(0), "SpecialCase: Sph3, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force3[2], real_t(0), "SpecialCase: Sph3, dir 2");

      WALBERLA_CHECK_FLOAT_EQUAL(force1[2], real_t(0), "SpecialCase:  Sph1, dir 2, sign check");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[2], -force4[2], "SpecialCase:  Sph1 - Sph4, dir 2");
      WALBERLA_CHECK_FLOAT_EQUAL(force4[0], real_t(0), "SpecialCase:  Sph4, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force4[1], real_t(0), "SpecialCase:  Sph4, dir 1");

      // clean up
      ps->clear();

   }

   ///////////////////
   // SPHERE - WALL //
   ///////////////////

   {
      Vector3<real_t> wallPosition(real_t(0), real_t(0), real_t(0));
      Vector3<real_t> wallNormal(real_t(0), real_t(0), real_t(1));

      auto planeShape = shapeStorage->create<mesa_pd::data::HalfSpace>( wallNormal.getNormalized() );

      real_t sphereRadius = real_t(2);
      auto sphereShape = shapeStorage->create<mesa_pd::data::Sphere>( sphereRadius );

      real_t relativeVelocity = real_t(0.7);

      real_t gapSize = real_t(0.3);
      real_t normalizationFactor = real_t(6) * math::pi * dynamicViscosity * sphereRadius * relativeVelocity;

      Vector3<real_t> position1(real_t(0), real_t(0), sphereRadius + gapSize);
      Vector3<real_t> position2(real_t( 3) * sphereRadius, real_t(0), sphereRadius + gapSize);
      Vector3<real_t> position3(real_t( 6) * sphereRadius, real_t(0), sphereRadius - real_t(0.1));
      Vector3<real_t> position4(real_t( 9) * sphereRadius, real_t(0), sphereRadius + real_t(0.3) * minimalGapSize);
      Vector3<real_t> position5(real_t(12) * sphereRadius, real_t(0), sphereRadius + real_t(1.1) * cutOffDistance);

      mesa_pd::data::Particle&& p1 = *ps->create();
      p1.setPosition(position1);
      p1.setShapeID(sphereShape);
      p1.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),relativeVelocity));
      auto idxSph1 = p1.getIdx();

      // mix up order to test Sph - HSp and HSp - Sph variants
      mesa_pd::data::Particle&& pW = *ps->create(true);
      pW.setPosition(wallPosition);
      pW.setShapeID(planeShape);
      auto idxWall = pW.getIdx();

      mesa_pd::data::Particle&& p2 = *ps->create();
      p2.setPosition(position2);
      p2.setShapeID(sphereShape);
      p2.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),-relativeVelocity));
      auto idxSph2 = p2.getIdx();

      mesa_pd::data::Particle&& p3 = *ps->create();
      p3.setPosition(position3);
      p3.setShapeID(sphereShape);
      p3.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),relativeVelocity));
      auto idxSph3 = p3.getIdx();

      mesa_pd::data::Particle&& p4 = *ps->create();
      p4.setPosition(position4);
      p4.setShapeID(sphereShape);
      p4.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),relativeVelocity));
      auto idxSph4 = p4.getIdx();

      mesa_pd::data::Particle&& p5 = *ps->create();
      p5.setPosition(position5);
      p5.setShapeID(sphereShape);
      p5.setLinearVelocity(Vector3<real_t>(real_t(0),real_t(0),relativeVelocity));
      auto idxSph5 = p5.getIdx();

      lbm_mesapd_coupling::LubricationCorrectionKernel lubCorrFctr(dynamicViscosity, [minimalGapSize](real_t /*r*/){return minimalGapSize;}, cutOffDistance, real_t(0), real_t(0));

      ps->forEachParticlePairHalf(false, mesa_pd::kernel::ExcludeInfiniteInfinite(), accessor,
                                  [&](const size_t idx1, const size_t idx2)
                                  {
                                     if (doubleCast(idx1, idx2, accessor, acd, accessor ))
                                     {
                                        //note: in parallel simulations, first use the contact filter to assign the pairs to a single CPU
                                        doubleCast(idx1, idx2, accessor, lubCorrFctr, accessor, acd.getContactNormal(), acd.getPenetrationDepth());
                                     }
                                  }
      );

      auto force1 = accessor.getForce(idxSph1);
      auto force2 = accessor.getForce(idxSph2);
      auto force3 = accessor.getForce(idxSph3);
      auto force4 = accessor.getForce(idxSph4);
      auto force5 = accessor.getForce(idxSph5);
      auto forceW = accessor.getForce(idxWall);

      WALBERLA_CHECK_FLOAT_EQUAL(forceW[0], real_t(0), "Wall, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(forceW[1], real_t(0), "Wall, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(forceW[2], real_t(0), "Wall, dir 2");

      real_t analyticalLubForceRegular = (analyticalNonDimForceSphWall(sphereRadius, gapSize) - analyticalNonDimForceSphWall(sphereRadius, cutOffDistance)) * normalizationFactor;
      WALBERLA_CHECK_FLOAT_EQUAL(force1[0], real_t(0), "Wall-Sph1, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force1[1], real_t(0), "Wall-Sph1, dir 1");
      WALBERLA_CHECK_LESS(force1[2], real_t(0), "Wall-Sph1, dir 2, sign check");
      WALBERLA_CHECK_LESS((std::abs(force1[2]) - std::abs(analyticalLubForceRegular))/ std::abs(analyticalLubForceRegular), real_t(0.2), "Wall-Sph1, dir 2, rel error check");

      WALBERLA_CHECK_FLOAT_EQUAL(force2[0], real_t(0), "Wall-Sph2, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force2[1], real_t(0), "Wall-Sph2, dir 1");
      WALBERLA_CHECK_GREATER(force2[2], real_t(0), "Wall-Sph2, dir 2, sign check");
      WALBERLA_CHECK_LESS((std::abs(force2[2]) - std::abs(analyticalLubForceRegular))/ std::abs(analyticalLubForceRegular), real_t(0.2), "Wall-Sph2, dir 2, rel error check");

      WALBERLA_CHECK_FLOAT_EQUAL(force3[0], real_t(0), "Wall-Sph3, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force3[1], real_t(0), "Wall-Sph3, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(force3[2], real_t(0), "Wall-Sph3, dir 2");

      real_t analyticalLubForceMin = (analyticalNonDimForceSphWall(sphereRadius, minimalGapSize) - analyticalNonDimForceSphWall(sphereRadius, cutOffDistance)) * normalizationFactor;
      WALBERLA_CHECK_FLOAT_EQUAL(force4[0], real_t(0), "Wall-Sph4, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force4[1], real_t(0), "Wall-Sph4, dir 1");
      WALBERLA_CHECK_LESS(force4[2], real_t(0), "Wall-Sph4, dir 2, sign check");
      WALBERLA_CHECK_LESS((std::abs(force4[2]) - std::abs(analyticalLubForceMin))/ std::abs(analyticalLubForceMin), real_t(0.2), "Wall-Sph4, dir 2, rel error check");

      WALBERLA_CHECK_FLOAT_EQUAL(force5[0], real_t(0), "Wall-Sph5, dir 0");
      WALBERLA_CHECK_FLOAT_EQUAL(force5[1], real_t(0), "Wall-Sph5, dir 1");
      WALBERLA_CHECK_FLOAT_EQUAL(force5[2], real_t(0), "Wall-Sph5, dir 2");

      // clean up
      ps->clear();

   }

   return 0;

}

} //namespace lubrication_correction_test

int main( int argc, char **argv ){
   lubrication_correction_test::main(argc, argv);
}
