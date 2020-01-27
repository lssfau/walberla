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
//! \file LubricationCorrectionKernel.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "mesa_pd/common/ParticleFunctions.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/IAccessor.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/data/shape/HalfSpace.h"

namespace walberla {
namespace lbm_mesapd_coupling {


/*! \brief Computes lubrication corrections between two spheres 1 and 2 acting on sphere 1.
 *
 * The formulas for the lubrication corrections for two spheres with different radii are taken from
 * Simeonov, Calantoni - Modeling mechanical contact and lubrication in Direct Numerical Simulations of colliding particles, IJMF, 2012
 *
 * interactionNormal12 is sphere2Center - sphere1Center, normalized
 *
 * radius1, radius2: radius of sphere 1 and 2
 * u1, u2: linear/translational velocity of sphere 1 and 2
 * omega1, omega2: angular velocity of sphere 1 and 2
 *
 * Returns force and torque via lubricationForce1 and lubricationTorque1
 *
 */
inline
void computeLubricationCorrectionSphereSphere( const Vector3<real_t> & interactionNormal12, const real_t gapSize, const real_t dynamicFluidViscosity,
                                               const real_t radius1, const real_t radius2,
                                               const Vector3<real_t> & u1, const Vector3<real_t> & u2,
                                               const Vector3<real_t> & omega1, const Vector3<real_t> & omega2,
                                               const real_t cutOffDistanceNormal, const real_t cutOffDistanceTangentialTranslational,
                                               const real_t cutOffDistanceTangentialRotational,
                                               Vector3<real_t> & lubricationForce1, Vector3<real_t> & lubricationTorque1 )
{
   WALBERLA_ASSERT_FLOAT_EQUAL(interactionNormal12.length(), real_t(1), "InteractionNormal has to be normalized!");

   real_t kappa = radius2 / radius1;
   real_t epsilon = gapSize / radius1;

   Vector3<real_t> e12 = -interactionNormal12;
   Vector3<real_t> u12 =  u2 - u1;
   Vector3<real_t> omega12 = omega1 + omega2;

   lubricationForce1 = Vector3<real_t>(real_t(0));
   lubricationTorque1 = Vector3<real_t>(real_t(0));

   if( gapSize < cutOffDistanceNormal )
   {
      // add lubrication force due to normal translation
      real_t epsres = cutOffDistanceNormal / radius1;
      Vector3<real_t> Fn1 = real_t(6) * walberla::math::pi * radius1 * dynamicFluidViscosity *
                            ( kappa * kappa / ( ( real_t(1) + kappa ) * ( real_t(1) + kappa ) ) ) * ( real_t(1) / epsilon - real_t(1) / epsres) *
                            math::dot(u12, e12) * e12;

      lubricationForce1 += Fn1;
   }

   if( gapSize < cutOffDistanceTangentialTranslational )
   {
      // add lubrication force and torque due to tangential translation
      real_t epsres = cutOffDistanceTangentialTranslational / radius1;
      real_t logEpsDivEpsres = std::log(epsilon/epsres);
      Vector3<real_t> Ftt1 = real_t(6) * walberla::math::pi * radius1 * dynamicFluidViscosity *
                             ( - real_t(4) * kappa * (real_t(2) + kappa + real_t(2) * kappa * kappa ) / ( real_t(15) * (real_t(1) + kappa) * (real_t(1) + kappa) * (real_t(1) + kappa) ) ) * logEpsDivEpsres *
                             (u12 - math::dot(u12,e12) * e12);

      lubricationForce1 += Ftt1;

      Vector3<real_t> Ttt1 = real_t(8) * walberla::math::pi * radius1 * radius1 * dynamicFluidViscosity *
                             ( - kappa * (real_t(4) + kappa) / (real_t(10) * (real_t(1)+kappa) * (real_t(1)+kappa) ) ) * logEpsDivEpsres *
                             math::cross(e12, u12);

      lubricationTorque1 += Ttt1;
   }

   if( gapSize < cutOffDistanceTangentialRotational )
   {
      // add lubrication force and torque due to tangential rotation
      real_t epsres = cutOffDistanceTangentialRotational / radius1;
      real_t logEpsDivEpsres = std::log(epsilon/epsres);
      Vector3<real_t> Ftr1 = real_t(6) * walberla::math::pi * radius1 * radius1 * dynamicFluidViscosity *
                             (real_t(2) * kappa * kappa / (real_t(15) * (real_t(1) + kappa) * (real_t(1) + kappa) ) ) * logEpsDivEpsres *
                             math::cross(omega12 + real_t(4) / kappa * omega1 + real_t(4) * kappa * omega2, e12);

      lubricationForce1 += Ftr1;

      Vector3<real_t> tempOmega = omega1 + kappa * omega2 / real_t(4);

      Vector3<real_t> Ttr1 = real_t(8) * walberla::math::pi * radius1 * radius1 * radius1 * dynamicFluidViscosity *
                             (real_t(2) * kappa / ( real_t(5) * (real_t(1) + kappa))) * logEpsDivEpsres *
                             (tempOmega - math::dot(tempOmega, e12) * e12);

      lubricationTorque1 += Ttr1;

   }

   // note: lubrication due to normal rotation is dropped here because of too low influence, see also Simeonov

}


/*! \brief Computes lubrication corrections between a sphere and a half space acting on the sphere.
 *
 * The formulas for the lubrication corrections are taken from
 * Simeonov, Calantoni - Modeling mechanical contact and lubrication in Direct Numerical Simulations of colliding particles, IJMF, 2012
 * with kappa -> infinity
 * *
 * radius1: radius of sphere
 * u1, u2: linear/translational velocity of sphere and half space
 * omega1: angular velocity of sphere 1, half space is assumed to have zero angular velocity
 *
 * Returns force and torque via lubricationForce1 and lubricationTorque1
 */
inline
void computeLubricationCorrectionSphereHalfSpace( const Vector3<real_t> & interactionNormal12, const real_t gapSize, const real_t dynamicFluidViscosity,
                                                  const real_t radius1,
                                                  const Vector3<real_t> & u1, const Vector3<real_t> & u2,
                                                  const Vector3<real_t> & omega1,
                                                  const real_t cutOffDistanceNormal, const real_t cutOffDistanceTangentialTranslational,
                                                  const real_t cutOffDistanceTangentialRotational,
                                                  Vector3<real_t> & lubricationForce1, Vector3<real_t> & lubricationTorque1 )
{
   WALBERLA_ASSERT_FLOAT_EQUAL(interactionNormal12.length(), real_t(1), "InteractionNormal has to be normalized!");

   real_t epsilon = gapSize / radius1;

   Vector3<real_t> e12 = -interactionNormal12;
   Vector3<real_t> u12 =  u2 - u1;
   Vector3<real_t> omega12 = omega1;

   lubricationForce1 = Vector3<real_t>(real_t(0));
   lubricationTorque1 = Vector3<real_t>(real_t(0));

   if( gapSize < cutOffDistanceNormal )
   {
      // add lubrication force due to normal translation
      real_t epsres = cutOffDistanceNormal / radius1;
      Vector3<real_t> Fn1 = real_t(6) * walberla::math::pi * radius1 * dynamicFluidViscosity *
                            ( real_t(1) / epsilon - real_t(1) / epsres) *
                            math::dot(u12, e12) * e12;

      lubricationForce1 += Fn1;

   }

   if( gapSize < cutOffDistanceTangentialTranslational )
   {
      // add lubrication force and torque due to tangential translation
      real_t epsres = cutOffDistanceTangentialTranslational / radius1;
      real_t logEpsDivEpsres = std::log(epsilon/epsres);
      Vector3<real_t> Ftt1 = real_t(6) * walberla::math::pi * radius1 * dynamicFluidViscosity *
                             ( - real_t(8) / real_t(15) ) * logEpsDivEpsres *
                             (u12 - math::dot(u12,e12) * e12);

      lubricationForce1 += Ftt1;

      Vector3<real_t> Ttt1 = real_t(8) * walberla::math::pi * radius1 * radius1 * dynamicFluidViscosity *
                             ( - real_t(1) / real_t(10) ) * logEpsDivEpsres *
                             math::cross(e12, u12);

      lubricationTorque1 += Ttt1;
   }

   if( gapSize < cutOffDistanceTangentialRotational )
   {
      // add lubrication force and torque due to tangential rotation
      real_t epsres = cutOffDistanceTangentialRotational / radius1;
      real_t logEpsDivEpsres = std::log(epsilon/epsres);
      Vector3<real_t> Ftr1 = real_t(6) * walberla::math::pi * radius1 * radius1 * dynamicFluidViscosity *
                             ( real_t(2) / real_t(15) ) * logEpsDivEpsres *
                             math::cross(omega12, e12);

      lubricationForce1 += Ftr1;

      Vector3<real_t> Ttr1 = real_t(8) * walberla::math::pi * radius1 * radius1 * radius1 * dynamicFluidViscosity *
                             ( real_t(2) / real_t(5) ) * logEpsDivEpsres *
                             (omega1 - math::dot(omega1, e12) * e12);

      lubricationTorque1 += Ttr1;

   }

   // note: lubrication due to normal rotation is dropped here because of too low influence, see also Simeonov

}


/**
 * Applies a correction for the unresolved lubrication forces and torques on the the interacting particles
 *
 * For gap sizes (in cells) below minimalGapSizeFunction(radius), no longer the actual gap size is used in the formulas.
 * This value/function has to be found by calibration with collision experiments.
 * For negative gap sizes, i.e. there is an overlap of the surfaces, no lubrication corrections are applied at any time.
 *
 * Three different cut off distance have to be specified: normal, tangential translational and tangential rotational.
 * These distances have to be determined by extra simulations to evaluate up to which resolution the forces can still be reliably predicted by the simulation itself.
 * Note that these cutoff distances are in cell size "units" and not non-dimensionalized as they are not physically motivated but instead numerical parameters.
 *
 * The formulas for the lubrication corrections for two spheres with different radii are taken from
 * Simeonov, Calantoni - Modeling mechanical contact and lubrication in Direct Numerical Simulations of colliding particles, IJMF, 2012
 *
 * By choosing the respective cutoff distances equal to 0, these contributions can be switched off.
 *
 * Should be used in combination with mesa_pd's contact detection algorithm to determine the interaction partners and the gap size.
 *
 */
class LubricationCorrectionKernel
{
public:

   explicit LubricationCorrectionKernel( real_t dynamicFluidViscosity,
                                         std::function<real_t(real_t)> minimalGapSizeFunction,
                                         real_t cutOffDistanceNormal = real_t(2) / real_t(3),
                                         real_t cutOffDistanceTangentialTranslational = real_t(0.5),
                                         real_t cutOffDistanceTangentialRotational = real_t(0.5) )
         : dynamicFluidViscosity_( dynamicFluidViscosity ), minimalGapSizeFunction_( minimalGapSizeFunction ),
           cutOffDistanceNormal_( cutOffDistanceNormal ), cutOffDistanceTangentialTranslational_( cutOffDistanceTangentialTranslational ),
           cutOffDistanceTangentialRotational_( cutOffDistanceTangentialRotational )
   {}

   template <typename Shape1_T, typename Shape2_T, typename ParticleAccessor_T>
   inline void operator()( const size_t /*idx1*/, const size_t /*idx2*/,
                           const Shape1_T& /*shape1*/, const Shape2_T& /*shape2*/,
                           ParticleAccessor_T& /*ac*/, const Vector3<real_t>& /*interactionNormal*/, const real_t& /*gapSize*/)
   {
      WALBERLA_ABORT("Lubrication correction not implemented!")
   }

   template <typename ParticleAccessor_T>
   inline void operator()( const size_t idx1, const size_t idx2,
                           const mesa_pd::data::Sphere& sphere1, const mesa_pd::data::Sphere& sphere2,
                           ParticleAccessor_T& ac, const Vector3<real_t>& interactionNormal, const real_t& gapSize)
   {
      WALBERLA_ASSERT_UNEQUAL(idx1, idx2, "interacting with itself!");

      // interaction normal is from sphere 2 to sphere 1

      if( gapSize > real_t(0))
      {
         real_t radius1 = sphere1.getRadius();
         real_t radius2 = sphere2.getRadius();
         Vector3<real_t> u1 = ac.getLinearVelocity(idx1);
         Vector3<real_t> u2 = ac.getLinearVelocity(idx2);
         Vector3<real_t> omega1 = ac.getAngularVelocity(idx1);
         Vector3<real_t> omega2 = ac.getAngularVelocity(idx2);

         // compute and add lubrication corrections on sphere 1
         auto gap1 = std::max(gapSize, minimalGapSizeFunction_(radius1) ); // TODO check this, maybe one should use an average radius here to assert symmetry of forces?
         Vector3<real_t> lubricationForce1(real_t(0));
         Vector3<real_t> lubricationTorque1(real_t(0));
         computeLubricationCorrectionSphereSphere(interactionNormal, gap1, dynamicFluidViscosity_, radius1, radius2, u1, u2, omega1, omega2,
                                                  cutOffDistanceNormal_, cutOffDistanceTangentialTranslational_, cutOffDistanceTangentialRotational_,
                                                  lubricationForce1, lubricationTorque1 );

         mesa_pd::addForceAtomic(idx1, ac, lubricationForce1);
         mesa_pd::addTorqueAtomic(idx1, ac, lubricationTorque1);

         // compute and add lubrication corrections on sphere 2
         auto gap2 = std::max(gapSize, minimalGapSizeFunction_(radius2) );
         Vector3<real_t> lubricationForce2(real_t(0));
         Vector3<real_t> lubricationTorque2(real_t(0));
         computeLubricationCorrectionSphereSphere( - interactionNormal, gap2, dynamicFluidViscosity_, radius2, radius1, u2, u1, omega2, omega1,
                                                  cutOffDistanceNormal_, cutOffDistanceTangentialTranslational_, cutOffDistanceTangentialRotational_,
                                                  lubricationForce2, lubricationTorque2 );

         mesa_pd::addForceAtomic(idx2, ac, lubricationForce2);
         mesa_pd::addTorqueAtomic(idx2, ac, lubricationTorque2);

      }
      // else: no lubrication correction

   }

   template <typename ParticleAccessor_T>
   inline void operator()( const size_t idx1, const size_t idx2,
                           const mesa_pd::data::Sphere& sphere, const mesa_pd::data::HalfSpace& /*halfSpace*/,
                           ParticleAccessor_T& ac, const Vector3<real_t>& interactionNormal, const real_t& gapSize)
   {
      // interaction normal is the normal of the half space

      if( gapSize > real_t(0))
      {
         real_t radius1 = sphere.getRadius();
         Vector3<real_t> u1 = ac.getLinearVelocity(idx1);
         Vector3<real_t> u2 = ac.getLinearVelocity(idx2);
         Vector3<real_t> omega1 = ac.getAngularVelocity(idx1);
         // angular velocity of half space is zero!

         // compute and add lubrication corrections on sphere 1
         auto gap1 = std::max(gapSize, minimalGapSizeFunction_(radius1) );
         Vector3<real_t> lubricationForce1(real_t(0));
         Vector3<real_t> lubricationTorque1(real_t(0));
         computeLubricationCorrectionSphereHalfSpace(interactionNormal, gap1, dynamicFluidViscosity_, radius1, u1, u2, omega1,
                                                     cutOffDistanceNormal_, cutOffDistanceTangentialTranslational_, cutOffDistanceTangentialRotational_,
                                                     lubricationForce1, lubricationTorque1 );

         mesa_pd::addForceAtomic(idx1, ac, lubricationForce1);
         mesa_pd::addTorqueAtomic(idx1, ac, lubricationTorque1);

         // NOTE: no lubrication corrections added to the half space!

      }
      // else: no lubrication correction
   }

   template <typename ParticleAccessor_T>
   inline void operator()( const size_t idx1, const size_t idx2,
                           const mesa_pd::data::HalfSpace& halfSpace, const mesa_pd::data::Sphere& sphere,
                           ParticleAccessor_T& ac, const Vector3<real_t>& interactionNormal, const real_t& gapSize)
   {
      operator()(idx2, idx1, sphere, halfSpace, ac, -interactionNormal, gapSize);
   }

   void setDynamicFluidViscosity( const real_t dynamicFluidViscosity){ dynamicFluidViscosity_ = dynamicFluidViscosity; }
   void setMinimalGapSizeFunction( std::function<real_t(real_t)> minimalGapSizeFunction ){ minimalGapSizeFunction_ = minimalGapSizeFunction; }
   void setNormalCutOffDistance( const real_t cutOffDistanceNormal){ cutOffDistanceNormal_ = cutOffDistanceNormal; }
   void setTangentialTranslationalCutOffDistance( const real_t cutOffDistanceTangentialTranslational){ cutOffDistanceTangentialTranslational_ = cutOffDistanceTangentialTranslational; }
   void setTangentialRotationalCutOffDistance( const real_t cutOffDistanceTangentialRotational){ cutOffDistanceTangentialRotational_ = cutOffDistanceTangentialRotational; }

   real_t getDynamicFluidViscosity() const { return dynamicFluidViscosity_; }
   real_t getMinimalGapSize(real_t radius) const { return minimalGapSizeFunction_(radius); }
   real_t getNormalCutOffDistance() const { return cutOffDistanceNormal_; }
   real_t getTangentialTranslationalCutOffDistance() const { return cutOffDistanceTangentialTranslational_; }
   real_t getTangentialRotationalCutOffDistance() const { return cutOffDistanceTangentialRotational_; }

private:

   real_t dynamicFluidViscosity_;
   std::function<real_t(real_t)> minimalGapSizeFunction_;
   real_t cutOffDistanceNormal_;
   real_t cutOffDistanceTangentialTranslational_;
   real_t cutOffDistanceTangentialRotational_;
};

} //namespace lbm_mesapd_coupling
} //namespace walberla