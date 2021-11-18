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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/common/ParticleFunctions.h>
#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/IAccessor.h>

#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Advanced DEM kernel
 *
 * This model is a linearized version of
 * Edward Biegert, Bernhard Vowinckel, Eckart Meiburg
 * A collision model for grain-resolving simulations of flows over dense, mobile, polydisperse granular sediment beds
 * https://doi.org/10.1016/j.jcp.2017.03.035
 *
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
class LinearSpringDashpot
{
public:
   LinearSpringDashpot(const uint_t numParticleTypes);
   LinearSpringDashpot(const LinearSpringDashpot& other) = default;
   LinearSpringDashpot(LinearSpringDashpot&& other) = default;
   LinearSpringDashpot& operator=(const LinearSpringDashpot& other) = default;
   LinearSpringDashpot& operator=(LinearSpringDashpot&& other) = default;

   template <typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor& ac,
                   const Vec3& contactPoint,
                   const Vec3& contactNormal,
                   const real_t& penetrationDepth,
                   const real_t& dt) const;

   real_t calcCoefficientOfRestitution(const size_t type1,
                                       const size_t type2,
                                       const real_t effectiveMass);

   real_t calcCollisionTime(const size_t type1,
                            const size_t type2,
                            const real_t effectiveMass);

   void setStiffnessAndDamping(const size_t type1,
                               const size_t type2,
                               const real_t coefficientOfRestitution,
                               const real_t collisionTime,
                               const real_t kappa,
                               const real_t effectiveMass);

   {% for param in parameters %}
   /// assumes this parameter is symmetric
   void set{{param | capFirst}}(const size_t type1, const size_t type2, const real_t& val);
   {%- endfor %}

   {% for param in parameters %}
   real_t get{{param | capFirst}}(const size_t type1, const size_t type2) const;
   {%- endfor %}
private:
   uint_t numParticleTypes_;
   {% for param in parameters %}
   std::vector<real_t> {{param}}_ {};
   {%- endfor %}
};

LinearSpringDashpot::LinearSpringDashpot(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   {% for param in parameters %}
   {{param}}_.resize(numParticleTypes * numParticleTypes, real_t(0));
   {%- endfor %}
}

{% for param in parameters %}
inline void LinearSpringDashpot::set{{param | capFirst}}(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   {{param}}_[numParticleTypes_*type1 + type2] = val;
   {{param}}_[numParticleTypes_*type2 + type1] = val;
}
{%- endfor %}

{% for param in parameters %}
inline real_t LinearSpringDashpot::get{{param | capFirst}}(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( {{param}}_[numParticleTypes_*type1 + type2],
                                {{param}}_[numParticleTypes_*type2 + type1],
                                "parameter matrix for {{param}} not symmetric!");
   return {{param}}_[numParticleTypes_*type1 + type2];
}
{%- endfor %}

template <typename Accessor>
inline void LinearSpringDashpot::operator()(const size_t p_idx1,
                                            const size_t p_idx2,
                                            Accessor& ac,
                                            const Vec3& contactPoint,
                                            const Vec3& contactNormal,
                                            const real_t& penetrationDepth,
                                            const real_t& dt) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (p_idx1 != p_idx2)
   {

      WALBERLA_ASSERT_FLOAT_EQUAL(math::sqrLength(contactNormal), real_t(1));

      real_t delta = -penetrationDepth;
      if (delta < real_t(0)) return;


      const Vec3 relVel ( -(getVelocityAtWFPoint(p_idx1, ac, contactPoint) - getVelocityAtWFPoint(p_idx2, ac, contactPoint)) );
      const Vec3 relVelN( math::dot(relVel, contactNormal) * contactNormal );
      const Vec3 relVelT( relVel - relVelN );


      // calculate the normal force based on a linear spring-dashpot force model
      Vec3 fN = getStiffnessN(ac.getType(p_idx1), ac.getType(p_idx2)) * delta * contactNormal +
                getDampingN(ac.getType(p_idx1), ac.getType(p_idx2)) * relVelN;


      // get further collision properties from the contact history or initialize them
      Vec3 tangentialSpringDisplacement( real_t(0) );
      real_t impactVelocityMagnitude(real_t(0));
      bool isSticking = false;
      auto contactHistory = ac.getOldContactHistoryRef(p_idx1).find(ac.getUid(p_idx2)); //TODO assert symmetry
      if(contactHistory != ac.getOldContactHistoryRef(p_idx1).end())
      {
         // get infos from the contact history
         tangentialSpringDisplacement = contactHistory->second.getTangentialSpringDisplacement();
         isSticking = contactHistory->second.getIsSticking();
         impactVelocityMagnitude = contactHistory->second.getImpactVelocityMagnitude();

      }else
      {
         // new contact: initialize values
         impactVelocityMagnitude = relVel.length();
      }

      //TODO: move to own tangential displacement integration kernel
      Vec3 rotatedTangentialDisplacement = tangentialSpringDisplacement - contactNormal * (contactNormal * tangentialSpringDisplacement);
      Vec3 newTangentialSpringDisplacement = rotatedTangentialDisplacement.sqrLength() <= real_t(0) ? // avoid division by zero
                                             Vec3(real_t(0)) :
                                             ( rotatedTangentialDisplacement * std::sqrt((tangentialSpringDisplacement.sqrLength() / rotatedTangentialDisplacement.sqrLength())));
      newTangentialSpringDisplacement = newTangentialSpringDisplacement + dt * relVelT;

      // calculate the tangential force based on a linear spring-dashpot force model
      const real_t stiffnessT = getStiffnessT(ac.getType(p_idx1), ac.getType(p_idx2));
      const real_t dampingT = getDampingT(ac.getType(p_idx1), ac.getType(p_idx2));
      Vec3 fTLS = stiffnessT * newTangentialSpringDisplacement +
                  dampingT * relVelT;

      const Vec3 t = fTLS.getNormalizedIfNotZero(); // tangential unit vector

      // calculate friction force
      const real_t fFrictionAbsStatic = getFrictionCoefficientStatic(ac.getType(p_idx1), ac.getType(p_idx2)) * fN.length(); // sticking, rolling
      const real_t fFrictionAbsDynamic = getFrictionCoefficientDynamic(ac.getType(p_idx1), ac.getType(p_idx2)) * fN.length(); // sliding

      const real_t tangentialVelocityThreshold = real_t(1e-8);

      real_t fFrictionAbs;
      if( isSticking && relVelT.length() < tangentialVelocityThreshold && fTLS.length() < fFrictionAbsStatic  )
      {
         fFrictionAbs = fFrictionAbsStatic;
      }
      else if( isSticking && fTLS.length() < fFrictionAbsDynamic )
      {
         // sticking
         fFrictionAbs = fFrictionAbsDynamic;
      }
      else
      {
         // slipping
         fFrictionAbs = fFrictionAbsDynamic;

         isSticking = false;

         // reset displacement vector
         if(stiffnessT > real_t(0) ) newTangentialSpringDisplacement = ( fFrictionAbs * t - dampingT * relVelT ) / stiffnessT;

         // if tangential force falls below coulomb limit, we are back in sticking
         if( fTLS.length() < fFrictionAbsDynamic )
         {
            isSticking = true;
         }
      }

      const real_t fTabs( std::min( fTLS.length(), fFrictionAbs) );
      const Vec3   fT   ( fTabs * t );

      //TODO check if tangential spring displacement is same for symmetric case
      auto& ch1 = ac.getNewContactHistoryRef(p_idx1)[ac.getUid(p_idx2)];
      ch1.setTangentialSpringDisplacement(newTangentialSpringDisplacement);
      ch1.setIsSticking(isSticking);
      ch1.setImpactVelocityMagnitude(impactVelocityMagnitude);

      auto& ch2 = ac.getNewContactHistoryRef(p_idx2)[ac.getUid(p_idx1)];
      ch2.setTangentialSpringDisplacement(newTangentialSpringDisplacement);
      ch2.setIsSticking(isSticking);
      ch2.setImpactVelocityMagnitude(impactVelocityMagnitude);

      // Add normal force at contact point
      addForceAtWFPosAtomic( p_idx1, ac,  fN, contactPoint );
      addForceAtWFPosAtomic( p_idx2, ac, -fN, contactPoint );

      // Add tangential force at contact point
      addForceAtWFPosAtomic( p_idx1, ac,  fT, contactPoint );
      addForceAtWFPosAtomic( p_idx2, ac, -fT, contactPoint );

   }
}

inline real_t LinearSpringDashpot::calcCoefficientOfRestitution(const size_t type1,
                                                                const size_t type2,
                                                                const real_t effectiveMass)
{
   auto a = real_t(0.5) * getDampingN(type1, type2) / effectiveMass;
   return std::exp(-a * math::pi / std::sqrt(getStiffnessN(type1, type2) / effectiveMass - a*a));
}

inline real_t LinearSpringDashpot::calcCollisionTime(const size_t type1,
                                                     const size_t type2,
                                                     const real_t effectiveMass)
{
   auto a = real_t(0.5) * getDampingN(type1, type2) / effectiveMass;
   return math::pi / std::sqrt( getStiffnessN(type1, type2)/effectiveMass - a*a);
}

inline void LinearSpringDashpot::setStiffnessAndDamping(const size_t type1,
                                                        const size_t type2,
                                                        const real_t coefficientOfRestitution,
                                                        const real_t collisionTime,
                                                        const real_t kappa,
                                                        const real_t effectiveMass)
{
   const real_t lnDryResCoeff = std::log(coefficientOfRestitution);

   setStiffnessN(type1, type2, effectiveMass * ( math::pi * math::pi + lnDryResCoeff * lnDryResCoeff ) / (collisionTime * collisionTime) );
   setDampingN  (type1, type2, - real_t(2) * effectiveMass * lnDryResCoeff / collisionTime );

   setStiffnessT(type1, type2, kappa * getStiffnessN(type1, type2) );
   setDampingT  (type1, type2, std::sqrt(kappa) * getDampingN(type1, type2) );
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla
