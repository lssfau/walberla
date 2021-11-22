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

#include <core/math/Constants.h>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
 * Basic DEM kernel
 *
 * DEM kernel with spring-dashpot interaction in normal direction and spring in tangential direction.
 *
 * \ingroup mesa_pd_kernel
 */
class SpringDashpotSpring
{
public:
   SpringDashpotSpring(const uint_t numParticleTypes);
   SpringDashpotSpring(const SpringDashpotSpring& other) = default;
   SpringDashpotSpring(SpringDashpotSpring&& other) = default;
   SpringDashpotSpring& operator=(const SpringDashpotSpring& other) = default;
   SpringDashpotSpring& operator=(SpringDashpotSpring&& other) = default;

   template <typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor& ac,
                   const Vec3& contactPoint,
                   const Vec3& contactNormal,
                   const real_t& penetrationDepth,
                   const real_t dt) const;

   
   /// assumes this parameter is symmetric
   void setStiffnessN(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setDampingN(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setStiffnessT(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setCoefficientOfFriction(const size_t type1, const size_t type2, const real_t& val);

   
   real_t getStiffnessN(const size_t type1, const size_t type2) const;
   real_t getDampingN(const size_t type1, const size_t type2) const;
   real_t getStiffnessT(const size_t type1, const size_t type2) const;
   real_t getCoefficientOfFriction(const size_t type1, const size_t type2) const;

   inline
   real_t calcCoefficientOfRestitution(const size_t type1,
                                       const size_t type2,
                                       const real_t meff)
   {
      auto a = real_t(0.5) * getDampingN(type1, type2) / meff;
      return std::exp(-a * math::pi / std::sqrt(getStiffnessN(type1, type2) / meff - a*a));
   }

   inline
   real_t calcCollisionTime(const size_t type1,
                            const size_t type2,
                            const real_t meff)
   {
      auto a = real_t(0.5) * getDampingN(type1, type2) / meff;
      return math::pi / std::sqrt( getStiffnessN(type1, type2)/meff - a*a);
   }

   inline
   void setParametersFromCOR(const size_t type1,
                             const size_t type2,
                             const real_t cor,
                             const real_t collisionTime,
                             const real_t meff)
   {
      const real_t lnDryResCoeff = std::log(cor);
      setStiffnessN(type1, type2, math::pi * math::pi * meff / ( collisionTime * collisionTime * ( real_t(1) - lnDryResCoeff * lnDryResCoeff / ( math::pi * math::pi + lnDryResCoeff* lnDryResCoeff ))  ));
      setDampingN( type1, type2, - real_t(2) * std::sqrt( meff * getStiffnessN(type1, type2) ) * ( lnDryResCoeff / std::sqrt( math::pi * math::pi + ( lnDryResCoeff * lnDryResCoeff ) ) ));
   }
private:
   uint_t numParticleTypes_;
   
   std::vector<real_t> stiffnessN_ {};
   std::vector<real_t> dampingN_ {};
   std::vector<real_t> stiffnessT_ {};
   std::vector<real_t> coefficientOfFriction_ {};
};

SpringDashpotSpring::SpringDashpotSpring(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   
   stiffnessN_.resize(numParticleTypes * numParticleTypes, real_t(0));
   dampingN_.resize(numParticleTypes * numParticleTypes, real_t(0));
   stiffnessT_.resize(numParticleTypes * numParticleTypes, real_t(0));
   coefficientOfFriction_.resize(numParticleTypes * numParticleTypes, real_t(0));
}


inline void SpringDashpotSpring::setStiffnessN(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   stiffnessN_[numParticleTypes_*type1 + type2] = val;
   stiffnessN_[numParticleTypes_*type2 + type1] = val;
}
inline void SpringDashpotSpring::setDampingN(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   dampingN_[numParticleTypes_*type1 + type2] = val;
   dampingN_[numParticleTypes_*type2 + type1] = val;
}
inline void SpringDashpotSpring::setStiffnessT(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   stiffnessT_[numParticleTypes_*type1 + type2] = val;
   stiffnessT_[numParticleTypes_*type2 + type1] = val;
}
inline void SpringDashpotSpring::setCoefficientOfFriction(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   coefficientOfFriction_[numParticleTypes_*type1 + type2] = val;
   coefficientOfFriction_[numParticleTypes_*type2 + type1] = val;
}


inline real_t SpringDashpotSpring::getStiffnessN(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( stiffnessN_[numParticleTypes_*type1 + type2],
                                stiffnessN_[numParticleTypes_*type2 + type1],
                                "parameter matrix for stiffnessN not symmetric!");
   return stiffnessN_[numParticleTypes_*type1 + type2];
}
inline real_t SpringDashpotSpring::getDampingN(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( dampingN_[numParticleTypes_*type1 + type2],
                                dampingN_[numParticleTypes_*type2 + type1],
                                "parameter matrix for dampingN not symmetric!");
   return dampingN_[numParticleTypes_*type1 + type2];
}
inline real_t SpringDashpotSpring::getStiffnessT(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( stiffnessT_[numParticleTypes_*type1 + type2],
                                stiffnessT_[numParticleTypes_*type2 + type1],
                                "parameter matrix for stiffnessT not symmetric!");
   return stiffnessT_[numParticleTypes_*type1 + type2];
}
inline real_t SpringDashpotSpring::getCoefficientOfFriction(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( coefficientOfFriction_[numParticleTypes_*type1 + type2],
                                coefficientOfFriction_[numParticleTypes_*type2 + type1],
                                "parameter matrix for coefficientOfFriction not symmetric!");
   return coefficientOfFriction_[numParticleTypes_*type1 + type2];
}

template <typename Accessor>
inline void SpringDashpotSpring::operator()(const size_t p_idx1,
                                      const size_t p_idx2,
                                      Accessor& ac,
                                      const Vec3& contactPoint,
                                      const Vec3& contactNormal,
                                      const real_t& penetrationDepth,
                                      const real_t dt) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (p_idx1 != p_idx2)
   {
      // skip if no penetration is present
      real_t delta = -penetrationDepth;
      if (delta < real_t(0)) return;

      // calculate relative velocities
      const Vec3   relVel ( -(getVelocityAtWFPoint(p_idx1, ac, contactPoint) - getVelocityAtWFPoint(p_idx2, ac, contactPoint)) );
      const real_t relVelN( math::dot(relVel, contactNormal) );
      const Vec3   relVelT( relVel - ( relVelN * contactNormal ) );
      const Vec3   contactTangent = relVelT.getNormalizedIfNotZero();

      // Calculating the normal force based on a linear spring-dashpot force model
      real_t fNabs = getStiffnessN(ac.getType(p_idx1), ac.getType(p_idx2)) * delta + getDampingN(ac.getType(p_idx1), ac.getType(p_idx2)) * relVelN;
      const Vec3 fN = fNabs * contactNormal;

      // get tangential displacement from contact history
      auto tangentialDisplacement = Vec3(real_t(0));
      auto contactHistory = ac.getOldContactHistoryRef(p_idx1).find(ac.getUid(p_idx2));
      if(contactHistory != ac.getOldContactHistoryRef(p_idx1).end())
      {
         // get infos from the contact history
         tangentialDisplacement = dot(contactHistory->second.getTangentialSpringDisplacement(), contactTangent) * contactTangent;
      }

      // accumulate tangential displacement
      tangentialDisplacement += relVelT * dt;

      // Calculating the tangential force
      const auto maxTangentialForce = fNabs * getCoefficientOfFriction(ac.getType(p_idx1), ac.getType(p_idx2));
      Vec3 fT = getStiffnessT(ac.getType(p_idx1), ac.getType(p_idx2)) * tangentialDisplacement;
      if (length(fT) > maxTangentialForce)
         fT = maxTangentialForce * fT.getNormalizedIfNotZero();

      // store new tangential displacements
      auto& ch1 = ac.getNewContactHistoryRef(p_idx1)[ac.getUid(p_idx2)];
      ch1.setTangentialSpringDisplacement(tangentialDisplacement);

      auto& ch2 = ac.getNewContactHistoryRef(p_idx2)[ac.getUid(p_idx1)];
      ch2.setTangentialSpringDisplacement(tangentialDisplacement);

      // Add normal force at contact point
      addForceAtWFPosAtomic( p_idx1, ac,  fN, contactPoint );
      addForceAtWFPosAtomic( p_idx2, ac, -fN, contactPoint );

      // Add tangential force at contact point
      addForceAtWFPosAtomic( p_idx1, ac,  fT, contactPoint );
      addForceAtWFPosAtomic( p_idx2, ac, -fT, contactPoint );
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla