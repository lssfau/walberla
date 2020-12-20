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
 * This DEM kernel supports spring&dashpot in normal direction as well as friction in tangential direction.
 *
 * \code
 * const walberla::mesa_pd::Vec3& getPosition(const size_t p_idx) const;
 *
 * const walberla::mesa_pd::Vec3& getLinearVelocity(const size_t p_idx) const;
 *
 * walberla::mesa_pd::Vec3& getForceRef(const size_t p_idx);
 *
 * const walberla::mesa_pd::Vec3& getAngularVelocity(const size_t p_idx) const;
 *
 * walberla::mesa_pd::Vec3& getTorqueRef(const size_t p_idx);
 *
 * const uint_t& getType(const size_t p_idx) const;
 *
 * const std::map<walberla::id_t, walberla::mesa_pd::Vec3>& getContactHistory(const size_t p_idx) const;
 * void setContactHistory(const size_t p_idx, const std::map<walberla::id_t, walberla::mesa_pd::Vec3>& v);
 *
 * \endcode
 * \ingroup mesa_pd_kernel
 */
class SpringDashpot
{
public:
   SpringDashpot(const uint_t numParticleTypes);
   SpringDashpot(const SpringDashpot& other) = default;
   SpringDashpot(SpringDashpot&& other) = default;
   SpringDashpot& operator=(const SpringDashpot& other) = default;
   SpringDashpot& operator=(SpringDashpot&& other) = default;

   template <typename Accessor>
   void operator()(const size_t p_idx1,
                   const size_t p_idx2,
                   Accessor& ac,
                   const Vec3& contactPoint,
                   const Vec3& contactNormal,
                   const real_t& penetrationDepth) const;

   
   /// assumes this parameter is symmetric
   void setStiffness(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setDampingN(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setDampingT(const size_t type1, const size_t type2, const real_t& val);
   /// assumes this parameter is symmetric
   void setFriction(const size_t type1, const size_t type2, const real_t& val);

   
   real_t getStiffness(const size_t type1, const size_t type2) const;
   real_t getDampingN(const size_t type1, const size_t type2) const;
   real_t getDampingT(const size_t type1, const size_t type2) const;
   real_t getFriction(const size_t type1, const size_t type2) const;

   inline
   real_t calcCoefficientOfRestitution(const size_t type1,
                                       const size_t type2,
                                       const real_t meff)
   {
      auto a = real_t(0.5) * getDampingN(type1, type2) / meff;
      return std::exp(-a * math::pi / std::sqrt(getStiffness(type1, type2) / meff - a*a));
   }

   inline
   real_t calcCollisionTime(const size_t type1,
                            const size_t type2,
                            const real_t meff)
   {
      auto a = real_t(0.5) * getDampingN(type1, type2) / meff;
      return math::pi / std::sqrt( getStiffness(type1, type2)/meff - a*a);
   }

   inline
   void setParametersFromCOR(const size_t type1,
                             const size_t type2,
                             const real_t cor,
                             const real_t collisionTime,
                             const real_t meff)
   {
      const real_t lnDryResCoeff = std::log(cor);
      setStiffness(type1, type2, math::pi * math::pi * meff / ( collisionTime * collisionTime * ( real_t(1) - lnDryResCoeff * lnDryResCoeff / ( math::pi * math::pi + lnDryResCoeff* lnDryResCoeff ))  ));
      setDampingN( type1, type2, - real_t(2) * std::sqrt( meff * getStiffness(type1, type2) ) * ( lnDryResCoeff / std::sqrt( math::pi * math::pi + ( lnDryResCoeff * lnDryResCoeff ) ) ));
   }
private:
   uint_t numParticleTypes_;
   
   std::vector<real_t> stiffness_ {};
   std::vector<real_t> dampingN_ {};
   std::vector<real_t> dampingT_ {};
   std::vector<real_t> friction_ {};
};

SpringDashpot::SpringDashpot(const uint_t numParticleTypes)
{
   numParticleTypes_ = numParticleTypes;
   
   stiffness_.resize(numParticleTypes * numParticleTypes, real_t(0));
   dampingN_.resize(numParticleTypes * numParticleTypes, real_t(0));
   dampingT_.resize(numParticleTypes * numParticleTypes, real_t(0));
   friction_.resize(numParticleTypes * numParticleTypes, real_t(0));
}


inline void SpringDashpot::setStiffness(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   stiffness_[numParticleTypes_*type1 + type2] = val;
   stiffness_[numParticleTypes_*type2 + type1] = val;
}
inline void SpringDashpot::setDampingN(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   dampingN_[numParticleTypes_*type1 + type2] = val;
   dampingN_[numParticleTypes_*type2 + type1] = val;
}
inline void SpringDashpot::setDampingT(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   dampingT_[numParticleTypes_*type1 + type2] = val;
   dampingT_[numParticleTypes_*type2 + type1] = val;
}
inline void SpringDashpot::setFriction(const size_t type1, const size_t type2, const real_t& val)
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   friction_[numParticleTypes_*type1 + type2] = val;
   friction_[numParticleTypes_*type2 + type1] = val;
}


inline real_t SpringDashpot::getStiffness(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( stiffness_[numParticleTypes_*type1 + type2],
                                stiffness_[numParticleTypes_*type2 + type1],
                                "parameter matrix for stiffness not symmetric!");
   return stiffness_[numParticleTypes_*type1 + type2];
}
inline real_t SpringDashpot::getDampingN(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( dampingN_[numParticleTypes_*type1 + type2],
                                dampingN_[numParticleTypes_*type2 + type1],
                                "parameter matrix for dampingN not symmetric!");
   return dampingN_[numParticleTypes_*type1 + type2];
}
inline real_t SpringDashpot::getDampingT(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( dampingT_[numParticleTypes_*type1 + type2],
                                dampingT_[numParticleTypes_*type2 + type1],
                                "parameter matrix for dampingT not symmetric!");
   return dampingT_[numParticleTypes_*type1 + type2];
}
inline real_t SpringDashpot::getFriction(const size_t type1, const size_t type2) const
{
   WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
   WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
   WALBERLA_ASSERT_FLOAT_EQUAL( friction_[numParticleTypes_*type1 + type2],
                                friction_[numParticleTypes_*type2 + type1],
                                "parameter matrix for friction not symmetric!");
   return friction_[numParticleTypes_*type1 + type2];
}

template <typename Accessor>
inline void SpringDashpot::operator()(const size_t p_idx1,
                                      const size_t p_idx2,
                                      Accessor& ac,
                                      const Vec3& contactPoint,
                                      const Vec3& contactNormal,
                                      const real_t& penetrationDepth) const
{
   static_assert(std::is_base_of<data::IAccessor, Accessor>::value, "please provide a valid accessor");

   if (p_idx1 != p_idx2)
   {
      // Global position of contact
      const Vec3& gpos( contactPoint );
      // The absolute value of the penetration length
      real_t delta = -penetrationDepth;
      if (delta < real_t(0)) return;

      const Vec3   relVel ( -(getVelocityAtWFPoint(p_idx1, ac, gpos) - getVelocityAtWFPoint(p_idx2, ac, gpos)) );
      const real_t relVelN( math::dot(relVel, contactNormal) );
      const Vec3   relVelT( relVel - ( relVelN * contactNormal ) );

      // Calculating the normal force based on a linear spring-dashpot force model
      real_t fNabs = getStiffness(ac.getType(p_idx1), ac.getType(p_idx2)) * delta + getDampingN(ac.getType(p_idx1), ac.getType(p_idx2)) * relVelN;
      const Vec3& fN = fNabs * contactNormal;

      // Calculating the tangential force based on the model by Haff and Werner
      const real_t fTabs( std::min( getDampingT(ac.getType(p_idx1), ac.getType(p_idx2)) * relVelT.length(), getFriction(ac.getType(p_idx1), ac.getType(p_idx2)) * fNabs ) );
      const Vec3   fT   ( fTabs * relVelT.getNormalizedOrZero() );

      // Add normal force at contact point
      addForceAtWFPosAtomic( p_idx1, ac,  fN, gpos );
      addForceAtWFPosAtomic( p_idx2, ac, -fN, gpos );

      // Add tangential force at contact point
      addForceAtWFPosAtomic( p_idx1, ac,  fT, gpos );
      addForceAtWFPosAtomic( p_idx2, ac, -fT, gpos );
   }
}

} //namespace kernel
} //namespace mesa_pd
} //namespace walberla