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
//! \file InitContactsForHCSITS.h
//! \author Tobias Leemann <tobias.leemann@fau.de>
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
#include <mesa_pd/data/IContactAccessor.h>
#include <mesa_pd/data/Flags.h>
#include <core/logging/Logging.h>

#include <vector>

namespace walberla {
namespace mesa_pd {
namespace kernel {

/**
* Init the datastructures for a contact for later use of the HCSITS-Solver.
* Call this kernel on all contacts that you want to treat with HCSITS before performing any relaxation timesteps.
* Use setErp() to set the error reduction parameter, which defines the share of the overlap that is resolved in one step
* (0 meaning after the relaxation the overlap is the same, 1 meaning the particles will have no overlap after relaxation).
* Use setFriction(a,b, cof) to define the coefficient of friction cof between materials a and b. It is assumed to be
 * symmetric w.r.t. the materials.
* \ingroup mesa_pd_kernel
*/
class InitContactsForHCSITS{
public:

   // Default constructor sets the default values
   explicit InitContactsForHCSITS(const uint_t numParticleTypes) :
   numParticleTypes_ (numParticleTypes),
   {%- for prop in properties %}
   {{prop.name}}_( {{prop.defValue}} ){{ "," if not loop.last}}
   {%- endfor %}
   {

      {% for param in material_parameters %}
      {{param}}_.resize(numParticleTypes * numParticleTypes, real_t(0));
      {%- endfor %}
   }

   // Getter and Setter Functions
   {%- for prop in properties %}
   inline const {{prop.type}}& get{{prop.name | capFirst}}() const {return {{prop.name}}_;}
   inline void set{{prop.name | capFirst}}({{prop.type}} v) { {{prop.name}}_ = v;}
   {% endfor %}

   // Getter and Setter Functions for material parameter combinations (they are assumed symmetric).
   {% for param in material_parameters %}
   inline void set{{param | capFirst}}(const size_t type1, const size_t type2, const real_t& val)
   {
      WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
      WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
      {{param}}_[numParticleTypes_*type1 + type2] = val;
      {{param}}_[numParticleTypes_*type2 + type1] = val;
   }
   {%- endfor %}

   {% for param in material_parameters %}
   inline real_t get{{param | capFirst}}(const size_t type1, const size_t type2) const
   {
      WALBERLA_ASSERT_LESS( type1, numParticleTypes_ );
      WALBERLA_ASSERT_LESS( type2, numParticleTypes_ );
      WALBERLA_ASSERT_FLOAT_EQUAL( {{param}}_[numParticleTypes_*type1 + type2],
              {{param}}_[numParticleTypes_*type2 + type1],
              "parameter matrix for {{param}} not symmetric!");
      return {{param}}_[numParticleTypes_*type1 + type2];
   }
   {%- endfor %}

   // List material parameters
   {% for param in material_parameters %}
   std::vector<real_t> {{param}}_ {};
   {%- endfor %}

   template<typename CAccessor, typename PAccessor>
   void operator()(size_t c, CAccessor &ca, PAccessor &pa);

private:
   // List properties
   const uint_t numParticleTypes_;

   {% for prop in properties %}
   {{prop.type}} {{prop.name }}_;
   {%- endfor %}

   {% for param in parameters %}
   std::vector<real_t> {{param}}_ {};
   {%- endfor %}

};


template<typename CAccessor, typename PAccessor>
inline void InitContactsForHCSITS::operator()(size_t c, CAccessor &ca, PAccessor &pa) {

      static_assert(std::is_base_of<data::IContactAccessor, CAccessor>::value, "please provide a valid contact accessor");
      static_assert(std::is_base_of<data::IAccessor, PAccessor>::value, "please provide a valid particle accessor");

   size_t bId1 = ca.getId1(c);
   size_t bId2 = ca.getId2(c);
   ca.setR1(c, ca.getPosition(c) - pa.getPosition(bId1));
   ca.setR2(c, ca.getPosition(c) - pa.getPosition(bId2));


   // construct vector perpendicular to the normal (cross product with cardinal basis vector where the 1 component is where the other vector has its component of smallest magnitude)
   if (std::fabs(ca.getNormalRef(c)[0]) < std::fabs(ca.getNormalRef(c)[1])) {
      if (std::fabs(ca.getNormalRef(c)[0]) < std::fabs(ca.getNormalRef(c)[2]))
         ca.setT(c, Vec3(0, ca.getNormalRef(c)[2], -ca.getNormalRef(c)[1]));   // = n x (1, 0, 0)^T
      else
         ca.setT(c, Vec3(ca.getNormalRef(c)[1], -ca.getNormalRef(c)[0], 0));   // = n x (0, 0, 1)^T
   } else {
      if (std::fabs(ca.getNormalRef(c)[1]) < std::fabs(ca.getNormalRef(c)[2]))
         ca.setT(c, Vec3(-ca.getNormalRef(c)[2], 0, ca.getNormalRef(c)[0]));   // = n x (0, 1, 0)^T
      else
         ca.setT(c, Vec3(ca.getNormalRef(c)[1], -ca.getNormalRef(c)[0], 0));   // = n x (0, 0, 1)^T
   }
   normalize(ca.getTRef(c));
   ca.setO(c, ca.getNormal(c) % ca.getT(c));

   Mat3 contactframe(ca.getNormal(c), ca.getT(c), ca.getO(c));

   // If the distance is negative then penetration is present. This is an error and should be corrected.
   // Correcting the whole error is not recommended since due to the linearization the errors cannot
   // completely fixed anyway and the error reduction will introduce artificial restitution.
   // However, if the distance is positive then it is not about error correction but the distance that
   // can still be overcome without penetration and thus no error correction parameter should be applied.
   if (ca.getDistance(c) < real_t(0.0)) {
      setMaximumPenetration(std::max(getMaximumPenetration(), -ca.getDistance(c)));
      ca.getDistanceRef(c) *= erp_;
   }

   ca.getMuRef(c) = getFriction(pa.getType(bId1), pa.getType(bId2));

   Mat3 diag = -(
           math::skewSymCrossProduct(ca.getR1(c),
                                     math::skewSymCrossProduct(pa.getInvInertia(bId1), ca.getR1(c)))
           + math::skewSymCrossProduct(ca.getR2(c),
                                       math::skewSymCrossProduct(pa.getInvInertia(bId2), ca.getR2(c))));
   diag[0] += pa.getInvMass(bId1) + pa.getInvMass(bId2);
   diag[4] += pa.getInvMass(bId1) + pa.getInvMass(bId2);
   diag[8] += pa.getInvMass(bId1) + pa.getInvMass(bId2);

   diag = contactframe.getTranspose() * diag * contactframe;

   // Diagonal block is known to be positive-definite and thus an inverse always exists.
   ca.getDiag_ntoRef(c) = diag;
   ca.getDiag_nto_invRef(c) = diag.getInverse();
   ca.getDiag_n_invRef(c) = math::inv(diag[0]);
   ca.getDiag_to_invRef(c) = Mat2(diag[4], diag[5], diag[7], diag[8]).getInverse();
}

}
}
}
