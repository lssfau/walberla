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
//! \file ShearRate.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"

#include <type_traits>



// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {

namespace internal {

    template<typename CollisionModel_T, class Enable = void>
    struct ShearRelaxationParameter {
        static_assert(never_true<CollisionModel_T>::value,
                      "For your current LB collision model, there is yet no implementation for calculating the shear rate!");
    };

    template<typename CollisionModel_T>
    struct ShearRelaxationParameter<CollisionModel_T, typename std::enable_if<std::is_same<typename CollisionModel_T::tag,
                                                                                                 collision_model::SRT_tag>::value>::type>
    {
        static inline real_t get(const CollisionModel_T &lm, cell_idx_t x, cell_idx_t y, cell_idx_t z)
        {
           return lm.omega(x, y, z);
        }
    };

    template<typename CollisionModel_T>
    struct ShearRelaxationParameter<CollisionModel_T, typename std::enable_if<std::is_same<typename CollisionModel_T::tag,
                                                                                                 collision_model::TRT_tag>::value>::type>
    {
        static inline real_t get(const CollisionModel_T &lm, cell_idx_t, cell_idx_t, cell_idx_t) {
           return lm.lambda_e();
        }
    };

    template<typename CollisionModel_T>
    struct ShearRelaxationParameter<CollisionModel_T, typename std::enable_if<std::is_same<typename CollisionModel_T::tag,
                                                                                               collision_model::MRT_tag>::value>::type>
    {
        static inline real_t get(const CollisionModel_T &lm, cell_idx_t, cell_idx_t, cell_idx_t) {
           return lm.omega();
        }
    };

    template<typename CollisionModel_T>
    struct ShearRelaxationParameter<CollisionModel_T, typename std::enable_if<std::is_same<typename CollisionModel_T::tag,
                                                                                                 collision_model::Cumulant_tag>::value>::type>
    {
        static inline real_t get(const CollisionModel_T &cm, cell_idx_t, cell_idx_t, cell_idx_t) {
           return cm.omega1();
        }
    };

} // namespace internal


template< typename LatticeModel_T >
struct ShearRate
{
   using Stencil = typename LatticeModel_T::Stencil;
   using ShearRelaxationParameter = typename internal::ShearRelaxationParameter<typename LatticeModel_T::CollisionModel>;

   template< typename FieldPtrOrIterator >
   static inline real_t get( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                             const Vector3< real_t > & velocity, const real_t rho )
   {
      const auto equilibrium = EquilibriumDistribution< LatticeModel_T >::get( velocity, rho );

      std::vector< real_t > nonEquilibrium( Stencil::Size );
      for( uint_t i = 0; i != Stencil::Size; ++i )
         nonEquilibrium[i] = it[i] - equilibrium[i];

      const real_t relaxationParam = ShearRelaxationParameter::get(latticeModel.collisionModel(),it.x(),it.y(),it.z());
      return get( nonEquilibrium, relaxationParam, rho );
   }

   template< typename PdfField_T >
   static inline real_t get( const LatticeModel_T & latticeModel,
                             const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                             const Vector3< real_t > & velocity, const real_t rho )
   {
      const auto equilibrium = EquilibriumDistribution< LatticeModel_T >::get( velocity, rho );

      const auto & xyz0 = pdf(x,y,z,0);

      std::vector< real_t > nonEquilibrium( Stencil::Size );
      for( uint_t i = 0; i != Stencil::Size; ++i )
         nonEquilibrium[i] = pdf.getF( &xyz0, i ) - equilibrium[i];

      const real_t relaxationParam = ShearRelaxationParameter::get(latticeModel.collisionModel(),x,y,z);
      return get( nonEquilibrium, relaxationParam, rho );
   }

   /// For incompressible LB you don't have to pass a value for 'rho' since for incompressible LB 'rho' is not used in this function!
   static inline real_t get( const std::vector< real_t > & nonEquilibrium, const real_t relaxationParam, const real_t rho = real_t(1) )
   {
      WALBERLA_ASSERT_EQUAL( nonEquilibrium.size(), Stencil::Size );

      real_t D2 = real_t(0);

      if( LatticeModel_T::compressible )
      {
         for( uint_t alpha = 0; alpha < 3; ++alpha )
            for( uint_t beta = 0; beta < 3; ++beta )
            {
               real_t curStrain = real_t(0);
               for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
                  curStrain += nonEquilibrium[ d.toIdx() ] * stencil::c[alpha][*d] * stencil::c[beta][*d];

               curStrain *= ( -real_t(3) / ( real_t(2) * rho ) ) * relaxationParam;
               D2 += curStrain * curStrain;
            }
      }
      else // incompressible LB
      {
         for( uint_t alpha = 0; alpha < 3; ++alpha )
            for( uint_t beta = 0; beta < 3; ++beta )
            {
               real_t curStrain = real_t(0);
               for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
                  curStrain += nonEquilibrium[ d.toIdx() ] * stencil::c[alpha][*d] * stencil::c[beta][*d];

               curStrain *= ( -real_t(3) / real_t(2) ) * relaxationParam;
               D2 += curStrain * curStrain;
            }
      }

      return real_t(2) * std::sqrt( D2 );
   }
};


} // namespace lbm
} // namespace walberla
