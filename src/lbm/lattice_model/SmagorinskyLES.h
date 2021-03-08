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
//! \file SmagorinskyLES.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "CollisionModel.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/iterators/IteratorMacros.h"
#include "field/EvaluationFilter.h"
#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"

namespace walberla {
namespace lbm {



/// Can be used together with the collision model "SRTField"
template< typename LatticeModel_T, typename Filter_T = field::DefaultEvaluationFilter >
class SmagorinskyLES
{
public:

   using PdfField_T = lbm::PdfField<LatticeModel_T>;
   using ScalarField_T = GhostLayerField<real_t, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   SmagorinskyLES( weak_ptr< StructuredBlockStorage > blocks, const ConstBlockDataID & pdfFieldId, const BlockDataID & omegaFieldId,
                   const real_t & kinematicViscosity, const real_t & smagorinskyConstant,
                   const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )

   : blocks_( blocks ), filter_( Filter_T() ), pdfFieldId_( pdfFieldId ), omegaFieldId_( omegaFieldId ),
     kinematicViscosity_( kinematicViscosity ), smagorinskyConstant_ ( smagorinskyConstant ),
     requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {}

   SmagorinskyLES( weak_ptr< StructuredBlockStorage > blocks, const Filter_T & filter,
                   const ConstBlockDataID & pdfFieldId, const BlockDataID & omegaFieldId,
                   const real_t & kinematicViscosity, const real_t & smagorinskyConstant,
                   const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                   const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )

   : blocks_( blocks ), filter_( filter ), pdfFieldId_( pdfFieldId ), omegaFieldId_( omegaFieldId ),
     kinematicViscosity_( kinematicViscosity ), smagorinskyConstant_ ( smagorinskyConstant ),
     requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {}

   // for LBM refinement time step callback (post stream function)
   void operator()( IBlock * block, const uint_t level, const uint_t = uint_t(0) );

   // can be used as a pre collide or post stream function for standard LBM without refinement
   void operator()()
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks );

      for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
         this->operator()( block.get(), blocks->getLevel(*block) );
   }

private:

   weak_ptr< StructuredBlockStorage > blocks_;

   Filter_T filter_;

   ConstBlockDataID pdfFieldId_;
   BlockDataID omegaFieldId_;

   real_t kinematicViscosity_;
   real_t smagorinskyConstant_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class SmagorinskyLES



template< typename LatticeModel_T, typename Filter_T >
void SmagorinskyLES< LatticeModel_T, Filter_T >::operator()( IBlock * block, const uint_t level, const uint_t )
{
   const PdfField_T *   pdfField    = block->getData< PdfField_T >( pdfFieldId_ );
         ScalarField_T * omegaField = block->getData< ScalarField_T >( omegaFieldId_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( omegaField );

   WALBERLA_ASSERT_EQUAL( pdfField->xyzSizeWithGhostLayer(), omegaField->xyzSizeWithGhostLayer() );

   {
      auto blocks = blocks_.lock(); // 'blocks' will be released at the end of this scope
      WALBERLA_CHECK_NOT_NULLPTR( blocks );
      WALBERLA_ASSERT_EQUAL( level, blocks->getLevel(*block) );
   }

   const real_t tau = real_t(1) / lbm::collision_model::levelDependentRelaxationParameter( level,
                                                                                           lbm::collision_model::omegaFromViscosity( kinematicViscosity_ ), uint_t(0) );

   const real_t factor = real_t(18) * std::sqrt( real_t(2) ) * smagorinskyConstant_ * smagorinskyConstant_;

   filter_( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( omegaField,

      if( filter_(x,y,z) )
      {
         Vector3< real_t > velocity;
         const real_t rho = pdfField->getDensityAndEquilibriumVelocity( velocity, x, y, z );

         const auto equilibrium = lbm::EquilibriumDistribution< LatticeModel_T >::get( velocity, rho );

         std::vector< real_t > nonEquilibrium( Stencil_T::Size );
         for( uint_t i = 0; i != Stencil_T::Size; ++i )
            nonEquilibrium[i] = pdfField->get(x,y,z,i) - equilibrium[i];

         real_t filteredMeanMomentum = real_t(0);
         for( uint_t alpha = 0; alpha < 3; ++alpha ) {
            for( uint_t beta = 0; beta < 3; ++beta )
            {
               real_t qij = real_t(0);
               for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
                  qij += nonEquilibrium[ d.toIdx() ] * real_c(stencil::c[alpha][*d]) * real_c(stencil::c[beta][*d]);
               filteredMeanMomentum += qij * qij;
            }
         }

         const real_t tauTurbulent = LatticeModel_T::compressible ? ( real_t(0.5) * ( std::sqrt( tau * tau + (factor / rho) * std::sqrt( real_t(2) * filteredMeanMomentum ) ) - tau ) ) :
                                                                    ( real_t(0.5) * ( std::sqrt( tau * tau + factor * std::sqrt( real_t(2) * filteredMeanMomentum ) ) - tau ) );

         omegaField->get(x,y,z) = real_t(1) / ( tau + tauTurbulent );
      }
   )
}



} // namespace lbm
} // namespace walberla
