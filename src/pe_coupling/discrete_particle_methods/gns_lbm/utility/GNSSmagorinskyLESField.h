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
//! \file GNSSmagorinskyLESField.h
//! \ingroup lbm
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/iterators/IteratorMacros.h"
#include "field/EvaluationFilter.h"
#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Adjusts locally the fluid viscosity based on a LES-Smagorinsky turbulence model when using GNS-LBM
 *
 * The original LES model is discussed in:
 * - S. Hou, J. Sterling, S. Chen, G. Doolen, A lattice Boltzmann subgrid model for high reynolds number flows, in: A. T.
 *   Lawniczak, R. Kapral (Eds.), Pattern formation and lattice gas automata, Vol. 6 of Fields Institute Communications,
 *   American Mathematical Society, Providence, RI, 1996, p. 151–166.
 * - H. Yu, S. S. Girimaji, L.-S. Luo, DNS and LES of decaying isotropic turbulence with and without frame rotation using
 *   lattice Boltzmann method, Journal of Computational Physics 209 (2) (2005) 599–616. doi:10.1016/j.jcp.2005.03.022.
 *
 * The second reference suggests to use a Smagorinsky constant of 0.1.
 *
 * An implementation of if for standard LBM is found in "lbm/lattice_model/SmagorinskyLES.h"
 *
 * This version is adapted to use the special GNS-LBM equilibrium distribution functions to calculate the
 * non-equilibrium PDFs. These are given in Z. Guo, T. S. Zhao - "Lattice Boltzmann model for incompressible flows
 * through porous media", Phys. Rev. E 66 (2002)036304. doi:10.1103/PhysRevE.66.036304.
 *
 * Combining LES with GNS-LBM was proposed in
 * - C. Rettinger, U. Ruede - "A Coupled Lattice Boltzmann Method and Discrete Element Method for Discrete Particle
 *   Simulations of Particulate Flows". arXiv preprint arXiv:1711.00336 (2017)
 *
 * Uses the value from the omegaField as tau0 and adds a calculated tauTurbulent to it.
 * Can be used together with the collision model "SRTField"
 */
template< typename LatticeModel_T, typename Filter_T = field::DefaultEvaluationFilter >
class GNSSmagorinskyLESField
{
public:

   using PdfField_T = lbm::PdfField<LatticeModel_T>;
   using ScalarField_T = GhostLayerField<real_t, 1>;
   using Stencil_T = typename LatticeModel_T::Stencil;

   static_assert( LatticeModel_T::CollisionModel::constant == false, "Only works with non-constant relaxation time fields!" );
   static_assert( LatticeModel_T::compressible == false,             "Only works with incompressible models!" );

   GNSSmagorinskyLESField( const shared_ptr< StructuredBlockStorage > & blocks,
                           const BlockDataID & omegaFieldID, const ConstBlockDataID & pdfFieldID,
                           const ConstBlockDataID & solidVolumeFractionFieldID,
                           const real_t & smagorinskyConstant,
                           const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                           const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )

   : blocks_( blocks ), filter_( Filter_T() ), omegaFieldID_( omegaFieldID ), pdfFieldID_( pdfFieldID ),
     solidVolumeFractionFieldID_( solidVolumeFractionFieldID ), smagorinskyConstant_ ( smagorinskyConstant ),
     requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {}

   GNSSmagorinskyLESField( const shared_ptr< StructuredBlockStorage > & blocks, const Filter_T & filter,
                           const BlockDataID & omegaFieldID, const ConstBlockDataID & pdfFieldID,
                           const ConstBlockDataID & solidVolumeFractionFieldID,
                           const real_t & smagorinskyConstant,
                           const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                           const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )

   : blocks_( blocks ), filter_( filter ), omegaFieldID_( omegaFieldID ), pdfFieldID_( pdfFieldID ),
     solidVolumeFractionFieldID_( solidVolumeFractionFieldID ), smagorinskyConstant_ ( smagorinskyConstant ),
     requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {}

   void operator()( IBlock * const block );

private:

   shared_ptr< StructuredBlockStorage > blocks_;

   Filter_T filter_;

   BlockDataID omegaFieldID_;
   ConstBlockDataID pdfFieldID_;
   ConstBlockDataID solidVolumeFractionFieldID_;

   real_t smagorinskyConstant_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class GNSSmagorinskyLESField



template< typename LatticeModel_T, typename Filter_T >
void GNSSmagorinskyLESField< LatticeModel_T, Filter_T >::operator()( IBlock * const block )
{
         ScalarField_T * omegaField = block->getData< ScalarField_T >( omegaFieldID_ );
   const PdfField_T *   pdfField    = block->getData< PdfField_T >( pdfFieldID_ );
   const ScalarField_T * svfField   = block->getData< ScalarField_T >( solidVolumeFractionFieldID_ );

   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( omegaField );


   WALBERLA_ASSERT_EQUAL( pdfField->xyzSizeWithGhostLayer(), omegaField->xyzSizeWithGhostLayer() );

   const real_t factor = real_t(18) * std::sqrt( real_t(2) ) * smagorinskyConstant_ * smagorinskyConstant_;

   filter_( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( omegaField,

      if( filter_(x,y,z) )
      {
         const real_t tau0 = real_t(1) / omegaField->get(x,y,z);

         Vector3< real_t > momentumDensity;
         const real_t rho = pdfField->getDensityAndEquilibriumMomentumDensity( momentumDensity, x, y, z );

         Vector3< real_t > velocity = momentumDensity / rho;

         const real_t porosity = real_t(1) - svfField->get(x,y,z);
         const real_t invPorosity = real_t(1) / porosity;

         // use special GNS equilibrium
         std::vector< real_t > equilibrium( Stencil_T::Size );
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const real_t wq = LatticeModel_T::w[ d.toIdx() ];
            const real_t vel = real_c(d.cx()) * velocity[0] + real_c(d.cy()) * velocity[1] + real_c(d.cz()) * velocity[2];
            real_t feq = wq * rho * ( real_t(1) - real_t(1.5) * invPorosity * velocity.sqrLength() +
                                      real_t(4.5) * invPorosity * vel * vel + real_t(3.0) * vel ); // modified feq
            feq -= wq; // walberla: center around 0 in incompressible case, attention if always needed!
            equilibrium[d.toIdx()] = feq;
         }

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

         const real_t tauTurbulent = LatticeModel_T::compressible ? ( real_t(0.5) * ( std::sqrt( tau0 * tau0 + (factor / rho) * std::sqrt( real_t(2) * filteredMeanMomentum ) ) - tau0 ) ) :
                                                                    ( real_t(0.5) * ( std::sqrt( tau0 * tau0 + factor * std::sqrt( real_t(2) * filteredMeanMomentum ) ) - tau0 ) );

         omegaField->get(x,y,z) = real_t(1) / ( tau0 + tauTurbulent );
      }
   )
}



} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
