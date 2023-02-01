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
//! \file GNSSweep.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"

#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"
#include "lbm/lattice_model/all.h"

#include "pe/Types.h"

#include "stencil/Directions.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {


/*!\brief Sweep for the LBM formulation of the generalized Navier Stokes equations
 *
 * Method is taken from:
 * Z. Guo, T. S. Zhao - "Lattice Boltzmann model for incompressible flows through porous media", Phys. Rev. E 66 (2002)036304. doi:10.1103/PhysRevE.66.036304.
 * Note: when an external forcing is present, use the GNSExternalForceToForceFieldAdder class to apply the correct external force
 * (i.e. fluid-volume-fraction * external-force )
 * Note: the calculated velocity is the volume averaged one!
 *
 * A LBM-forcing model that supports spatially (and temporally) varying forces has to be used, e.g. GuoField.
 *
 */
template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
class GNSSweep
   : public lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >
{
public:

   using PdfField_T = typename lbm::SweepBase<LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T>::PdfField_T;
   using Stencil_T = typename LatticeModel_T::Stencil;
   using ScalarField_T = GhostLayerField<real_t, 1>;

   static_assert( LatticeModel_T::ForceModel::constant == false, "Only works with non-constant force models!" );
   static_assert( LatticeModel_T::compressible == false,         "Only works with incompressible models!" );

   GNSSweep( const BlockDataID & pdfFieldID,
             const BlockDataID & solidVolumeFractionFieldID,
             const Filter_T & _filter = walberla::field::DefaultEvaluationFilter(),
             const DensityVelocityIn_T & _densityVelocityIn = lbm::DefaultDensityEquilibriumVelocityCalculation(),
             const DensityVelocityOut_T & _densityVelocityOut = lbm::DefaultDensityVelocityCallback() ) :
      lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >( pdfFieldID, _filter, _densityVelocityIn, _densityVelocityOut ),
      solidVolumeFractionFieldID_( solidVolumeFractionFieldID ) {}

   GNSSweep( const BlockDataID & src, const BlockDataID & dst,
             const BlockDataID & solidVolumeFractionFieldID,
             const Filter_T & _filter = walberla::field::DefaultEvaluationFilter(),
             const DensityVelocityIn_T & _densityVelocityIn = lbm::DefaultDensityEquilibriumVelocityCalculation(),
             const DensityVelocityOut_T & _densityVelocityOut = lbm::DefaultDensityVelocityCallback() ) :
      lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >( src, dst, _filter, _densityVelocityIn, _densityVelocityOut ),
      solidVolumeFractionFieldID_( solidVolumeFractionFieldID ) {}

   void operator()( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) )
   {
      streamCollide( block, numberOfGhostLayersToInclude );
   }

   void streamCollide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );
   void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) );

   inline ScalarField_T * getSolidVolumeFractionField( IBlock * const block ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      ScalarField_T * solidVolumeFractionField = block->getData<ScalarField_T>( solidVolumeFractionFieldID_ );
      WALBERLA_ASSERT_NOT_NULLPTR( solidVolumeFractionField );
      return solidVolumeFractionField;
   }

   inline void getFields( IBlock * const block, PdfField_T * & src, PdfField_T * & dst, ScalarField_T * & solidVolumeFractionField )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );

      src = this->getSrcField( block );
      dst = this->getDstField( block, src );
      solidVolumeFractionField = getSolidVolumeFractionField( block );

      WALBERLA_ASSERT_NOT_NULLPTR( src );
      WALBERLA_ASSERT_NOT_NULLPTR( dst );
      WALBERLA_ASSERT_NOT_NULLPTR( solidVolumeFractionField );

      WALBERLA_ASSERT_EQUAL( src->xyzSize(), dst->xyzSize() );
      WALBERLA_ASSERT_EQUAL( src->xyzSize(), solidVolumeFractionField->xyzSize() );
   }

private:
   const BlockDataID solidVolumeFractionFieldID_;
};


template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
void GNSSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T
>::streamCollide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src( NULL );
   PdfField_T * dst( NULL );
   ScalarField_T * solidVolumeFractionField( nullptr );

   getFields( block, src, dst, solidVolumeFractionField );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFractionField->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   const auto & lm = src->latticeModel();
   dst->resetLatticeModel( lm ); /* required so that member functions for getting density and equilibrium velocity can be called for dst! */

   this->filter( *block );
   this->densityVelocityIn( *block );
   this->densityVelocityOut( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,

      if( this->filter(x,y,z) )
      {
         const real_t fluidVolumeFraction = real_t(1) - solidVolumeFractionField->get(x,y,z); //i.e. porosity
         const real_t invFluidVolumeFraction = real_t(1)/fluidVolumeFraction;

         // stream pull
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
            dst->get( x, y, z, d.toIdx() ) = src->get( x-d.cx(), y-d.cy(), z-d.cz(), d.toIdx() );

         Vector3<real_t> velocity;
         real_t rho = this->densityVelocityIn( velocity, dst, x, y, z );

         this->densityVelocityOut( x, y, z, lm, velocity, rho );

         velocity /= rho;

         const real_t omega = lm.collisionModel().omega( x, y, z, velocity, rho );

         const Vector3<real_t> extForce = lm.forceModel().forceDensity(x,y,z);

         // collide
         for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
         {
            const real_t wq = LatticeModel_T::w[ d.toIdx() ];
            const Vector3<real_t> c( real_c(d.cx()), real_c(d.cy()), real_c(d.cz()) );
            const real_t forceTerm = real_t(3.0) * wq * ( real_t(1) - real_t(0.5) * omega ) *
                     ( ( c - invFluidVolumeFraction * velocity + ( real_t(3) * invFluidVolumeFraction * ( c * velocity ) * c ) ) * extForce ); // modified Guo forcing, multiplication by rho is wrong because of units

            const real_t vel = c * velocity;
            real_t feq = wq * rho * ( real_t(1) - real_t(1.5) * invFluidVolumeFraction * velocity.sqrLength() +
                                      real_t(4.5) * invFluidVolumeFraction * vel * vel + real_t(3.0) * vel ); // modified feq
            feq -= wq; // center PDFs around 0

            dst->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega ) * dst->get( x, y, z, d.toIdx() ) +
                                             omega * feq +
                                             forceTerm;
         }
      }

   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

   src->swapDataPointers( dst );

}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
void GNSSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T
>::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src( nullptr );
   PdfField_T * dst( nullptr );
   lbm::SweepBase<LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T>::getFields( block, src, dst );
   lbm::StreamPull< LatticeModel_T >::execute( src, dst, block, this->filter_, numberOfGhostLayersToInclude );
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
void GNSSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T
>::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude )
{
   PdfField_T * src = this->getSrcField( block );
   ScalarField_T * solidVolumeFractionField = getSolidVolumeFractionField( block );

   WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude );
   WALBERLA_ASSERT_GREATER_EQUAL( solidVolumeFractionField->nrOfGhostLayers(), numberOfGhostLayersToInclude );

   const auto & lm = src->latticeModel();

   this->filter( *block );
   this->densityVelocityIn( *block );
   this->densityVelocityOut( *block );

   WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ( src, numberOfGhostLayersToInclude,
   if( this->filter(x,y,z) )
   {
      const real_t fluidVolumeFraction = real_t(1) - solidVolumeFractionField->get(x,y,z); //i.e. porosity
      const real_t invFluidVolumeFraction = real_t(1)/fluidVolumeFraction;

      Vector3<real_t> velocity;
      real_t rho = this->densityVelocityIn( velocity, src, x, y, z );

      this->densityVelocityOut( x, y, z, lm, velocity, rho );

      velocity /= rho;

      const real_t omega = lm.collisionModel().omega( x, y, z, velocity, rho );

      const Vector3<real_t> extForce = lm.forceModel().forceDensity(x,y,z);

      // collide
      for( auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d )
      {
         const real_t wq = LatticeModel_T::w[ d.toIdx() ];
         const Vector3<real_t> c( real_c(d.cx()), real_c(d.cy()), real_c(d.cz()) );
         const real_t forceTerm = real_t(3.0) * wq * ( real_t(1) - real_t(0.5) * omega ) *
                  ( ( c - invFluidVolumeFraction * velocity + ( real_t(3) * invFluidVolumeFraction * ( c * velocity ) * c ) ) * extForce ); // modified Guo forcing

         const real_t vel = c * velocity;
         real_t feq = wq * rho * ( real_t(1) - real_t(1.5) * invFluidVolumeFraction * velocity.sqrLength() +
                                   real_t(4.5) * invFluidVolumeFraction * vel * vel + real_t(3.0) * vel ); // modified feq
         feq -= wq; // center PDFs around 0

         src->get( x, y, z, d.toIdx() ) = ( real_t(1.0) - omega ) * src->get( x, y, z, d.toIdx() ) +
                                          omega * feq +
                                          forceTerm;
      }
   }
   ) // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ

}

/////////////////////////////
// makeGNSSweep FUNCTIONS //
////////////////////////////

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
shared_ptr< GNSSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T > >
makeGNSSweep( const BlockDataID & pdfFieldID, const BlockDataID & solidVolumeFractionFieldID, const Filter_T & filter,
               const DensityVelocityIn_T & densityVelocityIn, const DensityVelocityOut_T & densityVelocityOut )
{
   using Sweep_T = GNSSweep<LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T>;
   return shared_ptr< Sweep_T >( new Sweep_T( pdfFieldID, solidVolumeFractionFieldID, filter, densityVelocityIn, densityVelocityOut ) );
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
shared_ptr< GNSSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T > >
makeGNSSweep( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & solidVolumeFractionFieldID, const Filter_T & filter,
               const DensityVelocityIn_T & densityVelocityIn, const DensityVelocityOut_T & densityVelocityOut )
{
   using Sweep_T = GNSSweep<LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T>;
   return shared_ptr< Sweep_T >( new Sweep_T( srcID, dstID, solidVolumeFractionFieldID, filter, densityVelocityIn, densityVelocityOut ) );
}

// only block data IDs of PDF data

template< typename LatticeModel_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter, lbm::DefaultDensityEquilibriumVelocityCalculation,
                       lbm::DefaultDensityVelocityCallback > >
makeGNSSweep( const BlockDataID & pdfFieldID, const BlockDataID & solidVolumeFractionFieldID )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::DefaultDensityVelocityCallback >
                       ( pdfFieldID, solidVolumeFractionFieldID,walberla::field::DefaultEvaluationFilter(),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback() );
}

template< typename LatticeModel_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter, lbm::DefaultDensityEquilibriumVelocityCalculation,
                       lbm::DefaultDensityVelocityCallback > >
makeGNSSweep( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & solidVolumeFractionFieldID )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::DefaultDensityVelocityCallback >
                       ( srcID, dstID, solidVolumeFractionFieldID, walberla::field::DefaultEvaluationFilter(),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback() );
}

// block data IDs of PDF data + flag field as filter

template< typename LatticeModel_T, typename FlagField_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                       lbm::DefaultDensityVelocityCallback > >
makeGNSSweep( const BlockDataID & pdfFieldID, const BlockDataID & solidVolumeFractionFieldID,
               const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToEvaluate )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::DefaultDensityVelocityCallback >
                       ( pdfFieldID, solidVolumeFractionFieldID, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldID, cellsToEvaluate ),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback() );
}

template< typename LatticeModel_T, typename FlagField_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                       lbm::DefaultDensityVelocityCallback > >
makeGNSSweep( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & solidVolumeFractionFieldID,
               const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToEvaluate )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::DefaultDensityVelocityCallback >
                       ( srcID, dstID, solidVolumeFractionFieldID, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldID, cellsToEvaluate ),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback() );
}

// block data IDs of PDF data + flag field as filter + block data ID of velocity field (out)

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
                       lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::VelocityCallback<VelocityField_T> > >
makeGNSSweep( const BlockDataID & pdfFieldID, const BlockDataID & solidVolumeFractionFieldID,
               const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToEvaluate, const BlockDataID & velocityFieldID )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::VelocityCallback<VelocityField_T> >
                       ( pdfFieldID, solidVolumeFractionFieldID, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldID, cellsToEvaluate ),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::VelocityCallback<VelocityField_T>( velocityFieldID ) );
}

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
                       lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::VelocityCallback<VelocityField_T> > >
makeGNSSweep( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & solidVolumeFractionFieldID,
               const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToEvaluate, const BlockDataID & velocityFieldID )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::VelocityCallback<VelocityField_T> >
                       ( srcID, dstID, solidVolumeFractionFieldID, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldID, cellsToEvaluate ),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::VelocityCallback<VelocityField_T>( velocityFieldID ) );
}

// block data IDs of PDF data + flag field as filter + block data IDs of velocity and density field (out)

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, typename DensityField_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
                       lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DensityVelocityCallback<VelocityField_T,DensityField_T> > >
makeGNSSweep( const BlockDataID & pdfFieldID, const BlockDataID & solidVolumeFractionFieldID,
               const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToEvaluate,
               const BlockDataID & velocityFieldID, const BlockDataID & densityFieldID )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::DensityVelocityCallback<VelocityField_T,DensityField_T> >
                       ( pdfFieldID, solidVolumeFractionFieldID, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldID, cellsToEvaluate ),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(),
                         lbm::DensityVelocityCallback<VelocityField_T,DensityField_T>( velocityFieldID, densityFieldID ) );
}

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, typename DensityField_T >
shared_ptr< GNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
                       lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DensityVelocityCallback<VelocityField_T,DensityField_T> > >
makeGNSSweep( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & solidVolumeFractionFieldID,
               const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToEvaluate,
               const BlockDataID & velocityFieldID, const BlockDataID & densityFieldID )
{
   return makeGNSSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>, lbm::DefaultDensityEquilibriumVelocityCalculation,
                         lbm::DensityVelocityCallback<VelocityField_T,DensityField_T> >
                       ( srcID, dstID, solidVolumeFractionFieldID, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldID, cellsToEvaluate ),
                         lbm::DefaultDensityEquilibriumVelocityCalculation(),
                         lbm::DensityVelocityCallback<VelocityField_T,DensityField_T>( velocityFieldID, densityFieldID ) );
}

} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
