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
//! \file Permeability.impl.h
//! \ingroup lbm
//! \author Tobias Schruff <schruff@iww.rwth-aachen.de>
//
//======================================================================================================================

namespace walberla {
namespace lbm {
namespace evaluations {


template< typename PdfField_T, typename BoundaryHandling_T >
Permeability<PdfField_T, BoundaryHandling_T>::Permeability( real_t viscosity, const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId, const FlagUID & fluid, const shared_ptr<blockforest::StructuredBlockForest> blocks )
   : nu_( viscosity ), pdfFieldId_( pdfFieldId ), boundaryHandlingId_( boundaryHandlingId ), fluid_( fluid ), blocks_( blocks ), time_( 0 ), lastU_( 0 )
{
   delta_               = std::numeric_limits<real_t>::max();
   k_                   = real_t(0);

   numSampleFluidNodes_ = uint_t(0);
   numC0FluidNodes_     = uint_t(0);
   numC1FluidNodes_     = uint_t(0); 
   flowAxis_            = uint_t(0);
   convCrit_            = real_t(0);
   interval_            = uint_t(0);

   initialized_  = false;
}


template< typename PdfField_T, typename BoundaryHandling_T >
void Permeability<PdfField_T, BoundaryHandling_T>::init( const Config::BlockHandle & config )
{
   sampleVolume_ = blocks_->getCellBBFromAABB( config.getParameter<AABB>( "sampleVolume", blocks_->getDomain() ) );
   flowAxis_     = config.getParameter<uint_t>( "flowAxis" );
   interval_     = config.getParameter<uint_t>( "calcFrequency" );
   convCrit_     = config.getParameter<real_t>( "convCriterion", real_t(1.0E-20) );

   initSampleVolume();
}


template< typename PdfField_T, typename BoundaryHandling_T >
void Permeability<PdfField_T, BoundaryHandling_T>::init( const AABB & sampleVolume, uint_t flowAxis, uint_t calcFrequency, real_t convCrit )
{
   sampleVolume_ = blocks_->getCellBBFromAABB( sampleVolume );
   flowAxis_     = flowAxis;
   interval_     = calcFrequency;
   convCrit_     = convCrit;

   initSampleVolume();
}


template< typename PdfField_T, typename BoundaryHandling_T >
void Permeability<PdfField_T, BoundaryHandling_T>::init( const CellInterval & sampleVolume, uint_t flowAxis, uint_t calcFrequency, real_t convCrit )
{
   sampleVolume_ = sampleVolume;
   flowAxis_     = flowAxis;
   interval_     = calcFrequency;
   convCrit_     = convCrit;

   initSampleVolume();
}


template< typename PdfField_T, typename BoundaryHandling_T >
void Permeability<PdfField_T, BoundaryHandling_T>::operator()()
{
   if( !initialized_ )
   {
      WALBERLA_LOG_WARNING( "You are trying to use an uninitialized Permeability class! Skipping evaluation ..." );
      return;
   }

   ++time_;
   if( ( time_ - 1u ) % interval_ != 0 || time_ <= 1 )
      return;

   // holds required local quantities
   std::vector<real_t> params( 6, 0.0 );
   auto & u   = params[0]; // flow velocity along flow axis
   auto & rho = params[1]; // fluid density
   auto & p0  = params[2]; // fluid pressure (density) in cross section 0
   auto & v0  = params[3]; // fluid velocity along flow axis in cross section 0
   auto & p1  = params[4]; // fluid pressure (density) in cross section 1
   auto & v1  = params[5]; // fluid velocity along flow axis in cross section 1

   // calculate required quantities on each block
   for( const auto &block : *blocks_ )
   {
      const auto * pdfField  = block.template getData<PdfField>( pdfFieldId_ );
      const auto * bHandling = block.template getData<BoundaryHandling>( boundaryHandlingId_ );
      const auto * flagField = bHandling->getFlagField();

      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
      WALBERLA_ASSERT_NOT_NULLPTR( bHandling );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      const flag_t domain = flagField->getFlag( fluid_ );

      // fluid velocity and density in sample volume
      CellInterval sampleVolume;
      blocks_->transformGlobalToBlockLocalCellInterval( sampleVolume, block, sampleVolume_ );
      sampleVolume.intersect( pdfField->xyzSize() );

      for( const auto &cell : sampleVolume )
      {
         if( flagField->isPartOfMaskSet( cell, domain ) )
         {
            u   += pdfField->getVelocity( cell )[flowAxis_];
            rho += pdfField->getDensity( cell );
         }
      }

      // fluid velocity and density in cross section c0
      CellInterval c0;
      blocks_->transformGlobalToBlockLocalCellInterval( c0, block, c0_ );
      c0.intersect( pdfField->xyzSize() );

      for( const auto &cell : c0 )
      {
         if( flagField->isPartOfMaskSet( cell, domain ) )
         {
            const real_t d = pdfField->getDensity( cell );
            p0 += d;
            v0 += pdfField->getVelocity( cell )[flowAxis_] * d;
         }
      }

      // fluid velocity and density in cross section c1
      CellInterval c1;
      blocks_->transformGlobalToBlockLocalCellInterval( c1, block, c1_ );
      c1.intersect( pdfField->xyzSize() );

      for( const auto &cell : c1 )
      {
         if( flagField->isPartOfMaskSet( cell, domain ) )
         {
            const real_t d = pdfField->getDensity( cell );
            p1 += d;
            v1 += pdfField->getVelocity( cell )[flowAxis_] * d;
         }
      }
   }

   // sum up parameters on root process
   mpi::reduceInplace( params, mpi::SUM );
 
   WALBERLA_ROOT_SECTION()
   {
     // calculate volume averaged quantities
     u   /= real_c( sampleVolume_.numCells() ); // volume averaged
     rho /= real_c( numSampleFluidNodes_ );     // average over fluid cells
     p0  /= real_c( numC0FluidNodes_ );         // average density in cross section c0
     p1  /= real_c( numC1FluidNodes_ );         // average density in cross section c1

     // convert density to pressure (P = rho / 3) and calculate average gradient
     const real_t pressureGradient = ( p0 - p1 ) / ( real_t(3) * real_c( sampleVolume_.size( flowAxis_ ) - 1 ) );

     if( math::isnan( u ) || math::isinf( u ) )
     {
        WALBERLA_LOG_WARNING( "Cannot determine permeability. Invalid mean fluid velocity " << u );

        delta_ = std::numeric_limits<real_t>::max();
        k_     = real_t(0);
     }

     else
     {
        // K = u * nu * rho / dP
        k_     = u * ( nu_ * rho ) / pressureGradient;
        // delta = dt / du
        delta_ = std::abs( u - lastU_ ) / real_c( interval_ );

        WALBERLA_LOG_DETAIL( "p0               = " << p0 );
        WALBERLA_LOG_DETAIL( "# c0 fluid nodes = " << numC0FluidNodes_ );
        WALBERLA_LOG_DETAIL( "mass flux Q0     = " << v0 );
        WALBERLA_LOG_DETAIL( "p1               = " << p1 );
        WALBERLA_LOG_DETAIL( "# c1 fluid nodes = " << numC1FluidNodes_ );
        WALBERLA_LOG_DETAIL( "mass flux Q1     = " << v1 );
        WALBERLA_LOG_DETAIL( "pressure grad.   = " << pressureGradient );
        WALBERLA_LOG_DETAIL( "mean density     = " << rho );
        WALBERLA_LOG_DETAIL( "mean velocity    = " << u );
        WALBERLA_LOG_DETAIL( "# sample nodes   = " << numSampleFluidNodes_ );
        WALBERLA_LOG_DETAIL( "nu               = " << nu_ );
        WALBERLA_LOG_DETAIL( "conv. crit.      = " << delta_ );

        WALBERLA_LOG_RESULT( "Permeability K   = " << k_ );

        lastU_ = u;

        if( hasConverged() )
        {
            WALBERLA_LOG_RESULT( "Permeability conv. criterion reached!" );
            WALBERLA_LOG_RESULT( "Delta          = " << currentDelta() );
            WALBERLA_LOG_RESULT( "Conv. crit.    = " << convCriterion() );
            WALBERLA_LOG_RESULT( "Rel. mass err. = " << std::abs( (v1 - v0) / v0 ) );
         }
      }
   }

   // broadcast delta and k to other processes to be accessible in member functions (currentDelta and currentValue)
   mpi::broadcastObject( delta_ );
   mpi::broadcastObject( k_     );
}


template< typename PdfField_T, typename BoundaryHandling_T >
void Permeability<PdfField_T, BoundaryHandling_T>::initSampleVolume()
{
   if( flowAxis_ == uint_t(0) ) // X AXIS
   {
      c0_ = CellInterval( sampleVolume_.xMin(), sampleVolume_.yMin(), sampleVolume_.zMin(),
                          sampleVolume_.xMin(), sampleVolume_.yMax(), sampleVolume_.zMax() );

      c1_ = CellInterval( sampleVolume_.xMax(), sampleVolume_.yMin(), sampleVolume_.zMin(),                           
                          sampleVolume_.xMax(), sampleVolume_.yMax(), sampleVolume_.zMax() );
   }
   else if( flowAxis_ == uint_t(1) ) // Y AXIS
   {
      c0_ = CellInterval( sampleVolume_.xMin(), sampleVolume_.yMin(), sampleVolume_.zMin(),
                          sampleVolume_.xMax(), sampleVolume_.yMin(), sampleVolume_.zMax() );

      c1_ = CellInterval( sampleVolume_.xMin(), sampleVolume_.yMax(), sampleVolume_.zMin(),
                          sampleVolume_.xMax(), sampleVolume_.yMax(), sampleVolume_.zMax() );
   }
   else if( flowAxis_ == uint_t(2) ) // Z AXIS
   {
      c0_ = CellInterval( sampleVolume_.xMin(), sampleVolume_.yMin(), sampleVolume_.zMin(),
                          sampleVolume_.xMax(), sampleVolume_.yMax(), sampleVolume_.zMin() );

      c1_ = CellInterval( sampleVolume_.xMin(), sampleVolume_.yMin(), sampleVolume_.zMax(),
                          sampleVolume_.xMax(), sampleVolume_.yMax(), sampleVolume_.zMax() );
   }
   else
      WALBERLA_ABORT( "Unrecognized value for flowAxis: " << flowAxis_ << ". Must be a value [0-2]! Aborting ..." );

   numSampleFluidNodes_ = 0;
   numC0FluidNodes_     = 0;
   numC1FluidNodes_     = 0;

   for( const auto &block : *blocks_ )
   {
      const auto * bHandling = block.template getData<BoundaryHandling>( boundaryHandlingId_ );
      const auto * flagField = bHandling->getFlagField();

      WALBERLA_ASSERT_NOT_NULLPTR( bHandling );
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      const flag_t domain = flagField->getFlag( fluid_ );

      CellInterval sampleVolume;
      blocks_->transformGlobalToBlockLocalCellInterval( sampleVolume, block, sampleVolume_ );
      sampleVolume.intersect( flagField->xyzSize() );

      CellInterval c0;
      blocks_->transformGlobalToBlockLocalCellInterval( c0, block, c0_ );
      c0.intersect( flagField->xyzSize() );

      CellInterval c1;
      blocks_->transformGlobalToBlockLocalCellInterval( c1, block, c1_ );
      c1.intersect( flagField->xyzSize() );

      for( const auto &cell : sampleVolume ) {
         if( flagField->isPartOfMaskSet( cell, domain ) )
            numSampleFluidNodes_++;
      }

      for( const auto &cell : c0 ) {
         if( flagField->isPartOfMaskSet( cell, domain ) )
            numC0FluidNodes_++;
      }

      for( const auto &cell : c1 ) {
         if( flagField->isPartOfMaskSet( cell, domain ) )
            numC1FluidNodes_++;
      }
   }

   mpi::allReduceInplace( numSampleFluidNodes_, mpi::SUM );
   mpi::allReduceInplace( numC0FluidNodes_    , mpi::SUM );
   mpi::allReduceInplace( numC1FluidNodes_    , mpi::SUM );

   initialized_ = true;
}


} // namespace evaluations
} // namespace lbm
} // namespace walberla