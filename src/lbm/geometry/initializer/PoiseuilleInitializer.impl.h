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
//! \file PoiseuilleInitializer.impl.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "PoiseuilleInitializer.h"

#include "core/Abort.h"
#include "geometry/initializer/Initializer.h"
#include "geometry/initializer/BoundaryFromDomainBorder.h"
#include "lbm/boundary/SimplePressure.h"
#include "lbm/boundary/UBB.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/ForceModel.h"


namespace walberla {
namespace lbm {
namespace initializer {


//===================================================================================================================
//
//  Helper function
//
//===================================================================================================================

template< typename BH_T, typename LM, typename SP, typename UBB > const uint_t Poiseuille<BH_T,LM,SP,UBB>::X_AXIS = 0;
template< typename BH_T, typename LM, typename SP, typename UBB > const uint_t Poiseuille<BH_T,LM,SP,UBB>::Y_AXIS = 1;
template< typename BH_T, typename LM, typename SP, typename UBB > const uint_t Poiseuille<BH_T,LM,SP,UBB>::Z_AXIS = 2;
template< typename BH_T, typename LM, typename SP, typename UBB > const uint_t Poiseuille<BH_T,LM,SP,UBB>::INVALID_AXIS = 3;



template< typename BH_T, typename LM, typename SP, typename UBB >
Poiseuille<BH_T,LM,SP,UBB>::Poiseuille( StructuredBlockStorage & blocks, BlockDataID & handlerID, BlockDataID  & pdfFieldID,
                                 field::FlagUID noSlipFlag, field::FlagUID ubbFlag,
                                 field::FlagUID pressureFlag1, field::FlagUID pressureFlag2 )
   : storage_( blocks ), handlerID_ ( handlerID ), pdfFieldID_( pdfFieldID ),
     noSlipFlag_( noSlipFlag ), ubbFlag_( ubbFlag ), pressureFlag1_( pressureFlag1 ), pressureFlag2_ ( pressureFlag2 )
{
   // No blocks here
   if ( storage_.begin() == storage_.end() )
      return;

   IBlock& firstBlock = * (storage_.begin() );
   lbm::PdfField<LM> * pdfField = firstBlock.getData< lbm::PdfField<LM> > ( pdfFieldID_ );
   latticeViscosity_ = pdfField->latticeModel().collisionModel().viscosity();

   maxPoint_ = Vector3<real_t> ( real_c( storage_.getNumberOfXCells() ) ,
                                 real_c( storage_.getNumberOfYCells() ) ,
                                 real_c( storage_.getNumberOfZCells() )  );

   midPoint_= maxPoint_ * real_c(0.5);
}


template< typename BH_T, typename LM, typename SP, typename UBB >
void Poiseuille<BH_T,LM,SP,UBB>::init( const Config::BlockHandle & blockHandle )
{
   // Parsing Scenario
   Scenario scenario;
   std::string scenarioStr = blockHandle.getParameter<std::string>( "scenario" );

   if      ( scenarioStr == "rect2D" ) scenario = RECT_2D;
   else if ( scenarioStr == "pipe"   ) scenario = PIPE;
   else {
      WALBERLA_ABORT("Error parsing Geometry block PoiseuilleInit: Unknown scenario " << scenarioStr
                      << "\nValid Scenarios are rect2D,rect3D and pipe" );
   }

   // Parsing boundary parameters
   BoundaryType boundarySetup;
   std::string boundaryStr  = blockHandle.getParameter<std::string>( "boundary" );

   if      ( boundaryStr == "pressureDriven" )  boundarySetup = PRESSURE_DRIVEN;
   else if ( boundaryStr == "velocityDriven" )  boundarySetup = VELOCITY_DRIVEN;
   else if ( boundaryStr == "forceDriven"    )  boundarySetup = FORCE_DRIVEN;
   else  {
      WALBERLA_ABORT("Error parsing Geometry block PoiseuilleInit: Unknown inflow_boundary: " << boundaryStr
                       << "\nValid values: pressureDriven, velocityDriven and forceDriven" );
   }

   uint_t flowAxisInput = blockHandle.getParameter<uint_t>("flowAxis");
   if (  flowAxisInput > 2 )
      WALBERLA_ABORT("Error parsing Geometry block PoiseuilleInit: flowAxis has to be either 0,1 or 2 ( for x,y and z )");
   Axis flowAxis = Axis( flowAxisInput );

   uint_t parabolaAxisRead = blockHandle.getParameter<uint_t>( "parabolaAxis", INVALID_AXIS );
   if ( parabolaAxisRead != X_AXIS && parabolaAxisRead != Y_AXIS && parabolaAxisRead != Z_AXIS )
      parabolaAxisRead = INVALID_AXIS;
   Axis parabolaAxis = Axis( parabolaAxisRead );
   if ( scenario == RECT_2D )
      parabolaAxis = checkParabolaAxis( parabolaAxis, flowAxis );


   real_t pressureDiff = 0;
   if ( blockHandle.isDefined( "pressureDiff") )
      pressureDiff = blockHandle.getParameter<real_t>("pressureDiff");
   else
      if ( blockHandle.isDefined ("velocity") )
         pressureDiff = pressureDiffFromVelocity( scenario, blockHandle.getParameter<real_t>("velocity"), flowAxis, parabolaAxis );
      else
         WALBERLA_ABORT("Error parsing Geometry block PoiseuilleInit: Neither pressureDiff nor velocity specified.");


   init ( scenario, boundarySetup, pressureDiff, flowAxis, parabolaAxis );
}



template< typename BH_T, typename LM, typename SP, typename UBB >
void Poiseuille<BH_T,LM,SP,UBB>::init( Scenario scenario, BoundaryType boundaryType, real_t pressureDiff, Axis flowAxis, Axis parabolaAxis )
{
   if ( scenario == RECT_2D )
      parabolaAxis = checkParabolaAxis( parabolaAxis, flowAxis );

   if ( scenario == PIPE ) {
      if ( storage_.isPeriodic( (flowAxis + 1 ) % 3 ) ||
           storage_.isPeriodic( (flowAxis + 2 ) % 3) )  {
         WALBERLA_ABORT("Error in Channel Setup:: In pipe scenario the axis normal to flowAxis must not be periodic.")
      }
   }

   // check that inflow-outflow axis is not periodic
   if ( (boundaryType != FORCE_DRIVEN) && storage_.isPeriodic( flowAxis ) ) {
      WALBERLA_ABORT("Error in Channel Setup:: Flow axis was marked as periodic and boundary is not set to forceDriven" );
   }
   if ( boundaryType == FORCE_DRIVEN && ! storage_.isPeriodic( flowAxis ) ) {
      WALBERLA_ABORT("Error in Channel Setup:: Flow axis was not marked as periodic and boundary is set to forceDriven" );
   }

   real_t velocity = velocityFromPressureDiff( scenario, pressureDiff, flowAxis, parabolaAxis  );

   switch ( boundaryType ) {
      case PRESSURE_DRIVEN:
         initPressureBoundary( flowAxis, pressureDiff );
         break;
      case VELOCITY_DRIVEN:
         initVelocityBoundary( scenario, flowAxis, velocity, parabolaAxis );
         break;
      case FORCE_DRIVEN:
         {
            Vector3<real_t> forceVec ( 0,0,0 );
            forceVec[ flowAxis ] = pressureDiff / real_c( storage_.getNumberOfCells( flowAxis ) );
            for( auto blockIt = storage_.begin(); blockIt != storage_.end(); ++blockIt ) {
               auto pdfField = blockIt->template getData< lbm::PdfField<LM> > ( pdfFieldID_ );
               typename LM::ForceModel & forceModel = pdfField->latticeModel().forceModel();
               bool res = forceModel.setConstantBodyForceIfPossible( forceVec );
               if ( ! res )
                  WALBERLA_ABORT( "Error in Channel Setup:: Trying to set up a force driven channel, " <<
                                  " but lattice model does not support force terms.")
            }
         }
         break;
      default:
         WALBERLA_ASSERT(false);
   }

   initNoSlipBoundaries( scenario, flowAxis, parabolaAxis );
   setInitVelocity( scenario, flowAxis, velocity, parabolaAxis );
}


template< typename BH_T, typename LM, typename SP, typename UBB >
typename Poiseuille<BH_T,LM,SP,UBB>::Axis
Poiseuille<BH_T,LM,SP,UBB>::checkParabolaAxis( Axis parabolaAxis, Axis flowAxis )
{
   if ( parabolaAxis != INVALID_AXIS )
   {
      //try to find a default parabola axis
      uint_t firstOption  = ( flowAxis + 1 ) % 3;
      uint_t secondOption = ( flowAxis + 2 ) % 3;

      if ( ! storage_.isPeriodic( firstOption ) )
         return firstOption;
      if ( ! storage_.isPeriodic( secondOption ) )
         return secondOption;

      WALBERLA_ABORT( "Error in Channel Setup:: No parabolaAxis given. \n" <<
                      "No axis (different from flowAxis) found which is not periodic and could be used as parabolaAxis" );
   }

   if ( parabolaAxis > 2 )  {
      WALBERLA_ABORT("Error in Channel Setup:: Invalid parabolaAxis: has to be between 0 and 2 ");
   }
   if( storage_.isPeriodic( parabolaAxis ) ) {
      WALBERLA_ABORT("Error in Channel Setup:: parabolaAxis is not allowed to be periodic");
   }
   return parabolaAxis;
}


template< typename BH_T, typename LM, typename SP, typename UBB >
void Poiseuille<BH_T,LM,SP,UBB>::initVelocityBoundary( Scenario scenario, Axis flowAxis, real_t maxVelocity, Axis parabolaAxis )
{
   for( auto blockIt = storage_.begin(); blockIt != storage_.end(); ++blockIt )
   {
      auto handling = blockIt->template getData<BH_T> ( handlerID_ );
      auto pdfField = blockIt->template getData< lbm::PdfField<LM> > ( pdfFieldID_ );

      SP & pressure = handling->template getBoundaryCondition< SP >( handling->getBoundaryUID( pressureFlag1_ ) );
      pressure.setLatticeDensity( real_t(1) );

      if ( storage_.atDomainMinBorder( flowAxis, *blockIt) )
      {
         // Set inflow
         CellInterval ci;
         auto direction = stencil::directionFromAxis( flowAxis, true );
         pdfField->getGhostRegion( direction, ci, 1, true );
         Cell offset( 0,0,0 );
         storage_.transformBlockLocalToGlobalCell( offset, *blockIt );

         for( auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt )
         {
            Cell  globalCell =  *cellIt + offset;

            Vector3<real_t> vel ( 0 );
            vel[ flowAxis ] = getVelocity( globalCell, scenario, flowAxis, maxVelocity, parabolaAxis );
            typename UBB::Velocity ubbVel ( vel );
            handling->forceBoundary( ubbFlag_, cellIt->x(), cellIt->y(), cellIt->z(), ubbVel  );
         }
      }

      if ( storage_.atDomainMaxBorder( flowAxis, *blockIt) )
      {
         // Set outflow
         CellInterval ci;
         auto direction = stencil::directionFromAxis( flowAxis, false );
         pdfField->getGhostRegion( direction, ci, 1, true );
         for( auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt )
            handling->forceBoundary( pressureFlag1_, cellIt->x(), cellIt->y(), cellIt->z()  );
      }
   }
}


template< typename BH_T, typename LM, typename SP, typename UBB >
void Poiseuille<BH_T,LM,SP,UBB>::initPressureBoundary( Axis flowAxis, real_t pressureDiff )
{
   using geometry::initializer::BoundaryFromDomainBorder;
   BoundaryFromDomainBorder<BH_T> borderInitializer( storage_, handlerID_ );

   stencil::Direction direction1 = stencil::directionFromAxis( flowAxis, true );
   stencil::Direction direction2 = stencil::directionFromAxis( flowAxis, false );

   borderInitializer.init( pressureFlag1_, direction1, -1 );
   borderInitializer.init( pressureFlag2_, direction2, -1 );

   for( auto blockIt = storage_.begin(); blockIt != storage_.end(); ++blockIt )
   {
      auto handling = blockIt->template getData<BH_T> ( handlerID_ );
      SP & pressure1 = handling->template getBoundaryCondition< SP >( handling->getBoundaryUID( pressureFlag1_ ) );
      SP & pressure2 = handling->template getBoundaryCondition< SP >( handling->getBoundaryUID( pressureFlag2_ ) );

      const real_t densityDiff = real_t(3) * pressureDiff;

      pressure1.setLatticeDensity ( real_t(1) + densityDiff / real_t(2) );
      pressure2.setLatticeDensity ( real_t(1) - densityDiff / real_t(2) );
   }
}

template< typename BH_T, typename LM, typename SP, typename UBB >
void Poiseuille<BH_T,LM,SP,UBB>::initNoSlipBoundaries( Scenario scenario, Axis flowAxis, Axis parabolaAxis )
{
   using geometry::initializer::BoundaryFromDomainBorder;
   BoundaryFromDomainBorder<BH_T> borderInitializer( storage_, handlerID_ );

   if ( scenario == RECT_2D )
   {
      borderInitializer.init( noSlipFlag_, stencil::directionFromAxis( parabolaAxis, false ), -1 );
      borderInitializer.init( noSlipFlag_, stencil::directionFromAxis( parabolaAxis, true ), -1 );
   }
   else if ( scenario == PIPE )
   {
      const real_t pipeRadius = getPipeRadius( scenario, flowAxis, parabolaAxis );
      const real_t pipeRadiusSq = pipeRadius * pipeRadius;

      for( auto blockIt = storage_.begin(); blockIt != storage_.end(); ++blockIt )
      {
         auto handling = blockIt->template getData<BH_T> ( handlerID_ );
         auto flagField = handling->getFlagField();


         Cell midCell (  storage_.getNumberOfXCells() / 2,
                         storage_.getNumberOfYCells() / 2,
                         storage_.getNumberOfZCells() / 2 );

         midCell[flowAxis] = 0;

         for( auto cellIt = flagField->beginWithGhostLayer(); cellIt != flagField->end(); ++cellIt )
         {
            Cell globalCell;
            storage_.transformBlockLocalToGlobalCell( globalCell, *blockIt, cellIt.cell() );
            globalCell[flowAxis] = 0;

            globalCell -= midCell;
            const real_t dist = real_c( globalCell[0] * globalCell[0] +
                                        globalCell[1] * globalCell[1] +
                                        globalCell[2] * globalCell[2] );
            if ( dist > pipeRadiusSq )
               handling->forceBoundary( noSlipFlag_, cellIt.x(), cellIt.y(), cellIt.z() );
         }
      }
   }
}

template< typename BH_T, typename LM, typename SP, typename UBB >
void Poiseuille<BH_T,LM,SP,UBB>::setInitVelocity( Scenario scenario, Axis flowAxis, real_t maxVelocity, Axis parabolaAxis )
{
   for( auto blockIt = storage_.begin(); blockIt != storage_.end(); ++blockIt )
   {
      auto handling = blockIt->template getData<BH_T> ( handlerID_ );
      auto pdfField = blockIt->template getData< lbm::PdfField<LM> > ( pdfFieldID_ );
      auto flagField = handling->getFlagField();
      // go over inner part - if velocity 0 -> set to noslip boundary
      for( auto cellIt = flagField->beginXYZ(); cellIt != flagField->end(); ++cellIt )
      {
         if ( !handling->isBoundary( cellIt ) )
         {
            Cell globalCell;
            storage_.transformBlockLocalToGlobalCell( globalCell, *blockIt, cellIt.cell() );

            real_t flowVel = getVelocity( globalCell, scenario, flowAxis, maxVelocity, parabolaAxis );
            Vector3<real_t> vel (0);
            vel[flowAxis] = flowVel;
            real_t density = pdfField->getDensity( cellIt.cell() ); // Preserve old density value
            pdfField->setDensityAndVelocity( cellIt.cell(), vel, density );
         }
      }
   }
}




template< typename BH_T, typename LM, typename SP, typename UBB >
real_t Poiseuille<BH_T,LM,SP,UBB>::getVelocity( const Cell & globalCell, Scenario scenario, Axis flowAxis, real_t maxVelocity, Axis parabolaAxis )
{
   Axis normalAxis1 = (flowAxis + 1) % 3;
   Axis normalAxis2 = (flowAxis + 2) % 3;

   if ( scenario == RECT_2D ) // rectangular setup with one axis noslip, one periodic
   {
      const real_t x   = real_c( globalCell[ parabolaAxis] );
      const real_t max = maxPoint_   [ parabolaAxis ];
      return - real_t(4) * x * (x-max) / ( max * max) * maxVelocity;
   }
   else if ( scenario == PIPE ) // pipe setup
   {
      const real_t xDiff = real_c( globalCell[normalAxis1] ) - midPoint_[ normalAxis1 ];
      const real_t yDiff = real_c( globalCell[normalAxis2] ) - midPoint_[ normalAxis2]  ;
      const real_t distSq = xDiff*xDiff + yDiff*yDiff;
      const real_t dist = std::sqrt( distSq );
      const real_t pipeRadius = getPipeRadius( scenario, flowAxis, parabolaAxis );
      if ( distSq < pipeRadius * pipeRadius )
         return - real_t(4) * dist * (dist - pipeRadius) / ( pipeRadius * pipeRadius) * maxVelocity;
      else  // outside pipe
         return 0;
   }
   return 0;
}


template< typename BH_T, typename LM, typename SP, typename UBB >
real_t Poiseuille<BH_T,LM,SP,UBB>::getPipeRadius( Scenario scenario, Axis flowAxis, Axis parabolaAxis ) const
{
   if ( scenario == PIPE )
   {
      const Axis normalAxis1 = ( flowAxis + 1) % 3;
      const Axis normalAxis2 = ( flowAxis + 2) % 3;
      return std::min( midPoint_[ normalAxis1], midPoint_[normalAxis2] );
   }
   else
      return midPoint_[ parabolaAxis ];
}

template< typename BH_T, typename LM, typename SP, typename UBB >
real_t Poiseuille<BH_T,LM,SP,UBB>::velocityFromPressureDiff( Scenario scenario, real_t pressureDiff, Axis flowAxis, Axis parabolaAxis )
{
   const real_t pipeRadius = getPipeRadius( scenario, flowAxis, parabolaAxis );
   real_t acceleration = pressureDiff / real_c( storage_.getNumberOfCells( flowAxis ) );
   real_t geometryFactor = (scenario == PIPE) ? real_t(4) : real_t(2);
   return acceleration / ( geometryFactor * latticeViscosity_ ) * pipeRadius * pipeRadius;
}

template< typename BH_T, typename LM, typename SP, typename UBB >
real_t Poiseuille<BH_T,LM,SP,UBB>::pressureDiffFromVelocity( Scenario scenario, real_t velocity, Axis flowAxis, Axis parabolaAxis )
{
   const real_t pipeRadius = getPipeRadius( scenario, flowAxis, parabolaAxis );
   real_t geometryFactor = (scenario == PIPE) ? real_t(4) : real_t(2);
   real_t acceleration = ( geometryFactor * velocity * latticeViscosity_ ) / ( pipeRadius * pipeRadius );
   return acceleration * real_c( storage_.getNumberOfCells( flowAxis ) );
}




}
} // namespace lbm
} // namespace walberla






