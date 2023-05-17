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
//! \file 01_BasicLBM.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"



namespace walberla {

using LatticeModel_T = lbm::D2Q9<lbm::collision_model::SRT>;
using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;

using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;



int main( int argc, char ** argv )
{
   walberla::Environment walberlaEnv( argc, argv );

   auto blocks = blockforest::createUniformBlockGridFromConfig( walberlaEnv.config() );

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock( "Parameters" );

   const real_t          omega           = parameters.getParameter< real_t >         ( "omega",           real_c( 1.4 ) );
   const Vector3<real_t> initialVelocity = parameters.getParameter< Vector3<real_t> >( "initialVelocity", Vector3<real_t>() );
   const uint_t          timesteps       = parameters.getParameter< uint_t >         ( "timesteps",       uint_c( 10 )  );

   const real_t remainingTimeLoggerFrequency = parameters.getParameter< real_t >( "remainingTimeLoggerFrequency", real_c(3.0) ); // in seconds

   // create fields
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ) );
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, initialVelocity, real_t(1) );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // create and initialize boundary handling
   const FlagUID fluidFlagUID( "Fluid" );

   auto boundariesConfig = walberlaEnv.config()->getOneBlock( "Boundaries" );

   using BHFactory = lbm::DefaultBoundaryHandlingFactory<LatticeModel_T, FlagField_T>;

   BlockDataID boundaryHandlingId = BHFactory::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldId, pdfFieldId, fluidFlagUID,
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity0", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< Vector3<real_t> >( "velocity1", Vector3<real_t>() ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure0", real_c( 1.0 ) ),
                                                                             boundariesConfig.getParameter< real_t > ( "pressure1", real_c( 1.0 ) ) );

   geometry::initBoundaryHandling<BHFactory::BoundaryHandling>( *blocks, boundaryHandlingId, boundariesConfig );
   geometry::setNonBoundaryCellsToDomain<BHFactory::BoundaryHandling> ( *blocks, boundaryHandlingId );

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( BHFactory::BoundaryHandling::getBlockSweep( boundaryHandlingId ), "boundary handling" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, fluidFlagUID ) ), "LB stream & collide" );

   // LBM stability check
   auto checkFunction = [](PdfField_T::value_type value) {return math::finite( value );};
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >( walberlaEnv.config(), blocks, pdfFieldId,
                                                                                                             flagFieldId, fluidFlagUID, checkFunction ) ),
                                  "LBM stability check" );

   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency ), "remaining time logger" );

   // add VTK output to time loop
   lbm::VTKOutput< LatticeModel_T, FlagField_T >::addToTimeloop( timeloop, blocks, walberlaEnv.config(), pdfFieldId, flagFieldId, fluidFlagUID );

   // create adaptors, so that the GUI also displays density and velocity
   // adaptors are like fields with the difference that they do not store values
   // but calculate the values based on other fields ( here the PdfField )
   field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::Density>       ( blocks, pdfFieldId, "DensityAdaptor" );
   field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::VelocityVector>( blocks, pdfFieldId, "VelocityAdaptor" );

   if( parameters.getParameter<bool>( "useGui", false ) )
   {
      GUI gui ( timeloop, blocks, argc, argv );
      lbm::connectToGui<LatticeModel_T> ( gui );
      gui.run();
   }
   else
      timeloop.run();

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
}