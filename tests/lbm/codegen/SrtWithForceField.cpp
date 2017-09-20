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
//! \file SrtWithForceField.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "SrtWithForceFieldModel.h"

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"



using namespace walberla;

typedef lbm::SrtWithForceFieldModel          LatticeModel_T;
typedef LatticeModel_T::Stencil              Stencil_T;
typedef LatticeModel_T::CommunicationStencil CommunicationStencil_T;
typedef lbm::PdfField< LatticeModel_T >      PdfField_T;

typedef GhostLayerField<real_t, LatticeModel_T::Stencil::D > ForceField_T;

typedef walberla::uint8_t    flag_t;
typedef FlagField< flag_t >  FlagField_T;



int main( int argc, char ** argv )
{
   walberla::Environment walberlaEnv( argc, argv );

   auto blocks = blockforest::createUniformBlockGridFromConfig( walberlaEnv.config() );

   // read parameters
   auto parameters = walberlaEnv.config()->getOneBlock( "Parameters" );

   const real_t          omega           = parameters.getParameter< real_t >         ( "omega",           real_c( 1.4 ) );
   const Vector3<real_t> initialVelocity = parameters.getParameter< Vector3<real_t> >( "initialVelocity", Vector3<real_t>() );
   const uint_t          timesteps       = parameters.getParameter< uint_t >         ( "timesteps",       uint_c( 10 )  );

   const double remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", 3.0 ); // in seconds

   // create force field

   // create fields
   BlockDataID forceFieldId = field::addToStorage<ForceField_T>(blocks, "Force", real_t(0.0) );

   LatticeModel_T latticeModel = LatticeModel_T( forceFieldId, omega );
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, initialVelocity, real_t(1) );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // create and initialize boundary handling
   const FlagUID fluidFlagUID( "Fluid" );

   auto boundariesConfig = walberlaEnv.config()->getOneBlock( "Boundaries" );

   typedef lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T > BHFactory;

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
   timeloop.add() << Sweep( LatticeModel_T::Sweep( pdfFieldId ), "LB stream & collide" );

   // LBM stability check
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >( walberlaEnv.config(), blocks, pdfFieldId,
                                                                                                             flagFieldId, fluidFlagUID ) ),
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
