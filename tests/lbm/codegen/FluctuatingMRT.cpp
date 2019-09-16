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
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "gui/all.h"
#include "timeloop/all.h"

#include "lbm/field/PdfField.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/gui/Connection.h"
#include "lbm/vtk/VTKOutput.h"

#include "FluctuatingMRT_LatticeModel.h"


using namespace walberla;

typedef lbm::FluctuatingMRT_LatticeModel          LatticeModel_T;
typedef LatticeModel_T::Stencil                   Stencil_T;
typedef LatticeModel_T::CommunicationStencil      CommunicationStencil_T;
typedef lbm::PdfField< LatticeModel_T >           PdfField_T;

typedef GhostLayerField< real_t, LatticeModel_T::Stencil::D > VectorField_T;
typedef GhostLayerField< real_t, 1 > ScalarField_T;

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
   const real_t          temperature     = real_t(0.01);
   const uint_t seed = uint_t(0);

   const double remainingTimeLoggerFrequency = parameters.getParameter< double >( "remainingTimeLoggerFrequency", 3.0 ); // in seconds

   // create fields
   BlockDataID forceFieldId = field::addToStorage<VectorField_T>( blocks, "Force", real_t( 0.0 ));

   LatticeModel_T latticeModel = LatticeModel_T( forceFieldId, omega, omega, seed, temperature, uint_t(0) );
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, initialVelocity, real_t(1) );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // create and initialize boundary handling
   const FlagUID fluidFlagUID( "Fluid" );

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction( [&](){ latticeModel.time_step = uint32_c(timeloop.getCurrentTimeStep()); }, "set RNG counter" )
                  << BeforeFunction( communication, "communication" )
                  << Sweep( LatticeModel_T::Sweep( pdfFieldId ), "LB stream & collide" );

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
