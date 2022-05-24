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

/* This tests momentum balance for a fluctuating MRT LB with a constant
   force applied
*/

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "timeloop/all.h"

#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"

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

   const real_t          omega           = parameters.getParameter< real_t >         ( "omega",           real_c(1.4) );
   const real_t          magic_number    = 3. / 16;
   const real_t          omega_2         = (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
   const Vector3<real_t> initialVelocity = parameters.getParameter< Vector3<real_t> >( "initialVelocity", Vector3<real_t>() );
   const uint_t          timesteps       = parameters.getParameter< uint_t >         ( "timesteps",       uint_c(10.0) );
   const real_t          temperature     = real_c(0.01);
   const uint_t          seed            = uint_t(0);

   const real_t remainingTimeLoggerFrequency = parameters.getParameter< real_t >( "remainingTimeLoggerFrequency", real_c(3.0) ); // in seconds

   // create fields
   real_t force             = real_c(2E-4); // Force to apply on each node on each axis
   BlockDataID forceFieldId = field::addToStorage<VectorField_T>( blocks, "Force", force, field::fzyx );

   LatticeModel_T latticeModel = LatticeModel_T( forceFieldId, omega, omega, omega_2, omega, seed, temperature, uint_t(0) );
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel, initialVelocity, real_c(1.0), uint_t(1), field::fzyx );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // create and initialize flag field
   const FlagUID fluidFlagUID( "Fluid" );
   geometry::setNonBoundaryCellsToDomain<FlagField_T>(*blocks, flagFieldId, fluidFlagUID);

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );

   // set the RNG counter to match the time step and propagate it to the fields' copies of the lattice model
   timeloop.add() << BeforeFunction( [&]() { latticeModel.time_step_ = uint32_c(timeloop.getCurrentTimeStep()); }, "set RNG counter" )
                  << Sweep( [&]( IBlock * block ) {
                        auto field = block->getData< PdfField_T >( pdfFieldId );
                        field->latticeModel().time_step_ = latticeModel.time_step_;
                     }, "set RNG counter" );
   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( LatticeModel_T::Sweep( pdfFieldId ), "LB stream & collide" );

   // LBM stability check
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >( walberlaEnv.config(), blocks, pdfFieldId,
                                                                                                             flagFieldId, fluidFlagUID ) ),
                                  "LBM stability check" );

   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency), "remaining time logger" );

   auto densityAdaptorId = field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::Density>        ( blocks, pdfFieldId, "DensityAdaptor" );
   auto velocityAdaptorId = field::addFieldAdaptor<lbm::Adaptor<LatticeModel_T>::VelocityVector>( blocks, pdfFieldId, "VelocityAdaptor" );

   timeloop.run();

   // Calculate momentum
   Vector3<real_t> momentum; // observed momentum
   int count = 0;      // count of lb nodes traversed
   for (auto block = blocks->begin(); block != blocks->end(); ++block)
   {
      auto v   = block->getData< lbm::Adaptor< LatticeModel_T >::VelocityVector >(velocityAdaptorId);
      auto rho = block->getData< lbm::Adaptor< LatticeModel_T >::Density >(densityAdaptorId);
      WALBERLA_FOR_ALL_CELLS_XYZ_OMP(v, omp critical, {
         momentum += rho->get(x, y, z) * v->get(x, y, z);
         count++;
      })
   }

   // check
   real_t expected_momentum = real_c(count) * force * real_c(timesteps);
   printf("%g %g %g | %g\n", momentum[0], momentum[1], momentum[2], expected_momentum);
   WALBERLA_CHECK_FLOAT_EQUAL(momentum[0], expected_momentum)
   WALBERLA_CHECK_FLOAT_EQUAL(momentum[1], expected_momentum)
   WALBERLA_CHECK_FLOAT_EQUAL(momentum[2], expected_momentum)

   return EXIT_SUCCESS;
}
