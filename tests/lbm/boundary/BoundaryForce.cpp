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
//! \file BoundaryForce.cpp
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//! \brief Calculates the drag force on a spherical no-slip boundary in a flow
//!        due to a velocity boundary condition and compares with Stokes' law.
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/math/Constants.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"

#include "geometry/initializer/BoundaryFromBody.h"
#include "geometry/initializer/BoundaryFromDomainBorder.h"
#include "geometry/InitBoundaryHandling.h"

#include "lbm/boundary/SimpleUBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "timeloop/SweepTimeloop.h"

using namespace walberla;



//////////////
// TYPEDEFS //
//////////////

using LatticeModel_T = lbm::D3Q19<lbm::collision_model::TRT>;
using Stencil_T = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

typedef lbm::SimpleUBB< LatticeModel_T, flag_t, false, true >  UBB_Sphere_T;
typedef lbm::SimpleUBB< LatticeModel_T, flag_t >  UBB_Wall_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, UBB_Sphere_T, UBB_Wall_T > BoundaryHandling_T;



///////////
// FLAGS //
///////////

const FlagUID        Fluid_Flag( "fluid" );
const FlagUID   UBB_Sphere_Flag( "sphere");
const FlagUID     UBB_Wall_Flag( "wall"  );



///////////////////////
// BOUNDARY HANDLING //
///////////////////////

class MyBoundaryHandling
{
public:
   
   MyBoundaryHandling( const BlockDataID & flagFieldId, const BlockDataID & pdfFieldId, const Vector3<real_t> & velocity ) :
   flagFieldId_( flagFieldId ), pdfFieldId_( pdfFieldId ), velocity_( velocity ) {}
   
   BoundaryHandling_T * operator()( IBlock* const block ) const;
   
private:
   
   const BlockDataID flagFieldId_;
   const BlockDataID  pdfFieldId_;
   
   const Vector3<real_t> velocity_;
   
}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block ) const
{
   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldId_ );
   
   const flag_t fluid = flagField->registerFlag( Fluid_Flag ); // register the fluid flag at the flag field
   
   return new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                  UBB_Sphere_T( "sphere", UBB_Sphere_Flag, pdfField, Vector3<real_t>() ),
                                  UBB_Wall_T(   "wall",   UBB_Wall_Flag,   pdfField, velocity_ ) );
}



//////////
// MAIN //
//////////

int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   
   int processes( MPIManager::instance()->numProcesses() );
   if( processes != 1 && processes != 2 && processes != 4 && processes != 8 )
      WALBERLA_ABORT( "The number of processes must be equal to either 1, 2, 4, or 8!" );
   
   const uint_t xBlocks = ( processes == 8 ) ? uint_c(2) : uint_c(1);
   const uint_t yBlocks = ( processes >= 4 ) ? uint_c(2) : uint_c(1);
   const uint_t zBlocks = ( processes >= 2 ) ? uint_c(2) : uint_c(1);

   const uint_t L = uint_t(64);
   
   auto blocks = blockforest::createUniformBlockGrid( xBlocks, yBlocks, zBlocks,
                                                     L / xBlocks,
                                                     L / yBlocks,
                                                     L / zBlocks,
                                                     real_c(1),
                                                     true,
                                                     false, false, true);
   
   const real_t omega = real_c(0.06);
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
   
   const Vector3<real_t> velocity(real_t(0), real_t(0), real_t(0.01));
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );
   BlockDataID boundaryHandlingId = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId, pdfFieldId, velocity ),
                                                                               "boundary handling" );
   
   // set boundary conditions
   geometry::initializer::BoundaryFromBody< BoundaryHandling_T > bodyInitializer( *blocks, boundaryHandlingId );
   const real_t R = real_t(5.4);
   geometry::Sphere sphere( Vector3<real_t>( real_t(L/2), real_t(L/2), real_t(L/2) ), R );
   bodyInitializer.template init< geometry::Sphere >( sphere, UBB_Sphere_Flag );
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > borderInitializer( *blocks, boundaryHandlingId );
   borderInitializer.init( UBB_Wall_Flag, stencil::N, cell_idx_t(-1) );
   borderInitializer.init( UBB_Wall_Flag, stencil::S, cell_idx_t(-1) );
   borderInitializer.init( UBB_Wall_Flag, stencil::E, cell_idx_t(-1) );
   borderInitializer.init( UBB_Wall_Flag, stencil::W, cell_idx_t(-1) );
   geometry::setNonBoundaryCellsToDomain< BoundaryHandling_T >( *blocks, boundaryHandlingId );
   
   uint_t timeSteps = uint_c(200);
   if( argc > 1 )
      timeSteps = uint_c(std::stoul( argv[1] ));
   SweepTimeloop timeloop( blocks->getBlockStorage(), timeSteps );
 
   blockforest::communication::UniformBufferedScheme< CommunicationStencil_T > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );
   
   timeloop.add() << BeforeFunction( communication, "communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "boundary handling" );

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag ) ), "LB stream & collide" );
   
   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );
  
   timeloop.run();
   
   // get the force acting on the sphere
   real_t fx = real_t(0), fy = real_t(0), fz = real_t(0);
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId );
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId );
      auto buid = boundaryHandling->getBoundaryUID( UBB_Sphere_Flag );
      auto & bc = boundaryHandling->getBoundaryCondition<UBB_Sphere_T>( buid );
      auto flag = flagField->getFlag(UBB_Sphere_Flag);
      
      // need to iterate over ghost cells too because they may contain part of the force
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(flagField, omp parallel for schedule(static) reduction(+:fx) reduction(+:fy) reduction(+:fz),
      {
         if( flagField->isFlagSet(x, y, z, flag))
         {
            const Vector3<real_t> & cellForce = bc.getForce(x, y, z);
            fx += cellForce[0];
            fy += cellForce[1];
            fz += cellForce[2];
         }
      });
   }
   Vector3<real_t> force(fx, fy, fz);
   
   mpi::allReduceInplace( force[0], mpi::SUM );
   mpi::allReduceInplace( force[1], mpi::SUM );
   mpi::allReduceInplace( force[2], mpi::SUM );
   
   real_t visc = lbm::collision_model::viscosityFromOmega( omega );
   Vector3<real_t> stokes = 6 * math::M_PI * visc * R * velocity;

   WALBERLA_LOG_RESULT_ON_ROOT("Expected force: " << stokes);
   WALBERLA_LOG_RESULT_ON_ROOT("Actual force: " << force);
   real_t err = (force-stokes).length() / stokes.length();
   WALBERLA_LOG_RESULT_ON_ROOT("Relative error: " << (real_t(100)*err) << " %");
   
   if( timeSteps > 100 )
   {
      WALBERLA_CHECK_LESS(err, 0.1);
      WALBERLA_LOG_PROGRESS_ON_ROOT("Test succeeded");
   }

   return EXIT_SUCCESS;
}
