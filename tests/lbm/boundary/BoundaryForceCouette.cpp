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
//! \file BoundaryForceCouette.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Calculates the force on top and bottom boundary in plane Couette flow
//!        and checks correctness against analytical force
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"

#include "geometry/initializer/BoundaryFromDomainBorder.h"
#include "geometry/InitBoundaryHandling.h"

#include "lbm/boundary/NoSlip.h"
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

typedef lbm::NoSlip< LatticeModel_T, flag_t, true >  BottomWall_T;
typedef lbm::SimpleUBB< LatticeModel_T, flag_t, false, true >  TopWall_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, BottomWall_T, TopWall_T > BoundaryHandling_T;



///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID BottomWall_Flag( "bottom wall" );
const FlagUID TopWall_Flag( "top wall" );



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
                                  BottomWall_T( "bottom", BottomWall_Flag, pdfField ),
                                  TopWall_T( "top", TopWall_Flag, pdfField, velocity_ ) );
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
   
   const uint_t xBlocks = ( processes == 8 ) ? uint_t(2) : uint_t(1);
   const uint_t yBlocks = ( processes >= 4 ) ? uint_t(2) : uint_t(1);
   const uint_t zBlocks = ( processes >= 2 ) ? uint_t(2) : uint_t(1);

   const uint_t L = uint_t(16);
   
   auto blocks = blockforest::createUniformBlockGrid( xBlocks, yBlocks, zBlocks,
                                                     L / xBlocks,
                                                     L / yBlocks,
                                                     L / zBlocks,
                                                     real_t(1),
                                                     true,
                                                     true, true, false);
   
   const real_t omega = real_c(1);
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
   
   const Vector3<real_t> velocity(real_t(0.01), real_t(0), real_t(0));
   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf field", latticeModel );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );
   BlockDataID boundaryHandlingId = blocks->addBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldId, pdfFieldId, velocity ),
                                                                               "boundary handling" );
   
   // set boundary conditions
   geometry::initializer::BoundaryFromBody< BoundaryHandling_T > bodyInitializer( *blocks, boundaryHandlingId );
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > borderInitializer( *blocks, boundaryHandlingId );
   borderInitializer.init( TopWall_Flag, stencil::T, cell_idx_t(-1) );
   borderInitializer.init( BottomWall_Flag, stencil::B, cell_idx_t(-1) );
   geometry::setNonBoundaryCellsToDomain< BoundaryHandling_T >( *blocks, boundaryHandlingId );
   
   uint_t timeSteps = uint_c(2000);
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
   
   // get the force acting on the walls
   real_t fXBottom = real_t(0), fXTop = real_t(0);
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      FlagField_T * flagField = block->getData< FlagField_T >( flagFieldId );
      BoundaryHandling_T * boundaryHandling = block->getData< BoundaryHandling_T >( boundaryHandlingId );
      auto topWallUid = boundaryHandling->getBoundaryUID( TopWall_Flag );
      auto bottomWallUid = boundaryHandling->getBoundaryUID( BottomWall_Flag );
      auto & topBc = boundaryHandling->getBoundaryCondition<TopWall_T>( topWallUid );
      auto & bottomBc = boundaryHandling->getBoundaryCondition<BottomWall_T>( bottomWallUid );
      auto topFlag = flagField->getFlag(TopWall_Flag);
      auto bottomFlag = flagField->getFlag(BottomWall_Flag);
      
      // need to iterate over ghost cells too because they may contain part of the force
      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ_OMP(flagField, omp parallel for schedule(static) reduction(+:fXTop) reduction(+:fXBottom),
      {
         if( flagField->isFlagSet(x, y, z, topFlag))
         {
            const Vector3<real_t> & cellForce = topBc.getForce(x, y, z);
            fXTop += cellForce[0];
         }

         if( flagField->isFlagSet(x, y, z, bottomFlag))
         {
            const Vector3<real_t> & cellForce = bottomBc.getForce(x, y, z);
            fXBottom += cellForce[0];
         }
      });
   }

   mpi::allReduceInplace( fXTop, mpi::SUM );
   mpi::allReduceInplace( fXBottom, mpi::SUM );
   
   real_t viscosity = lbm::collision_model::viscosityFromOmega( omega );
   real_t shearRate = velocity[0] / real_c(L);
   real_t wallArea = real_c(L) * real_c(L);
   real_t analyticalForce = shearRate * viscosity * wallArea;

   WALBERLA_LOG_RESULT_ON_ROOT("Expected force: " << analyticalForce);
   WALBERLA_LOG_RESULT_ON_ROOT("Actual force top: " << fXTop);
   WALBERLA_LOG_RESULT_ON_ROOT("Actual force bottom: " << fXBottom);

   real_t errorTop = std::fabs( ( std::fabs(fXTop) - analyticalForce ) / analyticalForce );
   real_t errorBottom = std::fabs( ( std::fabs(fXBottom) - analyticalForce ) / analyticalForce );

   WALBERLA_LOG_RESULT_ON_ROOT("Relative error top: " << (real_t(100)*errorTop) << " %");
   WALBERLA_LOG_RESULT_ON_ROOT("Relative error bottom: " << (real_t(100)*errorBottom) << " %");


   if( timeSteps > 100 )
   {
      WALBERLA_CHECK_LESS(errorTop, 0.01, "Error in force at top wall is too high! value = " <<  fXTop << " vs analytical = " << analyticalForce);
      WALBERLA_CHECK_LESS(errorBottom, 0.01, "Error in force at bottom wall is too high! value = " <<  fXBottom << " vs analytical = " << analyticalForce);
      WALBERLA_LOG_PROGRESS_ON_ROOT("Test succeeded");
   }

   return EXIT_SUCCESS;
}
