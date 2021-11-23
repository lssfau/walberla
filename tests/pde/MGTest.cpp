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
//! \file MGTest.cpp
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/math/Random.h"
#include "core/SharedFunctor.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/vtk/VTKWriter.h"

#include "pde/ResidualNorm.h"
#include "pde/iterations/VCycles.h"
#include "pde/sweeps/Multigrid.h"

#include "stencil/D3Q7.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cmath>

namespace walberla {



typedef GhostLayerField< real_t, 1 > PdeField_T;
using Stencil_T = stencil::D3Q7;
using StencilField_T = pde::VCycles<Stencil_T>::StencilField_T;



void initU( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & uId )
{
   real_t sum = 0;

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * u = block->getData< PdeField_T >( uId );
      CellInterval xyz = u->xyzSizeWithGhostLayer();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
         math::seedRandomGenerator( static_cast<unsigned int>( (p[0] * real_t(blocks->getNumberOfXCells()) + p[1]) * real_t(blocks->getNumberOfYCells()) + p[2] ) );
         u->get( *cell ) = math::realRandom( real_t(-10), real_t(10) );
         sum += u->get( *cell );
      }
   }
   WALBERLA_UNUSED(sum);
}



void copyWeightsToStencilField( const shared_ptr< StructuredBlockStorage > & blocks, const std::vector<real_t> & weights, const BlockDataID & stencilId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      StencilField_T * stencil = block->getData< StencilField_T >( stencilId );

      WALBERLA_FOR_ALL_CELLS_XYZ(stencil,
         for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
            stencil->get(x,y,z,dir.toIdx()) = weights[ dir.toIdx() ];
      );
   }
}



template <typename Field_T>
void clearField( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fieldId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      block->getData< Field_T >( fieldId )->set( typename Field_T::value_type() );
   }
}



void checkProlongateRestrict( const shared_ptr< StructuredBlockStorage > & blocks )
{
   auto getSize1 = [] ( const shared_ptr< StructuredBlockStorage > & sbs, IBlock * const b ) {
      return pde::VCycles<Stencil_T>::getSizeForLevel(1, sbs, b);
   };
   auto getSize2 = [] ( const shared_ptr< StructuredBlockStorage > & sbs, IBlock * const b ) {
      return pde::VCycles<Stencil_T>::getSizeForLevel(2, sbs, b);
   };
   
   BlockDataID originalId = field::addToStorage< PdeField_T >( blocks, "test_0", getSize2 );
   initU( blocks, originalId );
   BlockDataID coarseId = field::addToStorage< PdeField_T >( blocks, "test_1", getSize2 );
   BlockDataID fineId = field::addToStorage< PdeField_T >( blocks, "test_2", getSize1 );

   auto c2f = pde::ProlongateAndCorrect<Stencil_T>( blocks, coarseId, fineId );
   auto f2c = pde::Restrict<Stencil_T>( blocks, fineId, coarseId );

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * orig = block->getData< PdeField_T >( originalId );
      PdeField_T * coarse = block->getData< PdeField_T >( coarseId );
      
      coarse->set( orig );
      c2f( &*block );
      coarse->set( real_t(0.0) );
      f2c( &*block );
      
      CellInterval xyz = coarse->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         WALBERLA_CHECK_FLOAT_EQUAL(coarse->get(*cell), orig->get(*cell), "Restrict-after-Prolongate should be an identity operation");
      }
   }

   blocks->clearBlockData( originalId );
   blocks->clearBlockData( coarseId );
   blocks->clearBlockData( fineId );
   WALBERLA_LOG_RESULT_ON_ROOT("Confirmed that restrict-after-prolongate is identity operation")
}



int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   const uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   if( processes != uint_t(1) && processes != uint_t(8) )
      WALBERLA_ABORT( "The number of processes must be equal to 1 or 8!" );

   logging::Logging::printHeaderOnStream();
   WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }

   bool shortrun = false;
   for( int i = 1; i < argc; ++i )
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) shortrun = true;

   const uint_t numLvl = shortrun ? uint_t(3) : uint_t(5);

   const uint_t xBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t yBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t zBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t xCells = uint_t(128)/xBlocks;
   const uint_t yCells = uint_t(128)/yBlocks;
   const uint_t zCells = uint_t(128)/zBlocks;
   const real_t xSize = real_t(xCells * xBlocks);
   const real_t ySize = real_t(yCells * yBlocks);
   const real_t zSize = real_t(zCells * zBlocks);
   const real_t dx = xSize / real_c( xBlocks * xCells + uint_t(1) );
   const real_t dy = ySize / real_c( yBlocks * yCells + uint_t(1) );
   const real_t dz =  zSize / real_c( zBlocks * zCells + uint_t(1) );
   auto blocks = blockforest::createUniformBlockGrid( math::AABB( real_t(0.5) * dx, real_t(0.5) * dy, real_t(0.5) * dz,
                                                                  xSize - real_t(0.5) * dx, ySize - real_t(0.5) * dy, zSize - real_t(0.5) * dz ),
                                                      xBlocks, yBlocks, zBlocks,
                                                      xCells, yCells, zCells,
                                                      true,
                                                      true, true, true );

   // run other tests
   checkProlongateRestrict( blocks );

   // run the main test

   BlockDataID uId = field::addToStorage< PdeField_T >( blocks, "u", real_t(0), field::zyxf, uint_t(1) );

   initU( blocks, uId );

   BlockDataID fId = field::addToStorage< PdeField_T >( blocks, "f", real_t(0), field::zyxf, uint_t(1) );

   SweepTimeloop timeloop( blocks, uint_t(1) );

   std::vector< real_t > weights( Stencil_T::Size );
   weights[ Stencil_T::idx[ stencil::C ] ] = real_t(2) / ( blocks->dx() * blocks->dx() ) + real_t(2) / ( blocks->dy() * blocks->dy() ) + real_t(2) / ( blocks->dz() * blocks->dz() );
   weights[ Stencil_T::idx[ stencil::N ] ] = real_t(-1) / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::S ] ] = real_t(-1) / ( blocks->dy() * blocks->dy() );
   weights[ Stencil_T::idx[ stencil::E ] ] = real_t(-1) / ( blocks->dx() * blocks->dx() );
   weights[ Stencil_T::idx[ stencil::W ] ] = real_t(-1) / ( blocks->dx() * blocks->dx() );
   weights[ Stencil_T::idx[ stencil::T ] ] = real_t(-1) / ( blocks->dx() * blocks->dz() );
   weights[ Stencil_T::idx[ stencil::B ] ] = real_t(-1) / ( blocks->dx() * blocks->dz() );

   auto solverDCA = walberla::make_shared<pde::VCycles< Stencil_T > >( blocks, uId, fId, weights,
                                                                       shortrun ? uint_t(1) : uint_t(20),                                              // iterations
                                                                       numLvl,                                               							// levels
                                                                       3, 3, 10,                                                                       // pre-smoothing, post-smoothing, coarse-grid iterations
                                                                       pde::ResidualNorm< Stencil_T >( blocks->getBlockStorage(), uId, fId, weights ), // residual norm functor
                                                                       real_c(1e-12) );                                                                // target precision
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor(solverDCA), "Cell-centered multigrid V-cycles with uniformly constant stencil" );

   timeloop.run();

   if( !shortrun )
   {
      auto & convrate = solverDCA->convergenceRate();
      for (uint_t i = 1; i < convrate.size(); ++i)
      {
         WALBERLA_LOG_RESULT_ON_ROOT("Convergence rate in iteration " << i << ": " << convrate[i]);
         WALBERLA_CHECK_LESS(convrate[i], real_t(0.1));
      }

      vtk::writeDomainDecomposition( blocks );
      field::createVTKOutput< PdeField_T >( uId, *blocks, "solution" )();
   }

   // rerun the test with a stencil field and DCA

   clearField<PdeField_T>( blocks, uId);
   initU( blocks, uId );

   BlockDataID stencilId = field::addToStorage< StencilField_T >( blocks, "w" );

   SweepTimeloop timeloop2( blocks, uint_t(1) );

   copyWeightsToStencilField( blocks, weights, stencilId );

   pde::CoarsenStencilFieldsDCA<Stencil_T>  coarsenWithDCA( blocks, numLvl, uint_t(2));		// Set up DCA object with operator order 2 (Laplace)

   solverDCA = walberla::make_shared<pde::VCycles< Stencil_T, decltype(coarsenWithDCA) > >(
		   	   	   	   	   	   	   	   	   	   	   	   	   	  blocks, uId, fId, stencilId, coarsenWithDCA,
                                                              shortrun ? uint_t(1) : uint_t(20),                                              	// iterations
                                                              numLvl,                                               							// levels
                                                              3, 3, 10,                                                                       	// pre-smoothing, post-smoothing, coarse-grid iterations
                                                              pde::ResidualNormStencilField< Stencil_T >( blocks->getBlockStorage(), uId, fId, stencilId ), // residual norm functor
                                                              real_c(1e-12) );                                                                	// target precision
   timeloop2.addFuncBeforeTimeStep( makeSharedFunctor(solverDCA), "Cell-centered multigrid V-cycles with stencil field and direct coarsening " );

   timeloop2.run();

   if( !shortrun )
   {
      auto & convrate = solverDCA->convergenceRate();
      for (uint_t i = 1; i < convrate.size(); ++i)
      {
         WALBERLA_LOG_RESULT_ON_ROOT("Convergence rate in iteration " << i << ": " << convrate[i]);
         WALBERLA_CHECK_LESS(convrate[i], real_t(0.1));
      }
   }

   // rerun the test with a stencil field and GCA

   clearField<PdeField_T>( blocks, uId);
   initU( blocks, uId );

   SweepTimeloop timeloop3( blocks, uint_t(1) );

   pde::CoarsenStencilFieldsGCA<Stencil_T>  coarsenWithGCA( blocks, numLvl, real_t(2));		// Set up GCA object with overrelaxation factor 2 (only valid for Poisson equation)

   auto solverGCA = walberla::make_shared<pde::VCycles< Stencil_T, decltype(coarsenWithGCA) > >(
		                                                      blocks, uId, fId, stencilId, coarsenWithGCA,
                                                              shortrun ? uint_t(1) : uint_t(20),                                            // iterations
                                                              numLvl,                                               						// levels
                                                              3, 3, 10,                                                                       // pre-smoothing, post-smoothing, coarse-grid iterations
                                                              pde::ResidualNormStencilField< Stencil_T >( blocks->getBlockStorage(), uId, fId, stencilId ), // residual norm functor
                                                              real_c(1e-12) );                                                                // target precision
   timeloop3.addFuncBeforeTimeStep( makeSharedFunctor(solverGCA), "Cell-centered multigrid V-cycles with stencil field and Galerkin coarsening " );

   timeloop3.run();

   if( !shortrun )
   {
      auto & convrate = solverGCA->convergenceRate();
      for (uint_t i = 1; i < convrate.size(); ++i)
      {
         WALBERLA_LOG_RESULT_ON_ROOT("Convergence rate in iteration " << i << ": " << convrate[i]);
         WALBERLA_CHECK_LESS(convrate[i], real_t(0.1));
      }
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}