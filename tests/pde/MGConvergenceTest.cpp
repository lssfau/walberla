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
//! \file MGConvergenceTest.cpp
//! \ingroup pde
//! \author Dominik Bartuschat <dominik.bartuschat@fau.de>
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



//**********************************************************************************************************************
/*!\fn real_t initU( const shared_ptr< StructuredBlockStorage > & , const BlockDataID & , const real_t )
// \brief Initializes solution field with random values and ensures that mean value is zero,
//        such that a Laplace problem with periodic boundary conditions can be solved.
//
// \param blocks Shared pointer to structured block storage
// \param uId Block data ID for solution field
// \param cuboidValue Optional value for initial solution in cuboid (before setting average mean value to zero). Default value is 10
*/
void initU( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & uId )
{
   real_t sum = real_c(0);

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * u = block->getData< PdeField_T >( uId );
      CellInterval xyz = u->xyzSizeWithGhostLayer();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
         math::seedRandomGenerator( static_cast<unsigned int>( (p[0] * real_t(blocks->getNumberOfXCells()) + p[1]) * real_t(blocks->getNumberOfYCells()) + p[2]) );
         u->get( *cell ) = math::realRandom( real_t(-10), real_t(10) );
      }
   }
    
    
    // Set mean value to zero //
    uint_t numCells(0);
    
    // Initializing non-zero block with a given value in center of domain, relative to domain extension
    for( auto block = blocks->begin(); block != blocks->end(); ++block )
    {
        PdeField_T * u = block->getData< PdeField_T >( uId );
        CellInterval xyz = u->xyzSize();
        for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
        {
            sum += u->get( *cell );
            ++numCells;
        }
    }

    mpi::allReduceInplace( numCells, mpi::SUM );
    mpi::allReduceInplace( sum, mpi::SUM );
    
    real_t domainMeanVal(sum/real_c(numCells));
    // WALBERLA_LOG_RESULT_ON_ROOT("Mean value of all cell values: " << domainMeanVal );
    
    // Subtract mean value at each point to achieve zero mean
    for( auto block = blocks->begin(); block != blocks->end(); ++block )
    {
        PdeField_T * u = block->getData< PdeField_T >( uId );
        CellInterval xyz = u->xyzSizeWithGhostLayer();
        for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
        {
            u->get( *cell ) -= domainMeanVal;
        }
    }

    
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn real_t initURect( const shared_ptr< StructuredBlockStorage > & , const BlockDataID & , const real_t )
// \brief Initializes solution field with rectangular function, such that a Laplace problem with periodic boundary conditions can be solved.
//
// Sets a constant positive value in a cuboid at center of domain and reduces mean value in domain to zero afterwards.
//
// \param blocks Shared pointer to structured block storage
// \param uId Block data ID for solution field
// \param cuboidValue Optional value for initial solution in cuboid (before setting average mean value to zero). Default value is 10
*/
void initURect( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & uId, const real_t cuboidValue = 10 )
{
   real_t sumCuboidValues = real_c(0);
   uint_t numCells = 0;

   Vector3<uint_t> globalNumCells(blocks->getNumberOfXCells(), blocks->getNumberOfYCells(), blocks->getNumberOfZCells());
   WALBERLA_LOG_RESULT_ON_ROOT("Domain size: " << globalNumCells[0] << ", "  << globalNumCells[1] << ", "  << globalNumCells[2] << " cells");

   // extension of cuboid relative to domain size in each dimension
   real_t relCuboidExtension = real_c(0.5);
   Vector3<real_t> cuboidSize(globalNumCells);
   cuboidSize *= relCuboidExtension;

   WALBERLA_LOG_RESULT_ON_ROOT("Cuboid size: " << cuboidSize[0] << ", "  << cuboidSize[1] << ", "  << cuboidSize[2] );


   AABB cuboidAABB( real_t(0.5)*(real_c(globalNumCells[0]) - cuboidSize[0]), real_t(0.5)*(real_c(globalNumCells[1]) - cuboidSize[1]), real_t(0.5)*(real_c(globalNumCells[2]) - cuboidSize[2]),
                    real_t(0.5)*(real_c(globalNumCells[0]) + cuboidSize[0]), real_t(0.5)*(real_c(globalNumCells[1]) + cuboidSize[1]), real_t(0.5)*(real_c(globalNumCells[2]) + cuboidSize[2])
   );

   pde::Zeroize(blocks, uId);

   // Initializing non-zero block with a given value in center of domain, relative to domain extension
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * u = block->getData< PdeField_T >( uId );
      CellInterval xyz = u->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         // get global coordinate of current cell
         Cell globalCoordCell;
         blocks->transformBlockLocalToGlobalCell( globalCoordCell, *block, *cell );

         // set values in cuboid
         if(cuboidAABB.contains(real_c(globalCoordCell[0]),real_c(globalCoordCell[1]),real_c(globalCoordCell[2]))) {
            u->get( *cell ) = cuboidValue;
            sumCuboidValues += u->get( *cell );
         }
         ++numCells;
      }
   }

   // compute mean value of cuboid values and subtract it from all values in domain
   mpi::allReduceInplace( numCells, mpi::SUM );
   // WALBERLA_LOG_RESULT_ON_ROOT("Number of cells: " << numCells );

   mpi::allReduceInplace( sumCuboidValues, mpi::SUM );
   // WALBERLA_LOG_RESULT_ON_ROOT("Sum of overlapped cell values: " << sumCuboidValues );

   real_t domainMeanVal(sumCuboidValues/real_c(numCells));
   WALBERLA_LOG_RESULT_ON_ROOT("Mean value of all cell values: " << domainMeanVal );

   // Subtract mean value at each point to achieve zero mean
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdeField_T * u = block->getData< PdeField_T >( uId );
      CellInterval xyz = u->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         u->get( *cell ) -= domainMeanVal;
      }
   }

}
//**********************************************************************************************************************


//**********************************************************************************************************************
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
//**********************************************************************************************************************


//**********************************************************************************************************************
template <typename Field_T>
void clearField( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fieldId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      block->getData< Field_T >( fieldId )->set( typename Field_T::value_type() );
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn real_t runConvergenceConstStencil(const real_t , const real_t , const real_t ,
//                                   const uint_t , const uint_t , const uint_t ,
//                                   const real_t , const real_t , const real_t ,
//                                   const uint_t , const uint_t )
// \brief Runs convergence test with V-cycle with a fixed stencil for a given problem size
//        and returns average convergence rate of last few iterations.
//
// \param xDomainSize (Physical) domain size in x-direction
// \param yDomainSize (Physical) domain size in y-direction
// \param zDomainSize (Physical) domain size in z-direction
// \param xBlocks Number of blocks in x-direction
// \param yBlocks Number of blocks in y-direction
// \param zBlocks Number of blocks in z-direction
// \param xNumInnerCells  (Total) number of inner cells in x-direction
// \param yNumInnerCells  (Total) number of inner cells in y-direction
// \param zNumInnerCells  (Total) number of inner cells in z-direction
// \param numLvl Number of multigrid levels
// \param coarseIters Number of coarse-grid iterations
// \return Average convergence rate of last few V-cycles.
*/
real_t runConvergenceConstStencil(const real_t xDomainSize, const real_t yDomainSize, const real_t zDomainSize,
                                  const uint_t xBlocks, const uint_t yBlocks, const uint_t zBlocks,
                                  const real_t xNumInnerCells, const real_t yNumInnerCells, const real_t zNumInnerCells,
                                  const uint_t numLvl, const uint_t coarseIters)
{

   const uint_t xCells = uint_t(xNumInnerCells)/xBlocks;
   const uint_t yCells = uint_t(yNumInnerCells)/yBlocks;
   const uint_t zCells = uint_t(zNumInnerCells)/zBlocks;

   const real_t dx = xDomainSize / real_c( xBlocks * xCells );
   const real_t dy = yDomainSize / real_c( yBlocks * yCells );
   const real_t dz = zDomainSize / real_c( zBlocks * zCells );
   auto blocks = blockforest::createUniformBlockGrid( math::AABB( real_t(0), real_t(0), real_t(0),
                                                                  xDomainSize, yDomainSize, zDomainSize ),
                                                      xBlocks, yBlocks, zBlocks,    // number of blocks
                                                      xCells, yCells, zCells,       // number of cells per block
                                                      true,
                                                      true, true, true );


   WALBERLA_LOG_RESULT_ON_ROOT("Discretization dx: " << dx << ", " << dy << ", " << dz);

   BlockDataID uId = field::addToStorage< PdeField_T >( blocks, "u", real_t(0), field::zyxf, uint_t(1) );

   initU(blocks, uId);
   // initURect( blocks, uId );

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

   for (uint_t i = 0; i < Stencil_T::Size; ++i)
      WALBERLA_LOG_RESULT_ON_ROOT("Weights on finest level (" << i << ") = " << weights[i] );

   // compute initial residual before V-Cycle starts
   blockforest::communication::UniformBufferedScheme< Stencil_T > communication( blocks );
   communication.addPackInfo( make_shared< field::communication::PackInfo< PdeField_T > >( uId ) );
   communication();

   auto residualNorm = pde::ResidualNorm<Stencil_T>( blocks->getBlockStorage(), uId, fId, weights );
   real_t initialResidualNorm = residualNorm();
   WALBERLA_LOG_RESULT_ON_ROOT("Initial residual norm " << initialResidualNorm);

   auto solver = walberla::make_shared< pde::VCycles< Stencil_T > >( blocks, uId, fId, weights, uint_t( 20 ),                                         // iterations
                                                                     numLvl,                                                                          // levels
                                                                     3, 3, coarseIters,                                                               // pre-smoothing, post-smoothing, coarse-grid iterations
                                                                     pde::ResidualNorm< Stencil_T >( blocks->getBlockStorage(), uId, fId, weights ),  // residual norm functor
                                                                     initialResidualNorm*real_c( 1.5e-13 ) );                                         // target precision relative to initial residual
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( solver ), "Cell-centered multigrid V-cycles" );

   timeloop.run();

   const auto & convrate = solver->convergenceRate();

   for (uint_t i = 1; i < convrate.size(); ++i)
   {
      WALBERLA_LOG_RESULT_ON_ROOT("Convergence rate in iteration " << i << ": " << convrate[i]);
      WALBERLA_CHECK_LESS(convrate[i], real_t(0.1));
   }

   // computing average convergence rate of last few V-cycles
   real_t averagConvrate = real_c(0);
   uint_t numConsideredConvrates(3);

   if(convrate.size() < numConsideredConvrates)
      WALBERLA_LOG_RESULT_ON_ROOT("Too few V-cycles to compute average convergence rate");

   for (uint_t j = 0; j < numConsideredConvrates; ++j)
   {
      averagConvrate+=convrate[convrate.size()-1-j];
   }

   return averagConvrate/real_c(numConsideredConvrates);

}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn convergenceCheck(const std::vector<real_t> &finalConvRates, const real_t relDevCR_Limit)
// \brief Check that convergence rates (for different problem sizes) are nearly the same.
//        Throws error if difference higher than allowed threshold.
//
// \param finalConvRates Vector with average convergence rates of tests with different problem sizes.
// \param relDevCR_Limit Upper limit for allowed relative deviation of the convergence rates from first convergence rate.
*/
void convergenceCheck(const std::vector<real_t> &finalConvRates, const real_t relDevCR_Limit) {

   uint_t numTests = finalConvRates.size();

   WALBERLA_LOG_RESULT_ON_ROOT("----------------------------------------" );
   for (uint_t b = 0; b < numTests; ++b)
   {
      WALBERLA_LOG_RESULT_ON_ROOT("Final convergence rate of test: " << b+1 << ": " << finalConvRates[b]);
   }

   WALBERLA_MPI_BARRIER();
   WALBERLA_LOG_RESULT_ON_ROOT("----------------------------------------" );

   for (uint_t c = 1; c < numTests; ++c)
   {
      // relative deviation of current convergence rate to first convergence rate
      real_t relDevFirstConvRate = (finalConvRates[c]-finalConvRates[0]) / finalConvRates[0];
      WALBERLA_LOG_RESULT_ON_ROOT("Final convergence rate of test: " << c+1 << ": " << finalConvRates[c] << " has relative deviation of " << relDevFirstConvRate*100 << "% from convergence rate at smallest problem size" );

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_CHECK_LESS( fabs(relDevFirstConvRate), relDevCR_Limit, "Relative deviation of final convergence rate in test " << c+1 << " from first test is with abs(" << relDevFirstConvRate << ") higher than limit " << relDevCR_Limit*100 << "%" );
      }
   }
   WALBERLA_LOG_RESULT_ON_ROOT("----------------------------------------" );

}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\brief Checks if convergence rates of multigrid V-cycle with constant stencils are independent of problem size
//        and h-independent. Throws error if convergence rate variations higher than allowed threshold.
*/
int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   const uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   if( processes != uint_t(1) && processes != uint_t(8) )
      WALBERLA_ABORT( "The number of processes must be equal to 1 or 8!" );

   logging::Logging::printHeaderOnStream();
   WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }

   const uint_t xBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t yBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);
   const uint_t zBlocks = ( processes == uint_t(1) ) ? uint_t(1) : uint_t(2);

   // Limit of relative deviation of convergence rates
   const real_t relDevCR_Limit = real_c(0.16);
   const uint_t numTests(4);
   const uint_t initNumLvl(3);


   /////////////////////////////
   // Tests with cubic domain //
   /////////////////////////////

   WALBERLA_LOG_RESULT_ON_ROOT("Running convergence test with increasing cubic domain size and constant h:");

   // number of inner cells
   real_t xSize = real_c(32);
   real_t ySize = real_c(32);
   real_t zSize = real_c(32);

   // physical domain size
   real_t xDomainSize(xSize);
   real_t yDomainSize(ySize);
   real_t zDomainSize(zSize);

   uint_t numLvl = initNumLvl;
   uint_t coarseIters = 10;

   std::vector<real_t> finalConvRates;

   // run test with constant stencil and increasing physical domain size (s.th. h=1)
   for (uint_t a = 0; a < numTests; ++a)
   {

      // Increase domain size
      xDomainSize=xSize;
      yDomainSize=ySize;
      zDomainSize=zSize;

      const auto finalConvrate = runConvergenceConstStencil(xDomainSize, yDomainSize, zDomainSize, xBlocks, yBlocks, zBlocks, xSize, ySize, zSize, numLvl, coarseIters);
      finalConvRates.push_back(finalConvrate);

      xSize*=2;
      ySize*=2;
      zSize*=2;
      ++numLvl;

   }

   convergenceCheck(finalConvRates, relDevCR_Limit);


   WALBERLA_LOG_RESULT_ON_ROOT("Running convergence test with fixed cubic domain size and decreasing h:");

   xSize = 32;
   ySize = 32;
   zSize = 32;

   xDomainSize=xSize;
   yDomainSize=ySize;
   zDomainSize=zSize;

   numLvl = initNumLvl;
   finalConvRates.clear();

   // run test with constant stencil and fixed physical domain size (s.th. h is halved in each test)
   for (uint_t b = 0; b < numTests; ++b)
   {

      const auto finalConvrate = runConvergenceConstStencil(xDomainSize, yDomainSize, zDomainSize, xBlocks, yBlocks, zBlocks, xSize, ySize, zSize, numLvl, coarseIters);
      finalConvRates.push_back(finalConvrate);

      xSize*=2;
      ySize*=2;
      zSize*=2;
      ++numLvl;

   }

   convergenceCheck(finalConvRates, relDevCR_Limit);


   /////////////////////////////////
   // Tests with non-cubic domain //
   /////////////////////////////////

   // re-run with non-cubic domain
   WALBERLA_LOG_RESULT_ON_ROOT("Running convergence test with non-cubic domain:");

   // re-set initial problem size
   xSize = 64;
   ySize = 32;
   zSize = 32;
   numLvl = initNumLvl;

   coarseIters = 16;
   finalConvRates.clear();

   for (uint_t a = 0; a < numTests; ++a)
   {

      const auto finalConvrate = runConvergenceConstStencil(xSize, ySize, zSize, xBlocks, yBlocks, zBlocks, xSize, ySize, zSize, numLvl, coarseIters);
      finalConvRates.push_back(finalConvrate);

      xSize*=2;
      ySize*=2;
      zSize*=2;
      ++numLvl;

   }
   WALBERLA_LOG_RESULT_ON_ROOT( "Finished convergence rate check");

   convergenceCheck(finalConvRates, relDevCR_Limit);

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}
//**********************************************************************************************************************
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}