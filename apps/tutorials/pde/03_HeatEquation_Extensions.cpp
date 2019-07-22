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
//! \file 06_HeatEquation_Extensions.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/Reduce.h"

#include "field/Field.h"
#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "stencil/D2Q5.h"
#include "stencil/D2Q9.h"

#include "gui/Gui.h"
#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cmath>


namespace walberla {


typedef GhostLayerField<real_t,1> ScalarField;
using Stencil_T = stencil::D2Q5;



// function to initialise the field holding the unknown u
void initU( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & uID )
{
   // iterate all blocks with an iterator 'block'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      // get the field data out of the block
      auto u = block->getData< ScalarField > ( uID );

      // obtain a CellInterval object that holds information about the number of cells in x,y,z direction of the field excluding ghost layers
      CellInterval xyz = u->xyzSize();

      // iterate all (inner) cells in the field
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell ){

         // obtain the physical coordinate of the center of the current cell
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );

         // set the initial condition, given by the function u(x,y,0) = sin(PI*x)*sin(PI*y)
         u->get( *cell ) = std::sin( math::pi * p[0] ) * std::sin( math::pi * p[1] );
      }
   }
}



class JacobiIterationResidual
{
public:

   //constructor
   JacobiIterationResidual( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & rhsID, const std::vector< real_t > & weights,
                    const shared_ptr< StructuredBlockStorage > & blocks, blockforest::communication::UniformBufferedScheme<Stencil_T> & myCommScheme,
                    const uint_t & maxIterations, const real_t & residualThreshold )
   : srcID_( srcID ), dstID_( dstID ), rhsID_( rhsID ), weights_( weights ), blocks_( blocks ), myCommScheme_ ( myCommScheme ),
     maxIterations_( maxIterations ), residualThreshold_( residualThreshold )
   {

      // initialize the cells_ variable
      init();
   }

   // functor that carries out the full Jacobi iteration
   void operator()();

   // function to obtain the total number of cells in the domain and store it in cells_, needed for the calculation of the residual
   void init();

   // computed the weighted L2 norm of the residual, serving as a termination criterion for the Jacobi method
   real_t residualNorm();

private:
   const BlockDataID srcID_;
   const BlockDataID dstID_;
   const BlockDataID rhsID_;

   std::vector< real_t > weights_;
   const shared_ptr< StructuredBlockStorage > blocks_;
   blockforest::communication::UniformBufferedScheme< Stencil_T > myCommScheme_;

   const uint_t maxIterations_;

   real_t cells_;
   const real_t residualThreshold_;


};

void JacobiIterationResidual::operator()()
{
   // Jacobi iteration, at most maxIterations_ times
   for( uint_t i = 0; i < maxIterations_; ++i )
   {
      // communicate to update ghost layer cells
      myCommScheme_();

      // iterate all blocks
      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         // get the source, destination and rhs field data from the block
         auto src = block->getData< ScalarField >( srcID_ );
         auto dst = block->getData< ScalarField >( dstID_ );
         auto rhs = block->getData< ScalarField >( rhsID_ );

         // iterate all cells of the fields and carry out the Jacobi sweep
         WALBERLA_FOR_ALL_CELLS_XYZ(src,

            dst->get(x,y,z) =  rhs->get(x,y,z);

            // iterate the neighboring cells and multiply their value with the respective weight
            for( auto dir = Stencil_T::beginNoCenter(); dir != Stencil_T::end(); ++dir )
               dst->get(x,y,z) -= weights_[ dir.toIdx() ] * src->getNeighbor(x,y,z,*dir);

            dst->get(x,y,z) /= weights_[ Stencil_T::idx[ stencil::C ] ];
         )

         // swap source and destination fields
         src->swapDataPointers( dst );
      }

      // check if residual norm is below threshold, if so stop the iteration
      if( residualNorm() < residualThreshold_ )
      {
         // only the root process will output this, to avoid multiple output when running with more processes
         WALBERLA_ROOT_SECTION() { std::cout << "Terminated Jacobi iteration after " << i << " iterations." << std::endl; }
         break;
      }
   }
}

void JacobiIterationResidual::init()
{
   // temporal storage
   uint_t cells( uint_c(0) );

   // iterate all blocks
   for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
   {
      // get the number of cells in each block and sum it up
      auto u = block->getData< ScalarField >( srcID_ );
      cells += u->xyzSize().numCells();
   }

   cells_ = real_c( cells );

   // communicate with other processes and sum up their local cells_ values to get the global number of cells in the domain
   mpi::allReduceInplace( cells_, mpi::SUM );
}

real_t JacobiIterationResidual::residualNorm()
{
   real_t norm( real_c(0) );

   // iterate all blocks
   for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
   {
      // get the field of the unknown u and the right-hand side
      auto u   = block->getData< ScalarField >( srcID_ );
      auto rhs = block->getData< ScalarField >( rhsID_ );

      // temporal storage
      real_t residual( real_c(0) );

      // iterate all cells inside the block
      WALBERLA_FOR_ALL_CELLS_XYZ(u,

         // calculates the residual r = rhs - A*u
         // This means taking the right-hand side value and subtracting the stencil applied to the unknown u.
         residual =  rhs->get(x,y,z);

         // iterate the cells (including the center) and multiply their value with the respective weight
         // This corresponds to applying the stencil to the unknown u.
         for( auto dir = Stencil_T::begin(); dir != Stencil_T::end(); ++dir )
            residual -= weights_[ dir.toIdx() ] * u->getNeighbor(x,y,z,*dir);

         // store the squared residual value
         norm += residual * residual;
      )
   }

   // communicate with other processes and sum up their local squared residual values to get the global squared residual value for the entire domain
   mpi::allReduceInplace( norm, mpi::SUM );

   // divide by the number of cells and take the square root to finally arrive at the value for the weighted L2 norm of the residual
   norm = std::sqrt( norm / cells_ );

   return norm;
}



class Reinitialize
{
public:

   Reinitialize( const BlockDataID & srcID, const BlockDataID & rhsID )
      : srcID_( srcID ), rhsID_( rhsID ) {}

   void operator()( IBlock * block );

private:

   const BlockDataID srcID_;
   const BlockDataID rhsID_;
};

void Reinitialize::operator()( IBlock * block )
{
   // get the source and rhs field data from the block
   auto src = block->getData< ScalarField >( srcID_ );
   auto rhs = block->getData< ScalarField >( rhsID_ );

   // swap source and right-hand side fields as the old solution is the right-hand side when doing implicit time stepping
   src->swapDataPointers( rhs );
}




int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );

   // parameters
   // time step size
   real_t dt = real_c(0.01);
   // thermal diffusivity
   real_t kappa = real_c(1);

   // number of blocks in x,y,z direction
   const uint_t xBlocks = uint_c(1);
   const uint_t yBlocks = uint_c(1);
   const uint_t zBlocks = uint_c(1);

   // get the number of processes
   const uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   // check if this number is equal to the number of blocks, since there is one block per process
   if( processes != xBlocks * yBlocks * zBlocks )
      WALBERLA_ABORT( "The number of processes must be equal the total number of blocks" );


   // number of cells per block in x,y,z direction
   // two dimensional layout, so only one cell in z direction
   const uint_t xCells = uint_c(25);
   const uint_t yCells = uint_c(25);
   const uint_t zCells = uint_c(1);

   // domain size is [0,1]x[0,1]
   const real_t xMin = real_c(0);
   const real_t xMax = real_c(1);
   const real_t yMin = real_c(0);
   const real_t yMax = real_c(1);

   // mesh spacings in x and y direction
   const real_t dx = (xMax - xMin)/real_c( xBlocks*xCells + uint_c(1) );
   const real_t dy = (yMax - yMin)/real_c( yBlocks*yCells + uint_c(1) );

   // create axis-aligned bounding box to define domain
   // defines a rectangular domain by specifying (xMin, yMin, zMin, xMax, yMax, zMax)
   // care has to be taken regarding the cell-centered data layout in walberla
   // in z-direction only one cell layer with thickness dx is present
   auto aabb = math::AABB( xMin + real_c(0.5)*dx, yMin + real_c(0.5)*dy, real_c(0),
                           xMax - real_c(0.5)*dx, yMax - real_c(0.5)*dy, dx );

   // create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
             aabb,                               // axis-aligned bounding box
             xBlocks, yBlocks, zBlocks,          // number of blocks in x,y,z direction
             xCells,  yCells,  zCells,           // how many cells per block (x,y,z)
             true,                               // one block per process
             false, false, false );              // no periodicity


   // add fields with ghost layers to all blocks
   // source and destination fields for the unknowns u, required by the Jacobi method
   BlockDataID srcID = field::addToStorage< ScalarField >( blocks, "src", real_c(0), field::zyxf, uint_c(1));
   BlockDataID dstID = field::addToStorage< ScalarField >( blocks, "dst", real_c(0), field::zyxf, uint_c(1));
   // field to store the right-hand side of the equation
   BlockDataID rhsID = field::addToStorage< ScalarField >( blocks, "rhs", real_c(0), field::zyxf, uint_c(1));

   // set the field to the initial condition u(x,y,0)
   initU( blocks, srcID );

   // create the communication scheme for the D2Q5 stencil
   blockforest::communication::UniformBufferedScheme< Stencil_T > myCommScheme ( blocks );
   // add a PackInfo that packs/unpacks the source field
   // Since the RHS field is unchanged, there is no need to communicate it.
   myCommScheme.addPackInfo( make_shared< field::communication::PackInfo<ScalarField> >( srcID ) );

   // set up the stencil weights
   std::vector< real_t > weights( Stencil_T::Size );
   weights[ Stencil_T::idx[ stencil::C ] ] = real_c(2) * dt * kappa / ( dx * dx ) + real_c(2) * dt * kappa / ( dy * dy ) + real_c(1);
   weights[ Stencil_T::idx[ stencil::E ] ] = - dt * kappa / ( dx * dx );
   weights[ Stencil_T::idx[ stencil::W ] ] = - dt * kappa / ( dx * dx );
   weights[ Stencil_T::idx[ stencil::N ] ] = - dt * kappa / ( dy * dy );
   weights[ Stencil_T::idx[ stencil::S ] ] = - dt * kappa / ( dy * dy );

   // set up the timeloop object
   SweepTimeloop timeloop ( blocks, uint_c(20) );

   // add the routine that swaps the old solution to the right-hand side
   timeloop.add() << Sweep( Reinitialize( srcID, rhsID ), "Reinitialize" );

   // add the Jacobi iteration that computes the norm of the residual to stop the iteration
   timeloop.addFuncAfterTimeStep( JacobiIterationResidual( srcID, dstID, rhsID, weights, blocks, myCommScheme, uint_c(10000), real_c(1e-9) ), "JacobiIteration");

   // additional communication scheme that updates all ghost layers before writing the VTK output to yield correct visualization in ParaView
   //blockforest::UniformBufferedScheme< stencil::D2Q9 > vtkSyncScheme ( blocks );
   //vtkSyncScheme.addPackInfo( make_shared< field::communication::PackInfo<ScalarField>( srcID ) );
   //timeloop.addFuncAfterTimeStep( vtkSyncScheme , "SyncVTK" );

   // write VTK output, every timestep and with ghost layer
   timeloop.addFuncAfterTimeStep( field::createVTKOutput< ScalarField, float >( srcID, *blocks, "solution", uint_c(1), uint_c(1) ), "VTK" );

   // run the simulation for the specified number of time steps
   timeloop.run();

   return 0;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
