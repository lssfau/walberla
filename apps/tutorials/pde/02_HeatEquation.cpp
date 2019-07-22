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
//! \file 06_HeatEquation.cpp
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

#include "gui/Gui.h"
#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <cmath>


namespace walberla {


typedef GhostLayerField<real_t,1> ScalarField;
using Stencil_T = stencil::D2Q5;



// function to initialize the field holding the unknown u
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



class JacobiIteration
{
public:

   //constructor
   JacobiIteration( const BlockDataID & srcID, const BlockDataID & dstID, const BlockDataID & rhsID, const std::vector< real_t > & weights,
                    const shared_ptr< StructuredBlockStorage > & blocks, blockforest::communication::UniformBufferedScheme<Stencil_T> & myCommScheme,
                    const uint_t & maxIterations )
   : srcID_( srcID ), dstID_( dstID ), rhsID_( rhsID ), weights_( weights ), blocks_( blocks ),
     myCommScheme_ ( myCommScheme ), maxIterations_( maxIterations ) {}

   // functor that carries out the full Jacobi iteration
   void operator()();

private:
   const BlockDataID srcID_;
   const BlockDataID dstID_;
   const BlockDataID rhsID_;

   std::vector< real_t > weights_;
   const shared_ptr< StructuredBlockStorage > blocks_;
   blockforest::communication::UniformBufferedScheme< Stencil_T > myCommScheme_;

   const uint_t maxIterations_;
};

void JacobiIteration::operator()()
{
   // Jacobi iteration
   for( uint_t i = 0; i < maxIterations_; ++i )
   {
      // communicate to update ghost layer cells
      myCommScheme_();

      // iterate all blocks
      for( auto block = blocks_->begin(); block != blocks_->end(); ++block )
      {
         // get the source, destination, and rhs field data from the block
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
   }
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
             false,                              // one block per process - "false" means all blocks to one process
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

   // set up the time loop object
   SweepTimeloop timeloop ( blocks, uint_c(20) );

   // add the routine that swaps the old solution to the right-hand side
   timeloop.add() << Sweep( Reinitialize( srcID, rhsID ), "Reinitialize" );

   // add the complete Jacobi iteration with 10000 iterations
   // This can not be done as a sweep since it includes an interior iteration, independent of the time loop.
   timeloop.addFuncAfterTimeStep( JacobiIteration( srcID, dstID, rhsID, weights, blocks, myCommScheme, uint_c(10000) ), "JacobiIteration");

   // start the GUI and run the simulation
   GUI gui ( timeloop, blocks, argc, argv );
   gui.run();

   return 0;
}
}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
}
