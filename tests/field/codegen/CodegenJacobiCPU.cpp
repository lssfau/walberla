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
//! \file CodegenJacobiGPU.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "JacobiKernel2D.h"
#include "JacobiKernel3D.h"
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q7.h"

#include "timeloop/SweepTimeloop.h"


using namespace walberla;

using ScalarField = GhostLayerField<real_t, 1>;


void testJacobi2D()
{
   uint_t xSize = 20;
   uint_t ySize = 20;
   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
           uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
           xSize, ySize, uint_t(1),            // how many cells per block (x,y,z)
           real_c(1.0),                          // dx: length of one cell in physical coordinates
           false,                              // one block per process - "false" means all blocks to one process
           true, true, true );                 // full periodicity


   BlockDataID fieldID = field::addToStorage<ScalarField>(blocks, "Field", real_c(0.0));

   // Initialize a quarter of the field with ones, the rest remains 0
   // Jacobi averages the domain -> every cell should be at 0.25 at sufficiently many timesteps
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto f = blockIt->getData<ScalarField>( fieldID );
      for( cell_idx_t y = 0; y < cell_idx_c( f->ySize() / 2 ); ++y )
         for( cell_idx_t x = 0; x < cell_idx_c( f->xSize() / 2 ); ++x )
            f->get( x, y, 0 ) = real_c(1.0);
   }

   using CommScheme = blockforest::communication::UniformBufferedScheme<stencil::D2Q9>;
   using Packing = field::communication::PackInfo<ScalarField>;
   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(fieldID) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(800);
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( pystencils::JacobiKernel2D(fieldID), "Jacobi Kernel" );

   timeloop.run();

   auto firstBlock = blocks->begin();
   auto f = firstBlock->getData<ScalarField>( fieldID );

   WALBERLA_CHECK_FLOAT_EQUAL(f->get(0,0,0), real_c(1.0 / 4.0))
}


void testJacobi3D()
{
   uint_t xSize = 12;
   uint_t ySize = 12;
   uint_t zSize = 12;
   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
           uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
           xSize, ySize, zSize,                // how many cells per block (x,y,z)
           real_c(1.0),                          // dx: length of one cell in physical coordinates
           false,                              // one block per process - "false" means all blocks to one process
           true, true, true );                 // no periodicity


   BlockDataID fieldID = field::addToStorage<ScalarField>(blocks, "Field", real_c(0.0));

   // Initialize a quarter of the field with ones, the rest remains 0
   // Jacobi averages the domain -> every cell should be at 0.25 at sufficiently many timesteps
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto f = blockIt->getData<ScalarField>( fieldID );
      for( cell_idx_t z = 0; z < cell_idx_c( f->zSize() / 2); ++z )
         for( cell_idx_t y = 0; y < cell_idx_c( f->ySize() / 2 ); ++y )
            for( cell_idx_t x = 0; x < cell_idx_c( f->xSize() / 2 ); ++x )
               f->get( x, y, z ) = real_c(1.0);
   }

   using CommScheme = blockforest::communication::UniformBufferedScheme<stencil::D3Q7>;
   using Packing = field::communication::PackInfo<ScalarField>;
   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(fieldID) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(800); // number of timesteps
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( pystencils::JacobiKernel3D(fieldID), "Jacobi Kernel" );

   timeloop.run();

   auto firstBlock = blocks->begin();
   auto f = firstBlock->getData<ScalarField>( fieldID );
   WALBERLA_CHECK_FLOAT_EQUAL(f->get(0,0,0), real_c(1.0 / 8.0))
}


int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testJacobi2D();
   testJacobi3D();

   return EXIT_SUCCESS;
}
