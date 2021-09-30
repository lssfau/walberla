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
//! \file JacobiGpu.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "CudaJacobiKernel2D.h"
#include "CudaJacobiKernel3D.h"

#include "cuda/HostFieldAllocator.h"
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "cuda/FieldCopy.h"
#include "cuda/GPUField.h"
#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/communication/GPUPackInfo.h"
#include "cuda/FieldIndexing.h"

#include "field/AddToStorage.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/initializer/ScalarFieldFromGrayScaleImage.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q7.h"

#include "timeloop/SweepTimeloop.h"


using namespace walberla;

typedef GhostLayerField<double,1> ScalarField;
typedef cuda::GPUField<double> GPUField;


ScalarField * createField( IBlock* const block, StructuredBlockStorage* const storage )
{
   return new ScalarField (
            storage->getNumberOfXCells( *block ),   // number of cells in x direction per block
            storage->getNumberOfYCells( *block ),   // number of cells in y direction per block
            storage->getNumberOfZCells( *block ),   // number of cells in z direction per block
            1,                                      // one ghost layer
            double(0),                              // initial value
            field::fzyx,                            // layout
            make_shared<cuda::HostFieldAllocator<double> >()  // allocator for host pinned memory
            );
}

void testJacobi2D()
{
   uint_t xSize = 20;
   uint_t ySize = 20;

   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
           uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
           xSize, ySize, uint_t(1),            // how many cells per block (x,y,z)
           real_t(1),                          // dx: length of one cell in physical coordinates
           false,                              // one block per process - "false" means all blocks to one process
           true, true, true );                 // no periodicity


   BlockDataID cpuFieldID = blocks->addStructuredBlockData<ScalarField>( &createField, "CPU Field" );
   BlockDataID gpuField = cuda::addGPUFieldToStorage<ScalarField>( blocks, cpuFieldID, "GPU Field Src" );

   // Initialize a quarter of the field with ones, the rest remains 0
   // Jacobi averages the domain -> every cell should be at 0.25 at sufficiently many timesteps
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto f = blockIt->getData<ScalarField>( cpuFieldID );
      for( cell_idx_t y = 0; y < cell_idx_c( f->ySize() / 2 ); ++y )
         for( cell_idx_t x = 0; x < cell_idx_c( f->xSize() / 2 ); ++x )
            f->get( x, y, 0 ) = 1.0;
   }

   typedef blockforest::communication::UniformBufferedScheme<stencil::D2Q9> CommScheme;
   typedef cuda::communication::GPUPackInfo<GPUField> Packing;

   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(gpuField) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(800);
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( pystencils::CudaJacobiKernel2D(gpuField), "Jacobi Kernel" );


   cuda::fieldCpy<GPUField, ScalarField>( blocks, gpuField, cpuFieldID );
   timeloop.run();
   cuda::fieldCpy<ScalarField, GPUField>( blocks, cpuFieldID, gpuField );

   auto firstBlock = blocks->begin();
   auto f = firstBlock->getData<ScalarField>( cpuFieldID );
   WALBERLA_CHECK_FLOAT_EQUAL(f->get(0,0,0), real_t(1.0 / 4.0));
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
           real_t(1),                          // dx: length of one cell in physical coordinates
           false,                              // one block per process - "false" means all blocks to one process
           true, true, true );                 // no periodicity


   BlockDataID cpuFieldID = blocks->addStructuredBlockData<ScalarField>( &createField, "CPU Field" );
   BlockDataID gpuField = cuda::addGPUFieldToStorage<ScalarField>( blocks, cpuFieldID, "GPU Field Src" );

   // Initialize a quarter of the field with ones, the rest remains 0
   // Jacobi averages the domain -> every cell should be at 0.25 at sufficiently many timesteps
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto f = blockIt->getData<ScalarField>( cpuFieldID );
      for( cell_idx_t z = 0; z < cell_idx_c( f->zSize() / 2 ); ++z )
         for( cell_idx_t y = 0; y < cell_idx_c( f->ySize() / 2 ); ++y )
            for( cell_idx_t x = 0; x < cell_idx_c( f->xSize() / 2 ); ++x )
               f->get( x, y, z ) = 1.0;
   }

   typedef blockforest::communication::UniformBufferedScheme<stencil::D3Q7> CommScheme;
   typedef cuda::communication::GPUPackInfo<GPUField> Packing;

   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(gpuField) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(800);
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( pystencils::CudaJacobiKernel3D(gpuField), "Jacobi Kernel" );


   cuda::fieldCpy<GPUField, ScalarField>( blocks, gpuField, cpuFieldID );
   timeloop.run();
   cuda::fieldCpy<ScalarField, GPUField>( blocks, cpuFieldID, gpuField );

   auto firstBlock = blocks->begin();
   auto f = firstBlock->getData<ScalarField>( cpuFieldID );
   WALBERLA_CHECK_FLOAT_EQUAL(f->get(0,0,0), real_t(1.0 / 8.0));
}

int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testJacobi2D();
   testJacobi3D();

   return 0;
}
