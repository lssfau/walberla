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
//! \file FieldTransferTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "field/GhostLayerField.h"

#include "gpu/GPUField.h"
#include "gpu/FieldIndexing.h"
#include "gpu/FieldCopy.h"
#include "gpu/Kernel.h"

using namespace walberla;

namespace walberla{

void kernel_double(gpu::FieldAccessor<double> f );
}

GhostLayerField<double,1> * createCPUField( IBlock* const block, StructuredBlockStorage* const storage )
{
   return new GhostLayerField<double,1> (
      storage->getNumberOfXCells( *block ), // number of cells in x direction
      storage->getNumberOfYCells( *block ), // number of cells in y direction
      storage->getNumberOfZCells( *block ), // number of cells in z direction
      1,                                    // number of ghost layers
      double(1),                            // initial value
      field::fzyx);
}

gpu::GPUField<double> * createGPUField( IBlock* const block, StructuredBlockStorage* const storage )
{
   return new gpu::GPUField<double> (
      storage->getNumberOfXCells( *block ), // number of cells in x direction
      storage->getNumberOfYCells( *block ), // number of cells in y direction
      storage->getNumberOfZCells( *block ), // number of cells in z direction
      1,                                    // fSize
      1,                                    // number of ghost layers
      field::fzyx );
}


int main( int argc, char ** argv )
{
   walberla::Environment const env( argc, argv );
   debug::enterTestMode();

   shared_ptr< StructuredBlockForest > const blocks = blockforest::createUniformBlockGrid (
      uint_t(1),   uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
      uint_t(14),  uint_t(14), uint_t(14), // how many cells per block (x,y,z)
      real_c(0.5),                         // dx: length of one cell in physical coordinates
      false,                               // one block per process - "false" means all blocks to one process
      false, false, false );               // no periodicity



   BlockDataID const cpuFieldID = blocks->addStructuredBlockData< GhostLayerField<double,1> > ( &createCPUField, "CPUField" );


   BlockDataID const gpuFieldID = blocks->addStructuredBlockData< gpu::GPUField<double>    > ( &createGPUField, "GPUField" );

   for ( auto blockIterator = blocks->begin(); blockIterator != blocks->end(); ++blockIterator )
   {
      IBlock & currentBlock = *blockIterator;

      // get the field stored on the current block
      auto cpuField = currentBlock.getData< GhostLayerField<double,1> > ( cpuFieldID );
      auto gpuField = currentBlock.getData< gpu::GPUField<double>    > ( gpuFieldID );

      gpu::fieldCpy( *gpuField, *cpuField );

      auto myKernel = gpu::make_kernel( &kernel_double );
      auto indexing = gpu::FieldIndexing<double>::sliceBeforeGhostLayerXYZ( *gpuField, 1, stencil::W, true );
      myKernel.addFieldIndexingParam(indexing);
      myKernel();

      gpu::fieldCpy( *cpuField, *gpuField );

      WALBERLA_ASSERT_FLOAT_EQUAL( cpuField->get(0,0,0), real_t(2) )
   }

   return EXIT_SUCCESS;
}
