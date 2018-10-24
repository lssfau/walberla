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
//! \file MicroBenchmarkPdfCopy.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"
#include "blockforest/Initialization.h"

#include "field/Field.h"

#include "cuda/GPUField.h"
#include "cuda/FieldCopy.h"
#include "cuda/AddGPUFieldToStorage.h"

#include "MicroBenchmarkCopyKernel.h"
#include "MicroBenchmarkStreamKernel.h"


using namespace walberla;


int main( int argc, char **argv )
{
   debug::enterTestMode();
   mpi::Environment env( argc, argv );

   shared_ptr<StructuredBlockForest> blocks = blockforest::createUniformBlockGrid(1u, 1u, 1u,
           128u, 128u, 128u, 1.0, false, false, false, false);

   BlockDataID srcID = cuda::addGPUFieldToStorage<cuda::GPUField<double> >(blocks, "src", 19, field::fzyx, 1);
   BlockDataID dstID = cuda::addGPUFieldToStorage<cuda::GPUField<double> >(blocks, "dst", 19, field::fzyx, 1);

   int iterations = 3;

   pystencils::MicroBenchmarkCopyKernel copy(dstID, srcID);
   for( int i=0 ; i < iterations; ++i )
      for( auto &block: *blocks )
         copy( &block );


   pystencils::MicroBenchmarkStreamKernel stream(dstID, srcID);
   for( int i=0 ; i < iterations; ++i )
      for( auto &block: *blocks )
         stream( &block );

   WALBERLA_CUDA_CHECK(cudaDeviceSynchronize());

   return 0;
}
