//========================================================================================================================
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
//! \file GPUFieldPackInfoTest.cpp
//! \ingroup cuda
//! \author João Victor Tozatti Risso <jvtrisso@inf.ufpr.br>
//! \brief Short communication test to verify the equivalence of GPUPackInfo using a default stream and multiple streams.
//
//========================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/mpi/Environment.h"

#include "stencil/Directions.h"
#include "stencil/Iterator.h"
#include "stencil/D3Q27.h"

#include "domain_decomposition/BlockDataID.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "field/GhostLayerField.h"

#include "cuda/ErrorChecking.h"
#include "cuda/HostFieldAllocator.h"
#include "cuda/GPUField.h"
#include "cuda/FieldCopy.h"
#include "cuda/communication/GPUPackInfo.h"

#include <cuda_runtime.h>
#include <vector>

using namespace walberla;

using DataType = walberla::uint_t;
using StencilType = stencil::D3Q27;
using FieldType = field::GhostLayerField< DataType, StencilType::Size >;
using GPUFieldType = cuda::GPUField< DataType >;
using CommSchemeType = blockforest::communication::UniformBufferedScheme<StencilType>;
using GPUPackInfoType = cuda::communication::GPUPackInfo< GPUFieldType >;

static std::vector< cuda::Layout > fieldLayouts = { cuda::fzyx, cuda::zyxf };
static uint_t fieldLayoutIndex = 0;


FieldType * createField( IBlock* const block, StructuredBlockStorage* const storage )
{
   return new FieldType(
            storage->getNumberOfXCells( *block ),   // number of cells in x direction per block
            storage->getNumberOfYCells( *block ),   // number of cells in y direction per block
            storage->getNumberOfZCells( *block ),   // number of cells in z direction per block
            1,                                      // one ghost layer
            DataType(0),                            // initial value
            fieldLayouts[fieldLayoutIndex],         // layout
            make_shared<cuda::HostFieldAllocator< DataType > >() // allocator for host pinned memory
            );
}


GPUFieldType * createGPUField( IBlock* const block, StructuredBlockStorage* const storage )
{
   return new GPUFieldType(
            storage->getNumberOfXCells( *block ), // number of cells in x direction
            storage->getNumberOfYCells( *block ), // number of cells in y direction
            storage->getNumberOfZCells( *block ), // number of cells in z direction
            StencilType::Size,                    // number of cells for pdfs
            1,                                    // one ghost layer
            fieldLayouts[fieldLayoutIndex] );
}


void initFields( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fieldID )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      auto fieldPtr = block->getData< FieldType >( fieldID );

      for( auto fieldIt = fieldPtr->begin(); fieldIt != fieldPtr->end(); ++fieldIt )
         *fieldIt = math::intRandom< DataType >();
   }
}


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );


   const Vector3< uint_t > cells = Vector3< uint_t >( 4, 4, 4 );

   uint_t nProc = uint_c( MPIManager::instance()->numProcesses() );

   for(; fieldLayoutIndex < fieldLayouts.size(); ++fieldLayoutIndex )
   {
      auto blocks = blockforest::createUniformBlockGrid(nProc, 1, 1,                  // blocks
                                                        cells[0], cells[1], cells[2], // cells
                                                        1,                            // unit cell spacing
                                                        true,                        // one block per process
                                                        true, true, true);            // periodic in all directions

      BlockDataID sourceFieldId = blocks->addStructuredBlockData< FieldType >( &createField,
                                                                               "ScalarField" );

      BlockDataID syncGPUFieldId = blocks->addStructuredBlockData< GPUFieldType >( &createGPUField,
                                                                                   "syncGPUField" );

      BlockDataID asyncGPUFieldId = blocks->addStructuredBlockData< GPUFieldType >( &createGPUField,
                                                                                    "asyncGPUField" );

      math::seedRandomGenerator( numeric_cast<std::mt19937::result_type>( MPIManager::instance()->rank() ) );
      // Initialize CPU field with random values
      initFields( blocks, sourceFieldId );

      // Copy same CPU field to both GPU fields
      for( auto block = blocks->begin(); block != blocks->end(); ++block )
      {
         auto sourceFieldPtr = block->getData< FieldType >( sourceFieldId );

         auto syncGPUFieldPtr = block->getData< GPUFieldType >( syncGPUFieldId );
         cuda::fieldCpy( *syncGPUFieldPtr, *sourceFieldPtr );

         auto asyncGPUFieldPtr = block->getData< GPUFieldType >( asyncGPUFieldId );
         cuda::fieldCpy( *asyncGPUFieldPtr, *sourceFieldPtr );
      }

      // Setup communication schemes for synchronous GPUPackInfo
      CommSchemeType syncCommScheme(blocks);
      syncCommScheme.addPackInfo( make_shared< GPUPackInfoType >( syncGPUFieldId ) );

      // Setup communication scheme for asynchronous GPUPackInfo, which uses CUDA streams
      CommSchemeType asyncCommScheme(blocks);
      asyncCommScheme.addPackInfo( make_shared< GPUPackInfoType >( asyncGPUFieldId ) );

      // Perform one communication step for each scheme
      syncCommScheme();
      asyncCommScheme();

      // Check results
      FieldType syncFieldCpu( cells[0], cells[1], cells[2], 1, fieldLayouts[fieldLayoutIndex],
                              make_shared< cuda::HostFieldAllocator< DataType > >() );
      FieldType asyncFieldCpu( cells[0], cells[1], cells[2], 1, fieldLayouts[fieldLayoutIndex],
                               make_shared< cuda::HostFieldAllocator< DataType > >() );

      for( auto block = blocks->begin(); block != blocks->end(); ++block )
      {
         auto syncGPUFieldPtr = block->getData< GPUFieldType >( syncGPUFieldId );
         cuda::fieldCpy( syncFieldCpu, *syncGPUFieldPtr );

         auto asyncGPUFieldPtr = block->getData< GPUFieldType >( asyncGPUFieldId );
         cuda::fieldCpy( asyncFieldCpu, *asyncGPUFieldPtr );

         for( auto syncIt = syncFieldCpu.beginWithGhostLayerXYZ(), asyncIt = asyncFieldCpu.beginWithGhostLayerXYZ();
                  syncIt != syncFieldCpu.end();
                  ++syncIt, ++asyncIt )
            WALBERLA_CHECK_EQUAL( *syncIt, *asyncIt );
      }
   }


   return EXIT_SUCCESS;
}
