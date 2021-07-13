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
//! \file GeneratedFieldPackInfoTestGPU.cpp
//! \ingroup field
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//! \brief Tests if a GPU Field is correctly communicated using generated pack info
//
//======================================================================================================================

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "cuda/FieldCopy.h"
#include "cuda/communication/UniformGPUScheme.h"

#include "stencil/D3Q27.h"

// include generated files
#include "ScalarFieldCommunicationGPU.h"
#include "ScalarFieldPullReductionGPU.h"


namespace walberla {

using Stencil_T = stencil::D3Q27;

cuda::GPUField<int> * createGPUField( IBlock* const block, StructuredBlockStorage* const storage ) {

   return new cuda::GPUField<int> (
      storage->getNumberOfXCells( *block ), // number of cells in x direction
      storage->getNumberOfYCells( *block ), // number of cells in y direction
      storage->getNumberOfZCells( *block ), // number of cells in z direction
      1,                                    // fSize
      1,                                    // number of ghost layers
      field::fzyx );

}

cuda::GPUField<int> * createSmallGPUField( IBlock * const , StructuredBlockStorage * const ) {
   return new cuda::GPUField<int> (2, 2, 2, 1, 1, field::fzyx );
}


void testScalarField( std::shared_ptr<blockforest::StructuredBlockForest> & sbf, BlockDataID gpuFieldId ) {

   cuda::communication::UniformGPUScheme< Stencil_T > us{ sbf };
   us.addPackInfo(std::make_shared< pystencils::ScalarFieldCommunicationGPU >(gpuFieldId));

   for( auto & block : *sbf ) {

      auto & gpuField = *(block.getData< cuda::GPUField< int > >(gpuFieldId));

      field::GhostLayerField< int, 1 > cpuField(gpuField.xSize(), gpuField.ySize(), gpuField.zSize(), 1, 0,
                                                field::fzyx);
      cpuField.setWithGhostLayer(0);

      WALBERLA_CHECK_EQUAL(cpuField.xSize(), 2)
      WALBERLA_CHECK_EQUAL(cpuField.ySize(), 2)
      WALBERLA_CHECK_EQUAL(cpuField.zSize(), 2)

      // initialize the bottom boundary
      cpuField(0, 0, 0) = 1;
      cpuField(0, 1, 0) = 2;
      cpuField(1, 0, 0) = 3;
      cpuField(1, 1, 0) = 4;

      cuda::fieldCpy(gpuField, cpuField);

      // communicate
      us.communicate();

      cuda::fieldCpy(cpuField, gpuField);

      WALBERLA_CHECK_EQUAL(cpuField(0, 0, +2), 1)
      WALBERLA_CHECK_EQUAL(cpuField(0, 1, +2), 2)
      WALBERLA_CHECK_EQUAL(cpuField(1, 0, +2), 3)
      WALBERLA_CHECK_EQUAL(cpuField(1, 1, +2), 4)

   }
}

void testScalarFieldPullReduction( std::shared_ptr<blockforest::StructuredBlockForest> & sbf, BlockDataID gpuFieldId ) {

   cuda::communication::UniformGPUScheme< Stencil_T > us1{ sbf };
   us1.addPackInfo(std::make_shared< pystencils::ScalarFieldPullReductionGPU >(gpuFieldId));

   cuda::communication::UniformGPUScheme< Stencil_T > us2{ sbf };
   us2.addPackInfo(std::make_shared< pystencils::ScalarFieldCommunicationGPU >(gpuFieldId));

   for( auto & block : *sbf ) {

      auto& gpuField = *(block.getData< cuda::GPUField< int > >(gpuFieldId));

      field::GhostLayerField< int, 1 > cpuField(gpuField.xSize(), gpuField.ySize(), gpuField.zSize(), 1, 0,
                                                field::fzyx);
      cpuField.setWithGhostLayer(0);

      WALBERLA_CHECK_EQUAL(cpuField.xSize(), 2)
      WALBERLA_CHECK_EQUAL(cpuField.ySize(), 2)
      WALBERLA_CHECK_EQUAL(cpuField.zSize(), 2)

      // initialize the bottom ghost layer cells
      cpuField(0, 0, -1) = 1;
      cpuField(0, 1, -1) = 2;
      cpuField(1, 0, -1) = 3;
      cpuField(1, 1, -1) = 4;

      // initialize the top interior cells
      cpuField(0, 0, 1) = 1;
      cpuField(0, 1, 1) = 1;
      cpuField(1, 0, 1) = 1;
      cpuField(1, 1, 1) = 1;

      cuda::fieldCpy(gpuField, cpuField);

      // communicate pull += reduction
      us1.communicate();

      cuda::fieldCpy(cpuField, gpuField);

      // check values in top ghost layer
      WALBERLA_CHECK_EQUAL(cpuField(0, 0, 2), 0)
      WALBERLA_CHECK_EQUAL(cpuField(0, 1, 2), 0)
      WALBERLA_CHECK_EQUAL(cpuField(1, 0, 2), 0)
      WALBERLA_CHECK_EQUAL(cpuField(1, 1, 2), 0)

      // check values in top interior cells
      WALBERLA_CHECK_EQUAL(cpuField(0, 0, 1), 2)
      WALBERLA_CHECK_EQUAL(cpuField(0, 1, 1), 3)
      WALBERLA_CHECK_EQUAL(cpuField(1, 0, 1), 4)
      WALBERLA_CHECK_EQUAL(cpuField(1, 1, 1), 5)

      // communicate to sync ghost layers
      us2.communicate();

      cuda::fieldCpy(cpuField, gpuField);

      // check values in bottom ghost layer
      WALBERLA_CHECK_EQUAL(cpuField(0, 0, -1), 2)
      WALBERLA_CHECK_EQUAL(cpuField(0, 1, -1), 3)
      WALBERLA_CHECK_EQUAL(cpuField(1, 0, -1), 4)
      WALBERLA_CHECK_EQUAL(cpuField(1, 1, -1), 5)

      // check values in top interior cells
      WALBERLA_CHECK_EQUAL(cpuField(0, 0, 1), 2)
      WALBERLA_CHECK_EQUAL(cpuField(0, 1, 1), 3)
      WALBERLA_CHECK_EQUAL(cpuField(1, 0, 1), 4)
      WALBERLA_CHECK_EQUAL(cpuField(1, 1, 1), 5)

   }
}

int main(int argc, char **argv) {

   using blockforest::createUniformBlockGrid;

   debug::enterTestMode();
   Environment walberlaEnv(argc,argv);

   // Create a BlockForest with 2x2x2 cells per block
   uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
   auto blocks = createUniformBlockGrid(processes,1 ,1, //blocks
                                        2,2,2, //cells
                                        1, //dx
                                        true, //one block per process
                                        true,true,true); //periodicity

   // Create a Field with the same number of cells as the block
   BlockDataID scalarGPUFieldId = blocks->addStructuredBlockData<cuda::GPUField<int> > ( &createGPUField, "ScalarGPUField" );

   testScalarField( blocks, scalarGPUFieldId );

   // Create a BlockForest with 8x8x8 cells per block
   blocks = createUniformBlockGrid(processes,1 ,1, //blocks
                                   8,8,8,          //cells
                                   1,              //dx
                                   true,          //one block per process
                                   true,true,true);//periodicity

   // Create a Field with one quarter as many cells per dimension, i.e. a field with the same size as the one above
   scalarGPUFieldId = blocks->addStructuredBlockData<cuda::GPUField<int> > ( &createSmallGPUField, "ScalarGPUField" );

   testScalarField( blocks, scalarGPUFieldId );

   testScalarFieldPullReduction( blocks, scalarGPUFieldId );

   return 0;

}

} // namespace walberla

int main( int argc, char* argv[] ) {
   return walberla::main( argc, argv );
}