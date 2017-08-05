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
//! \file 03_GameOfLife.cpp
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "01_GameOfLife_kernels.h"
#include "cuda/HostFieldAllocator.h"
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"

#include "cuda/HostFieldAllocator.h"
#include "cuda/FieldCopy.h"
#include "cuda/GPUField.h"
#include "cuda/Kernel.h"
#include "cuda/AddGPUFieldToStorage.h"
#include "cuda/communication/GPUPackInfo.h"
#include "cuda/FieldIndexing.h"

#include "field/AddToStorage.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/initializer/ScalarFieldFromGrayScaleImage.h"
#include "geometry/structured/GrayScaleImage.h"

#include "gui/Gui.h"

#include "stencil/D2Q9.h"

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

class GameOfLifeSweepCUDA
{
   public:
      GameOfLifeSweepCUDA( BlockDataID gpuFieldSrcID, BlockDataID gpuFieldDstID )
         : gpuFieldSrcID_( gpuFieldSrcID ), gpuFieldDstID_( gpuFieldDstID )
      {
      }
      void operator() ( IBlock * block )
      {
         auto srcCudaField = block->getData< cuda::GPUField<double> > ( gpuFieldSrcID_ );
         auto dstCudaField = block->getData< cuda::GPUField<double> > ( gpuFieldDstID_ );

         auto myKernel = cuda::make_kernel( &gameOfLifeKernel );
         myKernel.addFieldIndexingParam( cuda::FieldIndexing<double>::xyz( *srcCudaField ) );
         myKernel.addFieldIndexingParam( cuda::FieldIndexing<double>::xyz( *dstCudaField ) );
         myKernel();

         srcCudaField->swapDataPointers( dstCudaField );
      }
   private:
      BlockDataID gpuFieldSrcID_;
      BlockDataID gpuFieldDstID_;
};


int main( int argc, char ** argv )
{
   walberla::Environment env( argc, argv );

   geometry::GrayScaleImage image ("GosperGliderGun.png");

   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
            uint_t(1) ,              uint_t(2),                           uint_t(1), // number of blocks in x,y,z direction
            image.size( uint_t(0) ), image.size( uint_t(1) ) / uint_t(2), uint_t(1), // how many cells per block (x,y,z)
            real_t(1),                                                               // dx: length of one cell in physical coordinates
            false,                                                                   // one block per process - "false" means all blocks to one process
            false, false, false );                                                   // no periodicity


   BlockDataID cpuFieldID = blocks->addStructuredBlockData<ScalarField>( &createField, "CPU Field" );

   // Initializing the field from an image
   using geometry::initializer::ScalarFieldFromGrayScaleImage;
   ScalarFieldFromGrayScaleImage fieldInitializer ( *blocks, cpuFieldID ) ;
   fieldInitializer.init( image, uint_t(2), false );

   BlockDataID gpuFieldSrcID = cuda::addGPUFieldToStorage<ScalarField>( blocks, cpuFieldID, "GPU Field Src" );
   BlockDataID gpuFieldDstID = cuda::addGPUFieldToStorage<ScalarField>( blocks, cpuFieldID, "GPU Field Dst" );


   typedef blockforest::communication::UniformBufferedScheme<stencil::D2Q9 > CommScheme;
   typedef cuda::communication::GPUPackInfo<GPUField> Packing;
   // Alternative, if CUDA enabled MPI is available
   //blockforest::communication::UniformDirectScheme<stencil::D2Q9 >
   //typedef field::communication::UniformMPIDatatypeInfo<GPUField> Packing

   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(gpuFieldSrcID) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(100); // number of timesteps for non-gui runs
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( GameOfLifeSweepCUDA(gpuFieldSrcID, gpuFieldDstID ), "GameOfLifeSweep" );

   timeloop.add() << Sweep( cuda::fieldCpyFunctor<ScalarField, GPUField >(cpuFieldID, gpuFieldDstID) );

   // Register VTK output
   timeloop.addFuncAfterTimeStep( field::createVTKOutput<ScalarField>( cpuFieldID, *blocks, "game_of_life" ) );
   
   // GUI output
   GUI gui ( timeloop, blocks, argc, argv );
   gui.run();

   return 0;
}
