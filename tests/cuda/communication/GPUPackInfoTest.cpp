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
//! \file GPUFieldPackInfoTest.cpp
//! \ingroup cuda
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \brief Tests if a GPUField is correctly packed into buffers
//
//======================================================================================================================

#include "field/GhostLayerField.h"

#include "cuda/GPUField.h"
#include "cuda/FieldCopy.h"
#include "cuda/communication/GPUPackInfo.h"

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "stencil/D3Q27.h"

#include <cstring>
#include <vector>
#include <cuda_runtime.h>

#define F_SIZE    19

using namespace walberla;

static std::vector< field::Layout > fieldLayouts = { field::fzyx, field::zyxf };
static uint_t fieldLayoutIndex = 0;

cuda::GPUField<int> * createGPUField( IBlock* const block, StructuredBlockStorage* const storage )
{
   return new cuda::GPUField<int> (
            storage->getNumberOfXCells( *block ), // number of cells in x direction
            storage->getNumberOfYCells( *block ), // number of cells in y direction
            storage->getNumberOfZCells( *block ), // number of cells in z direction
            F_SIZE,                               // fSize
            1,                                    // number of ghost layers
            fieldLayouts[fieldLayoutIndex] );
}

// Tester base class. The communicate() template method allows testing different communication methods.
class GPUPackInfoTester
{
public:

   typedef cuda::communication::GPUPackInfo< cuda::GPUField<int> > GPUPackInfoType;

   GPUPackInfoTester( IBlock* block, BlockDataID fieldId ) :
      block_( block ), fieldId_( fieldId ) {}

   virtual ~GPUPackInfoTester() {}

   void test( stencil::Direction dir )
   {
      cuda::GPUField<int> & gpuField = *(block_->getData<cuda::GPUField<int> >( fieldId_ ));

      field::GhostLayerField<int,F_SIZE> cpuField(
               gpuField.xSize(),       // number of cells in x direction
               gpuField.ySize(),       // number of cells in y direction
               gpuField.zSize(),       // number of cells in z direction
               1,                      // number of ghost layers
               0,                      // initial value
               fieldLayouts[fieldLayoutIndex]);
      cpuField.setWithGhostLayer( 0 );

      int val = 0;
      for ( auto it = cpuField.beginSliceBeforeGhostLayer( dir ); it != cpuField.end(); ++it )
      {
         *it = ++val;
      }
      cuda::fieldCpy( gpuField, cpuField );

      GPUPackInfoType gpuPackInfo( fieldId_ );

      communicate( gpuPackInfo, dir );
      cuda::fieldCpy( cpuField, gpuField );

      val = 0;
      for ( auto it = cpuField.beginGhostLayerOnly( stencil::inverseDir[dir] ); it != cpuField.end(); ++it )
      {
         WALBERLA_CHECK_EQUAL( *it, ++val );
      }

   }

protected:

   virtual void communicate( GPUPackInfoType& gpuPackInfo, stencil::Direction dir ) = 0;

   IBlock* block_;
   BlockDataID fieldId_;
};


// Tester for buffer communication
class GPUPackInfoBufferTester: public GPUPackInfoTester
{
public:
   GPUPackInfoBufferTester( IBlock* block, BlockDataID fieldId): GPUPackInfoTester( block, fieldId ) {}

protected:
   void communicate( GPUPackInfoType& gpuPackInfo, stencil::Direction dir )
   {
      mpi::GenericSendBuffer<> sendBuf;
      sendBuf.addDebugMarker( "Be" );
      gpuPackInfo.packData( block_, dir, sendBuf );
      sendBuf.addDebugMarker( "Af" );

      // Manually copy over the send to the receive buffer
      mpi::GenericRecvBuffer<> recvBuf;
      recvBuf.resize( sendBuf.size() );
      memcpy( recvBuf.ptr(), sendBuf.ptr(), sendBuf.size() * sizeof(mpi::GenericSendBuffer<>::ElementType) );

      recvBuf.readDebugMarker( "Be" );
      gpuPackInfo.unpackData( block_,  stencil::inverseDir[dir], recvBuf );
      recvBuf.readDebugMarker( "Af" );
   }
};


// Tester for local communication
class GPUPackInfoLocalTester: public GPUPackInfoTester
{
public:
   GPUPackInfoLocalTester( IBlock* block, BlockDataID fieldId ): GPUPackInfoTester( block, fieldId ) {}

protected:
   void communicate( GPUPackInfoType& gpuPackInfo, stencil::Direction dir )
   {
      gpuPackInfo.communicateLocal( block_, block_, dir );
   }
};


int main(int argc, char **argv)
{
   using blockforest::createUniformBlockGrid;

   debug::enterTestMode();
   MPIManager::instance()->initializeMPI(&argc,&argv);

   for(; fieldLayoutIndex < fieldLayouts.size(); ++fieldLayoutIndex )
   {
      // Create BlockForest
      uint_t processes = uint_c( MPIManager::instance()->numProcesses() );
      auto blocks = createUniformBlockGrid(processes,1,1,  //blocks
                                           2,2,2,          //cells
                                           1,              //dx
                                           false,          //one block per process
                                           true,true,true);//periodicity

      BlockDataID scalarGPUFieldId = blocks->addStructuredBlockData<cuda::GPUField<int> >(
              &createGPUField, "ScalarGPUField" );

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         GPUPackInfoBufferTester bufferTester( &(*blockIt), scalarGPUFieldId );
         GPUPackInfoLocalTester localTester( &(*blockIt), scalarGPUFieldId );

         for( auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir )
         {
            localTester.test( *dir );
            bufferTester.test( *dir );
         }
      }
   }

   return 0;
}
