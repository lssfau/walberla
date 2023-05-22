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
//! \ingroup gpu
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//! \brief Tests if a GPUField is correctly packed into buffers
//
//======================================================================================================================

#include "gpu/communication/GPUPackInfo.h"

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "field/GhostLayerField.h"

#include "stencil/D3Q27.h"

#include <cstring>
#include <vector>

#include "gpu/FieldCopy.h"
#include "gpu/GPUField.h"

#define F_SIZE 19

using namespace walberla;

static std::vector< field::Layout > fieldLayouts = { field::fzyx, field::zyxf };
static uint_t fieldLayoutIndex                   = 0;

gpu::GPUField< int >* createGPUField(IBlock* const block, StructuredBlockStorage* const storage)
{
   return new gpu::GPUField< int >(storage->getNumberOfXCells(*block), // number of cells in x direction
                                   storage->getNumberOfYCells(*block), // number of cells in y direction
                                   storage->getNumberOfZCells(*block), // number of cells in z direction
                                   F_SIZE,                             // fSize
                                   1,                                  // number of ghost layers
                                   fieldLayouts[fieldLayoutIndex]);
}

// Tester base class. The communicate() template method allows testing different communication methods.
class GPUPackInfoTester
{
 public:
   using GPUPackInfoType = gpu::communication::GPUPackInfo< gpu::GPUField< int > >;

   GPUPackInfoTester(IBlock* block, BlockDataID fieldId) : block_(block), fieldId_(fieldId) {}

   virtual ~GPUPackInfoTester() = default;

   void test(stencil::Direction dir)
   {
      gpu::GPUField< int >& gpuField = *(block_->getData< gpu::GPUField< int > >(fieldId_));

      field::GhostLayerField< int, F_SIZE > cpuField(gpuField.xSize(), // number of cells in x direction
                                                     gpuField.ySize(), // number of cells in y direction
                                                     gpuField.zSize(), // number of cells in z direction
                                                     1,                // number of ghost layers
                                                     0,                // initial value
                                                     fieldLayouts[fieldLayoutIndex]);
      cpuField.setWithGhostLayer(0);

      int val = 0;
      for (auto it = cpuField.beginSliceBeforeGhostLayer(dir); it != cpuField.end(); ++it)
      {
         *it = ++val;
      }
      gpu::fieldCpy(gpuField, cpuField);

      GPUPackInfoType gpuPackInfo(fieldId_);

      communicate(gpuPackInfo, dir);
      gpu::fieldCpy(cpuField, gpuField);

      val = 0;
      for (auto it = cpuField.beginGhostLayerOnly(stencil::inverseDir[dir]); it != cpuField.end(); ++it)
      {
         WALBERLA_CHECK_EQUAL(*it, ++val)
      }
   }

 protected:
   virtual void communicate(GPUPackInfoType& gpuPackInfo, stencil::Direction dir) = 0;

   IBlock* block_;
   BlockDataID fieldId_;
};

// Tester for buffer communication
class GPUPackInfoBufferTester : public GPUPackInfoTester
{
 public:
   GPUPackInfoBufferTester(IBlock* block, BlockDataID fieldId) : GPUPackInfoTester(block, fieldId) {}

 protected:
   void communicate(GPUPackInfoType& gpuPackInfo, stencil::Direction dir) override
   {
      mpi::GenericSendBuffer<> sendBuf;
      sendBuf.addDebugMarker("Be");
      gpuPackInfo.packData(block_, dir, sendBuf);
      sendBuf.addDebugMarker("Af");

      // Manually copy over the send to the receive buffer
      mpi::GenericRecvBuffer<> recvBuf;
      recvBuf.resize(sendBuf.size());
      memcpy(recvBuf.ptr(), sendBuf.ptr(), sendBuf.size() * sizeof(mpi::GenericSendBuffer<>::ElementType));

      recvBuf.readDebugMarker("Be");
      gpuPackInfo.unpackData(block_, stencil::inverseDir[dir], recvBuf);
      recvBuf.readDebugMarker("Af");
   }
};

// Tester for local communication
class GPUPackInfoLocalTester : public GPUPackInfoTester
{
 public:
   GPUPackInfoLocalTester(IBlock* block, BlockDataID fieldId) : GPUPackInfoTester(block, fieldId) {}

 protected:
   void communicate(GPUPackInfoType& gpuPackInfo, stencil::Direction dir) override
   {
      gpuPackInfo.communicateLocal(block_, block_, dir);
   }
};

int main(int argc, char** argv)
{
   using blockforest::createUniformBlockGrid;

   debug::enterTestMode();
   MPIManager::instance()->initializeMPI(&argc, &argv);

   for (; fieldLayoutIndex < fieldLayouts.size(); ++fieldLayoutIndex)
   {
      // Create BlockForest
      uint_t const processes = uint_c(MPIManager::instance()->numProcesses());
      auto blocks            = createUniformBlockGrid(processes, 1, 1,   // blocks
                                                      2, 2, 2,           // cells
                                                      1,                 // dx
                                                      false,             // one block per process
                                                      true, true, true); // periodicity

      BlockDataID const scalarGPUFieldId =
         blocks->addStructuredBlockData< gpu::GPUField< int > >(&createGPUField, "ScalarGPUField");

      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      {
         GPUPackInfoBufferTester bufferTester(&(*blockIt), scalarGPUFieldId);
         GPUPackInfoLocalTester localTester(&(*blockIt), scalarGPUFieldId);

         for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
         {
            localTester.test(*dir);
            bufferTester.test(*dir);
         }
      }
   }

   return EXIT_SUCCESS;
}
