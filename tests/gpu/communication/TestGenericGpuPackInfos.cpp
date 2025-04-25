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
//! \file TestGenericGpuPackInfos.cpp
//! \ingroup gpu
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \brief Test the behaviour of generic GPU field pack infos for all available communication stencils
//
//======================================================================================================================

#include "blockforest/all.h"

#include "core/all.h"

#include "field/all.h"

#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/GPUField.h"
#include "gpu/communication/MemcpyPackInfo.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "stencil/all.h"

/**
 * This test suite tests the available generic GPU field pack infos available in waLBerla.
 * Each pack info is tested with all available communication stencils
 * and different numbers of ghost layers.
 */

namespace TestGenericPackInfos
{

using namespace walberla;

using CpuField_T = field::GhostLayerField< real_t, 3 >;
using GpuField_T = gpu::GPUField< real_t >;

template< typename Stencil, typename PackInfo = gpu::communication::MemcpyPackInfo< GpuField_T > >
class PackInfoTest
{
 public:
   using CommScheme = gpu::communication::UniformGPUScheme< Stencil >;

   PackInfoTest(const uint_t numGhostLayers)
   {
      auto mpiManager             = mpi::MPIManager::instance();
      uint_t numProcs             = uint_c(mpiManager->numProcesses());
      Vector3< uint_t > numBlocks = math::getFactors3D(numProcs);

      blocks_ = blockforest::createUniformBlockGrid( //
         numBlocks[0], numBlocks[1], numBlocks[2],   //
         4, 4, 4,                                    // cellsPerBlock
         1.0,                                        // dx
         true,                                       // oneBlockPerProcess
         true, true, true                             // periodic
      );

      cpuFieldId_ = field::addToStorage< CpuField_T >(blocks_, "cpuField", real_c(0.0), field::fzyx, numGhostLayers);
      gpuFieldId_ = gpu::addGPUFieldToStorage< CpuField_T >(blocks_, cpuFieldId_, "gpuField");

      commScheme_   = std::make_unique< CommScheme >(blocks_);
      auto packInfo = std::make_shared< PackInfo >(gpuFieldId_, numGhostLayers);
      commScheme_->addPackInfo(packInfo);
   }

   void run()
   {
      CellInterval allCells{ { 0, 0, 0 },
                             {
                                blocks_->getNumberOfXCellsPerBlock() - 1,
                                blocks_->getNumberOfYCellsPerBlock() - 1,
                                blocks_->getNumberOfZCellsPerBlock() - 1,
                             } };

      for (auto& b : *blocks_)
      {
         CpuField_T& cpuField = *b.getData< CpuField_T >(cpuFieldId_);
         for (Cell c : allCells)
         {
            Vector3< real_t > cc{ blocks_->getBlockLocalCellCenter(b, c) };
            cpuField.get(c, 0) = cc[0];
            cpuField.get(c, 1) = cc[1];
            cpuField.get(c, 2) = cc[2];
         }
      }

      gpu::fieldCpy< GpuField_T, CpuField_T >(blocks_, gpuFieldId_, cpuFieldId_);

      (*commScheme_)();

      gpu::fieldCpy< CpuField_T, GpuField_T >(blocks_, cpuFieldId_, gpuFieldId_);

      for (auto& b : *blocks_)
      {
         CpuField_T& cpuField = *b.getData< CpuField_T >(cpuFieldId_);
         for (auto dIt = Stencil::beginNoCenter(); dIt != Stencil::end(); ++dIt)
         {
            CellInterval glInterval;
            cpuField.getGhostRegion(*dIt, glInterval, 1);
            for (Cell c : glInterval)
            {
               Vector3< real_t > cc{ blocks_->getBlockLocalCellCenter(b, c) };
               blocks_->mapToPeriodicDomain( cc );
               WALBERLA_CHECK_FLOAT_EQUAL(cpuField.get(c, 0), cc[0]);
               WALBERLA_CHECK_FLOAT_EQUAL(cpuField.get(c, 1), cc[1]);
               WALBERLA_CHECK_FLOAT_EQUAL(cpuField.get(c, 2), cc[2]);
            }
         }
      }
   }

 private:
   std::shared_ptr< StructuredBlockForest > blocks_;
   BlockDataID cpuFieldId_;
   BlockDataID gpuFieldId_;
   std::unique_ptr< CommScheme > commScheme_;
};

void run()
{
   using namespace stencil;

   for(uint_t numGhostLayers = 1; numGhostLayers <= 3; ++numGhostLayers){
      PackInfoTest< D2Q4 >(numGhostLayers).run();
      PackInfoTest< D2Q5 >(numGhostLayers).run();
      PackInfoTest< D2Q9 >(numGhostLayers).run();

      PackInfoTest< D3Q6 >(numGhostLayers).run();
      PackInfoTest< D3Q7 >(numGhostLayers).run();
      PackInfoTest< D3Q15 >(numGhostLayers).run();
      PackInfoTest< D3Q19 >(numGhostLayers).run();
      PackInfoTest< D3Q27 >(numGhostLayers).run();
   }
}

} // namespace TestGenericPackInfos

int main(int argc, char** argv)
{
   walberla::mpi::Environment env{ argc, argv };

   TestGenericPackInfos::run();
}
