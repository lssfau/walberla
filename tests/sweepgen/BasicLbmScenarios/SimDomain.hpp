#pragma once

#include "blockforest/all.h"

#include "core/all.h"

#include "field/all.h"
#include "field/communication/StencilRestrictedPackInfo.h"

#include "stencil/Directions.h"

#include <memory>

#include "gen/LbmAlgorithms.hpp"
#include "walberla/experimental/Sweep.hpp"

#if defined(LBM_SCENARIOS_GPU_BUILD)
#   include "gpu/AddGPUFieldToStorage.h"
#   include "gpu/FieldCopy.h"
#   include "gpu/GPUField.h"
#   include "gpu/GPUWrapper.h"
#   include "gpu/HostFieldAllocator.h"
#   include "gpu/communication/GPUPackInfo.h"
#   include "gpu/communication/MemcpyPackInfo.h"
#   include "gpu/communication/UniformGPUScheme.h"
#endif

namespace BasicLbmScenarios
{

using namespace walberla;
using namespace walberla::experimental;

using PdfField_T    = field::GhostLayerField< real_t, gen::LbStencil::Q >;
using ScalarField_T = field::GhostLayerField< real_t, 1 >;
using VectorField_T = field::GhostLayerField< real_t, gen::LbStencil::D >;
using FlagField_T   = FlagField< uint8_t >;

using CpuCommScheme   = blockforest::communication::UniformBufferedScheme< gen::LbStencil >;
using CpuPdfsPackInfo = field::communication::StencilRestrictedPackInfo< PdfField_T, gen::LbStencil >;

#if defined(LBM_SCENARIOS_GPU_BUILD)
using CommonGpuField = gpu::GPUField< PdfField_T::value_type >;

using GpuCommScheme = gpu::communication::UniformGPUScheme< gen::LbStencil >;
// using GpuPdfsPackInfo = gpu::communication::MemcpyPackInfo< CommonGpuField >;
using GpuPdfsPackInfo = gen::comm::GpuPdfsPackInfo;
#endif

struct SimDomain
{
   const std::shared_ptr< StructuredBlockForest > blocks;

   struct
   {
      const BlockDataID pdfsId;
      const BlockDataID rhoId;
      const BlockDataID uId;
      const BlockDataID flagFieldId;
   } cpuFields;

   CpuCommScheme commCpu;

#if defined(LBM_SCENARIOS_GPU_BUILD)
   struct
   {
      const BlockDataID pdfsId;
      const BlockDataID rhoId;
      const BlockDataID uId;
   } gpuFields;

   std::unique_ptr< GpuCommScheme > commGpu;

   void initFromFields(const Vector3< real_t > force)
   {
      gen::bulk::LbInitFromFields initialize{ gpuFields.pdfsId, gpuFields.rhoId, gpuFields.uId, force };

      for (auto& b : *blocks)
      {
         initialize(&b);
      }

      wait();
   }

   void initConstant(const real_t rho, const Vector3< real_t > u, const Vector3< real_t > force)
   {
      gen::bulk::LbInitConstant initialize{ gpuFields.pdfsId, force, rho, u };

      for (auto& b : *blocks)
      {
         initialize(&b);
      }

      wait();
   }

   gen::bulk::LbStreamCollide streamCollideSweep(const real_t omega, const Vector3< real_t > force)
   {
      return { gpuFields.pdfsId, gpuFields.rhoId, gpuFields.uId, force, omega };
   }

   auto freeSlipBottom()
   {
      return sweep::SweepFactory{ blocks }.atDomainBorder< stencil::Direction::B >(
         gen::bc_grid_aligned::FreeSlipBottom{ gpuFields.pdfsId });
   }

   auto noSlipTop()
   {
      return sweep::SweepFactory{ blocks }.atDomainBorder< stencil::Direction::T >(gen::bc_grid_aligned::NoSlipTop{ gpuFields.pdfsId });
   }

   auto irregularFreeSlipFactory() { return gen::bc_sparse::FreeSlipIrregularFactory(blocks, gpuFields.pdfsId); }

   void wait() { WALBERLA_GPU_CHECK(gpuDeviceSynchronize()); }

   void syncGhostLayers()
   {
      // WALBERLA_GPU_CHECK(gpuPeekAtLastError());
      (*commGpu)();
   }

   void fields2host()
   {
      wait();
      gpu::fieldCpy< PdfField_T, CommonGpuField >(blocks, cpuFields.pdfsId, gpuFields.pdfsId);
      gpu::fieldCpy< ScalarField_T, CommonGpuField >(blocks, cpuFields.rhoId, gpuFields.rhoId);
      gpu::fieldCpy< VectorField_T, CommonGpuField >(blocks, cpuFields.uId, gpuFields.uId);
   }

   void fields2device()
   {
      wait();
      gpu::fieldCpy< CommonGpuField, PdfField_T >(blocks, gpuFields.pdfsId, cpuFields.pdfsId);
      gpu::fieldCpy< CommonGpuField, ScalarField_T >(blocks, gpuFields.rhoId, cpuFields.rhoId);
      gpu::fieldCpy< CommonGpuField, VectorField_T >(blocks, gpuFields.uId, cpuFields.uId);
   }

#else

   void initFromFields(const Vector3< real_t > force)
   {
      gen::bulk::LbInitFromFields initialize{ cpuFields.pdfsId, cpuFields.rhoId, cpuFields.uId, force };

      for (auto& b : *blocks)
      {
         initialize(&b);
      }
   }

   void initConstant(const real_t rho, const Vector3< real_t > u, const Vector3< real_t > force)
   {
      gen::bulk::LbInitConstant initialize{ cpuFields.pdfsId, force, rho, u };

      for (auto& b : *blocks)
      {
         initialize(&b);
      }
   }

   gen::bulk::LbStreamCollide streamCollideSweep(const real_t omega, const Vector3< real_t > force)
   {
      return { cpuFields.pdfsId, cpuFields.rhoId, cpuFields.uId, force, omega };
   }

   auto freeSlipTop()
   {
      return sweep::BorderSweep< stencil::Direction::T, gen::bc_grid_aligned::FreeSlipTop >{
         blocks, gen::bc_grid_aligned::FreeSlipTop{ cpuFields.pdfsId }
      };
   }

   auto freeSlipBottom()
   {
      return sweep::BorderSweep< stencil::Direction::B, gen::bc_grid_aligned::FreeSlipBottom >{
         blocks, gen::bc_grid_aligned::FreeSlipBottom{ cpuFields.pdfsId }
      };
   }

   auto noSlipTop()
   {
      return sweep::BorderSweep< stencil::Direction::T, gen::bc_grid_aligned::NoSlipTop >{
         blocks, gen::bc_grid_aligned::NoSlipTop{ cpuFields.pdfsId }
      };
   }

   auto irregularFreeSlipFactory() { return gen::bc_sparse::FreeSlipIrregularFactory(blocks, cpuFields.pdfsId); }

   void syncGhostLayers() { commCpu(); }

   void wait() { /* NOP */ }

   void fields2host() { /* NOP */ }
   void fields2device() { /* NOP */ }

#endif
};

struct SimDomainBuilder
{
   std::array< uint_t, 3 > blocks;
   std::array< uint_t, 3 > cellsPerBlock;
   std::array< bool, 3 > periodic;

   SimDomain build()
   {
      auto sbfs =
         blockforest::createUniformBlockGrid(blocks[0], blocks[1], blocks[2], cellsPerBlock[0], cellsPerBlock[1],
                                             cellsPerBlock[2], 1.0, true, periodic[0], periodic[1], periodic[2]);

#if defined(LBM_SCENARIOS_GPU_BUILD)
      auto hostAlloc = make_shared< gpu::HostFieldAllocator< PdfField_T::value_type > >();
#else
      auto hostAlloc = make_shared< field::StdFieldAlloc< PdfField_T::value_type > >();
#endif

      const BlockDataID pdfsId = field::addToStorage< PdfField_T >(sbfs, "f", real_c(0.0), field::fzyx, 1, hostAlloc);
      const BlockDataID rhoId =
         field::addToStorage< ScalarField_T >(sbfs, "rho", real_c(0.0), field::fzyx, 1, hostAlloc);
      const BlockDataID uId = field::addToStorage< VectorField_T >(sbfs, "u", real_c(0.0), field::fzyx, 1, hostAlloc);
      const BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(sbfs, "flagField");

      CpuCommScheme commCpu{ sbfs };
      auto pdfsPackInfo = std::make_shared< CpuPdfsPackInfo >(pdfsId);
      commCpu.addPackInfo(pdfsPackInfo);

#if defined(LBM_SCENARIOS_GPU_BUILD)
      const BlockDataID pdfsIdGpu = gpu::addGPUFieldToStorage< PdfField_T >(sbfs, pdfsId, "f_gpu");
      const BlockDataID rhoIdGpu  = gpu::addGPUFieldToStorage< ScalarField_T >(sbfs, rhoId, "rho_gpu");
      const BlockDataID uIdGpu    = gpu::addGPUFieldToStorage< VectorField_T >(sbfs, uId, "u_gpu");

      auto commGpu         = std::make_unique< GpuCommScheme >(sbfs);
      auto gpuPdfsPackInfo = std::make_shared< GpuPdfsPackInfo >(pdfsIdGpu);
      commGpu->addPackInfo(gpuPdfsPackInfo);
      // commCpu.addPackInfo(gpuPdfsPackInfo);
#endif

      return {
         .blocks = sbfs, //
         .cpuFields = { //
            .pdfsId = pdfsId,
            .rhoId = rhoId,
            .uId = uId,
            .flagFieldId = flagFieldId
         }, //
         .commCpu = commCpu, //
#if defined(LBM_SCENARIOS_GPU_BUILD)
         .gpuFields = { //
            .pdfsId = pdfsIdGpu,
            .rhoId = rhoIdGpu,
            .uId = uIdGpu
         },
         .commGpu = std::move(commGpu)
#endif
      };
   }
};

} // namespace BasicLbmScenarios
