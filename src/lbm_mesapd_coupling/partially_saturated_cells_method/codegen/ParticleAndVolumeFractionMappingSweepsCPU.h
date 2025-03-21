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
//! \file ParticleAndVolumeFractionMappingSweepsCPU.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/mapping/ParticleBoundingBox.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"

#include "mesa_pd/common/AABBConversion.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/kernel/SingleCast.h"

#include <cassert>
#include <functional>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/Sphere.h>

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{
// The CPU version is on purpose in the gpu namespace to avoid changes in the application codes
namespace gpu
{

template< int Weighting_T >
inline void calculateWeighting(real_t* const  /*weighting*/, const real_t& /*epsilon*/, const real_t& /*tau*/)
{
   WALBERLA_STATIC_ASSERT(Weighting_T == 1 || Weighting_T == 2);
}

template<>
inline void calculateWeighting< 1 >(real_t* const weighting, const real_t& epsilon, const real_t& /*tau*/)
{
   *weighting = epsilon;
}
template<>
inline void calculateWeighting< 2 >(real_t* const weighting, const real_t& epsilon, const real_t& tau)
{
   *weighting = epsilon * (tau - real_t(0.5)) / ((real_t(1) - epsilon) + (tau - real_t(0.5)));
}

template< int Weighting_T >
void mapParticles(IBlock& blockIt, const ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA,
                  const real_t* const spherePositions, const real_t* const sphereRadii, const real_t* const f_rs,
                  const size_t* const numParticlesSubBlocks, const size_t* const particleIDsSubBlocks,
                  const Vector3< uint_t > subBlocksPerDim)
{
   auto nOverlappingParticlesField =
      blockIt.getData< nOverlappingParticlesField_T >(particleAndVolumeFractionSoA.nOverlappingParticlesFieldID);
   auto BsField  = blockIt.getData< BsField_T >(particleAndVolumeFractionSoA.BsFieldID);
   auto idxField = blockIt.getData< idxField_T >(particleAndVolumeFractionSoA.idxFieldID);
   auto BField   = blockIt.getData< BField_T >(particleAndVolumeFractionSoA.BFieldID);

   real_t dx = blockIt.getAABB().xSize() / real_t(nOverlappingParticlesField->xSize());

   WALBERLA_FOR_ALL_CELLS_XYZ(
      BField,
      for (size_t i = 0; i < MaxParticlesPerCell; i++) {
         BsField->get(x, y, z, i)  = real_t(0.0);
         idxField->get(x, y, z, i) = size_t(0);
      } nOverlappingParticlesField->get(x, y, z) = uint_t(0);
      BField->get(x, y, z)                       = real_t(0.0);
      const Vector3< real_t > cellCenter =
         Vector3< real_t >(real_t(x) + real_t(0.5) * dx, real_t(y) + real_t(0.5) * dx, real_t(z) + real_t(0.5) * dx) +
         blockIt.getAABB().minCorner();
      const Vector3< size_t > subBlockIndex(size_t(real_t(x) / blockIt.getAABB().xSize() * real_t(subBlocksPerDim[0])),
                                            size_t(real_t(y) / blockIt.getAABB().ySize() * real_t(subBlocksPerDim[1])),
                                            size_t(real_t(z) / blockIt.getAABB().zSize() * real_t(subBlocksPerDim[2])));
      const size_t linearizedSubBlockIndex = subBlockIndex[2] * subBlocksPerDim[0] * subBlocksPerDim[1] +
                                             subBlockIndex[1] * subBlocksPerDim[0] + subBlockIndex[0];

      for (size_t i = 0; i < numParticlesSubBlocks[linearizedSubBlockIndex]; i++) {
         size_t idxMapped = particleIDsSubBlocks[linearizedSubBlockIndex +
                                                 i * subBlocksPerDim[0] * subBlocksPerDim[1] * subBlocksPerDim[2]];
         const Vector3< real_t > minCornerSphere(spherePositions[idxMapped * 3] - sphereRadii[idxMapped],
                                                 spherePositions[idxMapped * 3 + 1] - sphereRadii[idxMapped],
                                                 spherePositions[idxMapped * 3 + 2] - sphereRadii[idxMapped]);
         const Vector3< real_t > maxCornerSphere(spherePositions[idxMapped * 3] + sphereRadii[idxMapped],
                                                 spherePositions[idxMapped * 3 + 1] + sphereRadii[idxMapped],
                                                 spherePositions[idxMapped * 3 + 2] + sphereRadii[idxMapped]);
         if (cellCenter[0] + dx > minCornerSphere[0] && cellCenter[0] - dx < maxCornerSphere[0] &&
             cellCenter[1] + dx > minCornerSphere[1] && cellCenter[1] - dx < maxCornerSphere[1] &&
             cellCenter[2] + dx > minCornerSphere[2] && cellCenter[2] - dx < maxCornerSphere[2])
         {
            const Vector3< real_t > cellSphereVector(spherePositions[idxMapped * 3] - cellCenter[0],
                                                     spherePositions[idxMapped * 3 + 1] - cellCenter[1],
                                                     spherePositions[idxMapped * 3 + 2] - cellCenter[2]);

            const real_t D =
               real_t(sqrt(cellSphereVector[0] * cellSphereVector[0] + cellSphereVector[1] * cellSphereVector[1] +
                           cellSphereVector[2] * cellSphereVector[2])) -
               sphereRadii[idxMapped];

            real_t epsilon = -D + f_rs[idxMapped];
            epsilon        = std::max(epsilon, real_t(0));
            epsilon        = std::min(epsilon, real_t(1));

            // Store overlap fraction only if there is an intersection
            if (epsilon > 0.0)
            {
               // Check that the maximum number of overlapping particles has not yet been reached
               assert(nOverlappingParticlesField->get(x, y, z) < MaxParticlesPerCell);
               BsField->get(x, y, z, nOverlappingParticlesField->get(x, y, z)) = epsilon;
               calculateWeighting< Weighting_T >(&BsField->get(x, y, z, nOverlappingParticlesField->get(x, y, z)),
                                                 BsField->get(x, y, z, nOverlappingParticlesField->get(x, y, z)),
                                                 real_t(1.0) / particleAndVolumeFractionSoA.omega_);
               idxField->get(x, y, z, nOverlappingParticlesField->get(x, y, z)) = idxMapped;
               BField->get(x, y, z) += BsField->get(x, y, z, nOverlappingParticlesField->get(x, y, z));
               nOverlappingParticlesField->get(x, y, z) += 1;
            }
         }
      }

      // Normalize fraction field (Bs) if sum over all fractions (B) > 1
      if (BField->get(x, y, z) > 1) {
         for (size_t i = 0; i < nOverlappingParticlesField->get(x, y, z); i++)
         {
            BsField->get(x, y, z, i) /= BField->get(x, y, z);
         }
         BField->get(x, y, z) = 1.0;
      })
}

template< typename ParticleAccessor_T, typename ParticleSelector_T, int Weighting_T >
class SphereFractionMappingSweep
{
 public:
   SphereFractionMappingSweep(const shared_ptr< StructuredBlockStorage >& blockStorage,
                              const shared_ptr< ParticleAccessor_T >& ac,
                              const ParticleSelector_T& mappingParticleSelector,
                              ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA,
                              const Vector3< uint_t > subBlockSize)
      : blockStorage_(blockStorage), ac_(ac), mappingParticleSelector_(mappingParticleSelector),
        particleAndVolumeFractionSoA_(particleAndVolumeFractionSoA), subBlockSize_(subBlockSize)
   {
      static_assert(std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value,
                    "Provide a valid accessor as template");
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         auto aabb = blockIt->getAABB();
         if (size_t(aabb.xSize()) % subBlockSize_[0] != 0 || size_t(aabb.ySize()) % subBlockSize_[1] != 0 ||
             size_t(aabb.zSize()) % subBlockSize_[2] != 0)
         {
            WALBERLA_ABORT("Number of cells per block (" << aabb << ") is not divisible by subBlockSize ("
                                                         << subBlockSize_ << ").")
         }
      }
   }

   void operator()(IBlock* block)
   {
      size_t numMappedParticles = 0;
      particleAndVolumeFractionSoA_.mappingUIDs.clear();
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            numMappedParticles++;
            // Store UIDs to make sure that the particles have not changed between the mapping and the PSM sweep
            particleAndVolumeFractionSoA_.mappingUIDs.push_back(ac_->getUid(idx));
         }
      }

      if (numMappedParticles == uint_t(0)) return;

      // Allocate memory storing the particle information needed for the overlap fraction computations
      const size_t scalarArraySize = numMappedParticles * sizeof(real_t);

      if (particleAndVolumeFractionSoA_.positions != nullptr) { free(particleAndVolumeFractionSoA_.positions); }
      particleAndVolumeFractionSoA_.positions = (real_t*) malloc(3 * scalarArraySize);
      real_t* radii                           = (real_t*) malloc(scalarArraySize);
      real_t* f_r = (real_t*) malloc(scalarArraySize); // f_r is described in https://doi.org/10.1108/EC-02-2016-0052

      // Store particle information inside the memory
      size_t idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            for (size_t d = 0; d < 3; ++d)
            {
               particleAndVolumeFractionSoA_.positions[idxMapped * 3 + d] = ac_->getPosition(idx)[d];
            }
            // If other shapes than spheres are mapped, ignore them here
            if (ac_->getShape(idx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE)
            {
               const real_t radius = static_cast< mesa_pd::data::Sphere* >(ac_->getShape(idx))->getRadius();
               radii[idxMapped]    = radius;
               real_t Va           = real_t(
                  (1.0 / 12.0 - radius * radius) * atan((0.5 * sqrt(radius * radius - 0.5)) / (0.5 - radius * radius)) +
                  1.0 / 3.0 * sqrt(radius * radius - 0.5) +
                  (radius * radius - 1.0 / 12.0) * atan(0.5 / sqrt(radius * radius - 0.5)) -
                  4.0 / 3.0 * radius * radius * radius * atan(0.25 / (radius * sqrt(radius * radius - 0.5))));
               f_r[idxMapped] = Va - radius + real_t(0.5);
            }
            idxMapped++;
         }
      }

      // Update fraction mapping
      // Split the block into sub-blocks and sort the particle indices into each overlapping sub-block. This way, in
      // the particle mapping, each iteration only has to check the potentially overlapping particles.
      auto blockAABB = block->getAABB();
      const Vector3< uint_t > subBlocksPerDim =
         Vector3< uint_t >(uint_t(blockAABB.xSize()) / subBlockSize_[0], uint_t(blockAABB.ySize()) / subBlockSize_[1],
                           uint_t(blockAABB.zSize()) / subBlockSize_[2]);
      const size_t numSubBlocks = subBlocksPerDim[0] * subBlocksPerDim[1] * subBlocksPerDim[2];
      std::vector< std::vector< size_t > > subBlocks(numSubBlocks);

      idxMapped = 0;
      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_))
         {
            if (ac_->getShape(idx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE)
            {
               auto sphereAABB = mesa_pd::getParticleAABB(idx, *ac_);
               if (blockAABB.intersects(sphereAABB))
               {
                  auto intersectionAABB = blockAABB.getIntersection(sphereAABB);
                  intersectionAABB.translate(-blockAABB.minCorner());
                  mesa_pd::Vec3 blockScaling = mesa_pd::Vec3(real_t(subBlocksPerDim[0]) / blockAABB.sizes()[0],
                                                             real_t(subBlocksPerDim[1]) / blockAABB.sizes()[1],
                                                             real_t(subBlocksPerDim[2]) / blockAABB.sizes()[2]);

                  for (size_t z = size_t(intersectionAABB.zMin() * blockScaling[2]);
                       z < size_t(ceil(intersectionAABB.zMax() * blockScaling[2])); ++z)
                  {
                     for (size_t y = size_t(intersectionAABB.yMin() * blockScaling[1]);
                          y < size_t(ceil(intersectionAABB.yMax() * blockScaling[1])); ++y)
                     {
                        for (size_t x = size_t(intersectionAABB.xMin() * blockScaling[0]);
                             x < size_t(ceil(intersectionAABB.xMax() * blockScaling[0])); ++x)
                        {
                           size_t index = z * subBlocksPerDim[0] * subBlocksPerDim[1] + y * subBlocksPerDim[0] + x;
                           subBlocks[index].push_back(idxMapped);
                        }
                     }
                  }
               }
            }
            idxMapped++;
         }
      }

      size_t maxParticlesPerSubBlock = 0;
      std::for_each(subBlocks.begin(), subBlocks.end(), [&maxParticlesPerSubBlock](std::vector< size_t >& subBlock) {
         maxParticlesPerSubBlock = std::max(maxParticlesPerSubBlock, subBlock.size());
      });

      size_t* numParticlesPerSubBlock = (size_t*) malloc(numSubBlocks * sizeof(size_t));
      size_t* particleIDsSubBlocks    = (size_t*) malloc(numSubBlocks * maxParticlesPerSubBlock * sizeof(size_t));

      // Copy data from std::vector to memory
      for (size_t z = 0; z < subBlocksPerDim[2]; ++z)
      {
         for (size_t y = 0; y < subBlocksPerDim[1]; ++y)
         {
            for (size_t x = 0; x < subBlocksPerDim[0]; ++x)
            {
               size_t index = z * subBlocksPerDim[0] * subBlocksPerDim[1] + y * subBlocksPerDim[0] + x;
               numParticlesPerSubBlock[index] = subBlocks[index].size();
               for (size_t k = 0; k < subBlocks[index].size(); k++)
               {
                  particleIDsSubBlocks[index + k * numSubBlocks] = subBlocks[index][k];
               }
            }
         }
      }

      mapParticles(*block, particleAndVolumeFractionSoA_, particleAndVolumeFractionSoA_.positions, radii, f_r,
                   numParticlesPerSubBlock, particleIDsSubBlocks, subBlocksPerDim);

      free(numParticlesPerSubBlock);
      free(particleIDsSubBlocks);

      free(radii);
      free(f_r);
   }

   shared_ptr< StructuredBlockStorage > blockStorage_;
   const shared_ptr< ParticleAccessor_T > ac_;
   const ParticleSelector_T& mappingParticleSelector_;
   ParticleAndVolumeFractionSoA_T< Weighting_T >& particleAndVolumeFractionSoA_;
   const Vector3< uint_t > subBlockSize_;
};

} // namespace gpu
} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
