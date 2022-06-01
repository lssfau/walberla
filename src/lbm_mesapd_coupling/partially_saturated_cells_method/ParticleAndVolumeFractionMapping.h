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
//! \file ParticleAndVolumeFractionMapping.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/mapping/ParticleBoundingBox.h"
#include "lbm_mesapd_coupling/overlapping/OverlapFraction.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/kernel/SingleCast.h"

#include <functional>
#include <mesa_pd/data/ParticleStorage.h>

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

template< typename ParticleAccessor_T, typename ParticleSelector_T >
class ParticleAndVolumeFractionMapping
{
 public:
   ParticleAndVolumeFractionMapping(const shared_ptr< StructuredBlockStorage >& blockStorage,
                                    const shared_ptr< ParticleAccessor_T >& ac,
                                    const ParticleSelector_T& mappingParticleSelector,
                                    const BlockDataID& particleAndVolumeFractionFieldID,
                                    const uint_t superSamplingDepth = uint_t(4))
      : blockStorage_(blockStorage), ac_(ac), mappingParticleSelector_(mappingParticleSelector),
        particleAndVolumeFractionFieldID_(particleAndVolumeFractionFieldID), superSamplingDepth_(superSamplingDepth)
   {
      static_assert(std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value,
                    "Provide a valid accessor as template");
   }

   void operator()()
   {
      // clear the field
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         ParticleAndVolumeFractionField_T* particleAndVolumeFractionField =
            blockIt->getData< ParticleAndVolumeFractionField_T >(particleAndVolumeFractionFieldID_);

         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(particleAndVolumeFractionField,
                                                          particleAndVolumeFractionField->nrOfGhostLayers(),
                                                          (particleAndVolumeFractionField->get(x, y, z)).clear();)
      }

      for (size_t idx = 0; idx < ac_->size(); ++idx)
      {
         if (mappingParticleSelector_(idx, *ac_)) { update(idx); }
      }
   }

 private:
   void update(const size_t idx)
   {
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         ParticleAndVolumeFractionField_T* particleAndVolumeFractionField =
            blockIt->getData< ParticleAndVolumeFractionField_T >(particleAndVolumeFractionFieldID_);

         CellInterval cellBB =
            getParticleCellBB(idx, *ac_, *blockIt, *blockStorage_, particleAndVolumeFractionField->nrOfGhostLayers());

         uint_t level = blockStorage_->getLevel(*blockIt);
         Vector3< real_t > dxVec(blockStorage_->dx(level), blockStorage_->dy(level), blockStorage_->dz(level));

         // compute the overlap fraction for each cell that has to be considered and save the information
         for (auto cellIt = cellBB.begin(); cellIt != cellBB.end(); ++cellIt)
         {
            Cell cell(*cellIt);

            Vector3< real_t > cellCenter;
            cellCenter = blockStorage_->getBlockLocalCellCenter(*blockIt, cell);

            real_t fraction = singleCast_(idx, *ac_, overlapFractionFctr_, ac_, cellCenter, dxVec, superSamplingDepth_);

            id_t particleUid = ac_->getUid(idx);
            if (fraction > real_t(0)) { particleAndVolumeFractionField->get(cell).emplace_back(particleUid, fraction); }
         }
      }
   }

   shared_ptr< StructuredBlockStorage > blockStorage_;
   const shared_ptr< ParticleAccessor_T > ac_;
   ParticleSelector_T mappingParticleSelector_;
   const BlockDataID particleAndVolumeFractionFieldID_;
   const uint_t superSamplingDepth_;

   mesa_pd::kernel::SingleCast singleCast_;
   OverlapFractionFunctor overlapFractionFctr_;
};

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
