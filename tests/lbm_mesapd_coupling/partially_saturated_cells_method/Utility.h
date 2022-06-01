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
//! \file Utility.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{

//*******************************************************************************************************************
/*!\brief Calculating the sum over all fraction values. This can be used as a sanity check since it has to be roughly
 * equal to the volume of all particles.
 *
 */
//*******************************************************************************************************************
class FractionFieldSum
{
 public:
   FractionFieldSum(const shared_ptr< StructuredBlockStorage >& blockStorage,
                    BlockDataID particleAndVolumeFractionFieldID)
      : blockStorage_(blockStorage), particleAndVolumeFractionFieldID_(particleAndVolumeFractionFieldID)
   {}

   real_t operator()()
   {
      real_t sum = 0.0;

      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt)
      {
         psm::ParticleAndVolumeFractionField_T* particleAndVolumeFractionField =
            blockIt->getData< psm::ParticleAndVolumeFractionField_T >(particleAndVolumeFractionFieldID_);

         const cell_idx_t xSize = cell_idx_c(particleAndVolumeFractionField->xSize());
         const cell_idx_t ySize = cell_idx_c(particleAndVolumeFractionField->ySize());
         const cell_idx_t zSize = cell_idx_c(particleAndVolumeFractionField->zSize());

         for (cell_idx_t z = 0; z < zSize; ++z)
         {
            for (cell_idx_t y = 0; y < ySize; ++y)
            {
               for (cell_idx_t x = 0; x < xSize; ++x)
               {
                  for (auto particleAndVolumeFraction : particleAndVolumeFractionField->get(x, y, z))
                  {
                     sum += particleAndVolumeFraction.second;
                  }
               }
            }
         }
      }

      WALBERLA_MPI_SECTION() { mpi::allReduceInplace(sum, mpi::SUM); }

      return sum;
   }

 private:
   shared_ptr< StructuredBlockStorage > blockStorage_;
   BlockDataID particleAndVolumeFractionFieldID_;
};

} // namespace lbm_mesapd_coupling
} // namespace walberla
