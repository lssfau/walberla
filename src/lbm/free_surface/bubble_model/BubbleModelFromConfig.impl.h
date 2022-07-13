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
//! \file BubbleModelFromConfig.impl.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Create and initialize BubbleModel with information from a config object.
//
//======================================================================================================================

#include "BubbleModelFromConfig.h"
#include "DisjoiningPressureBubbleModel.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
template< typename Stencil_T >
std::shared_ptr< BubbleModelBase > createFromConfig(const std::shared_ptr< StructuredBlockForest >& blockForest,
                                                    ConstBlockDataID fillFieldID,
                                                    const Config::BlockHandle& configBlock)
{
   if (!configBlock)
   {
      auto bubbleModel = std::make_shared< BubbleModel< Stencil_T > >(blockForest);
      bubbleModel->initFromFillLevelField(fillFieldID);
      return bubbleModel;
   }

   bool enableBubbleSplits = configBlock.getParameter< bool >("enableBubbleSplits", false);

   std::shared_ptr< BubbleModel< Stencil_T > > bubbleModel;

   real_t constantLatticeDensity = configBlock.getParameter< real_t >("constantLatticeDensity", real_c(-1));
   if (constantLatticeDensity > 0) { return std::make_shared< BubbleModelConstantPressure >(constantLatticeDensity); }

   auto disjoiningPressureCfg = configBlock.getBlock("DisjoiningPressure");
   if (disjoiningPressureCfg)
   {
      const real_t maxDistance = disjoiningPressureCfg.getParameter< real_t >(
         "maxBubbleDistance"); // d_range in the paper from Koerner et al., 2005
      const real_t disjoiningPressureConstant = disjoiningPressureCfg.getParameter< real_t >(
         "disjoiningPressureConstant"); // c_pi in the paper from Koerner et al., 2005
      const uint_t distFieldUpdateInter = disjoiningPressureCfg.getParameter< uint_t >("distanceFieldUpdateInterval");

      bubbleModel = std::make_shared< DisjoiningPressureBubbleModel< Stencil_T > >(
         blockForest, maxDistance, disjoiningPressureConstant, enableBubbleSplits, distFieldUpdateInter);

      WALBERLA_LOG_PROGRESS_ON_ROOT("Using Disjoining Pressure BubbleModel");
   }
   else
   {
      bubbleModel = std::make_shared< bubble_model::BubbleModel< Stencil_T > >(blockForest, enableBubbleSplits);
      WALBERLA_LOG_PROGRESS_ON_ROOT("Using Standard BubbleModel");
   }

   bubbleModel->initFromFillLevelField(fillFieldID);

   // get all atmosphere blocks
   std::vector< Config::BlockHandle > atmosphereBlocks;
   configBlock.getBlocks("atmosphere", atmosphereBlocks);
   for (auto it = atmosphereBlocks.begin(); it != atmosphereBlocks.end(); ++it)
   {
      Vector3< real_t > position = it->getParameter< Vector3< real_t > >("position");
      real_t density             = it->getParameter< real_t >("density", 1.0);
      Cell cell;
      blockForest->getCell(cell, position[0], position[1], position[2]);
      bubbleModel->setAtmosphere(cell, density);
   }

   return bubbleModel;
}

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
