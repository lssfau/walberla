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
//! \file BubbleModelFromConfig.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Create and initialize BubbleModel with information from a config object.
//
//======================================================================================================================

#pragma once

#include "core/config/Config.h"

#include "BubbleModel.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
* Create and initialize BubbleModel with information from a config object.
*
* Example configuration block:
* \verbatim
  {
      // For scenarios where no splits can occur - expensive split detection can be switched off - defaults to False
      enableBubbleSplits False;

      // optional atmosphere block(s) : Atmosphere bubbles are bubbles with constant pressure
      atmosphere {
         position  < 1.5, 175.5, 1 >;   // a point inside the bubble that should becomes atmosphere
         density     1.1;           // the value of constant pressure. Default value: 1.0
      }

      // optional disjoining pressure
      // Disjoining pressure model holds bubbles apart that are close to each other
      DisjoiningPressure
      {
         maxDisjoiningPressure 0.1; // maximum value of disjoining pressure - defaults to 0.2
         maxDistance 4;             // for bubbles with distance greater 'maxDistance' cells
                                    // there is no disjoining pressure - defaults to 10 cells
      }
   }
   \endverbatim
*
*
* The fill level field must be fully initialized.
* If the handle of the config block returns nullptr, the bubble model is created with default values.
***********************************************************************************************************************/
template< typename Stencil_T >
std::shared_ptr< BubbleModelBase > createFromConfig(const std::shared_ptr< StructuredBlockForest >& blockStorage,
                                                    ConstBlockDataID fillFieldID,
                                                    const Config::BlockHandle& configBlock = Config::BlockHandle());

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

#include "BubbleModelFromConfig.impl.h"