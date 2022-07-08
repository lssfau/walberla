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
//! \file BubbleDefinitions.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Type definitions for the bubble_model.
//
//======================================================================================================================

#pragma once

#include "field/GhostLayerField.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
using BubbleID                   = uint32_t;
const uint32_t INVALID_BUBBLE_ID = uint32_t(-1);

using BubbleField_T = GhostLayerField< BubbleID, 1 >;
using ScalarField_T = GhostLayerField< real_t, 1 >;

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla
