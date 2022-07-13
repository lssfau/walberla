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
//! \file Geometry.h
//! \ingroup bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Geometrical helper functions for the bubble model.
//
//======================================================================================================================

#pragma once

#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/FlagField.h"

#include "geometry/bodies/BodyOverlapFunctions.h"

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
/***********************************************************************************************************************
 * Set fill levels in a scalar ghost layer field (GhostLayerField<real_t,1>) using a geometric body object. Specifying
 * isGas determines whether the body is initialized as gas bubble or liquid drop.
 *
 * The function overlapVolume(body, cellmidpoint, dx) has to be defined for the body object.
 *
 * The volume fraction of the body is _SUBTRACTED_ from the each cells' current fill level, including ghost layer cells.
 * Therefore, the field should be initialized with 1 in each cell ( "everywhere fluid").
 * Then:
 *  - if a cell is inside a body, its fill level is 0
 *  - if a cell is completely outside the body, its fill level is not changed and remains 1
 *  - if a cell is partially inside the sphere, the amount of overlap is subtracted from the current fill level; it can
 *    not become negative and is limited to 0
 **********************************************************************************************************************/
template< typename Body_T >
void addBodyToFillLevelField(StructuredBlockStorage& blockStorage, BlockDataID fillFieldID, const Body_T& body,
                             bool isGas);

/***********************************************************************************************************************
 * Set flag field according to given fill level field.
 *
 * FillLevel <=0   : gas flag (and both other flags are removed if set)
 * FillLevel >=1   : fluid flag (and both other flags are removed if set)
 * otherwise       : interface flag (and both other flags are removed if set)
 *
 * The three FlagUIDs for liquid, gas and interface must already be registered at the flag field.
 **********************************************************************************************************************/
template< typename flag_t >
void setFlagFieldFromFillLevels(FlagField< flag_t >* flagField, const GhostLayerField< real_t, 1 >* fillField,
                                const FlagUID& liquid, const FlagUID& gas, const FlagUID& interFace);

/***********************************************************************************************************************
 * Check if interface layer is valid and correctly separates gas and fluid cells:
 * - gas cell must not have a fluid cell in its (Stencil_T) neighborhood
 * - fluid cell must not have a gas cell in its (Stencil_T) neighborhood
 *
 * The three FlagUIDs for liquid, gas and interface must already be registered at the flag field.
 **********************************************************************************************************************/
template< typename flag_t, typename Stencil_T >
bool checkForValidInterfaceLayer(const FlagField< flag_t >* flagField, const FlagUID& liquid, const FlagUID& gas,
                                 const FlagUID& interFace, bool printWarnings);

} // namespace bubble_model
} // namespace free_surface
} // namespace walberla

#include "Geometry.impl.h"