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
//! \file InitializerFunctions.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

#include "python_coupling/DictWrapper.h"
#pragma once

namespace walberla
{
void initPhaseField_sphere(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID, real_t R,
                           Vector3< real_t > bubbleMidPoint, bool bubble = true, real_t W = 5);

void init_Taylor_bubble(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID, real_t D = 5,
                        real_t H = 2, real_t DT = 20, real_t Donut_x0 = 40);

void init_bubble_field(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID, real_t R,
                       real_t W = 5);

void initPhaseField_RTI(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID phaseFieldID, real_t W = 5,
                        const bool pipe = true);

void initTubeWithCylinder(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID flagFieldID,
                          field::FlagUID boundaryFlagUID, real_t R_in, real_t eccentricity, real_t start_transition,
                          real_t length_transition, bool const eccentricity_or_pipe_ratio);

} // namespace walberla
