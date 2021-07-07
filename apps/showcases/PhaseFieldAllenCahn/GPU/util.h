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
//! \file util.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "python_coupling/DictWrapper.h"

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"

#include "field/communication/PackInfo.h"
#include "field/FlagField.h"
#include "field/vtk/VTKWriter.h"
#include "GenDefines.h"
#pragma once

namespace walberla {

    void calc_total_velocity(const shared_ptr <StructuredBlockStorage> &blocks, std::array<real_t, 5> &total_velocity,
                             BlockDataID phaseFieldID, BlockDataID velocityFieldID, ConstBlockDataID flagFieldID, FlagUID fluidFlagUID);

    void flood_fill(PhaseField_T &phaseField, VelocityField_T &velocityField, CellInterval boundingBox,
                    real_t &volume, uint_t &nrOfCells,
                    std::array<real_t, 3> &center_of_mass, std::array<real_t, 4> &total_velocity);

} // namespace walberla

