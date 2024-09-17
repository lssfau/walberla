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
//! \author Ravi Ayyala Somayajula <ravi.k.ayyala@fau.de>
//
//======================================================================================================================

#include "core/Environment.h"
#include "core/logging/Initialization.h"
#include "core/math/Constants.h"

#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/VTKWriter.h"

//#include "python_coupling/DictWrapper.h"
//#include "GeneralInfoHeader.h"
#pragma once

namespace walberla
{
void initConcentrationField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& ConcentrationFieldID,const math::AABB& domainAABB,Vector3< uint_t > domainSize);
void initConcentrationFieldGaussian(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& ConcentrationFieldID,const math::AABB& domainAABB,Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,const Vector3<real_t> uInflow,const Vector3< real_t > x_0);
void initFluidField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,const Vector3<real_t> uInflow);
void initFluidFieldPoiseuille(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,
                              const Vector3< real_t > uInflow,Vector3< uint_t > domainSize);
void analyticalSolGaussian(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& AnalyticalConcentrationFieldID,const math::AABB& domainAABB,Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,const Vector3<real_t> uInflow,const Vector3< real_t > x_0, const real_t time, uint_t advection_period);
std::vector<real_t> computeErrorL2(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& NumericalSolFieldID, BlockDataID& AnalyticalSolFieldID,BlockDataID& ErrorFieldID,
                                                const math::AABB& domainAABB);
} // namespace walberla