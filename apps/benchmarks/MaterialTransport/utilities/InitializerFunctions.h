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
void initFluidField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,const Vector3<real_t> uInflow,Vector3< uint_t > domainSize);
void initFluidFieldPoiseuille(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& FluidFieldID,
                              const Vector3< real_t > uInflow,Vector3< uint_t > domainSize);
void analyticalSolGaussian(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& AnalyticalConcentrationFieldID,const math::AABB& domainAABB,Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,const Vector3<real_t> uInflow,const Vector3< real_t > x_0, const real_t time, uint_t advection_period);
std::vector<real_t> computeErrorL2(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID& NumericalSolFieldID, BlockDataID& AnalyticalSolFieldID,BlockDataID& ErrorFieldID,
                                                const math::AABB& domainAABB);
void initConcentrationFieldSinusoidal(const shared_ptr< StructuredBlockStorage >& blocks,
                                      BlockDataID& ConcentrationFieldID, const math::AABB& domainAABB,
                                      Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,
                                      const Vector3< real_t > uInflow, const Vector3< real_t > x_0,const real_t dx,const real_t dt);
void analyticalSolSinusoidal(const shared_ptr< StructuredBlockStorage >& blocks,
                             BlockDataID& AnalyticalConcentrationFieldID, const math::AABB& domainAABB,
                             Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,
                             const Vector3< real_t > uInflow, const Vector3< real_t > x_0, const real_t time,const real_t dx,const real_t dt);

void initConcentrationFieldPacket(const shared_ptr< StructuredBlockStorage >& blocks,
                                  BlockDataID& ConcentrationFieldID, const math::AABB& domainAABB,
                                  Vector3< uint_t > domainSize, const real_t sigma_0, const real_t sigma_D,
                                  const Vector3< real_t > uInflow, const Vector3< real_t > x_0,const real_t dx,const real_t dt,const real_t diffusivity);

void analyticalSolPacket(const shared_ptr< StructuredBlockStorage >& blocks,
                         BlockDataID& AnalyticalConcentrationFieldID, const math::AABB& domainAABB,
                         Vector3< uint_t > domainSize, const real_t sigma_0, const real_t diffusivity,
                         const Vector3< real_t > uInflow, const Vector3< real_t > x_0, const real_t time,const real_t dx,const real_t dt);

std::vector< real_t > NusseltNumbers(const shared_ptr< StructuredBlockStorage >& blocks,
               BlockDataID& ConcentrationFieldID,BlockDataID& VelocityFieldID, const math::AABB& domainAABB,
               Vector3< uint_t > domainSize,const real_t delta_theta,const real_t ratio);

real_t Simpsons (real_t bb, real_t aa, BlockDataID& ConcentrationFieldID,const shared_ptr< StructuredBlockStorage >& blocks);
} // namespace walberla