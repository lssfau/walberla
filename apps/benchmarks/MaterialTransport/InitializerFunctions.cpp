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
using DensityField_concentration_T = walberla::field::GhostLayerField<double, 1>;
using VelocityField_fluid_T = walberla::field::GhostLayerField<double, 3>;

namespace walberla{

void initConcentrationField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID ConcentrationFieldID,const math::AABB& domainAABB,Vector3< uint_t > domainSize){
   const real_t radius = real_c(domainSize[1]/10);
   for (auto& block : *blocks)
   {
      auto ConcentrationField = block.getData< DensityField_concentration_T >(ConcentrationFieldID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(ConcentrationField,{
         Cell globalCell;
         const auto cellAABB = blocks->getBlockLocalCellAABB(block, Cell(x,y,z));
         auto cellCenter = cellAABB.center();
         real_t posX = cellCenter[0];
         real_t posY = cellCenter[1];
         real_t posZ = cellCenter[2];

         real_t distance = real_c(
            sqrt(pow((domainAABB.center()[0] - posX), 2) +
                 pow((domainAABB.center()[1] - posY), 2) +
                 pow((domainAABB.center()[2] - posZ), 2)));

         if(distance <= radius){
            ConcentrationField->get(x,y,z) = real_t(1.05);
         }
         else{ConcentrationField->get(x,y,z) = real_t(0.5);}

      })  // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }

}  // initConcentrationField


void initFluidField(const shared_ptr< StructuredBlockStorage >& blocks, BlockDataID FluidFieldID,const Vector3<real_t> uInflow){

   for (auto& block : *blocks)
   {
      auto FluidVelocityField = block.getData< VelocityField_fluid_T >(FluidFieldID);

      WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(FluidVelocityField,{
         Cell globalCell;
         const auto cellAABB = blocks->getBlockLocalCellAABB(block, Cell(x,y,z));
         FluidVelocityField->get(x,y,z,0) = uInflow[0];
         FluidVelocityField->get(x,y,z,1) = uInflow[1];
         FluidVelocityField->get(x,y,z,2) = uInflow[2];

      })  // WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ
   }
}



}  // namespace walberla