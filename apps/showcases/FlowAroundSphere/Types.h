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
//! \file Types.h
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
# pragma once

#include "domain_decomposition/BlockDataID.h"
#include "lbm_generated/field/PdfField.h"
#include "FlowAroundSphereInfoHeader.h"

using namespace walberla;

using StorageSpecification_T = lbm::FlowAroundSphereStorageSpecification;
using Stencil_T              = StorageSpecification_T::Stencil;
using CommunicationStencil_T = StorageSpecification_T::CommunicationStencil;

using PdfField_T           = lbm_generated::PdfField< StorageSpecification_T >;
using FlagField_T          = FlagField< uint8_t >;

struct IDs
{
   BlockDataID pdfField;
   BlockDataID velocityField;
   BlockDataID densityField;
   BlockDataID omegaField;
   BlockDataID flagField;

   BlockDataID pdfFieldGPU;
   BlockDataID velocityFieldGPU;
   BlockDataID densityFieldGPU;
   BlockDataID omegaFieldGPU;
};
