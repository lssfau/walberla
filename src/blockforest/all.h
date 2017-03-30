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
//! \file all.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Collective header file for module blockforest
//
//======================================================================================================================

#pragma once

#include "AABBRefinementSelection.h"
#include "Block.h"
#include "BlockDataHandling.h"
#include "BlockForest.h"
#include "BlockForestEvaluation.h"
#include "BlockID.h"
#include "BlockNeighborhoodConstruction.h"
#include "BlockNeighborhoodSection.h"
#include "BlockReconstruction.h"
#include "HilbertCurveConstruction.h"
#include "Initialization.h"
#include "PhantomBlock.h"
#include "PhantomBlockForest.h"
#include "SetupBlock.h"
#include "SetupBlockForest.h"
#include "StructuredBlockForest.h"
#include "Types.h"
#include "Utility.h"

#include "communication/all.h"
#include "loadbalancing/all.h"

#include "blockforest/CMakeDefs.h"
