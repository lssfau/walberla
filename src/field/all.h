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
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Collective header file for module field
//
//======================================================================================================================

#pragma once

#include "AccuracyEvaluation.h"
#include "AccuracyEvaluationLinePlot.h"
#include "AddToStorage.h"
#include "CellCounter.h"
#include "EvaluationFilter.h"
#include "Field.h"
#include "FieldClone.h"
#include "FileIO.h"
#include "FlagField.h"
#include "FlagFunctions.h"
#include "Gather.h"
#include "GhostLayerField.h"
#include "MassEvaluation.h"
#include "Printers.h"
#include "StabilityChecker.h"
#include "SwapableCompare.h"
#include "SymmetryCheck.h"
#include "VolumetricFlowRateEvaluation.h"

#include "adaptors/all.h"
#include "allocation/all.h"
#include "blockforest/all.h"
#include "communication/all.h"
#include "interpolators/all.h"
#include "iterators/all.h"
#include "refinement/all.h"
#include "vtk/all.h"

#include "field/CMakeDefs.h"
