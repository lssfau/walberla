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
//! \file PdfFieldPackInfo.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \brief Optimized packing for PDF fields, by communicating only required components.
//
//======================================================================================================================

#pragma once

#include "lbm/field/PdfField.h"
#include "field/communication/StencilRestrictedPackInfo.h"
#include "communication/UniformPackInfo.h"
#include "core/debug/Debug.h"
#include "stencil/Directions.h"


namespace walberla {
namespace lbm {



/**
 * \brief Optimized PackInfo for PDF fields (communicates only components pointing to neighbor)
 *
 * In principle a PDF field can be communicated using a FieldPackInfo which
 * copies all entries (=values for f) of the ghost layer.
 *
 * For PDF fields, however, it is sufficient to communicate only those discrete velocities
 * which point into the direction of the neighboring block
 *
 * This PackInfo sends only these required components and therefore the sizes of the
 * messages decrease.
 *
 * see also documentation for FieldPackInfo
 *
 * \warning For fields with nrOfGhostLayers > 1:  only the components pointing towards
 *          the boundary are communicated, which may not be the desired behavior
 *          for the 'inner' ghost layers
 *
 * \ingroup lbm
 */
template< typename LatticeModel_T >
using PdfFieldPackInfo = field::communication::StencilRestrictedPackInfo< PdfField<LatticeModel_T>,
                                                                          typename LatticeModel_T::Stencil>;


} // namespace lbm
} // namespace walberla
