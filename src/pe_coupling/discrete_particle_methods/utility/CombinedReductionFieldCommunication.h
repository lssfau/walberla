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
//! \file CombinedReductionFieldCommunication.h
//! \ingroup pe_coupling
//! \author Tobias Schruff <schruff@iww.rwth-aachen.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/StructuredBlockForest.h"

#include "field/communication/PackInfo.h"
#include "field/communication/UniformPullReductionPackInfo.h"

#include "stencil/D3Q27.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

/*!\brief Special communication routine to += the values in the ghost layers and update the ghost layers.
 *
 * Reduces the field with += to obtain the values from the block neighbors' ghost layers and add them onto the corresponding values inside the domain.
 * This is followed by a regular ghost layer communication to synchronize the ghost layer values.
 *
 * This functionality is typically used to synchronize the solid volume fraction field and the interaction force fields,
 * after they have been filled by a distributor.

 * For more infos on distributors, see src/field/distributors.
 *
 */
template< typename GhostLayerField_T >
class CombinedReductionFieldCommunication
{
public:

   CombinedReductionFieldCommunication( const shared_ptr<StructuredBlockForest> & bf, const BlockDataID & glFieldID )
      : pullReductionScheme_( bf ), commScheme_( bf )
   {
      pullReductionScheme_.addPackInfo( make_shared< field::communication::UniformPullReductionPackInfo<std::plus, GhostLayerField_T> >( glFieldID ) );
      commScheme_.addPackInfo( make_shared< field::communication::PackInfo<GhostLayerField_T> >( glFieldID ) );
      
   }

   void operator()()
   {
      pullReductionScheme_();
      commScheme_();
   }

private:
   blockforest::communication::UniformBufferedScheme<stencil::D3Q27> pullReductionScheme_;
   blockforest::communication::UniformBufferedScheme<stencil::D3Q27> commScheme_; 
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
