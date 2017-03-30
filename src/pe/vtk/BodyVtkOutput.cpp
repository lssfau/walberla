

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
//! \file BodyVtkOutput.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "BodyVtkOutput.h"
#include "pe/rigidbody/BodyStorage.h"

#include "core/selectable/IsSetSelected.h"
#include "core/uid/GlobalState.h"


namespace walberla {
namespace pe {

std::vector< DefaultBodyVTKOutput::Attributes > DefaultBodyVTKOutput::getAttributes() const
{
   std::vector< Attributes > attributes;
   attributes.push_back( Attributes( vtk::typeToString< float >(), "Velocity", uint_c(3) ) );
   attributes.push_back( Attributes( vtk::typeToString< int >(), "rank", uint_c(1) ) );
   attributes.push_back( Attributes( vtk::typeToString< int >(), "shadow", uint_c(1) ) );

   return attributes;
}

void DefaultBodyVTKOutput::configure()
{
   bodies_.clear();
   for( auto blockIt = blockStorage_.begin(); blockIt != blockStorage_.end(); ++blockIt )
   {

      const Storage& bs = *(blockIt->getData<const Storage>( storageID_ ));

      for( auto it = bs[0].begin(); it != bs[0].end(); ++it )
      {
         bodies_.push_back( *it );
      }
   }
}

std::vector< math::Vector3< real_t > > DefaultBodyVTKOutput::getPoints()
{
   std::vector< math::Vector3< real_t > > result( bodies_.size() );
   auto resultIt = result.begin();
   for( auto it = bodies_.begin(); it != bodies_.end(); ++it, ++resultIt )
   {
      *resultIt = (*it)->getPosition();
   }
   return result;
}

} // namespace pe
} // namespace walberla

