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
//! \file ForceOnBodiesAdder.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "ForceOnBodiesAdder.h"

#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "pe/rigidbody/BodyIterators.h"


namespace walberla {
namespace pe_coupling {

void ForceOnBodiesAdder::operator()()
{
   for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
   {
      for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
      {
         if( !bodySelectorFct_(bodyIt.getBodyID()) ) continue;
         bodyIt->addForce ( force_ );
      }
   }
}

void ForceOnBodiesAdder::updateForce( const Vector3<real_t> & newForce )
{
   force_ = newForce;
}

} // namespace pe_coupling
} // namespace walberla
