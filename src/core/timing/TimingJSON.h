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
//! \file Timer.h
//! \ingroup core
//!
//
//======================================================================================================================

#include "core/extern/json.hpp"
#include "core/timing/Timer.h"
#include "core/timing/TimingNode.h"
#include "core/timing/TimingTree.h"


namespace walberla {
namespace timing {


/*! Converts timer to json The signature is required by the json library
// \relates Timer
*/
template < typename TP > // Timing policy
void to_json( nlohmann::json& j, const Timer< TP >& timer )
{
   j = nlohmann::json{{"total", timer.total()},
                      {"average", timer.average()},
                      {"count", timer.getCounter()},
                      {"min", timer.min()},
                      {"max", timer.max()},
                      {"variance", timer.variance()}};
}

/// Converts a TimingNode to json. The signature is required by the json library
/// \relates TimingNode
template < typename TP > // Timing policy
void to_json( nlohmann::json& j, const TimingNode< TP >& tn )
{
   /// ignore the first timer in the timing node since it is empty
   if( tn.last_ == nullptr )
   {
      j = nlohmann::json( tn.tree_ );
   } else
   {
      j           = nlohmann::json( tn.timer_ );
      j["childs"] = nlohmann::json( tn.tree_ );
   }
}

/// Converts a TimingTree to json. The signature is required by the json library
/// \relates TimingTree
template < typename TP > // Timing policy
void to_json( nlohmann::json& j, const TimingTree< TP >& tt )
{
   j = nlohmann::json( tt.getRawData() );
}


}
}