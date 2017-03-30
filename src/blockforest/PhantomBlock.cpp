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
//! \file PhantomBlock.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockForest.h"
#include "PhantomBlock.h"
#include "PhantomBlockForest.h"



namespace walberla {
namespace blockforest {



PhantomBlock::NeighborBlock::NeighborBlock( const PhantomBlockForest & phantomForest, const BlockID & id, const uint_t process, const Set<SUID> & state ) :
   id_( id ), process_( process ), state_( state ), aabb_( phantomForest.getBlockForest().getAABBFromBlockId(id) )
{}



uint_t PhantomBlock::getProcess() const {

   return phantomForest_.getBlockForest().getProcess();
}



} // namespace blockforest
} // namespace walberla
