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
//! \file InfoCollection.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "BlockInfo.h"

#include "blockforest/BlockForest.h"

#include <map>

namespace walberla {
namespace pe {

typedef std::map<blockforest::BlockID, BlockInfo>  InfoCollection;
typedef std::pair<blockforest::BlockID, BlockInfo> InfoCollectionPair;

/**
 * @brief Fills \a InfoCollection with up to date information.
 *
 * The number of local particles is used as the computational weight.
 * The number of shadow particles is used as the communication weight.
 *
 */
void createWithNeighborhoodLocalShadow( const BlockForest& bf,
                                        const BlockDataID storageID,
                                        InfoCollection& ic );

/**
 * @brief Fills \a InfoCollection with up to date information.
 *
 * The number of contacts is used as the computational weight.
 * The number of shadow particles is used as the communication weight.
 *
 */
void createWithNeighborhoodContactsShadow( BlockForest& bf,
                                           const BlockDataID storageID,
                                           const BlockDataID fcdID,
                                           InfoCollection& ic );

}
}
