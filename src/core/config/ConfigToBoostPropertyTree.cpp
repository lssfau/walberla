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
//! \file ConfigToBoostPropertyTree.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"
#ifdef WALBERLA_BUILD_WITH_BOOST

#include "ConfigToBoostPropertyTree.h"

#include <boost/property_tree/ptree.hpp>


namespace walberla {
namespace config {

boost::property_tree::iptree configBlockHandleToBoostPropertyTree( const Config::BlockHandle & blockHandle )
{
	boost::property_tree::iptree propTree;

	//typedef std::pair<const std::string &, const std::string &> PairType;
	for( Config::const_iterator it = blockHandle.begin(); it != blockHandle.end(); ++it )
		propTree.put( it->first, it->second );

	Config::Blocks blocks;
	blockHandle.getBlocks(blocks);

	for( const Config::BlockHandle & handle : blocks )
		propTree.add_child( handle.getKey(), configBlockHandleToBoostPropertyTree( handle ) );

	return propTree;
}

boost::property_tree::iptree configToBoostPropertyTree( const Config & config )
{
	return configBlockHandleToBoostPropertyTree( config.getGlobalBlock() );
}

} // namespace config
} // namespace walberla

#endif
