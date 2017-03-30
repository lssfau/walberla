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
//! \file PropertyTreeTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/config/Config.h"
#include "core/config/ConfigToBoostPropertyTree.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include <boost/property_tree/info_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>


void write_property_tree( std::ostringstream & ss, const boost::property_tree::iptree & ptree, const std::string & path)
{
	if(ptree.empty())
		ss << path << " = " << ptree.data() << "\n";
	else
		for( boost::property_tree::iptree::const_iterator it = ptree.begin(); it != ptree.end(); ++it )
			write_property_tree(ss, it->second, path.empty() ? it->first : path + "." + it->first);
}

int main( int argc, char** argv )
{
	walberla::debug::enterTestMode();
	walberla::MPIManager::instance()->initializeMPI(&argc, &argv);

	if( argc != 2 )
	   WALBERLA_ABORT( "Wrong number of Arguments!\nUsage: PropertyTreeTest <InputFileName>" );

	using walberla::Config;
	using boost::property_tree::iptree;

	Config config;

	config.readParameterFile( argv[1] );

	iptree propertyTree = walberla::config::configToBoostPropertyTree( config );

	std::ostringstream ss;

	//write_xml( ss, propertyTree );
	//write_info( ss, propertyTree );

	write_property_tree(ss, propertyTree, "");

	std::cout << ss.str() << std::endl;

	return 0;
}
