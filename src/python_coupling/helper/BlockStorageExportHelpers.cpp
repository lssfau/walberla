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
//! \file BlockStorageExportHelpers.cpp
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "BlockStorageExportHelpers.h"

#include "python_coupling/PythonWrapper.h"

namespace walberla {
namespace python_coupling {
namespace py = pybind11;

#ifdef WALBERLA_BUILD_WITH_PYTHON

void NoSuchBlockData::translate(  const NoSuchBlockData & e ) {
   throw py::cast_error(e.what());
}

void BlockDataNotConvertible::translate(  const BlockDataNotConvertible & e ) {
   throw py::cast_error(e.what());
}
#else

void NoSuchBlockData::translate(  const NoSuchBlockData &  ) {}

void BlockDataNotConvertible::translate(  const BlockDataNotConvertible &  ) {}

#endif



BlockDataID blockDataIDFromString( BlockStorage & bs, const std::string & stringID )
{
   auto ids = bs.getBlockDataIdentifiers();

   for ( uint_t i =0; i < ids.size(); ++i )
      if ( ids[i] == stringID )
         return BlockDataID( i );

   throw NoSuchBlockData();

}


BlockDataID blockDataIDFromString( IBlock & block, const std::string & stringID )
{
   return blockDataIDFromString ( block.getBlockStorage(), stringID );
}

BlockDataID blockDataIDFromString( StructuredBlockStorage & bs, const std::string & stringID )
{
   return blockDataIDFromString( bs.getBlockStorage(), stringID );
}


} // namespace python_coupling
} // namespace walberla


