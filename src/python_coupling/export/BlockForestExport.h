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
//! \file BlockForestExport.h
//! \ingroup blockforest
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"


#ifdef WALBERLA_BUILD_WITH_PYTHON

#   include <pybind11/pybind11.h>

#   include "BlockForestCommunicationExport.h"

namespace walberla {
namespace blockforest {

   struct NoSuchBlockData : public std::out_of_range {
      explicit NoSuchBlockData ( ) : std::out_of_range( "No blockdata with the given name found" ) {}
      explicit NoSuchBlockData ( const std::string & w ) : std::out_of_range(w) {}
      static void translate( const NoSuchBlockData & e );
   };
   struct BlockDataNotConvertible : public std::runtime_error {
      explicit BlockDataNotConvertible (  ) : std::runtime_error( "This blockdata is not accessible from Python" ) {}
      explicit BlockDataNotConvertible ( const std::string & w ) : std::runtime_error(w) {}
      static void translate( const BlockDataNotConvertible &  e );
   };

   BlockDataID blockDataIDFromString( IBlock & block, const std::string & stringID );
   BlockDataID blockDataIDFromString( BlockStorage & bs, const std::string & stringID );
   BlockDataID blockDataIDFromString( StructuredBlockStorage & bs, const std::string & stringID );

   namespace py = pybind11;

   void exportBlockForest(py::module_ &m);


   template<typename... Stencils>
   void exportModuleToPython(py::module_ &m)
   {
      exportBlockForest(m);
      exportUniformBufferedScheme<Stencils...>(m);
      exportUniformDirectScheme<Stencils...>(m);
   }

} // namespace blockforest
} // namespace walberla


#endif // WALBERLA_BUILD_WITH_PYTHON
