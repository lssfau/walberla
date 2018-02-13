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
//! \file PythonExports.cpp
//! \ingroup domain_decomposition
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/Initialization.h"
#include "blockforest/SetupBlock.h"
#include "core/logging/Logging.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "python_coupling/Manager.h"
#include "python_coupling/helper/ConfigFromDict.h"

#include "stencil/D3Q7.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include <boost/algorithm/string.hpp>

#ifdef _MSC_VER
#  pragma warning(push)
// disable warning boost/python/raw_function.hpp(55): warning C4267: 'argument' : conversion from 'size_t' to 'int', possible loss of data
#  pragma warning( disable : 4267 )
#endif //_MSC_VER
#include <boost/python/raw_function.hpp>
#ifdef _MSC_VER
#  pragma warning(pop)
#endif //_MSC_VER


using namespace boost::python;


namespace walberla {
namespace blockforest {

using walberla::blockforest::communication::UniformBufferedScheme;

bool checkForThreeTuple( object obj )
{
   if( ! extract<tuple> ( obj ).check() )
      return false;

   tuple t = extract<tuple> ( obj );
   return len(t) == 3;
}


object python_createUniformBlockGrid(tuple args, dict kw)
{
   if( len(args) > 0 ) {
      PyErr_SetString( PyExc_ValueError, "This function takes only keyword arguments" );
      throw boost::python::error_already_set();
   }

   using boost::python::stl_input_iterator;

   boost::python::list keys = kw.keys();
   for( auto it = stl_input_iterator<std::string>( keys ); it != stl_input_iterator<std::string>(); ++it )
   {
      if ( *it != "cells" &&
           *it != "cellsPerBlock" &&
           *it != "blocks" &&
           *it != "periodic" &&
           *it != "dx" &&
           *it != "oneBlockPerProcess"  )
      {
         PyErr_SetString( PyExc_ValueError, (std::string("Unknown Parameter: ") + (*it) ).c_str() );
         throw boost::python::error_already_set();
      }
   }

   if( kw.has_key("cells ") && ! checkForThreeTuple( kw["cells"] ) ) {
      PyErr_SetString( PyExc_ValueError, "Parameter 'cells' has to be tuple of length 3, indicating cells in x,y,z direction" );
      throw boost::python::error_already_set();
   }
   if( kw.has_key("cellsPerBlock ") && ! checkForThreeTuple( kw["cellsPerBlock"] ) ) {
      PyErr_SetString( PyExc_ValueError, "Parameter 'cellsPerBlock' has to be tuple of length 3, indicating cells in x,y,z direction" );
      throw boost::python::error_already_set();
   }
   if( kw.has_key("blocks ") && ! checkForThreeTuple( kw["blocks"] ) ) {
      PyErr_SetString( PyExc_ValueError, "Parameter 'blocks' has to be tuple of length 3, indicating cells in x,y,z direction" );
      throw boost::python::error_already_set();
   }

   bool keepGlobalBlockInformation = false;
   if ( kw.has_key("keepGlobalBlockInformation") )
   {
      if ( extract<bool>( kw["keepGlobalBlockInformation"] ).check() )
         keepGlobalBlockInformation = extract<bool>( kw["keepGlobalBlockInformation"] );
      else
      {
         PyErr_SetString( PyExc_ValueError, "Parameter 'keepGlobalBlockInformation' has to be a boolean" );
         throw boost::python::error_already_set();
      }
   }

   shared_ptr<Config> cfg = python_coupling::configFromPythonDict( kw );

   try {
      shared_ptr< StructuredBlockForest > blocks = createUniformBlockGridFromConfig( cfg->getGlobalBlock(), NULL, keepGlobalBlockInformation );
      return object(blocks);
   }
   catch( std::exception & e)
   {
      PyErr_SetString( PyExc_ValueError, e.what() );
      throw boost::python::error_already_set();
   }

}


object createUniformNeighborScheme(  const shared_ptr<StructuredBlockForest> & bf,
                                     const std::string & stencil )
{
   if ( boost::iequals(stencil, "D3Q7") )
      return object ( make_shared< UniformBufferedScheme<stencil::D3Q7> > ( bf ) );
   if ( boost::iequals(stencil, "D3Q19") )
      return object ( make_shared< UniformBufferedScheme<stencil::D3Q19> > ( bf ) );
   if ( boost::iequals(stencil, "D3Q27") )
      return object ( make_shared< UniformBufferedScheme<stencil::D3Q27> > ( bf ) );
   else {
      PyErr_SetString( PyExc_RuntimeError, "Unknown stencil. Allowed values 'D3Q27', 'D3Q19', 'D3Q7'");
      throw error_already_set();
      return object();
   }
}

template<typename Stencil>
void exportUniformBufferedScheme()
{
   typedef UniformBufferedScheme<Stencil> UNS;

   class_< UNS, shared_ptr<UNS>, boost::noncopyable >( "UniformBufferedScheme", no_init )
            .def( "__call__",             &UNS::operator()             )
            .def( "communicate",          &UNS::communicate            )
            .def( "startCommunication",   &UNS::startCommunication     )
            .def( "wait",                 &UNS::wait                   )
            .def( "addPackInfo",          &UNS::addPackInfo            )
            .def( "addDataToCommunicate", &UNS::addDataToCommunicate   )
    ;

}




void exportBlockForest()
{
   class_< StructuredBlockForest,
           shared_ptr<StructuredBlockForest>,
           bases<StructuredBlockStorage>, boost::noncopyable > ( "StructuredBlockForest", no_init );

   class_< SetupBlock, boost::noncopyable > ( "SetupBlock", no_init )
            .def( "getLevel",    &SetupBlock::getLevel    )
            .def( "getWorkload", &SetupBlock::getWorkload )
            .def( "setWorkload", &SetupBlock::setWorkload )
            .def( "getMemory",   &SetupBlock::getMemory   )
            .def( "setMemory",   &SetupBlock::setMemory   )
   ;

#ifdef _MSC_VER
#  pragma warning(push)
// disable warning boost/python/raw_function.hpp(55): warning C4267: 'argument' : conversion from 'size_t' to 'int', possible loss of data
#  pragma warning( disable : 4267 )
#endif //_MSC_VER
   def( "createUniformBlockGrid", raw_function(python_createUniformBlockGrid) );
#ifdef _MSC_VER
#  pragma warning(pop)
#endif //_MSC_VER

}

} // namespace domain_decomposition
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON
