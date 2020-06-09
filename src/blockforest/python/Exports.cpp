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
#include "blockforest/SetupBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "core/logging/Logging.h"
#include "core/StringUtility.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "python_coupling/Manager.h"
#include "python_coupling/helper/ConfigFromDict.h"

#include "stencil/D3Q7.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include <memory>

#include <sstream>

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

bool checkForThreeTuple( object obj ) //NOLINT
{
   if( ! extract<tuple> ( obj ).check() )
      return false;

   tuple t = extract<tuple> ( obj );
   return len(t) == 3;
}


object python_createUniformBlockGrid(tuple args, dict kw) //NOLINT
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
      shared_ptr< StructuredBlockForest > blocks = createUniformBlockGridFromConfig( cfg->getGlobalBlock(), nullptr, keepGlobalBlockInformation );
      return object(blocks);
   }
   catch( std::exception & e)
   {
      PyErr_SetString( PyExc_ValueError, e.what() );
      throw boost::python::error_already_set();
   }

}

shared_ptr<StructuredBlockForest> createStructuredBlockForest( Vector3<uint_t> blocks,
                                                               Vector3<uint_t> cellsPerBlock,
                                                               Vector3<bool> periodic,
                                                               object blockExclusionCallback = object(),
                                                               object workloadMemoryCallback = object(),
                                                               object refinementCallback = object(),
                                                               const real_t dx = 1.0,
                                                               memory_t processMemoryLimit = std::numeric_limits<memory_t>::max(),
                                                               const bool keepGlobalBlockInformation = false)
{
   using namespace blockforest;
   Vector3<real_t> bbMax;
   for( uint_t i=0; i < 3; ++i )
      bbMax[i] = real_c( blocks[i] * cellsPerBlock[i] ) * dx;
   AABB domainAABB ( Vector3<real_t>(0),  bbMax );

   SetupBlockForest sforest;

   auto blockExclusionFunc = [&blockExclusionCallback] ( std::vector<walberla::uint8_t>& excludeBlock, const SetupBlockForest::RootBlockAABB& aabb ) -> void
   {
      for( uint_t i = 0; i != excludeBlock.size(); ++i )
      {
         AABB bb = aabb(i);
         auto pythonReturnVal = blockExclusionCallback(bb);
         if( ! extract<bool>( pythonReturnVal ).check() ) {
            PyErr_SetString( PyExc_ValueError, "blockExclusionCallback has to return a boolean");
            throw boost::python::error_already_set();
         }

         bool returnVal = extract<bool>(pythonReturnVal);
         if ( returnVal )
            excludeBlock[i] = 1;
      }
   };

   auto workloadMemoryFunc = [&workloadMemoryCallback] ( SetupBlockForest & forest )-> void
   {
      std::vector< SetupBlock* > blockVector;
      forest.getBlocks( blockVector );

      for( uint_t i = 0; i != blockVector.size(); ++i ) {
         blockVector[i]->setMemory( memory_t(1) );
         blockVector[i]->setWorkload( workload_t(1) );
         workloadMemoryCallback( boost::python::ptr(blockVector[i]) );
      }
   };

   auto refinementFunc = [&refinementCallback] ( SetupBlockForest & forest )-> void
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         SetupBlock * sb = &(*block);
         auto pythonRes = refinementCallback( boost::python::ptr(sb) );
         if( ! extract<bool>( pythonRes ).check() ) {
            PyErr_SetString( PyExc_ValueError, "refinementCallback has to return a boolean");
            throw boost::python::error_already_set();
         }
         bool returnVal = extract<bool>( pythonRes );
         if( returnVal )
            block->setMarker( true );
      }
   };

   if ( blockExclusionCallback ) {
      if( !PyCallable_Check( blockExclusionCallback.ptr() ) ) {
         PyErr_SetString( PyExc_ValueError, "blockExclusionCallback has to be callable");
         throw boost::python::error_already_set();
      }
      sforest.addRootBlockExclusionFunction( blockExclusionFunc );
   }

   if ( workloadMemoryCallback ) {
      if( !PyCallable_Check( workloadMemoryCallback.ptr() ) ) {
         PyErr_SetString( PyExc_ValueError, "workloadMemoryCallback has to be callable");
         throw boost::python::error_already_set();
      }
      sforest.addWorkloadMemorySUIDAssignmentFunction( workloadMemoryFunc );
   }
   else
      sforest.addWorkloadMemorySUIDAssignmentFunction( uniformWorkloadAndMemoryAssignment );

   if ( refinementCallback ) {
      if( !PyCallable_Check( refinementCallback.ptr() ) ) {
         PyErr_SetString( PyExc_ValueError, "refinementCallback has to be callable");
         throw boost::python::error_already_set();
      }
      sforest.addRefinementSelectionFunction( refinementFunc );
   }

   sforest.init( domainAABB, blocks[0], blocks[1], blocks[2], periodic[0], periodic[1], periodic[2] );

   // calculate process distribution
   sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalanceWeighted(),
                        uint_c( MPIManager::instance()->numProcesses() ),
                        real_t(0), processMemoryLimit );

   if( !MPIManager::instance()->rankValid() )
      MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation );

   auto sbf = std::make_shared< StructuredBlockForest >( bf, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] );
   sbf->createCellBoundingBoxes();

   return sbf;
}

object createUniformNeighborScheme(  const shared_ptr<StructuredBlockForest> & bf,
                                     const std::string & stencil )
{
   if ( string_icompare(stencil, "D3Q7") == 0 )
      return object ( make_shared< UniformBufferedScheme<stencil::D3Q7> > ( bf ) );
   if ( string_icompare(stencil, "D3Q19") == 0 )
      return object ( make_shared< UniformBufferedScheme<stencil::D3Q19> > ( bf ) );
   if ( string_icompare(stencil, "D3Q27") == 0 )
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

std::string printSetupBlock(const SetupBlock & b )
{
   std::stringstream out;
   out <<  "SetupBlock at " << b.getAABB();
   return out.str();
}



void exportBlockForest()
{
   class_< StructuredBlockForest, //NOLINT
           shared_ptr<StructuredBlockForest>,
           bases<StructuredBlockStorage>, boost::noncopyable > ( "StructuredBlockForest", no_init );

   class_< SetupBlock, boost::noncopyable > ( "SetupBlock", no_init )
            .add_property("level", &SetupBlock::getLevel)
            .add_property("workload", &SetupBlock::getWorkload, &SetupBlock::setWorkload)
            .add_property("memory",  &SetupBlock::getMemory, &SetupBlock::setMemory)
            .add_property("aabb", make_function(&SetupBlock::getAABB, return_value_policy<copy_const_reference>()))
            .def("__repr__", &printSetupBlock)
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

   def( "createCustomBlockGrid", createStructuredBlockForest,
               (arg("blocks"), arg("cellsPerBlock"), arg("periodic"),
                arg("blockExclusionCallback") = object(),
                arg("workloadMemoryCallback") = object(),
                arg("refinementCallback") = object() ,
                arg("dx") = 1.0,
                arg("processMemoryLimit") = std::numeric_limits<memory_t>::max(),
                arg("keepGlobalBlockInformation") = false ) );
}

} // namespace domain_decomposition
} // namespace walberla


#endif //WALBERLA_BUILD_WITH_PYTHON
