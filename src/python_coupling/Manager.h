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
//! \file Manager.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//! \brief Singleton managing Python coupling
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"
#include "PythonWrapper.h"
#include "python_coupling/helper/MplHelpers.h"
#include "core/singleton/Singleton.h"
#include "domain_decomposition/IBlock.h"
#include <pybind11/pybind11.h>

void init_module_walberla_cpp();

namespace walberla {
namespace python_coupling {
namespace py = pybind11;


   class BlockDataToObjectTester
   {
   public:
      BlockDataToObjectTester( IBlock * block, BlockDataID blockDataID )
         : block_( block), blockDataID_( blockDataID )
      {}

      template<typename TypeToTest>
      void operator() ( NonCopyableWrap<TypeToTest> )
      {
         using py::object;
         if ( result_.is(py::object()) &&  block_->isDataClassOrSubclassOf<TypeToTest>( blockDataID_ ) ) {
            result_ = py::cast( block_->getData<TypeToTest>( blockDataID_ ) );
         }
      }

      py::object getResult() { return result_; }
   private:
      IBlock * block_;
      BlockDataID blockDataID_;
      py::object result_;
   };

   template<typename... Types>
   py::object testBlockData( IBlock & block, BlockDataID blockDataID )
   {
      BlockDataToObjectTester tester( &block, blockDataID );
      for_each_noncopyable_type< Types... > ( std::ref(tester) );
      return tester.getResult();
   }


   class Manager : public singleton::Singleton< Manager >
   {
   public:
      WALBERLA_BEFRIEND_SINGLETON;

      typedef std::function<void(py::module_&)> ExporterFunction;
      typedef std::function< py::object ( IBlock&, BlockDataID ) > BlockDataToObjectFunction;


      ~Manager();

      void triggerInitialization();

      void addEntryToPythonPath( const std::string & path );


      void addExporterFunction( const ExporterFunction & f)
      {
         exporterFunctions_.push_back( f );
      }

      template<typename... Types>
      void addBlockDataConversion() { blockDataToObjectFunctions_.push_back( &testBlockData<Types...>  ); }

      void exportAll(py::module_ &m);
    protected:
      void addPath( const std::string & path );
      friend void ::init_module_walberla_cpp();

      friend py::object IBlock_getData( py::object, const std::string &  );
      py::object pythonObjectFromBlockData( IBlock & block, BlockDataID  id );

      Manager();


      std::vector< ExporterFunction >          exporterFunctions_;
      std::vector< BlockDataToObjectFunction > blockDataToObjectFunctions_;

      /// Adding entries to python path is possible after initialization
      /// calls to addEntryToPythonPath() are store here and added after initialization
      std::vector<std::string> entriesForPythonPath_;

      bool initialized_;
   };




} // namespace python_coupling
} // namespace walberla


