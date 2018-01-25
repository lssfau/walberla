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
//! \brief Singleton managing Python coupling
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"
#include "PythonWrapper.h"
#include "python_coupling/helper/MplHelpers.h"
#include "core/singleton/Singleton.h"
#include "domain_decomposition/IBlock.h"

void init_module_walberla_cpp();

namespace walberla {
namespace python_coupling {


   class BlockDataToObjectTester
   {
   public:
      BlockDataToObjectTester( IBlock * block, BlockDataID blockDataID )
         : block_( block), blockDataID_( blockDataID )
      {}

      template<typename TypeToTest>
      void operator() ( NonCopyableWrap<TypeToTest> )
      {
         using boost::python::object;
         if ( result_ == object() &&  block_->isDataClassOrSubclassOf<TypeToTest>( blockDataID_ ) ) {
            result_ = object( boost::python::ptr( block_->getData<TypeToTest>( blockDataID_ ) ) );
         }
      }

      boost::python::object getResult() { return result_; }
   private:
      IBlock * block_;
      BlockDataID blockDataID_;
      boost::python::object result_;
   };

   template<typename TypeList>
   boost::python::object testBlockData( IBlock & block, BlockDataID blockDataID )
   {
      BlockDataToObjectTester tester( &block, blockDataID );
      for_each_noncopyable_type< TypeList > ( boost::ref(tester) );
      return tester.getResult();
   }


   class Manager : public singleton::Singleton< Manager >
   {
   public:
      WALBERLA_BEFRIEND_SINGLETON;

      typedef std::function<void()> ExporterFunction;
      typedef std::function< boost::python::object ( IBlock&, BlockDataID ) > BlockDataToObjectFunction;


      ~Manager();

      void triggerInitialization();

      void addEntryToPythonPath( const std::string & path );


      void addExporterFunction( const ExporterFunction & f ) { exporterFunctions_.push_back( f ); }

      template<typename TypeList>
      void addBlockDataConversion() { blockDataToObjectFunctions_.push_back( &testBlockData<TypeList>  ); }


   protected:
      void addPath( const std::string & path );
      friend void ::init_module_walberla_cpp();
      void exportAll();

      friend boost::python::object IBlock_getData( boost::python::object, const std::string &  );
      boost::python::object pythonObjectFromBlockData( IBlock & block, BlockDataID  id );

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


