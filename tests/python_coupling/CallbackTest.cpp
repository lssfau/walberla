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
//! \file CallbackTest.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//!
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/Field.h"

#include "python_coupling/DictWrapper.h"
#include "python_coupling/Manager.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/export/FieldExports.h"

using namespace walberla;


int main( int argc, char ** argv )
{
   auto pythonManager = python_coupling::Manager::instance();
   pythonManager->addExporterFunction( field::exportModuleToPython<Field<int,1>> );

   if ( argc != 2 ) {
      WALBERLA_ABORT_NO_DEBUG_INFO("Wrong parameter count: \nUsage: \n ./CallbackTest CallbackTest.py");
   }
   std::string pythonFile ( argv[1] );

   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );
   MPIManager::instance()->useWorldComm();


   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
             3,   2,  4,         // number of blocks in x,y,z direction
             10,  8, 12,         // how many cells per block (x,y,z)
             0.5,                // dx: length of one cell in physical coordinates
             false,              // one block per process - "false" means all blocks to one process
             false,false,false); // no periodicity


   // This callback should sum up the two given integers
   python_coupling::PythonCallback cb1 ( pythonFile, "cb1" );
   WALBERLA_ASSERT( cb1.isCallable() );
   cb1.data().exposeValue("input1", 5);
   cb1.data().exposeValue("input2", 10);
   cb1();
   int result = cb1.data().get<int>( "returnValue" );
   WALBERLA_CHECK_EQUAL( result, 15 );

   typedef GhostLayerField<int,1> ScalarField;
   ScalarField f ( 3,2,2, 1, 25 );
   python_coupling::PythonCallback cb2 ( pythonFile, "cb2" );
   WALBERLA_ASSERT( cb2.isCallable() );
   cb2.data().exposePtr("field", &f );
   cb2();
   WALBERLA_CHECK_EQUAL( f(0,0,0), 42 );
   WALBERLA_CHECK_EQUAL( f(-1,-1,-1), 5 );


   //python_coupling::Shell shell("MyGreatInputShell");
   //shell();

   return 0;
}
