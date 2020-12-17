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
//! \file FieldExportTest.cpp
//! //! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Random.h"
#include "core/math/Vector2.h"
#include "core/math/Vector3.h"

#include "python_coupling/DictWrapper.h"
#include "python_coupling/Manager.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/PythonWrapper.h"
#include "python_coupling/export/BlockForestExport.h"
#include "python_coupling/export/FieldExports.h"

#include "stencil/D2Q9.h"

using namespace walberla;



int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );

   auto pythonManager = python_coupling::Manager::instance();

   pythonManager->addExporterFunction( field::exportModuleToPython<Field<int, 3>, Field<real_t, 3>> );
   pythonManager->addBlockDataConversion< Field<int, 3>, Field<real_t, 3> >() ;
   pythonManager->addExporterFunction( blockforest::exportModuleToPython<stencil::D2Q9> );


   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid( 1,1,1, 20,20,1, real_t(1.0), false, true,true,true );

   auto srcIntFieldID = field::addToStorage< GhostLayerField<int, 3> >( blocks, "srcIntFieldID", int(0), field::fzyx, 1 );
   auto dstIntFieldID = field::addToStorage< GhostLayerField<int, 3> >( blocks, "dstIntFieldID", int(0), field::fzyx, 1 );

   auto srcDoubleFieldID = field::addToStorage< GhostLayerField<real_t, 3> >( blocks, "srcDoubleFieldID", real_t(0.0), field::fzyx, 1 );
   auto dstDoubleFieldID = field::addToStorage< GhostLayerField<real_t, 3> >( blocks, "dstDoubleFieldID", real_t(0.0), field::fzyx, 1 );

   // random init
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto srcIntField = blockIt->getData<GhostLayerField<int, 3> >( srcIntFieldID );
      auto srcDoubleField = blockIt->getData<GhostLayerField<real_t, 3> >( srcDoubleFieldID );
      for( auto cellIt = srcIntField->begin(); cellIt != srcIntField->end(); ++cellIt )
         *cellIt = math::intRandom( int(0), int(42) );
      for( auto cellIt = srcDoubleField->begin(); cellIt != srcDoubleField->end(); ++cellIt )
         *cellIt = math::realRandom( real_t(0.0), real_t(42.0) );
   }

   // call python function which should copy over the values to the Vector fields
   std::string pythonFile ( argv[1] );
   python_coupling::PythonCallback cb ( pythonFile, "theCallback" );
   WALBERLA_ASSERT( cb.isCallable() );
   cb.data().exposeValue("blocks", blocks);
   cb();


   // check for equivalence
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto srcIntField = blockIt->getData<GhostLayerField<int, 3> >( srcIntFieldID );
      auto dstIntField = blockIt->getData<GhostLayerField<int, 3> >( dstIntFieldID );
      auto srcDoubleField = blockIt->getData<GhostLayerField<real_t, 3> >( srcDoubleFieldID );
      auto dstDoubleField = blockIt->getData<GhostLayerField<real_t, 3> >( dstDoubleFieldID );

      {
         for(cell_idx_t z = 0; z < cell_idx_c(srcIntField->zSize()); ++z)
            for(cell_idx_t y = 0; y < cell_idx_c(srcIntField->zSize()); ++y)
               for(cell_idx_t x = 0; x < cell_idx_c(srcIntField->zSize()); ++x)
               {
                  WALBERLA_CHECK_EQUAL( srcIntField->get(x,y,z, 0), dstIntField->get(x,y,z, 0) );
                  WALBERLA_CHECK_EQUAL( srcIntField->get(x,y,z, 1), dstIntField->get(x,y,z, 1) );
               }
         for(cell_idx_t z = 0; z < cell_idx_c(srcDoubleField->zSize()); ++z)
            for(cell_idx_t y = 0; y < cell_idx_c(srcDoubleField->zSize()); ++y)
               for(cell_idx_t x = 0; x < cell_idx_c(srcDoubleField->zSize()); ++x)
               {
                  WALBERLA_CHECK_FLOAT_EQUAL( srcDoubleField->get(x,y,z, 0), dstDoubleField->get(x,y,z, 0) );
                  WALBERLA_CHECK_FLOAT_EQUAL( srcDoubleField->get(x,y,z, 1), dstDoubleField->get(x,y,z, 1) );
                  WALBERLA_CHECK_FLOAT_EQUAL( srcDoubleField->get(x,y,z, 2), dstDoubleField->get(x,y,z, 2) );
               }
      }
   }


   return 0;
}
