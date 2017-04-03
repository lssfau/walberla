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
//! \file FieldExportTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "python_coupling/PythonWrapper.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Vector2.h"
#include "core/math/Vector3.h"
#include "core/math/Random.h"

#include "blockforest/Initialization.h"
#include "blockforest/python/Exports.h"
#include "field/python/Exports.h"
#include "field/iterators/IteratorMacros.h"
#include "python_coupling/Manager.h"
#include "python_coupling/PythonCallback.h"
#include "python_coupling/DictWrapper.h"
#include "stencil/D2Q9.h"


#include <boost/mpl/vector.hpp>

using namespace walberla;


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment mpiEnv( argc, argv );

   auto pythonManager = python_coupling::Manager::instance();
   typedef boost::mpl::vector<
           Field<Vector2<int>,1>,
           Field<Vector3<int>,1>,
           Field<int,2>,
           Field<int,3> > FieldTypes;

   typedef boost::mpl::vector< stencil::D2Q9> Stencils;

   pythonManager->addExporterFunction( field::exportModuleToPython<FieldTypes> );
   pythonManager->addBlockDataConversion< FieldTypes >() ;
   pythonManager->addExporterFunction( blockforest::exportModuleToPython<Stencils> );


   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid( 1,1,1, 20,20,1, real_t(1.0), false, true,true,true );

   auto sca2FieldID = field::addToStorage< GhostLayerField<int,2> >( blocks, "sca2Field", int(0), field::fzyx, 1 );
   auto sca3FieldID = field::addToStorage< GhostLayerField<int,3> >( blocks, "sca3Field", int(0), field::fzyx, 1 );

   auto vec2FieldID = field::addToStorage< GhostLayerField<Vector2<int>,1> >( blocks, "vec2Field", Vector2<int>(), field::zyxf, 1 );
   auto vec3FieldID = field::addToStorage< GhostLayerField<Vector3<int>,1> >( blocks, "vec3Field", Vector3<int>(), field::zyxf, 1 );

   // random init
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto sca2Field = blockIt->getData<GhostLayerField<int,2> >( sca2FieldID );
      auto sca3Field = blockIt->getData<GhostLayerField<int,3> >( sca3FieldID );
      for( auto cellIt = sca2Field->begin(); cellIt != sca2Field->end(); ++cellIt )
         *cellIt = math::intRandom( int(0), int(42) );
      for( auto cellIt = sca3Field->begin(); cellIt != sca3Field->end(); ++cellIt )
         *cellIt = math::intRandom( int(0), int(42) );
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
      auto sca2Field = blockIt->getData<GhostLayerField<int,2> >( sca2FieldID );
      auto sca3Field = blockIt->getData<GhostLayerField<int,3> >( sca3FieldID );
      auto vec2Field = blockIt->getData<GhostLayerField<Vector2<int>,1 > >( vec2FieldID );
      auto vec3Field = blockIt->getData<GhostLayerField<Vector3<int>,1 > >( vec3FieldID );

      {
         WALBERLA_FOR_ALL_CELLS(scaIt, sca2Field, vecIt, vec2Field,
            WALBERLA_ASSERT_EQUAL(scaIt[0], (*vecIt)[0]);
            WALBERLA_ASSERT_EQUAL(scaIt[1], (*vecIt)[1]);
			)
      }

      {
		  WALBERLA_FOR_ALL_CELLS(scaIt, sca3Field, vecIt, vec3Field,
            WALBERLA_ASSERT_EQUAL(scaIt[0], (*vecIt)[0]);
            WALBERLA_ASSERT_EQUAL(scaIt[1], (*vecIt)[1]);
            WALBERLA_ASSERT_EQUAL(scaIt[2], (*vecIt)[2]);
			)
      }
   }

   return 0;
}
