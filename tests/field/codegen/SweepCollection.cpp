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
//! \file SweepCollection.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "timeloop/SweepTimeloop.h"
#include "SweepCollection.h"

using namespace walberla;

typedef GhostLayerField<real_t, 1> ScalarField;
using SweepCollection_T = pystencils::SweepCollection;

void testSweepCollection()
{
   uint_t xSize = 20;
   uint_t ySize = 20;
   uint_t zSize = 20;
   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
           uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
           xSize, ySize, zSize,            // how many cells per block (x,y,z)
           real_c(1.0),                          // dx: length of one cell in physical coordinates
           false,                              // one block per process - "false" means all blocks to one process
           true, true, true );                 // full periodicity


   const real_t initField1 = real_c(1.0);
   const real_t initField2 = real_c(0.0);
   const real_t initField3 = real_c(0.0);
   const real_t a = real_c(2.0);

   const BlockDataID field1ID = field::addToStorage<ScalarField>(blocks, "Field1", initField1);
   const BlockDataID field2ID = field::addToStorage<ScalarField>(blocks, "Field2", initField2);
   const BlockDataID field3ID = field::addToStorage<ScalarField>(blocks, "Field3", initField3);

   SweepCollection_T sweepCollection(blocks, field1ID, field2ID, field3ID, a);

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(100);
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << Sweep( sweepCollection.fct1(SweepCollection_T::ALL), "fc1" );
   timeloop.add() << Sweep( sweepCollection.fct2(SweepCollection_T::ALL), "fc2" );

   timeloop.run();

   auto firstBlock = blocks->begin();
   auto field1 = firstBlock->getData<ScalarField>( field1ID );
   auto field2 = firstBlock->getData<ScalarField>( field2ID );
   auto field3 = firstBlock->getData<ScalarField>( field3ID );

   WALBERLA_CHECK_FLOAT_EQUAL(field1->get(0,0,0), initField1)
   WALBERLA_CHECK_FLOAT_EQUAL(field2->get(0,0,0), initField1 * real_c(2.0) * a)
   WALBERLA_CHECK_FLOAT_EQUAL(field3->get(0,0,0), initField1 * real_c(2.0) * a * real_c(2.0) * a)
}


int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testSweepCollection();
   return EXIT_SUCCESS;
}
