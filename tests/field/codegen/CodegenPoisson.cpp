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
//! \file CodegenPoisson.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "Poisson.h"
#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "gui/Gui.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q7.h"

#include "timeloop/SweepTimeloop.h"


using namespace walberla;

typedef GhostLayerField<double,1> ScalarField;


void testPoisson()
{
   uint_t L0 = 50;
   uint_t xSize = L0;
   uint_t ySize = L0;

   double h = 1 / (double(L0) + 1);
   double rhs = 1.0;
   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
           uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
           xSize, ySize, uint_t(1),            // how many cells per block (x,y,z)
           real_t(1),                          // dx: length of one cell in physical coordinates
           false,                              // one block per process - "false" means all blocks to one process
           false, false, false );              // no periodicity


   BlockDataID fieldID = field::addToStorage<ScalarField>(blocks, "Field", real_t(0.0));

   typedef blockforest::communication::UniformBufferedScheme<stencil::D2Q9> CommScheme;
   typedef field::communication::PackInfo<ScalarField> Packing;
   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(fieldID) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(2500);
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( pystencils::Poisson(fieldID, h, rhs), "Poisson Kernel" );

   timeloop.run();

   auto firstBlock = blocks->begin();
   auto f = firstBlock->getData<ScalarField>( fieldID );
   WALBERLA_CHECK_FLOAT_EQUAL(f->get(int(L0/2), int(L0/2), 0), real_t(7.28886456002774547e-02));
}


int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testPoisson();

   return 0;
}
