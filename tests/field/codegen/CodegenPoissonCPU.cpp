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
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "stencil/D2Q9.h"
#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"


using namespace walberla;

typedef GhostLayerField<real_t, 1> ScalarField_T;

// U with Dirichlet Boundary
void initU( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & srcId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      if( blocks->atDomainYMaxBorder( *block ) )
      {
         ScalarField_T * src = block->getData< ScalarField_T >( srcId );
         CellInterval xyz = src->xyzSizeWithGhostLayer();
         xyz.yMin() = xyz.yMax();
         for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
         {
            const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );
            src->get( *cell ) = std::sin( math::pi * p[0] ) * std::sinh( math::pi * p[1] );
         }
      }
   }
}

// right hand side
void initF( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & fId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      ScalarField_T * f = block->getData< ScalarField_T >( fId );
      CellInterval xyz = f->xyzSize();
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell )
      {
         f->get( *cell ) = 0.0;
      }
   }
}

void testPoisson()
{

   const uint_t xCells = uint_t(100);
   const uint_t yCells = uint_t(100);
   const real_t xSize = real_t(1);
   const real_t ySize = real_t(1);
   const real_t dx = xSize / real_c( xCells + uint_t(1) );
   const real_t dy = ySize / real_c( yCells + uint_t(1) );

   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
           math::AABB( real_t(0.5) * dx, real_t(0.5) * dy, real_t(0),
                       xSize - real_t(0.5) * dx, ySize - real_t(0.5) * dy, dx ),
           uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
           xCells, yCells, uint_t(1),          // how many cells per block (x,y,z)
           false,                              // one block per process - "false" means all blocks to one process
           false, false, false );              // no periodicity


   BlockDataID fieldID = field::addToStorage<ScalarField_T>(blocks, "Field", real_t(0.0));
   initU( blocks, fieldID );

   BlockDataID fId = field::addToStorage< ScalarField_T >( blocks, "f", real_t(0.0));
   initF( blocks, fId );

   typedef blockforest::communication::UniformBufferedScheme<stencil::D2Q9> CommScheme;
   typedef field::communication::PackInfo<ScalarField_T> Packing;
   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(fieldID) );

   // Create Timeloop
   const uint_t numberOfTimesteps = uint_t(10000);
   SweepTimeloop timeloop ( blocks, numberOfTimesteps );

   // Registering the sweep
   timeloop.add() << BeforeFunction(  commScheme, "Communication" )
                  << Sweep( pystencils::Poisson(fId, fieldID, dx, dy), "Poisson Kernel" );

   timeloop.run();

   auto firstBlock = blocks->begin();
   auto f = firstBlock->getData<ScalarField_T>( fieldID );
   WALBERLA_CHECK_LESS(f->get(50,99,0) - std::sin( math::pi  * 0.5 ) * std::sinh( math::pi * 0.99 ), 0.01);
}


int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testPoisson();

   return EXIT_SUCCESS;
}
