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
//! \file EK.cpp
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "EKContinuity.h"
#include "EKFlux.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Constants.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "stencil/D2Q9.h"

#include "timeloop/SweepTimeloop.h"


using namespace walberla;

typedef GhostLayerField<real_t, 1> DensityField_T;
typedef GhostLayerField<real_t, 2> FluxField_T;

void initC(const shared_ptr< StructuredBlockStorage > & blocks, BlockDataID cID)
{
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto c = blockIt->getData<DensityField_T>( cID );
      uint_t size = c->ySize() * c->xSize();
      c->get( cell_idx_c(c->ySize()/2), cell_idx_c(c->ySize()/2), 0 ) = real_c(size)/real_c(4);
      c->get( cell_idx_c(c->ySize()/2), cell_idx_c(c->ySize()/2)-1, 0 ) = real_c(size)/real_c(4);
      c->get( cell_idx_c(c->ySize()/2)-1, cell_idx_c(c->ySize()/2), 0 ) = real_c(size)/real_c(4);
      c->get( cell_idx_c(c->ySize()/2)-1, cell_idx_c(c->ySize()/2)-1, 0 ) = real_c(size)/real_c(4);
   }
}

void checkC(const shared_ptr< StructuredBlockStorage > & blocks, BlockDataID cID, real_t D, uint_t t)
{
   real_t sum(0.0);
   uint_t size(0);
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto c = blockIt->getData<DensityField_T>( cID );
      size += c->ySize() * c->xSize();
   }
   
   std::cout << "#cell pos actual expected deviation" << std::endl;
   for(auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto c = blockIt->getData<DensityField_T>( cID );
      
      for( cell_idx_t y = 0; y < cell_idx_c( c->ySize() ); ++y )
      {
         for( cell_idx_t x = 0; x < cell_idx_c( c->xSize() ); ++x )
         {
            const real_t val = c->get(x, y, 0);
            sum += val;
            
            const real_t x0 = real_c(x) - real_c(c->xSize()/2) + real_c(0.5);
            const real_t y0 = real_c(y) - real_c(c->ySize()/2) + real_c(0.5);
            const real_t r2 = x0*x0 + y0*y0;
            
            // solution to the diffusion equation in 2D
            const real_t ref = real_c(size)/(real_c(4.0)*math::pi*D*real_c(t)) * std::exp(-r2/real_c(4*D*real_c(t)));
            
            real_t rel = std::abs(real_c(1) - val/ref);
            if(ref >= 1) // we only expect good agreement where the values are larger than sum/size
            {
               if(x == cell_idx_c(c->ySize()/2)) // print out the values in the middle
               {
                  std::cout << y << " " << y0 << " " << val << " " << ref << " " << (rel*100) << "%" << std::endl;
               }
               
               WALBERLA_CHECK_LESS(rel, 0.03) // 3% deviation is okay
            }
         }
      }
   }
   WALBERLA_CHECK_FLOAT_EQUAL(sum, real_c(size))
}

void testEK()
{
   uint_t xSize = 100;
   uint_t ySize = 100;
   real_t D = real_c(0.1);
   
   // Create blocks
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid (
      uint_t(1) , uint_t(1),  uint_t(1),  // number of blocks in x,y,z direction
      xSize, ySize, uint_t(1),            // how many cells per block (x,y,z)
      real_c(1.0),                          // dx: length of one cell in physical coordinates
      false,                              // one block per process - "false" means all blocks to one process
      true, true, true );                 // full periodicity

   BlockDataID c = field::addToStorage<DensityField_T>(blocks, "c", real_c(0.0), field::fzyx);
   initC(blocks, c);
   BlockDataID j = field::addToStorage<FluxField_T>(blocks, "j", real_c(0.0), field::fzyx);
   
   typedef blockforest::communication::UniformBufferedScheme<stencil::D2Q9> CommScheme;
   typedef field::communication::PackInfo<DensityField_T> Packing;
   CommScheme commScheme(blocks);
   commScheme.addDataToCommunicate( make_shared<Packing>(c) );
   
   const uint_t numberOfTimesteps = uint_t(200);
   SweepTimeloop timeloop(blocks, numberOfTimesteps);

   // Registering the sweep
   timeloop.add() << BeforeFunction( commScheme, "Communication" )
                  << Sweep( pystencils::EKFlux(c, j, D), "EK flux Kernel" );
   timeloop.add() << Sweep( pystencils::EKContinuity(c, j), "EK continuity Kernel" );

   timeloop.run();
   
   checkC(blocks, c, D, numberOfTimesteps);
}



int main( int argc, char ** argv )
{
   mpi::Environment env( argc, argv );
   debug::enterTestMode();

   testEK();

   return EXIT_SUCCESS;
}
