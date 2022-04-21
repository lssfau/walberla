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
//! \file DirectionBasedReduceCommTest.cpp
//! \ingroup comm
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Compares distributed calculation (with comm) to one big local computation
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/DirectionBasedReduceScheme.h"

#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "domain_decomposition/MakeBlockDataInitFunction.h"

#include "field/AddToStorage.h"
#include "field/communication/ReducePackInfo.h"

#include "stencil/Directions.h"

#include "timeloop/SweepTimeloop.h"


namespace walberla {

typedef GhostLayerField<real_t,1> ScalarField;

template<bool init_>
class SumSweep
{
   public:
      SumSweep(BlockDataID valFieldID, BlockDataID sumFieldID)
         : valFieldID_(valFieldID), sumFieldID_(sumFieldID){}

      void operator() (IBlock * block) const
      {
         ScalarField * valField = block->getData<ScalarField>(valFieldID_);
         ScalarField * sumField = block->getData<ScalarField>(sumFieldID_);
         WALBERLA_ASSERT_NOT_NULLPTR( valField )
         WALBERLA_ASSERT_NOT_NULLPTR( sumField )

         if( init_ ){
            auto itSum0 = sumField->beginGhostLayerOnly(stencil::T);
            for( ; itSum0 != sumField->end(); ++itSum0 )
               *itSum0 = real_t(0);
         }

         auto itVal = valField->rbegin();
         auto itSum = sumField->rbegin();
         for( ; itVal != valField->rend(); ++itVal, ++itSum )
            *itSum = *itVal + itSum.neighbor( cell_idx_t(0), cell_idx_t(0), cell_idx_t(1) );
      }

   protected:
      BlockDataID valFieldID_;
      BlockDataID sumFieldID_;
};


class CompareSweep
{
   public:
      CompareSweep(BlockDataID smlFieldID, BlockDataID bigFieldID)
         : smlFieldID_(smlFieldID), bigFieldID_(bigFieldID)
      {}

      void operator() (IBlock * block)
      {
         ScalarField * sf = block->getData<ScalarField>(smlFieldID_);
         ScalarField * bf = block->getData<ScalarField>(bigFieldID_);
         WALBERLA_ASSERT_NOT_NULLPTR( sf )
         WALBERLA_ASSERT_NOT_NULLPTR( bf )

         const AABB & bb = block->getAABB();
         const cell_idx_t offset [3] = { cell_idx_c(bb.min(uint_t(0u))),
                                         cell_idx_c(bb.min(uint_t(1u))),
                                         cell_idx_c(bb.min(uint_t(2u))) };

         for(ScalarField::iterator i = sf->begin(); i != sf->end(); ++i )
         {
            real_t globalValue = bf->get(i.x()+offset[0],i.y()+offset[1], i.z()+offset[2], i.f());
            real_t localValue  = *i;
            WALBERLA_CHECK_FLOAT_EQUAL(localValue,globalValue)
         }
      }

   protected:
      BlockDataID smlFieldID_;
      BlockDataID bigFieldID_;
};

} // namespace walberla

int main(int argc, char **argv)
{
   using namespace walberla;
   debug::enterTestMode();
   mpi::Environment env( argc, argv );

   const uint_t cells [] = { uint_t(5u), uint_t(2u), uint_t(7u) };
   const uint_t blockCount [] = { uint_t(1u), uint_t(1u), uint_c( MPIManager::instance()->numProcesses() ) };
   const uint_t nrOfTimeSteps = uint_t(3u);
   bool periodic = false;
   const field::Layout layout = field::fzyx;

   // Create BlockForest
   uint_t glbCellsX = cells[0] * blockCount[0];
   uint_t glbCellsY = cells[1] * blockCount[1];
   uint_t glbCellsZ = cells[2] * blockCount[2];

   using blockforest::createUniformBlockGrid;

   auto blocks = createUniformBlockGrid(blockCount[0],blockCount[1],blockCount[2],  //blocks
                                        cells[0],cells[1],cells[2],                 //cells
                                        1,                                          //dx
                                        true,                                       //one block per process
                                        periodic,periodic,periodic, true);          //periodicity & global information


   // In addition to the normal GhostLayerField's  we allocated additionally a field containing the whole global simulation domain for each block
   // we can then check if the GhostLayer communication is correct, by comparing the small field to the corresponding part of the big field
   BlockDataID valFieldLoc = field::addToStorage<ScalarField>(blocks, "Val", real_t(1) );
   BlockDataID sumFieldLoc = field::addToStorage<ScalarField>(blocks, "Sum", real_t(0) );

   BlockDataID valFieldGlb = blocks->addBlockData<ScalarField>( makeBlockDataInitFunction<ScalarField>( glbCellsX, glbCellsY, glbCellsZ, uint_t(1), real_t(1), layout), "Global Src" );
   BlockDataID sumFieldGlb = blocks->addBlockData<ScalarField>( makeBlockDataInitFunction<ScalarField>( glbCellsX, glbCellsY, glbCellsZ, uint_t(1), real_t(0), layout), "Global Dst" );

   // small local fields
   blockforest::DirectionBasedReduceScheme<stencil::B> schemeB( blocks );
   schemeB.addPackInfo( make_shared<field::communication::ReducePackInfo< std::plus, ScalarField > >( sumFieldLoc, real_t(0) ) );

   // Create TimeLoop
   SweepTimeloop timeLoop (blocks->getBlockStorage(), nrOfTimeSteps);

   timeLoop.add() << Sweep ( CompareSweep( sumFieldLoc, sumFieldGlb ) );
   timeLoop.add() << Sweep ( SumSweep< true  >( valFieldLoc, sumFieldLoc ) );
   timeLoop.add() << BeforeFunction( schemeB )
                  << Sweep ( SumSweep< false >( valFieldLoc, sumFieldLoc ) );
   timeLoop.add() << Sweep ( SumSweep< true  >( valFieldGlb, sumFieldGlb ) );

   timeLoop.run();
}
