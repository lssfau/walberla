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
//! \file StabilityCheckerTest.cpp
//! \ingroup field
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "timeloop/all.h"



namespace walberla {

using Field_T = GhostLayerField<real_t, 1>;

class TestSweep
{
 public:
   TestSweep(BlockDataID fieldID) : fieldID_(fieldID) {}

   void operator()(IBlock* const block)
   {
      Field_T* field = block->getData< Field_T >(fieldID_);

      WALBERLA_FOR_ALL_CELLS(fieldIt, field, { *fieldIt += Field_T::value_type(1); }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   BlockDataID fieldID_;
};

int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   auto blocks = blockforest::createUniformBlockGrid( 1, 1, 1,
                                                      4, 4, 4,
                                                      1.0);

   BlockDataID fieldID = field::addToStorage<Field_T>( blocks, "Field", Field_T::value_type(0));
   SweepTimeloop timeloop(blocks->getBlockStorage(), uint_c(2));

   timeloop.add() << Sweep(TestSweep(fieldID), "Test Sweep");

   // LBM stability check
   auto checkFunction = [](Field_T::value_type value) {return value < math::abs(Field_T::value_type(5));};
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< Field_T >( blocks, fieldID, uint_c(1), checkFunction) ),"Stability check" );
   timeloop.run();

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   walberla::main(argc, argv);
}
