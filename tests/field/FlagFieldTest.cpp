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
//! \file FlagFieldTest.cpp
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "field/FlagField.h"
#include "field/Printers.h"
#include "core/debug/TestSubsystem.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include <iostream>
#include <string>
#include <vector>


using namespace walberla;
namespace wlb = walberla;

using std::vector;
using std::string;
using std::cout;
using std::endl;

void registerTest()
{
   FlagField<walberla::uint8_t> ff (1,1,1,1);

   ff.registerFlag("Flag1");

   ff.registerFlag("Flag5",4);
   ff.registerFlag("Flag3",2);

   ff.registerFlag("Flag2");
   ff.registerFlag("Flag4");
   ff.registerFlag("Flag6");
   ff.registerFlag("Flag7");
   ff.registerFlag("Flag8");

   // Test that registration of more flags
   // than possible results in an error
   bool overFlow=false;
   try {
      ff.registerFlag("invalid");
   }
   catch( std::runtime_error & /*err*/) {
      overFlow = true;
   }
   WALBERLA_CHECK(overFlow);

   vector<string> names;
   names.emplace_back("Flag1");
   names.emplace_back("Flag2");
   names.emplace_back("Flag3");
   names.emplace_back("Flag4");
   names.emplace_back("Flag5");
   names.emplace_back("Flag6");
   names.emplace_back("Flag7");
   names.emplace_back("Flag8");

   for(size_t i=0; i<names.size(); ++i)
   {
      WALBERLA_CHECK(ff.flagExists(names[i]));
      WALBERLA_CHECK(ff.flagExists(i));
      WALBERLA_CHECK_EQUAL(ff.getFlag(names[i]), 1 << i );
      WALBERLA_CHECK_EQUAL( ff.getFlagUID( numeric_cast<FlagField<walberla::uint8_t>::flag_t>( 1 << i ) ), FlagUID(names[i]) );
   }

   //ff.printRegistered(cout);
}

void accessTest()
{
   FlagField<wlb::uint8_t> ff (3,3,3,1);
   WALBERLA_CHECK_EQUAL( ff.xSize(), 3 );

   wlb::uint8_t flag = ff.registerFlag("MyFlag");

   WALBERLA_CHECK_EQUAL( ff.get(0,0,0), 0);
   ff.addFlag(0,0,0,flag);

   WALBERLA_CHECK (ff.isFlagSet(0,0,0,flag));

   ff.removeFlag(0,0,0,flag);
   WALBERLA_CHECK (! ff.isFlagSet(0,0,0,flag) );

   WALBERLA_CHECK( ff.get(0,0,0) == 0);
}

void iteratorTest()
{
   FlagField<wlb::uint8_t> ff (3,3,3,1);
   unsigned char flag = ff.registerFlag("MyFlag");

   for( auto i = ff.begin(); i != ff.end(); ++i )
      addFlag(i, flag);

   for( auto i = ff.begin(); i != ff.end(); ++i )
      WALBERLA_CHECK( isFlagSet(i, flag) );
}

void shallowCopyTest()
{

   // Test shallow copy
   using FField = FlagField<wlb::uint8_t>;
   FField ff ( 3,3,3,1 );
   ff.registerFlag("FirstFlag");

   shared_ptr< FField > sliced = shared_ptr< FField > ( ff.getSlicedField(CellInterval(0,1,0,2,2,2)) );
   WALBERLA_CHECK_NOT_NULLPTR( sliced );

   sliced->registerFlag("SecondFlag");

   WALBERLA_CHECK( ff.flagExists("SecondFlag")  );


   // Test deep copy
   FField * f1 = new FField ( 3,3,3,1 );
   auto flag1 = f1->registerFlag( "Flag1" );
   auto flag2 = f1->registerFlag( "Flag2" );
   f1->addFlag( 0,0,0,flag1 );
   f1->addFlag( 1,1,1,flag2 );

   FField * f2 = f1->clone();
   delete f1;

   WALBERLA_CHECK( f2->flagExists( "Flag1" ) );
   WALBERLA_CHECK( f2->flagExists( "Flag2" ) );
   WALBERLA_CHECK( f2->isFlagSet( 0,0,0,flag1 ) );
   WALBERLA_CHECK( f2->isFlagSet( 1,1,1,flag2 ) );

   delete f2;

}

void printingTest()
{
   using FField = FlagField<wlb::uint8_t>;
   FField ff ( 3,3,3,1 );
   auto ns = ff.registerFlag("NoSlip");
   auto fs = ff.registerFlag("FreeSlip");
   auto ob = ff.registerFlag("Obstacle");

   for(cell_idx_t i=0; i < cell_idx_c( ff.xSize() ); ++i)
      ff.addFlag(i,0,0,ns);

   for(cell_idx_t i=0; i < cell_idx_c( ff.ySize() ); ++i)
      ff.addFlag(0,i,0,fs);

   ff.addFlag(0,1,0,ob);
   ff.addFlag(1,0,0,ob);

   //printSlice( cout,ff, 2,0 );
}

void neighborhoodTest()
{
   using FField = FlagField<wlb::uint8_t>;
   FField ff ( 3,3,3,1 );
   auto i = ff.registerFlag("Interface");
   auto l = ff.registerFlag("Liquid");
   auto g = ff.registerFlag("Gas");

   for( auto it = ff.begin(); it != ff.end(); ++it )
      addFlag(it, l);

   ff.addFlag(1,0,0,i);
   ff.addFlag(1,2,2,g);
   for( auto it = ff.begin(); it != ff.end(); ++it )
   {
      if (it.x() == 1 && it.y() == 1 && it.z() == 1) {
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q19>(it, i) );
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q19>(it, g) );
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q19>(it, l) );
      }

      if (it.x() == 1 && it.y() == 2 && it.z() == 2) {
         WALBERLA_CHECK( !isFlagInNeighborhood<stencil::D3Q19>(it, i) );
         WALBERLA_CHECK( !isFlagInNeighborhood<stencil::D3Q19>(it, g) );
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q19>(it, l) );
      }

   }

   ff.removeFlag(1,0,0,i);
   ff.addFlag(0,0,0,i);
   ff.removeFlag(1,2,2,g);
   ff.addFlag(2,2,2,g);
   for( auto it = ff.begin(); it != ff.end(); ++it )
   {
      if (it.x() == 1 && it.y() == 1 && it.z() == 1) {
         WALBERLA_CHECK( field::isFlagInNeighborhood<stencil::D3Q27>(it, i) );
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q27>(it, g) );
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q27>(it, l) );
      }

      if (it.x() == 2 && it.y() == 2 && it.z() == 2) {
         WALBERLA_CHECK( !isFlagInNeighborhood<stencil::D3Q27>(it, i) );
         WALBERLA_CHECK( !isFlagInNeighborhood<stencil::D3Q27>(it, g) );
         WALBERLA_CHECK( isFlagInNeighborhood<stencil::D3Q27>(it, l) );
      }

   }
}

int main()
{
   debug::enterTestMode();

   registerTest();
   accessTest();
   iteratorTest();
   shallowCopyTest();
   printingTest();
   neighborhoodTest();
   return 0;
}
