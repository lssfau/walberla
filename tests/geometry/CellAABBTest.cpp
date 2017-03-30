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
//! \file CellAABBTest.cpp
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "geometry/structured/BasicVoxelFileReader.h"
#include "core/debug/TestSubsystem.h"


int main()
{
   walberla::debug::enterTestMode();

   using walberla::geometry::CellAABB;

   {
      CellAABB ca;

      WALBERLA_CHECK_EQUAL( ca.xBegin, 0 );
      WALBERLA_CHECK_EQUAL( ca.yBegin, 0 );
      WALBERLA_CHECK_EQUAL( ca.zBegin, 0 );
      WALBERLA_CHECK_EQUAL( ca.xEnd, 0 );
      WALBERLA_CHECK_EQUAL( ca.yEnd, 0 );
      WALBERLA_CHECK_EQUAL( ca.zEnd, 0 );

      WALBERLA_CHECK_EQUAL( ca.xSize(), 1 );
      WALBERLA_CHECK_EQUAL( ca.ySize(), 1 );
      WALBERLA_CHECK_EQUAL( ca.zSize(), 1 );

      WALBERLA_CHECK_EQUAL( ca.numCells(), 1 );
   }
   {
      CellAABB ca(1, 2, 3, 4, 5, 6);

      WALBERLA_CHECK_EQUAL( ca.xBegin, 1 );
      WALBERLA_CHECK_EQUAL( ca.yBegin, 2 );
      WALBERLA_CHECK_EQUAL( ca.zBegin, 3 );
      WALBERLA_CHECK_EQUAL( ca.xEnd, 4 );
      WALBERLA_CHECK_EQUAL( ca.yEnd, 5 );
      WALBERLA_CHECK_EQUAL( ca.zEnd, 6 );

      WALBERLA_CHECK_EQUAL( ca.xSize(), 4 );
      WALBERLA_CHECK_EQUAL( ca.ySize(), 4 );
      WALBERLA_CHECK_EQUAL( ca.zSize(), 4 );

      WALBERLA_CHECK_EQUAL( ca.numCells(), 4*4*4 );
   }
   {
      CellAABB ca(0, 0, 0, 0, 0, 0);
      WALBERLA_CHECK_EQUAL( ca.numCells(), 1 );
   }
   {
      CellAABB ca(0, 0, 0, 1, 1, 1);
      WALBERLA_CHECK_EQUAL( ca.numCells(), 8 );
   }
   {
      CellAABB ca(0, 0, 0, 9, 9, 9);
      WALBERLA_CHECK_EQUAL( ca.numCells(), 1000 );
   }
   {
      CellAABB ca(1, 1, 1, 10, 10, 10);
      WALBERLA_CHECK_EQUAL( ca.numCells(), 1000 );
   }
}
