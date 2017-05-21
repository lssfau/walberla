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
//! \file
//! \author Paulo Carvalho <prcjunior@inf.ufpr.br>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"

#include "field/GhostLayerField.h"


#ifdef DBG_PRINT_ON
#define DBG_PRINT_FIELD( f )     printField( f )
#define DBG_PRINT( fmt, ... )    printf( fmt, ##__VA_ARGS__ )
#else
#define DBG_PRINT_FIELD( f )
#define DBG_PRINT(fmt, ...)
#endif

template<typename Field_T>
void printField( Field_T& field )
{
   using namespace walberla;
   cell_idx_t fs = 0;
   cell_idx_t zs = -(cell_idx_t)field.nrOfGhostLayers();
   cell_idx_t ys = -(cell_idx_t)field.nrOfGhostLayers();
   cell_idx_t xs = -(cell_idx_t)field.nrOfGhostLayers();
   cell_idx_t nf = (cell_idx_t)field.fSize();
   cell_idx_t nz = (cell_idx_t)(field.zSize() + field.nrOfGhostLayers());
   cell_idx_t ny = (cell_idx_t)(field.ySize() + field.nrOfGhostLayers());
   cell_idx_t nx = (cell_idx_t)(field.xSize() + field.nrOfGhostLayers());
   for ( cell_idx_t f = fs; f < nf; ++f ) {
      std::cout << "{";
      for ( cell_idx_t z = zs; z < nz; ++z ) {
         std::cout << ( z == zs ? "[" : " [" );
         for ( cell_idx_t y = ys; y < ny; ++y ) {
            std::cout << "(";
            for ( cell_idx_t x = xs; x < nx; ++x ) {
               std::cout << field( x, y, z, f ) << ( x == nx-1 ? "" : " " );
            }
            std::cout << ( y == ny-1 ? ")" : ") " );
         }
         std::cout << "]\n";
      }
      std::cout << "}\n";
   }
}


#define CHECK_FIELD_EQUAL( f1, f2 ) WALBERLA_CHECK( checkFieldEqual( f1, f2 ), "Field differ" )

template< typename Field_T >
bool checkFieldEqual( Field_T& field1, Field_T& field2 )
{
   using namespace walberla;
   WALBERLA_ASSERT( field1.xSize() == field2.xSize() &&
                    field1.ySize() == field2.ySize() &&
                    field1.zSize() == field2.zSize() &&
                    field1.fSize() == field2.fSize() );

   WALBERLA_FOR_ALL_CELLS_XYZ( &field2,
      for ( uint_t f = 0; f < field1.fSize(); ++f )
      {
         if ( field1.get( x, y, z, f ) != field2.get( x, y, z, f ) )
         {
            return false;
         }
      }
   )
   return true;
}
