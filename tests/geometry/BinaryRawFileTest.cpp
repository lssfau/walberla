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
//! \file BinaryRawFileTest.cpp
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/DataTypes.h"
#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"

#include "field/AddToStorage.h"
#include "field/vtk/VTKWriter.h"

#include "geometry/structured/BinaryRawFile.h"

#include "vtk/Initialization.h"

#include <string>
#include <vector>

namespace walberla {
namespace geometry {

void test( const std::string & filename, const Vector3<uint_t> & size, const std::string & datatype, const bool vtkOut )
{
   WALBERLA_LOG_INFO( "Testing unscaled" );
   BinaryRawFile brf( filename, size, datatype );

   auto blocks = blockforest::createUniformBlockGrid( uint_t(1), uint_t( 1 ), uint_t( 1 ), 
      size[0], size[1], size[2], 
      real_t( 1 ), 
      uint_t( 1 ), uint_t( 1 ), uint_t( 1 ) );

   typedef GhostLayerField< uint8_t, uint_t( 1 ) > ScalarField;

   BlockDataID scalarFieldID = field::addToStorage<ScalarField>( blocks, "BinaryRawFile" );

   for ( auto & block : *blocks)
   {
      auto field = block.getData<ScalarField>( scalarFieldID );
      
      CellInterval ci( 0, 0, 0, cell_idx_c( size[0] ) - 1, cell_idx_c( size[1] ) - 1, cell_idx_c( size[2] ) - 1 );

      for (const Cell c : ci)
      {
         field->get( c ) = brf.get( uint_c( c[0] ), uint_c( c[1] ), uint_c( c[2] ) );
      }
   }

   if (vtkOut)
   {
      WALBERLA_LOG_INFO( "Writing unscaled" );
      auto vtkOutput = vtk::createVTKOutput_BlockData( blocks, "BinaryRawFile" );
      vtkOutput->addCellDataWriter( make_shared< field::VTKWriter< ScalarField, uint8_t > >( scalarFieldID, "BinaryRawFile" ) );
      writeFiles( vtkOutput, true )();
   }
}

void testScaled( const std::string & filename, const Vector3<uint_t> & size, const std::string & datatype, const bool vtkOut )
{
   WALBERLA_LOG_INFO( "Testing scaled" );
   BinaryRawFile brf( filename, size, datatype );

   Vector3<uint_t> scaledSize( std::max( uint_t( 1 ), size[0] / uint_t( 2 ) ),
                               std::max( uint_t( 1 ), size[1] / uint_t( 3 ) ),
                               std::max( uint_t( 1 ), size[2] / uint_t( 5 ) ) );

   auto blocks = blockforest::createUniformBlockGrid( uint_t( 1 ), uint_t( 1 ), uint_t( 1 ),
      scaledSize[0], scaledSize[1], scaledSize[2],
      real_t( 1 ),
      uint_t( 1 ), uint_t( 1 ), uint_t( 1 ) );

   BinaryRawFileInterpolator brfi( blocks->getDomain(), brf, BinaryRawFileInterpolator::NEAREST_NEIGHBOR );

   typedef GhostLayerField< uint8_t, uint_t( 1 ) > ScalarField;

   BlockDataID scalarFieldID = field::addToStorage<ScalarField>( blocks, "BinaryRawFile" );

   for (auto & block : *blocks)
   {
      auto field = block.getData<ScalarField>( scalarFieldID );

      CellInterval ci( 0, 0, 0, cell_idx_c( scaledSize[0] ) - 1, cell_idx_c( scaledSize[1] ) - 1, cell_idx_c( scaledSize[2] ) - 1 );

      for (const Cell c : ci)
      {
         auto pos = blocks->getBlockLocalCellCenter( block, c );
         field->get( c ) = brfi.get( pos );
      }
   }

   if (vtkOut)
   {
      WALBERLA_LOG_INFO( "Writing scaled" );
      auto vtkOutput = vtk::createVTKOutput_BlockData( blocks, "BinaryRawFileScaled" );
      vtkOutput->addCellDataWriter( make_shared< field::VTKWriter< ScalarField, uint8_t > >( scalarFieldID, "BinaryRawFile" ) );
      writeFiles( vtkOutput, true )();
   }
}

} // namespace geometry
} // namespace walberla

int main( int argc, char * argv[] )
{
   walberla::mpi::Environment env( argc, argv );

   std::vector< std::string > args( argv, argv + argc );

   auto it = std::find( args.begin(), args.end(), std::string( "--vtk" ) );
   const bool vtk = it != args.end();
   if (it != args.end())
      args.erase( it );

   if (args.size() != 6u)
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args.front() << " [--vtk] FILENAME XSIZE YSIZE ZSIZE DATATYPE" );

   const std::string filename = args[1];

   walberla::Vector3< walberla::uint_t > size;
   try
   {
      size.set( std::stoull( args[2] ), 
                std::stoull( args[3] ),
                std::stoull( args[4] ) );
   }
   catch (const std::invalid_argument &)
   {
      WALBERLA_ABORT_NO_DEBUG_INFO( "USAGE: " << args.front() << " FILENAME XSIZE YSIZE ZSIZE DATATYPE" );
   }

   const std::string datatype = args[5];

   walberla::geometry::test( filename, size, datatype, vtk );
   walberla::geometry::testScaled( filename, size, datatype, vtk );
}

