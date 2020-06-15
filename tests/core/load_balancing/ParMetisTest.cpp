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
//! \file ParMetisTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Test for MetisWrapper
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/load_balancing/ParMetisWrapper.h"
#include "core/math/Vector2.h"
#include "core/math/IntegerFactorization.h"
#include "core/logging/Logging.h"
#include "core/stringToNum.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"
#include "field/vtk/VTKWriter.h"
#include "field/communication/PackInfo.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"

#include "vtk/VTKOutput.h"
#include "vtk/BlockCellDataWriter.h"

int main( int argc, char * argv[] )
{
   using namespace walberla;

   mpi::Environment env( argc, argv );

   debug::enterTestMode();

   std::vector< std::string > args( argv, argv + argc );
   Vector2<uint_t> fieldSize;
   uint_t          partitions;
   bool            vtk = true;

   try {
      fieldSize.set( stringToNum< uint_t >( args.at(1) ), stringToNum< uint_t >( args.at(2) ) );
      partitions = stringToNum< uint_t >( args.at(3) );

      auto it = std::find( args.begin(), args.end(), "--no-vtk" );
      if(it != args.end())
      {
         vtk = false;
         args.erase( it );
      }
   }
   catch( std::exception & e )
   {
      WALBERLA_ABORT( "USAGE: " << args[0] << " SizeX SizeY #Partitions\n\n" << e.what() );
   }

   int numProcesses = MPIManager::instance()->numProcesses();
   auto gridSize = math::getFactors( uint_c( numProcesses ), uint_t(2) );

   auto blocks = blockforest::createUniformBlockGrid( gridSize[0],  gridSize[1], uint_t(1),
                                                      fieldSize[0], fieldSize[1], uint_t(1),
                                                      real_t(1),
                                                      true,
                                                      true, true, false );

   typedef field::GhostLayerField< int64_t, 1 > FieldType;

   auto domainId    = field::addToStorage< FieldType >( blocks, "domain", int64_t(-1), field::zyxf, uint_t(1) );
   auto partFieldId = field::addToStorage< FieldType >( blocks, "partitions", int64_t(-1), field::zyxf, uint_t(1) );

   auto & domain    = *( blocks->begin()->getData< FieldType >( domainId    ) );
   auto & partField = *( blocks->begin()->getData< FieldType >( partFieldId ) );

   WALBERLA_CHECK_EQUAL( std::distance( blocks->begin(), blocks->end() ), 1 );

   int64_t ctr = int64_c( fieldSize[0] * fieldSize[1] ) * int64_c( MPIManager::instance()->rank() );
   Cell c( cell_idx_t(0), cell_idx_t(0), cell_idx_t(0) );
   for( c[0] = 0; c[0] < cell_idx_c( domain.xSize() ); ++c[0] )
      for( c[1] = 0; c[1] < cell_idx_c( domain.ySize() ); ++c[1] )
      {
         domain.get( c ) = ctr++;
      }

   blockforest::communication::UniformBufferedScheme< stencil::D2Q9 > scheme( blocks );
   scheme.addPackInfo( make_shared< field::communication::PackInfo< FieldType > >( domainId ) );
   scheme();

   std::vector< int64_t > vtxdist;
   for( int i = 0; i < MPIManager::instance()->numProcesses() + 1; ++i )
   {
      vtxdist.push_back( int64_c(i) * int64_c( fieldSize[0] * fieldSize[1] ) );
   }

   int64_t ncon  = int64_t(1);
   std::vector< int64_t > xadj, adjncy;
   xadj.push_back( int64_t(0) );

   c = Cell( cell_idx_t(0), cell_idx_t(0), cell_idx_t(0) );
   for( c[0] = 0; c[0] < cell_idx_c( domain.xSize() ); ++c[0] )
      for( c[1] = 0; c[1] < cell_idx_c( domain.ySize() ); ++c[1] )
      {
         for( auto dirIt = stencil::D2Q9::beginNoCenter(); dirIt != stencil::D2Q9::end(); ++dirIt )
            adjncy.push_back( domain.get( c + *dirIt ) );
       
         xadj.push_back( int64_c( adjncy.size() ) );         
      }

   WALBERLA_CHECK_EQUAL( xadj.size(), fieldSize[0] * fieldSize[1] + uint_t(1) );
   WALBERLA_CHECK_EQUAL( adjncy.size(), fieldSize[0] * fieldSize[1] * uint_t(8) );

   int64_t wgtflag = 0;
   int64_t numflag = 0;
   int64_t nparts = int64_c( partitions );
   std::vector< double > tpwgts( partitions, 1.0 / numeric_cast<double>( partitions ) );
   std::vector< double > ubvec( numeric_cast<size_t>(ncon), 1.05 );
   int64_t options[] = {0,0,0};
   int64_t edgecut;
   std::vector< int64_t > part( fieldSize[0] * fieldSize[1] );
   MPI_Comm comm = MPIManager::instance()->comm();

   int result = core::ParMETIS_V3_PartKway( &(vtxdist.front()), &(xadj.front()), &(adjncy.front()), nullptr, nullptr, &wgtflag, &numflag, &ncon, &nparts,
                                            &(tpwgts.front()), &(ubvec.front()), options, &edgecut, &(part.front()), &comm );


   WALBERLA_CHECK_EQUAL( result, core::METIS_OK );  

   c = Cell( cell_idx_t(0), cell_idx_t(0), cell_idx_t(0) );
   auto it = part.begin();
   for( c[0] = 0; c[0] < cell_idx_c( domain.xSize() ); ++c[0] )
      for( c[1] = 0; c[1] < cell_idx_c( domain.ySize() ); ++c[1] )
      {
         if( domain.get( c ) == int64_t( -1 ) )
            continue;

         WALBERLA_CHECK_UNEQUAL( it, part.end() );

         partField.get( c ) = *it++;

      }
   WALBERLA_CHECK_EQUAL( it, part.end() );

   if( vtk )
   {
      auto vtkOutput = vtk::createVTKOutput_BlockData( blocks, "Metis", 1, 0 );
      vtkOutput->addCellDataWriter( make_shared< field::VTKWriter< FieldType, int64_t > >( domainId, "domain" ) );
      vtkOutput->addCellDataWriter( make_shared< field::VTKWriter< FieldType, int64_t > >( partFieldId, "partitions" ) );
      vtk::writeFiles( vtkOutput )();
   }

   return EXIT_SUCCESS;
}
