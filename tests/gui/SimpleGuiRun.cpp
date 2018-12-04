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
//! \file SimpleGuiRun.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "gui/Gui.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Vector3.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/gui/Connection.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"


namespace walberla {

typedef GhostLayerField<real_t,19> PdfField;
typedef GhostLayerField<real_t,1>  ScalarField;
typedef GhostLayerField<Vector3<real_t>,1 > VectorField;
using FField = FlagField<walberla::uint32_t>;

using LatticeModel = lbm::D3Q19<lbm::collision_model::SRT>;


int main(int argc, char **argv )
{
   walberla::Environment env( argc, argv );

   const uint_t cells [] = { 6,5,3 };
   const uint_t blockCount [] = { 4, 3,2 };
   const uint_t nrOfTimeSteps = 20;

   // Create BlockForest
   auto blocks = blockforest::createUniformBlockGrid(blockCount[0],blockCount[1],blockCount[2],  //blocks
                                        cells[0],cells[1],cells[2], //cells
                                        1,                          //dx
                                        false,                      //one block per process
                                        true,true,true);            //periodicity

   LatticeModel latticeModel( lbm::collision_model::SRT(1.5 ) );

   // In addition to the normal GhostLayerField's  we allocated additionally a field containing the whole global simulation domain for each block
   // we can then check if the GhostLayer communication is correct, by comparing the small field to the corresponding part of the big field

   BlockDataID pdfField     = lbm::addPdfFieldToStorage( blocks, "PdfField", latticeModel );

   BlockDataID scalarField1 = field::addToStorage<ScalarField>( blocks, "ScalarFieldOneGl", real_t(0), field::zyxf,  1 );
   BlockDataID scalarField2 = field::addToStorage<ScalarField>( blocks, "ScalarFieldTwoGl", real_t(0), field::zyxf,  2 );
   BlockDataID vectorField  = field::addToStorage<VectorField>( blocks, "VectorField", Vector3<real_t>(0,0,0) );
   BlockDataID flagField    = field::addFlagFieldToStorage<FField>( blocks, "FlagField" );


   // Init src field with some values
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) // block loop
   {
      // Init PDF field
      PdfField * src = blockIt->getData<PdfField>(pdfField);
      for( auto cellIt = src->begin(); cellIt != src->end(); ++cellIt ) // over all x,y,z,f
      {
         Cell cell ( cellIt.x(), cellIt.y(), cellIt.z() );
         blocks->transformBlockLocalToGlobalCell(cell, *blockIt);
         *cellIt = real_c( ( cell[0] + cell[1] + cell[2] + cellIt.f() ) % cell_idx_t(42) );
      }

      // Init scalarField1
      ScalarField * sf = blockIt->getData<ScalarField> ( scalarField1 );
      for( auto cellIt = sf->beginWithGhostLayer(); cellIt != sf->end(); ++cellIt ) // over all x,y,z
      {
         Cell cell ( cellIt.x(), cellIt.y(), cellIt.z() );
         blocks->transformBlockLocalToGlobalCell(cell, *blockIt);
         *cellIt = real_c( ( cell[0] + cell[1] + cell[2] ) % cell_idx_t(42) );
      }

      // Init scalarField2
      sf = blockIt->getData<ScalarField> ( scalarField2 );
      for( auto cellIt = sf->beginWithGhostLayer(); cellIt != sf->end(); ++cellIt ) // over all x,y,z
      {
         Cell cell ( cellIt.x(), cellIt.y(), cellIt.z() );
         blocks->transformBlockLocalToGlobalCell(cell, *blockIt);
         *cellIt = real_c( ( cell[0] + cell[1] + cell[2] ) % cell_idx_t(42) );
      }

      // Init vector field
      VectorField * vf = blockIt->getData<VectorField> ( vectorField );
      for ( auto cellIt = vf->beginWithGhostLayer(); cellIt != vf->end(); ++cellIt )
      {
         Cell cell ( cellIt.x(), cellIt.y(), cellIt.z() );
         blocks->transformBlockLocalToGlobalCell(cell, *blockIt);
         *cellIt = Vector3<real_t>( real_c(cell[0]), real_c(cell[1]), real_c(cell[2]) );
      }

      // Init Flag field
      FField * ff = blockIt->getData<FField> ( flagField );
      auto flag1 = ff->registerFlag( "AFlag 1" );
      auto flag2 = ff->registerFlag( "BFlag 2" );
      for ( auto cellIt = ff->beginWithGhostLayer(); cellIt != ff->end(); ++cellIt )
      {
         Cell cell ( cellIt.x(), cellIt.y(), cellIt.z() );
         blocks->transformBlockLocalToGlobalCell( cell, *blockIt );
         if ( ( cell[0] + cell[1] + cell[2] ) % 2 )
            addFlag( cellIt, flag1);
         else
            addFlag( cellIt, flag2);
      }

   }

   // Create TimeLoop
   SweepTimeloop timeloop (blocks, nrOfTimeSteps );

   GUI gui (timeloop, blocks, argc, argv);
   lbm::connectToGui<LatticeModel>( gui );
   gui.run();
   //timeloop.singleStep();
   return EXIT_SUCCESS;
}
}// namespace walberla


int main(int argc, char **argv){
  return walberla::main( argc, argv );
}




