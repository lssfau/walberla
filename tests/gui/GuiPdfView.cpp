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
//! \file GuiPdfView.cpp
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "gui/Gui.h"

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Vector3.h"

#include "field/AddToStorage.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"


namespace walberla {


typedef GhostLayerField<real_t,19> PdfField;
typedef GhostLayerField<real_t,1>  ScalarField;
typedef GhostLayerField<Vector3<real_t>,1 > VectorField;
using FField = FlagField<walberla::uint32_t>;


int main(int argc, char **argv)
{
   walberla::Environment env( argc, argv );

   const uint_t cells [] = { 1,1,1 };
   const uint_t blockCount [] = { 1, 1,1 };
   const uint_t nrOfTimeSteps = 20;

   // Create BlockForest
   auto blocks = blockforest::createUniformBlockGrid(blockCount[0],blockCount[1],blockCount[2],  //blocks
                                                     cells[0],cells[1],cells[2], //cells
                                                     1,                          //dx
                                                     false,                      //one block per process
                                                     true,true,true);            //periodicity


   // In addition to the normal GhostLayerField's  we allocated additionally a field containing the whole global simulation domain for each block
   // we can then check if the GhostLayer communication is correct, by comparing the small field to the corresponding part of the big field
   BlockDataID pdfField     = field::addToStorage<PdfField>( blocks, "Src" );

   // Init src field with some values
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) // block loop
   {
      // Init PDF field
      PdfField * src = blockIt->getData<PdfField>(pdfField);
      for( auto cellIt = src->beginWithGhostLayerXYZ(); cellIt != src->end(); ++cellIt ) // over all x,y,z,f
      {
         for( auto d = stencil::D3Q19::begin(); d != stencil::D3Q19::end(); ++d )
            cellIt.getF( d.toIdx() ) = real_c( d.toIdx() );
      }

   }

   // Create TimeLoop
   SweepTimeloop timeloop (blocks, nrOfTimeSteps );

   GUI gui (timeloop, blocks, argc, argv);
   gui.run();
   //timeloop.singleStep();
   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}