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
//! \file CurveGatherTest.cpp
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Tests the MPIGatherScheme together with CurveGatherPackInfo
//
//======================================================================================================================

#include "gather/CurveGatherPackInfo.h"
#include "gather/GnuPlotGraphWriter.h"
#include "gather/MPIGatherScheme.h"

#include "blockforest/Initialization.h"

#include "core/Environment.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/interpolators/TrilinearInterpolator.h"

#include "gui/Gui.h"

#include "timeloop/SweepTimeloop.h"


namespace walberla {

typedef GhostLayerField<real_t,1> GlField;



GlField * createField( IBlock* const block, StructuredBlockStorage* const storage )
{
   auto f = new  GlField ( storage->getNumberOfXCells( *block ),
                           storage->getNumberOfYCells( *block ),
                           storage->getNumberOfZCells( *block ),
                           1);

   for ( auto i = f->beginWithGhostLayer(); i != f->end(); ++i )
   {
      Cell globalCell;

      storage->transformBlockLocalToGlobalCell( globalCell, *block, i.cell() );
      *i = real_c( globalCell.z() ) * storage->dz();
   }
   return f;
}


int main( int argc, char ** argv )
{
   walberla::Environment walberlaEnv( argc, argv );

   const real_t dx = real_c(0.1);
   const walberla::uint_t cells = uint_c( real_t(1) / dx );
   const walberla::uint_t nrBlocks[] = { 2,2,1 };

   using blockforest::createUniformBlockGrid;
   shared_ptr<StructuredBlockForest>
   blocks = createUniformBlockGrid( nrBlocks[0], nrBlocks[1], nrBlocks[2],
                                    cells / nrBlocks[0], cells / nrBlocks[1], cells / nrBlocks[2],
                                    dx,
                                    true,                // one block per process
                                    false, false, false, // periodicity
                                    false );             // do NOT keep global information

   BlockDataID fieldId = blocks->addStructuredBlockData<GlField >( &createField, "MyField" );


   using namespace gather;

   auto gnuplotWriter = make_shared<GnuPlotGraphWriter>( "gnuplot" );

   const char * curveX = "0.5 * cos(t) + 0.5";
   const char * curveY = "0.2";
   const char * curveZ = "0.499 * sin(t) + 0.5";
   auto curvePackInfo = make_shared<CurveGatherPackInfo<GlField> >( blocks, fieldId,
                                                         curveX, curveY, curveZ,
                                                         real_t(0), real_c(2 * 3.141),
                                                         uint_t(100), gnuplotWriter );

   MPIGatherScheme gatherScheme( blocks->getBlockStorage() );
   gatherScheme.addPackInfo( curvePackInfo );
   gatherScheme();

   //SweepTimeloop timeloop( blocks, 10 );
   //GUI gui ( timeloop, blocks, argc, argv );
   //gui.run();

   return 0;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}