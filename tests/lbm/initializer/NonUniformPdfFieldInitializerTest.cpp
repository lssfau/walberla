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
//! \file NonUniformPdfFieldInitializerTest.cpp
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "lbm/all.h"

#include <vector>
#include <string>
#include <sstream>


namespace walberla {


using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using PdfFieldInitializer_T = lbm::initializer::PdfFieldInitializer<LatticeModel_T>;

static void refinementSelectionFunction( SetupBlockForest& forest ) {

   SetupBlock* block = forest.getRootBlock( 0 );

   if( !block->hasChildren() )
      block->setMarker( true );
}

static void workloadMemorySUIDAssignmentFunction( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      blocks[i]->setMemory( 1.0 );
      blocks[i]->setWorkload( 1.0 );
   }
}


struct LinearDensityCellProfile
{
   real_t operator()( const Cell & globalCell ) const
   {
      return real_t(1.0) + real_t(0.5) * real_c( globalCell.x() );
   }
};

struct LinearDensityCenterProfile
{
   real_t operator()( const Vector3<real_t> & center ) const
   {
      return real_t(1.0) + real_t(0.5) * center[0];
   }
};


template<typename InitFunc>
void testCellInit( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks )
{
   lbm::initializer::PdfFieldInitializer< LatticeModel_T > initializer( pdfFieldId, blocks );

   InitFunc linearDensityProfile;
   initializer.initDensity( linearDensityProfile );

   /// check for cells in ghost layers between coarse and fine blocks

   // get neighbouring coarse and fine block
   auto coarseBlock = blocks->getBlock(real_t(15), real_t(5), real_t(5));
   auto   fineBlock = blocks->getBlock(real_t(7.5), real_t(2.5), real_t(2.5));

   // retrieve fields
   auto coarsePdfField = coarseBlock->getData<PdfField_T>(pdfFieldId);
   auto finePdfField = fineBlock->getData<PdfField_T>(pdfFieldId);
   WALBERLA_ASSERT_NOT_NULLPTR(coarseBlock)
   WALBERLA_ASSERT_NOT_NULLPTR(fineBlock)

   // ghost layers in coarse block
   CellInterval coarseGL;
   coarsePdfField->getGhostRegion(stencil::directionFromAxis(0, true), coarseGL, 1);
   CellInterval fineInner;
   finePdfField->getSliceBeforeGhostLayer(stencil::directionFromAxis(0, false), fineInner, 2);

   Cell fineCellXMin {fineInner.xMin(), fineInner.yMin(), fineInner.zMin()};
   auto fineValXMin = finePdfField->getDensity(fineCellXMin);

   Cell fineCellXMax {fineInner.xMax(), fineInner.yMin(), fineInner.zMin()};
   auto fineValXMax = finePdfField->getDensity(fineCellXMax);

   Cell coarseCell {coarseGL.xMin(), coarseGL.yMin(), coarseGL.zMin()};
   auto coarseVal = coarsePdfField->getDensity(coarseCell);

   WALBERLA_CHECK_FLOAT_EQUAL(coarseVal, real_t(0.5) * (fineValXMin + fineValXMax))
}

int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment env( argc, argv );

   if( !walberla::MPIManager::instance()->rankValid() )
      walberla::MPIManager::instance()->useWorldComm();

   SetupBlockForest sforest;

   sforest.addRefinementSelectionFunction( refinementSelectionFunction );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadMemorySUIDAssignmentFunction );

   AABB domain( 0, 0, 0, 20, 10, 10 );
   sforest.init( domain, 2, 1, 1, true, false, false );
   sforest.assignAllBlocksToRootProcess();

   sforest.writeVTKOutput("initialDomainDecomposition");

   // setup StructuredBlockForest

   auto blocks = std::make_shared<StructuredBlockForest>( std::make_shared< BlockForest >( 0, sforest, true ), 10, 10, 10 );
   blocks->createCellBoundingBoxes();


   LatticeModel_T lm( real_c(1) );

   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf_field", lm );

   testCellInit<LinearDensityCellProfile>(pdfFieldId, blocks);
   testCellInit<LinearDensityCenterProfile>(pdfFieldId, blocks);

   return 0;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}