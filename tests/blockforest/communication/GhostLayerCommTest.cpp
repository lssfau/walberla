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
//! \file GhostLayerCommTest.cpp
//! \ingroup comm
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Compares distributed calculation (with comm) to one big local computation
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/MakeBlockDataInitFunction.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include <iostream>


namespace walberla {
namespace wlb = walberla;

typedef GhostLayerField<real_t,19> PdfField;


template<typename Sten>
class StreamingSweep
{
   public:
      StreamingSweep(BlockDataID src, BlockDataID dst)
         : src_(src), dst_(dst)
      {}

      void operator() (IBlock * block) const
      {
         using namespace stencil;

         PdfField * src = block->getData<PdfField>(src_);
         PdfField * dst = block->getData<PdfField>(dst_);
         WALBERLA_ASSERT_NOT_NULLPTR( src )
         WALBERLA_ASSERT_NOT_NULLPTR( dst )

         for(PdfField::iterator i = dst->begin(); i != dst->end(); ++i )
         {
            stencil::Direction d = Sten::dir[i.f()];
            *i = src->get( i.x() - cx[d], i.y() - cy[d], i.z() - cz[d], i.f() );
         }
         src->swapDataPointers( dst );
      }

   protected:
      BlockDataID src_;
      BlockDataID dst_;
};

template<typename Sten>
class ManualPeriodicBoundarySweep
{
   public:
      ManualPeriodicBoundarySweep(BlockDataID field)
         : fieldId_(field), packInfo_(field)
      {}

      void operator() (IBlock * block)
      {
         /*
         PdfField * field = block->getData<PdfField>(fieldId_);
         for( cell_idx_t y = -1; y <= cell_idx_c( field->ySize() ); ++y )
         {
            for( cell_idx_t x = -1; x <= cell_idx_c( field->xSize() ); ++x )
               std::cout << field->get(x,y,0,0) << "\t";

            std::cout << std::endl;
         }*/

         for(typename Sten::iterator d = Sten::beginNoCenter(); d != Sten::end(); ++d) {
            packInfo_.communicateLocal(block,block,*d); //copy to itself
         }

         /*
         for( cell_idx_t y = -1; y <= cell_idx_c( field->ySize() ); ++y )
         {
            for( cell_idx_t x = -1; x <= cell_idx_c( field->xSize() ); ++x )
               std::cout << field->get(x,y,0,0) << "\t";

            std::cout << std::endl;
         }*/
      }


   protected:
      BlockDataID fieldId_;
      field::communication::PackInfo<PdfField> packInfo_;
};



class CompareSweep
{
   public:
      CompareSweep(BlockDataID smallField, BlockDataID bigField)
         : smallField_(smallField), bigField_(bigField)
      {}

      void operator() (IBlock * block)
      {
         PdfField * sf = block->getData<PdfField>(smallField_);
         PdfField * bf = block->getData<PdfField>(bigField_);

         const AABB & bb = block->getAABB();
         const cell_idx_t offset [3] = { wlb::cell_idx_c(bb.min(0)),
                                         wlb::cell_idx_c(bb.min(1)),
                                         wlb::cell_idx_c(bb.min(2)) };


#if 0
         using namespace stencil;
         using namespace std;

         WALBERLA_CRITICAL_SECTION_START

         cout << "Big Field " << MPIManager::instance()->rank() << endl;
         for( cell_idx_t y = -1; y <= cell_idx_c( bf->ySize() ); ++y )
         {
            for( cell_idx_t x = -1; x <= cell_idx_c( bf->xSize() ); ++x )
               cout << bf->get(x,y,0,D3Q19::idx[W]) << "\t";

            cout << endl;
         }

         cout << "Small Field " << MPIManager::instance()->rank() <<endl;
         for( cell_idx_t y = -1; y <= cell_idx_c( sf->ySize() ); ++y )
         {
            for( cell_idx_t x = -1; x <= cell_idx_c( sf->xSize() ); ++x )
               cout << sf->get(x,y,0,D3Q19::idx[W]) << "\t";

            cout << endl;
         }

         WALBERLA_CRITICAL_SECTION_END
#endif

         for(PdfField::iterator i = sf->begin(); i != sf->end(); ++i )
         {
            real_t globalValue = bf->get(i.x()+offset[0],i.y()+offset[1], i.z()+offset[2], i.f());
            real_t localValue  = *i;
            WALBERLA_CHECK_FLOAT_EQUAL(localValue,globalValue)
         }
      }


   protected:
      BlockDataID smallField_;
      BlockDataID bigField_;
};



int main(int argc, char **argv)
{
   debug::enterTestMode();

   //walberla::init(argc,argv);
   auto mpiManager = MPIManager::instance();
   mpiManager->initializeMPI(&argc,&argv);

   const uint_t cells [] = { 5,2,7 };
   const uint_t blockCount [] = { uint_c( mpiManager->numProcesses() ), 1, 1 };
   const uint_t nrOfTimeSteps = 30;
   bool periodic = true;
   const field::Layout layout = field::fzyx;

   // Create BlockForest
   uint_t glbCellsX = cells[0] * blockCount[0];
   uint_t glbCellsY = cells[1] * blockCount[1];
   uint_t glbCellsZ = cells[2] * blockCount[2];

   using blockforest::createUniformBlockGrid;

   auto blocks = createUniformBlockGrid(blockCount[0],blockCount[1],blockCount[2],  //blocks
                                        cells[0],cells[1],cells[2],                 //cells
                                        1,                                          //dx
                                        true,                                       //one block per process
                                        periodic,periodic,periodic);                //periodicity


   // In addition to the normal GhostLayerField's  we allocated additionally a field containing the whole global simulation domain for each block
   // we can then check if the GhostLayer communication is correct, by comparing the small field to the corresponding part of the big field
   BlockDataID srcField = field::addToStorage<PdfField>( blocks,"Src", 0.0, layout );
   BlockDataID dstField = field::addToStorage<PdfField>( blocks,"Dst", 0.0, layout );

   BlockDataID srcFieldGlb = blocks->addBlockData<PdfField>( makeBlockDataInitFunction<GhostLayerField<real_t,19> >(glbCellsX,glbCellsY,glbCellsZ, uint_t(1), layout),
                                                             "Global Src" );
   BlockDataID dstFieldGlb = blocks->addBlockData<PdfField>( makeBlockDataInitFunction<GhostLayerField<real_t,19> >(glbCellsX,glbCellsY,glbCellsZ, uint_t(1), layout),
                                                             "Global Dst" );

   // Init src field with some values
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt ) // block loop
   {
      // Init local field
      PdfField * src = blockIt->getData<PdfField>(srcField);
      PdfField * dst = blockIt->getData<PdfField>(dstField);

      for(auto cellIt = src->beginWithGhostLayer(); cellIt != src->end(); ++cellIt ) // over all x,y,z,f
      {
         Cell cell = cellIt.cell();
         blocks->transformBlockLocalToGlobalCell(cell, *blockIt);
         *cellIt = real_c( ( cell[0] + cell[1] + cell[2] + cellIt.f() ) % cell_idx_t(42) );
         dst->get(cellIt) = *cellIt;
      }

      // Init global field
      PdfField * glbSrc = blockIt->getData<PdfField>(srcFieldGlb);
      PdfField * glbDst = blockIt->getData<PdfField>(dstFieldGlb);

      for(auto cellIt = glbSrc->beginWithGhostLayer(); cellIt != src->end(); ++cellIt )  // over all x,y,z,f
      {
         const Cell & cell = cellIt.cell();
         *cellIt = real_c( ( cell[0] + cell[1] + cell[2] + cellIt.f() ) % cell_idx_t(42) );
         glbDst->get(cellIt) = *cellIt;
         WALBERLA_CHECK_FLOAT_EQUAL ( (*glbDst)(cell[0],cell[1],cell[2],cellIt.f()), (*glbSrc)(cell[0],cell[1],cell[2],cellIt.f()) )
      }

      WALBERLA_LOG_INFO("Initialization complete")
   }

   // Create TimeLoop
   SweepTimeloop timeLoop (blocks->getBlockStorage(), nrOfTimeSteps);

   timeLoop.add() << Sweep ( CompareSweep(srcField, srcFieldGlb) );


   // small local fields
   blockforest::communication::UniformBufferedScheme<stencil::D3Q19> scheme(blocks);
   scheme.addPackInfo( make_shared<field::communication::PackInfo< PdfField > >(srcField) );

   timeLoop.add() << BeforeFunction( scheme )
                  << Sweep ( StreamingSweep<stencil::D3Q19>( srcField, dstField ) );

   // big global fields
   if( periodic )
      timeLoop.add() << Sweep ( ManualPeriodicBoundarySweep<stencil::D3Q19>(srcFieldGlb) );


   timeLoop.add() << Sweep ( StreamingSweep<stencil::D3Q19>(srcFieldGlb,dstFieldGlb) );


   timeLoop.run();

   return EXIT_SUCCESS;
}
}

int main(int argc, char **argv) {
  return walberla::main(argc, argv);
}



