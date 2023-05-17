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
//! \file CudaMPI.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "blockforest/Initialization.h"


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/mpi/Datatype.h"

#include "gpu/GPUField.h"

#include "field/communication/MPIDatatypes.h"
#include "field/AddToStorage.h"
#include "timeloop/SweepTimeloop.h"

#include "gui/Gui.h"


using namespace walberla;


void fullFieldTransfer()
{
   Field<double,4>  h_f1 ( 3, 4, 2, 42.0, field::fzyx );
   Field<double,4>  h_f2 ( 3, 4, 2, 27.0, field::fzyx );

   gpu::GPUField<double> d_f ( 3, 4, 2, 4, 0, field::fzyx );


   // Transfer h_f1 from CPU to GPU d_f

   auto h_f1_datatype = mpi::Datatype ( field::communication::mpiDatatype( h_f1 ) );
   auto h_f2_datatype = mpi::Datatype ( field::communication::mpiDatatype( h_f2 ) );
   auto d_f_datatype  = mpi::Datatype ( field::communication::mpiDatatype( d_f )  );

   WALBERLA_LOG_DEVEL("ISend");
   MPI_Request request1;
   MPI_Isend( h_f1.data(), 1, h_f1_datatype, 0, 0, MPI_COMM_WORLD, &request1 );

   WALBERLA_LOG_DEVEL("IRecv");
   MPI_Request request2;
   MPI_Irecv( d_f.data(), 1, d_f_datatype, 0, 0, MPI_COMM_WORLD, &request2 );

   MPI_Wait( &request1, MPI_STATUS_IGNORE );
   MPI_Wait( &request2, MPI_STATUS_IGNORE );

   // Transfer GPU field d_f back to CPU into h_f2

   MPI_Request request3;
   WALBERLA_LOG_DEVEL("ISend");
   MPI_Isend( d_f.data(), 1, d_f_datatype, 0, 0, MPI_COMM_WORLD , &request3 );

   MPI_Request request4;
   WALBERLA_LOG_DEVEL("IRecv");
   MPI_Irecv( h_f2.data(), 1, h_f2_datatype, 0, 0, MPI_COMM_WORLD, &request4 );

   MPI_Wait( &request3, MPI_STATUS_IGNORE );
   MPI_Wait( &request4, MPI_STATUS_IGNORE );

   WALBERLA_CHECK_EQUAL( h_f1, h_f2 );
}


void blockStorageAndGui( int argc, char ** argv )
{
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid(
            uint_c(1) , uint_c(1), uint_c(1), // number of blocks in x,y,z direction
            uint_c(5) , uint_c(7), uint_c(3), // number of blocks in x,y,z direction
            real_c(1),                        // dx: length of one cell in physical coordinates
            false,                            // one block per process? - "false" means all blocks to one process
            true, true, true );               // no periodicity

   typedef GhostLayerField<real_t,1> ScalarField;
   BlockDataID cpuFieldID1 = field::addToStorage<ScalarField>( blocks, "CPUField 1", real_c(42), field::fzyx, uint_c(1) );
   BlockDataID cpuFieldID2 = field::addToStorage<ScalarField>( blocks, "CPUField 2", real_c(0),  field::fzyx, uint_c(1) );

   typedef gpu::GPUField<real_t> GPUField;
   BlockDataID gpuFieldID = blocks->addStructuredBlockData< GPUField >(
            [&] ( IBlock * block, StructuredBlockStorage * const s ) {
               return new GPUField( s->getNumberOfXCells(*block),
                                    s->getNumberOfYCells(*block),
                                    s->getNumberOfZCells(*block),
                                    1 , 1);
             },
             "GPU Field" );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      // get the field stored on the current block
      ScalarField * h_f1 = blockIt->getData<ScalarField>( cpuFieldID1 );
      ScalarField * h_f2 = blockIt->getData<ScalarField>( cpuFieldID2 );
      GPUField    * d_f  = blockIt->getData<GPUField>   ( gpuFieldID );

      auto h_f1_datatype = mpi::Datatype ( field::communication::mpiDatatypeSliceBeforeGhostlayer( *h_f1, stencil::W, 1, true ) );
      auto h_f2_datatype = mpi::Datatype ( field::communication::mpiDatatypeSliceBeforeGhostlayer( *h_f2, stencil::W, 1, true ) );
      auto d_f_datatype  = mpi::Datatype ( field::communication::mpiDatatypeSliceBeforeGhostlayer( *d_f , stencil::W, 1, true ) );

      MPI_Sendrecv( h_f1->data() , 1, h_f1_datatype, 0, 0,
                    d_f->data(),   1, d_f_datatype , 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

      MPI_Sendrecv( d_f->data(), 1, d_f_datatype, 0, 0,
                    h_f2->data(),  1, h_f2_datatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );


   }

   SweepTimeloop timeloop( blocks, 4 );
   GUI gui( timeloop, blocks, argc, argv );
   gui.run();

}


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   fullFieldTransfer();
   //blockStorageAndGui(argc, argv);


   return 0;
}
