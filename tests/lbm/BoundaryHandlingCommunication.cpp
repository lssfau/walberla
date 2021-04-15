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
//! \file BoundaryHandlingCommunication.cpp
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/UBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"
#include "lbm/sweeps/CellwiseSweep.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"
#include "boundary/communication/HandlingPackInfo.h"

#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"

#include "timeloop/SweepTimeloop.h"


//#define TEST_USES_VTK_OUTPUT
#ifdef TEST_USES_VTK_OUTPUT
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"
#include "vtk/VTKOutput.h"
#endif

//#define DEVEL_OUTPUT
#ifdef DEVEL_OUTPUT
#include "core/math/Sample.h"
#endif


namespace walberla{

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID    UBB_Flag( "velocity bounce back" );
const FlagUID NoSlip_Flag( "no slip" );

const uint_t FieldGhostLayers  = uint_t(1);
const real_t GlobalOmega       = real_t(1.4);



template< typename LatticeModel_T >
class MyBoundaryHandling
{
public:

   typedef lbm::NoSlip< LatticeModel_T, flag_t >  NoSlip_T;
   typedef lbm::UBB< LatticeModel_T, flag_t >     UBB_T;

   typedef BoundaryHandling< FlagField_T, typename LatticeModel_T::Stencil, NoSlip_T, UBB_T > BoundaryHandling_T;

   MyBoundaryHandling( const std::string & id, const BlockDataID & flagField, const BlockDataID & pdfField, const real_t velocity ) :
      id_( id ), flagField_( flagField ), pdfField_( pdfField ), velocity_( velocity ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const std::string id_;

   const BlockDataID flagField_;
   const BlockDataID  pdfField_;

   const real_t velocity_;
};

template< typename LatticeModel_T >
typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T *
MyBoundaryHandling<LatticeModel_T>::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagField_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfField_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * boundaryHandling = new BoundaryHandling_T( std::string("boundary handling ")+id_, flagField, fluid,
                                                                   NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                                                   UBB_T( "velocity bounce back", UBB_Flag, pdfField ) );

   // closed no slip channel, periodic in x-direction

   CellInterval domainBB = storage->getDomainCellBB();
   storage->transformGlobalToBlockLocalCellInterval( domainBB, *block );

   domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.xMax() += cell_idx_c( FieldGhostLayers );

   // no slip SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   boundaryHandling->forceBoundary( NoSlip_Flag, south );

   // no slip NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   boundaryHandling->forceBoundary( NoSlip_Flag, north );

   // no slip BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   boundaryHandling->forceBoundary( NoSlip_Flag, bottom );

   // no slip TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   boundaryHandling->forceBoundary( NoSlip_Flag, top );

   // velocity bounce back ball :-)

   const Vector3<real_t> center( real_t(36), real_t(27), real_t(33) );
   const real_t sqrRadius( real_t(7*7) );

   for( auto cell = flagField->beginWithGhostLayer(); cell != flagField->end(); ++cell )
   {
      const cell_idx_t x = cell.x();
      const cell_idx_t y = cell.y();
      const cell_idx_t z = cell.z();

      Vector3<real_t> cellCenter;
      storage->getBlockLocalCellCenter( *block, Cell(x,y,z), cellCenter[0], cellCenter[1], cellCenter[2] );

      Vector3<real_t> distance = center - cellCenter;
      if( distance.sqrLength() <= sqrRadius )
         boundaryHandling->forceBoundary( UBB_Flag, x, y, z, typename UBB_T::Velocity( velocity_, real_t(0), real_t(0) ) );
   }

   boundaryHandling->fillWithDomain( domainBB );

   return boundaryHandling;
}



template< typename LatticeModel_T >
void removeFlagsInGhostLayer( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      auto boundaryHandling = block->template getData< typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T >( boundaryHandlingId );
      auto flagField = boundaryHandling->getFlagField();

      const CellInterval & innerCells = flagField->xyzSize();
      for( auto cell = flagField->beginWithGhostLayer(); cell != flagField->end(); ++cell )
      {
         if( !innerCells.contains( cell.x(), cell.y(), cell.z() ) )
         {
            boundaryHandling->clear( cell.x(), cell.y(), cell.z() );
         }
      }
   }
}



template< uint_t StencilSize >
void check( const shared_ptr< StructuredBlockForest > & blocks, const ConstBlockDataID & flagFieldId,
            const ConstBlockDataID & pdf1, const ConstBlockDataID & pdf2 )
{
   typedef GhostLayerField< real_t, StencilSize > Field_T;

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      const FlagField_T * flagField = block->template getData< FlagField_T >( flagFieldId );
      const flag_t fluid = flagField->getFlag( Fluid_Flag );

      const Field_T * field1 = block->template getData< Field_T >( pdf1 );
      const Field_T * field2 = block->template getData< Field_T >( pdf2 );

      const cell_idx_t xSize = cell_idx_c( blocks->getNumberOfXCells( *block ) );
      const cell_idx_t ySize = cell_idx_c( blocks->getNumberOfYCells( *block ) );
      const cell_idx_t zSize = cell_idx_c( blocks->getNumberOfZCells( *block ) );

      #ifdef DEVEL_OUTPUT
      math::Sample samples[StencilSize];
      #endif

      for( cell_idx_t z = cell_idx_c(0); z < zSize; ++z )
         for( cell_idx_t y = cell_idx_c(0); y < ySize; ++y )
            for( cell_idx_t x = cell_idx_c(0); x < xSize; ++x )
               if( flagField->isFlagSet(x,y,z,fluid) )
                  for( cell_idx_t f = 0; f < cell_idx_c( StencilSize ); ++f )
                  {
                     #ifdef DEVEL_OUTPUT
                     samples[f].insert( std::fabs( field1->get(x,y,z,f) - field2->get(x,y,z,f) ) );
                     #else
                     WALBERLA_CHECK_FLOAT_EQUAL( field1->get(x,y,z,f), field2->get(x,y,z,f) );
                     #endif
                  }

      #ifdef DEVEL_OUTPUT
      for( cell_idx_t f = 0; f < StencilSize; ++f )
      {
         WALBERLA_LOG_DEVEL( "Direction " << f << ": " << samples[f].format( "[%min, %max], %mean, %med" ) );
         WALBERLA_CHECK_FLOAT_EQUAL( samples[f].range(), 0.0 );
      }
      #endif
   }
}



int main( int argc, char ** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   if( MPIManager::instance()->numProcesses() != 8 )
      WALBERLA_ABORT( "The number of processes must be equal to 8!" );

   auto blocks = blockforest::createUniformBlockGrid( uint_t(6),  uint_t(6),  uint_t(6),
                                                      uint_t(12), uint_t(9), uint_t(11),
                                                      real_c(1.0),
                                                      uint_t(2),  uint_t(2),  uint_t(2),
                                                      true, false, false ); // periodicity

   const real_t velocity = real_t(0.0005);

#ifdef TEST_USES_VTK_OUTPUT
   SweepTimeloop timeloop( blocks->getBlockStorage(), uint_t(201) );
#else
   SweepTimeloop timeloop( blocks->getBlockStorage(), uint_t(3) );
#endif

   ///////////
   // D3Q19 //
   ///////////

   using D3Q19_TRT_COMP = lbm::D3Q19<lbm::collision_model::TRT>;
   D3Q19_TRT_COMP d3q19comp = D3Q19_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( GlobalOmega ) );

   BlockDataID pdfFieldId1  = lbm::addPdfFieldToStorage( blocks, "pdf field (D3Q19 with ghosts set)", d3q19comp );
   BlockDataID flagFieldId1 = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field (D3Q19 with ghosts set)" );

   BlockDataID boundaryHandlingId1 = blocks->addStructuredBlockData< MyBoundaryHandling< D3Q19_TRT_COMP >::BoundaryHandling_T >(
            MyBoundaryHandling< D3Q19_TRT_COMP >( " flag field (D3Q19 with ghosts set)", flagFieldId1, pdfFieldId1, velocity ),
                                                  "boundary handling flag field (D3Q19 with ghosts set)" );

   blockforest::communication::UniformBufferedScheme< D3Q19_TRT_COMP::CommunicationStencil > communication1( blocks );
   communication1.addPackInfo( make_shared< lbm::PdfFieldPackInfo< D3Q19_TRT_COMP > >( pdfFieldId1 ) );
   communication1.setLocalMode( blockforest::START );

   timeloop.add()
      << BeforeFunction( communication1, "communication (D3Q19 with ghosts set)" )
      << Sweep( MyBoundaryHandling<D3Q19_TRT_COMP>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId1 ), "boundary sweep (D3Q19 with ghosts set)" );

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_COMP, FlagField_T >( pdfFieldId1, flagFieldId1, Fluid_Flag ) ),
                            "LB stream & collide (D3Q19 with ghosts set)" );

   // ---

   BlockDataID pdfFieldId2  = lbm::addPdfFieldToStorage( blocks, "pdf field (D3Q19 with communication)", d3q19comp );
   BlockDataID flagFieldId2 = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field (D3Q19 with communication)" );

   BlockDataID boundaryHandlingId2 = blocks->addStructuredBlockData< MyBoundaryHandling< D3Q19_TRT_COMP >::BoundaryHandling_T >(
            MyBoundaryHandling< D3Q19_TRT_COMP >( " flag field (D3Q19 with communication)", flagFieldId2, pdfFieldId2, velocity ),
                                                  "boundary handling flag field (D3Q19 with communication)" );

   removeFlagsInGhostLayer<D3Q19_TRT_COMP>( blocks, boundaryHandlingId2 );

   blockforest::communication::UniformBufferedScheme< D3Q19_TRT_COMP::CommunicationStencil > communication2( blocks );
   communication2.addPackInfo( make_shared< lbm::PdfFieldPackInfo< D3Q19_TRT_COMP > >( pdfFieldId2 ) );
   communication2.addPackInfo( make_shared< boundary::HandlingPackInfo< MyBoundaryHandling< D3Q19_TRT_COMP >::BoundaryHandling_T > >( boundaryHandlingId2 ) );
   communication2.setLocalMode( blockforest::BUFFER );

   timeloop.add()
      << BeforeFunction( communication2, "communication (D3Q19 with communication)" )
      << Sweep( MyBoundaryHandling<D3Q19_TRT_COMP>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId2 ), "boundary sweep (D3Q19 with communication)" );

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_COMP, FlagField_T >( pdfFieldId2, flagFieldId2, Fluid_Flag ) ),
                            "LB stream & collide (D3Q19 with communication)" );

   ///////////
   // D3Q27 //
   ///////////

   using D3Q27_TRT_COMP = lbm::D3Q27<lbm::collision_model::TRT>;
   D3Q27_TRT_COMP d3q27comp = D3Q27_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( GlobalOmega ) );

   BlockDataID pdfFieldId3  = lbm::addPdfFieldToStorage( blocks, "pdf field (D3Q27 with ghosts set)", d3q27comp );
   BlockDataID flagFieldId3 = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field (D3Q27 with ghosts set)" );


   BlockDataID boundaryHandlingId3 = blocks->addStructuredBlockData< MyBoundaryHandling< D3Q27_TRT_COMP >::BoundaryHandling_T >(
            MyBoundaryHandling< D3Q27_TRT_COMP >( " flag field (D3Q27 with ghosts set)", flagFieldId3, pdfFieldId3, velocity ),
                                                  "boundary handling flag field (D3Q27 with ghosts set)" );

   blockforest::communication::UniformBufferedScheme< D3Q27_TRT_COMP::CommunicationStencil > communication3( blocks );
   communication3.addPackInfo( make_shared< lbm::PdfFieldPackInfo< D3Q27_TRT_COMP > >( pdfFieldId3 ) );
   communication3.setLocalMode( blockforest::WAIT );

   timeloop.add()
      << BeforeFunction( communication3, "communication (D3Q27 with ghosts set)" )
      << Sweep( MyBoundaryHandling<D3Q27_TRT_COMP>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId3 ), "boundary sweep (D3Q27 with ghosts set)" );

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_TRT_COMP, FlagField_T >( pdfFieldId3, flagFieldId3, Fluid_Flag ) ),
                            "LB stream & collide (D3Q27 with ghosts set)" );

   // ---

   BlockDataID pdfFieldId4  = lbm::addPdfFieldToStorage( blocks, "pdf field (D3Q27 with communication)", d3q27comp );
   BlockDataID flagFieldId4 = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field (D3Q27 with communication)" );

   BlockDataID boundaryHandlingId4 = blocks->addStructuredBlockData< MyBoundaryHandling< D3Q27_TRT_COMP >::BoundaryHandling_T >(
            MyBoundaryHandling< D3Q27_TRT_COMP >( " flag field (D3Q27 with communication)", flagFieldId4, pdfFieldId4, velocity ),
                                                  "boundary handling flag field (D3Q27 with communication)" );

   removeFlagsInGhostLayer<D3Q27_TRT_COMP>( blocks, boundaryHandlingId4 );

   blockforest::communication::UniformBufferedScheme< D3Q27_TRT_COMP::CommunicationStencil > communication4( blocks );
   communication4.addPackInfo( make_shared< lbm::PdfFieldPackInfo< D3Q27_TRT_COMP > >( pdfFieldId4 ) );
   communication4.addPackInfo( make_shared< boundary::HandlingPackInfo< MyBoundaryHandling< D3Q27_TRT_COMP >::BoundaryHandling_T > >( boundaryHandlingId4 ) );
   communication4.setLocalMode( blockforest::BUFFER );

   timeloop.add()
      << BeforeFunction( communication4, "communication (D3Q27 with communication)" )
      << Sweep( MyBoundaryHandling<D3Q27_TRT_COMP>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId4 ), "boundary sweep (D3Q27 with communication)" );

   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_TRT_COMP, FlagField_T >( pdfFieldId4, flagFieldId4, Fluid_Flag ) ),
                            "LB stream & collide (D3Q27 with communication)" );

#ifdef TEST_USES_VTK_OUTPUT

   field::createVTKOutput< FlagField_T >( flagFieldId1, *blocks, "flag_field_d3q19_no_comm", uint_t(1), uint_t(1) )();
   field::createVTKOutput< FlagField_T >( flagFieldId2, *blocks, "flag_field_d3q19_comm", uint_t(1), uint_t(1) )();

   auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, "fluid_field", uint_t(10) );

   field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldId1 );
   fluidFilter.addFlag( Fluid_Flag );
   pdfFieldVTKWriter->addCellInclusionFilter( fluidFilter );

   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_TRT_COMP, float > >( pdfFieldId1, "velocity (D3Q19 with ghosts set)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_TRT_COMP, float > >( pdfFieldId1, "density (D3Q19 with ghosts set)" ) );

   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_TRT_COMP, float > >( pdfFieldId2, "velocity (D3Q19 with communication)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_TRT_COMP, float > >( pdfFieldId2, "density (D3Q19 with communication)" ) );

   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q27_TRT_COMP, float > >( pdfFieldId3, "velocity (D3Q27 with ghosts set)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q27_TRT_COMP, float > >( pdfFieldId3, "density (D3Q27 with ghosts set)" ) );

   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q27_TRT_COMP, float > >( pdfFieldId4, "velocity (D3Q27 with communication)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q27_TRT_COMP, float > >( pdfFieldId4, "density (D3Q27 with communication)" ) );

   timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK" );

#endif

   WcTimingPool timeloopTiming;
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();
   
   WALBERLA_MPI_WORLD_BARRIER();

   check<19>( blocks, flagFieldId1, pdfFieldId1, pdfFieldId2 );
   check<27>( blocks, flagFieldId3, pdfFieldId3, pdfFieldId4 );

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
