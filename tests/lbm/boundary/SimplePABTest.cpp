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
//! \file SimplePABTest.cpp
//! \ingroup lbm
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimplePAB.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/IntegerFactorization.h"
#include "core/mpi/MPIManager.h"
#include "core/stringToNum.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <stdexcept>


namespace walberla {

   using flag_t = uint8_t;
   typedef lbm::D3Q19< lbm::collision_model::SRT, false > LatticeModel;
   using Stencil = LatticeModel::Stencil;
   using CommunicationStencil = LatticeModel::CommunicationStencil;
   using PDFField = lbm::PdfField<LatticeModel>;
   using MyFlagField = FlagField<flag_t>;


shared_ptr< StructuredBlockForest > makeStructuredBlockStorage( uint_t channelWidth, uint_t channelLength )
{
   const uint_t numProcesses = uint_c( MPIManager::instance()->numProcesses() );

   std::vector< real_t > weighting;
   weighting.push_back( real_c(1) / real_c(channelWidth)  );
   weighting.push_back( real_c(1) / real_c(channelWidth)  );
   weighting.push_back( real_c(1) / real_c(channelLength) );

   std::vector<uint_t> processes = math::getFactors( numProcesses, 3, weighting );

   uint_t blockSize[] = { uint_c( std::ceil( real_c(channelWidth)  / real_c(processes[0]) ) ),
                          uint_c( std::ceil( real_c(channelWidth)  / real_c(processes[1]) ) ),
                          uint_c( std::ceil( real_c(channelLength) / real_c(processes[2]) ) ) };

   WALBERLA_ASSERT_GREATER_EQUAL( blockSize[0] * processes[0], channelWidth  );
   WALBERLA_ASSERT_GREATER_EQUAL( blockSize[1] * processes[1], channelWidth  );
   WALBERLA_ASSERT_GREATER_EQUAL( blockSize[2] * processes[2], channelLength );

   return blockforest::createUniformBlockGrid(
      processes[0], processes[1], processes[2],  // blocks/processes in x/y/z direction
      blockSize[0], blockSize[1], blockSize[2], // cells per block in x/y/z direction
      real_t(1),                                 // cell size
      true,                                      // one block per process
      false,        false,        false,         // periodicity
      false );
}


const FlagUID& getFluidFlag()          { static FlagUID uid( "Fluid" );          return uid; }
const FlagUID& getNoSlipFlag()         { static FlagUID uid( "NoSlip" );         return uid; }
const FlagUID& getSimplePABLeftFlag()  { static FlagUID uid( "SimplePABLeft" );  return uid; }
const FlagUID& getSimplePABRightFlag() { static FlagUID uid( "SimplePABRight" ); return uid; }

class BoundaryHandlingCreator
{
public:
   typedef lbm::NoSlip<LatticeModel,flag_t>                              MyNoSlip;
   typedef lbm::SimplePAB<LatticeModel,MyFlagField>                      MySimplePAB;
   typedef BoundaryHandling< MyFlagField, Stencil, MyNoSlip, MySimplePAB, MySimplePAB > BoundaryHandlingT;

   BoundaryHandlingCreator( const BlockDataID& flagField, const BlockDataID& pdfField,
                            const real_t leftLatticeDensity, const real_t rightLatticeDensity,
                            uint_t channelLength, uint_t channelWidth, real_t omega )
      : flagField_( flagField ), pdfField_( pdfField ),
        channelLength_( channelLength ), channelWidth_( channelWidth ),
        omega_( omega ),
        leftLatticeDensity_( leftLatticeDensity ), rightLatticeDensity_( rightLatticeDensity )
   { }

   BoundaryHandlingT * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagField_;
   const BlockDataID pdfField_;

   uint_t channelLength_, channelWidth_;

   real_t omega_;

   const real_t leftLatticeDensity_, rightLatticeDensity_;
};

BoundaryHandlingCreator::BoundaryHandlingT * BoundaryHandlingCreator::operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   MyFlagField* flagField = block->getData<MyFlagField>( flagField_ );

   const auto fluidFlag = flagField->registerFlag( getFluidFlag() );

   /*
      FlagFieldT * const flagField, const real_t latticeDensity, const real_t omega, FlagUID domain, FlagUID noSlip )
   */

   BoundaryHandlingT * handling = new BoundaryHandlingT( "Boundary Handling", flagField, fluidFlag,
         MyNoSlip( "NoSlip", getNoSlipFlag(), block->getData< PDFField >( pdfField_ ) ),
         MySimplePAB( "SimplePABLeft",  getSimplePABLeftFlag(),  block->getData< PDFField >( pdfField_ ), block->getData< MyFlagField >( flagField_ ), leftLatticeDensity_,  omega_, getFluidFlag(), getNoSlipFlag() ),
         MySimplePAB( "SimplePABRight", getSimplePABRightFlag(), block->getData< PDFField >( pdfField_ ), block->getData< MyFlagField >( flagField_ ), rightLatticeDensity_,  omega_, getFluidFlag(), getNoSlipFlag() )
      );

   CellInterval channelBB(0, 0, 0,
      cell_idx_c(channelWidth_) - 1,  cell_idx_c(channelWidth_) - 1,  cell_idx_c(channelLength_) - 1 );

   // PAB WEST
   CellInterval west( channelBB.xMin() - 1, channelBB.yMin() - 1, channelBB.zMin() - 1,
      channelBB.xMin() - 1, channelBB.yMax() + 1, channelBB.zMax() + 1 );
   storage->transformGlobalToBlockLocalCellInterval( west, *block );
   handling->forceBoundary( getNoSlipFlag(), west );

   // PAB EAST
   CellInterval east( channelBB.xMax() + 1, channelBB.yMin() - 1, channelBB.zMin() - 1,
      channelBB.xMax() + 1, channelBB.yMax() + 1, channelBB.zMax() + 1 );
   storage->transformGlobalToBlockLocalCellInterval( east, *block );
   handling->forceBoundary( getNoSlipFlag(), east );

   // no slip SOUTH
   CellInterval south( channelBB.xMin() - 1, channelBB.yMin() - 1, channelBB.zMin() - 1,
      channelBB.xMax() + 1, channelBB.yMin() - 1, channelBB.zMax() + 1 );
   storage->transformGlobalToBlockLocalCellInterval( south, *block );
   handling->forceBoundary( getNoSlipFlag(), south );

   // no slip NORTH
   CellInterval north( channelBB.xMin() - 1, channelBB.yMax() + 1, channelBB.zMin() - 1,
      channelBB.xMax() + 1, channelBB.yMax() + 1, channelBB.zMax() + 1 );
   storage->transformGlobalToBlockLocalCellInterval( north, *block );
   handling->forceBoundary( getNoSlipFlag(), north );

   // no slip BOTTOM
   CellInterval bottom( channelBB.xMin() - 1, channelBB.yMin() - 1, channelBB.zMin() - 1,
      channelBB.xMax() + 1, channelBB.yMax() + 1, channelBB.zMin() - 1 );
   storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
   handling->forceBoundary( getSimplePABLeftFlag(), bottom );

   // no slip TOP
   CellInterval top( channelBB.xMin() - 1, channelBB.yMin() - 1, channelBB.zMax() + 1,
      channelBB.xMax() + 1, channelBB.yMax() + 1, channelBB.zMax() + 1 );
   storage->transformGlobalToBlockLocalCellInterval( top, *block );
   handling->forceBoundary( getSimplePABRightFlag(), top );

   CellInterval localChannelBB;
   storage->transformGlobalToBlockLocalCellInterval( localChannelBB, *block, channelBB );

   handling->fillWithDomain( localChannelBB );

   WALBERLA_ASSERT( handling->checkConsistency() );

   return handling;
}


int main( int argc, char **argv )
{
   MPIManager::instance()->initializeMPI( &argc, &argv );

   uint_t channelLength, channelWidth, numTimesteps;
   real_t deltaDensity, omega;

   try {
      std::vector<std::string> args( argv, argv + argc );

      if( args.size() != 6 )
         throw std::invalid_argument( "Wrong number of command line arguments!" );

      channelLength = stringToNum<uint_t>( args[1] );
      channelWidth  = stringToNum<uint_t>( args[2] );
      omega         = stringToNum<real_t>( args[3] );
      deltaDensity  = stringToNum<real_t>( args[4] );
      numTimesteps  = stringToNum<uint_t>( args[5] );
   }
   catch( std::exception & )
   {
      WALBERLA_ABORT( "Wrong Parameters!\n\nUsage:\nSimplePABTest CHANNEL_LENGTH CHANNEL_WIDTH OMEGA DELTA_DENSITY NUM_TIMESTEPS" );
   }

   auto blockStorage = makeStructuredBlockStorage( channelWidth, channelLength );

   LatticeModel latticeModel = LatticeModel( lbm::collision_model::SRT( omega ) );
   BlockDataID pdfField = lbm::addPdfFieldToStorage( blockStorage, "PDF field", latticeModel );

   BlockDataID flagField = field::addFlagFieldToStorage<MyFlagField>( blockStorage, "Flag field" );

   BlockDataID boundaryHandling = blockStorage->addStructuredBlockData( "Boundary Handling" )
      << StructuredBlockDataCreator<BoundaryHandlingCreator::BoundaryHandlingT>(
      BoundaryHandlingCreator(flagField, pdfField, real_t(1), real_t(1) + deltaDensity, channelLength, channelWidth, omega ) );

   SweepTimeloop timeloop( blockStorage->getBlockStorage(), numTimesteps );

   blockforest::communication::UniformBufferedScheme< CommunicationStencil > scheme( blockStorage );
   scheme.addPackInfo( make_shared< field::communication::PackInfo< PDFField > >( pdfField ) );
   timeloop.addFuncBeforeTimeStep( scheme, "Communication" );

   timeloop.add() << Sweep( BoundaryHandlingCreator::BoundaryHandlingT::getBlockSweep(boundaryHandling) );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel, MyFlagField >( pdfField, flagField, getFluidFlag() ) ), "LBM_SRT" );

   timeloop.run();

   auto vtkOut = vtk::createVTKOutput_BlockData( *blockStorage );
   auto densityWriter  = make_shared<lbm::DensityVTKWriter<LatticeModel>>( pdfField, "Density" );
   auto velocityWriter = make_shared<lbm::VelocityVTKWriter<LatticeModel>>( pdfField, "Velocity" );
   //auto flagsWriter    = make_shared<field::DumpField<flag_t, 1>>( flagField, "Flags" );
   vtkOut->addCellDataWriter( densityWriter );
   vtkOut->addCellDataWriter( velocityWriter );
   //vtkOut->addCellDataWriter( flagsWriter );
   vtkOut->write();

   return 0;
}

} // namespace walberla



int main(int argc, char **argv) {

   walberla::debug::enterTestMode();

   return walberla::main(argc,argv);
}
