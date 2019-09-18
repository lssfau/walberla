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
//! \file UniformGrid.cpp
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/communication/UniformDirectScheme.h"
#include "blockforest/loadbalancing/Cartesian.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/config/Config.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/IntegerFactorization.h"
#include "core/math/Vector3.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/Reduce.h"
#include "core/timing/TimingPool.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/FlagUID.h"
#include "field/communication/PackInfo.h"
#include "field/communication/UniformMPIDatatypeInfo.h"
#include "field/iterators/FieldIterator.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/BlockForestEvaluation.h"
#include "lbm/PerformanceEvaluation.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/communication/PdfFieldMPIDatatypeInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/sweeps/SplitSweep.h"
#include "lbm/sweeps/SweepWrappers.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "sqlite/SQLite.h"

#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/BlockCellDataWriter.h"
#include "vtk/Initialization.h"
#include "vtk/VTKOutput.h"

#include <cstdlib>
#include <iostream>
#include <memory>



namespace uniform_grid {

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;
using walberla::real_t;

//////////////
// TYPEDEFS //
//////////////

typedef lbm::D3Q19< lbm::collision_model::SRT,           false > D3Q19_SRT_INCOMP;
typedef lbm::D3Q19< lbm::collision_model::SRT,           true  > D3Q19_SRT_COMP;
typedef lbm::D3Q19< lbm::collision_model::TRT,           false > D3Q19_TRT_INCOMP;
typedef lbm::D3Q19< lbm::collision_model::TRT,           true  > D3Q19_TRT_COMP;
typedef lbm::D3Q19< lbm::collision_model::D3Q19MRT,      false > D3Q19_MRT_INCOMP;
typedef lbm::D3Q27< lbm::collision_model::D3Q27Cumulant, true  > D3Q27_CUMULANT_COMP;


template< typename LatticeModel_T >
struct Types
{
   using Stencil_T = typename LatticeModel_T::Stencil;
   using CommunicationStencil_T = typename LatticeModel_T::CommunicationStencil;
   using PdfField_T = lbm::PdfField< LatticeModel_T >;
};

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID    UBB_Flag( "velocity bounce back" );
const FlagUID NoSlip_Flag( "no slip" );






/////////////////////
// OUTPUT HELPERS  //
/////////////////////

template< typename LatticeModel_T, class Enable = void >
struct StencilString;

template< typename LatticeModel_T >
struct StencilString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q19 >::value >::type >
{
   static const char * str() { return "D3Q19"; }
};


template< typename LatticeModel_T >
struct StencilString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::Stencil, stencil::D3Q27 >::value >::type >
{
  static const char * str() { return "D3Q27"; }
};

template< typename LatticeModel_T, class Enable = void >
struct CollisionModelString;

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::SRT_tag >::value >::type >
{
   static const char * str() { return "SRT"; }
};

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::TRT_tag >::value >::type >
{
   static const char * str() { return "TRT"; }
};

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::MRT_tag >::value >::type >
{
   static const char * str() { return "MRT"; }
};

template< typename LatticeModel_T >
struct CollisionModelString< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                          lbm::collision_model::Cumulant_tag >::value >::type >
{
  static const char * str() { return "Cumulant"; }
};





/////////////
// CONFIG  //
/////////////

static inline void getCells( const Config::BlockHandle & configBlock, uint_t & xCells, uint_t & yCells, uint_t & zCells )
{
   xCells = configBlock.getParameter< uint_t >( "xCells", 10 );
   yCells = configBlock.getParameter< uint_t >( "yCells", 10 );
   zCells = configBlock.getParameter< uint_t >( "zCells", 10 );
}

static inline void getCellsAndProcesses( const Config::BlockHandle & configBlock,
                                         const uint_t numberOfProcesses, const uint_t blocksPerProcess,
                                         uint_t & xCells, uint_t & yCells, uint_t & zCells,
                                         uint_t & xProcesses, uint_t & yProcesses, uint_t & zProcesses,
                                         uint_t & xBlocks, uint_t & yBlocks, uint_t & zBlocks )
{
   getCells( configBlock, xCells, yCells, zCells );

   std::vector< real_t > weighting;
   weighting.push_back( configBlock.getParameter< real_t >( "xWeight", real_t(1) ) / real_c(xCells) );
   weighting.push_back( configBlock.getParameter< real_t >( "yWeight", real_t(1) ) / real_c(yCells) );
   weighting.push_back( configBlock.getParameter< real_t >( "zWeight", real_t(1) ) / real_c(zCells) );

   std::vector< uint_t > processes = math::getFactors( numberOfProcesses, 3, weighting );

   WALBERLA_CHECK_EQUAL( processes.size(), 3 );
   WALBERLA_CHECK_EQUAL( processes[0] * processes[1] * processes[2], numberOfProcesses );

   xProcesses = processes[0];
   yProcesses = processes[1];
   zProcesses = processes[2];

   Vector3< uint_t > blocksPerProcess3D = math::getFactors3D( blocksPerProcess );

   xBlocks = processes[0] * blocksPerProcess3D[0];
   yBlocks = processes[1] * blocksPerProcess3D[1];
   zBlocks = processes[2] * blocksPerProcess3D[2];
}






///////////////////////////
// BLOCK STRUCTURE SETUP //
///////////////////////////

void createSetupBlockForest( blockforest::SetupBlockForest & sforest, const Config::BlockHandle & configBlock, const uint_t numberOfProcesses, const uint_t blocksPerProcess )
{
   uint_t numberOfXCellsPerBlock;
   uint_t numberOfYCellsPerBlock;
   uint_t numberOfZCellsPerBlock;
   uint_t numberOfXProcesses;
   uint_t numberOfYProcesses;
   uint_t numberOfZProcesses;
   uint_t numberOfXBlocks;
   uint_t numberOfYBlocks;
   uint_t numberOfZBlocks;

   getCellsAndProcesses( configBlock, numberOfProcesses, blocksPerProcess,
                         numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock,
                         numberOfXProcesses, numberOfYProcesses, numberOfZProcesses,
                         numberOfXBlocks, numberOfYBlocks, numberOfZBlocks );

   sforest.addWorkloadMemorySUIDAssignmentFunction( blockforest::uniformWorkloadAndMemoryAssignment );

   sforest.init( AABB( real_t(0), real_t(0), real_t(0), real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                        real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                        real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                                                        numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, false, false, false );

   if( MPIManager::instance()->numProcesses() > 1 )
   {
      shared_ptr< std::vector<uint_t> > processIdMap;

      WALBERLA_MPI_SECTION()
      {
         if ( MPIManager::instance()->isCartesianCommValid() )
         {
            MPIManager::instance()->createCartesianComm(numberOfXProcesses, numberOfYProcesses, numberOfZProcesses, false, false, false);

            processIdMap = make_shared<std::vector<uint_t> >(numberOfXProcesses * numberOfYProcesses * numberOfZProcesses);

            for (uint_t z = 0; z != numberOfZProcesses; ++z) {
               for (uint_t y = 0; y != numberOfYProcesses; ++y) {
                  for (uint_t x = 0; x != numberOfXProcesses; ++x)
                  {
                     (*processIdMap)[z * numberOfXProcesses * numberOfYProcesses + y * numberOfXProcesses + x] =
                           uint_c(MPIManager::instance()->cartesianRank(x, y, z));
                  }
               }
            }
         }
         else {
            WALBERLA_LOG_WARNING_ON_ROOT( "Your version of OpenMPI contains a bug. See waLBerla issue #73 for more "
                                          "information. As a workaround, MPI_COMM_WORLD instead of a "
                                          "Cartesian MPI communicator is used." );
            MPIManager::instance()->useWorldComm();
         }
      }
      else
         MPIManager::instance()->useWorldComm();

      sforest.balanceLoad( blockforest::CartesianDistribution( numberOfXProcesses, numberOfYProcesses, numberOfZProcesses, processIdMap.get() ),
                           numberOfXProcesses * numberOfYProcesses * numberOfZProcesses );
   }
   else
   {
      MPIManager::instance()->useWorldComm();

      sforest.balanceLoad( blockforest::CartesianDistribution( numberOfXProcesses, numberOfYProcesses, numberOfZProcesses, nullptr ),
                           numberOfXProcesses * numberOfYProcesses * numberOfZProcesses, real_t(0), 0, true );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "SetupBlockForest created successfully:\n" << sforest );
}



shared_ptr< blockforest::StructuredBlockForest > createStructuredBlockForest( const Config::BlockHandle & configBlock )
{
   uint_t numberOfXCellsPerBlock;
   uint_t numberOfYCellsPerBlock;
   uint_t numberOfZCellsPerBlock;
   getCells( configBlock, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );

   const uint_t blocksPerProcess = configBlock.getParameter< uint_t >( "blocksPerProcess", uint_t( 1 ) );

   if( configBlock.isDefined( "sbffile" ) )
   {
      std::string sbffile = configBlock.getParameter< std::string >( "sbffile" );

      WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure: loading from file \'" << sbffile << "\' ..." );

      return blockforest::createUniformBlockGrid( sbffile, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock, false );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure ..." );

   blockforest::SetupBlockForest sforest;
   createSetupBlockForest( sforest, configBlock, uint_c( MPIManager::instance()->numProcesses() ), blocksPerProcess );

   auto bf = std::make_shared< blockforest::BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, false );

   auto sbf = std::make_shared< blockforest::StructuredBlockForest >( bf, numberOfXCellsPerBlock,
                                                                                                            numberOfYCellsPerBlock,
                                                                                                            numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();

   return sbf;
}






///////////////////////
// BOUNDARY HANDLING //
///////////////////////

template< typename LatticeModel_T >
class MyBoundaryHandling
{
public:

   typedef lbm::NoSlip< LatticeModel_T, flag_t >    NoSlip_T;
   typedef lbm::SimpleUBB< LatticeModel_T, flag_t > UBB_T;

   typedef BoundaryHandling< FlagField_T, typename Types<LatticeModel_T>::Stencil_T, NoSlip_T, UBB_T > BoundaryHandling_T;



   MyBoundaryHandling( const BlockDataID & flagField, const BlockDataID & pdfField, const real_t velocity ) :
      flagField_( flagField ), pdfField_( pdfField ), velocity_( velocity ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagField_;
   const BlockDataID  pdfField_;

   const real_t velocity_;

}; // class MyBoundaryHandling

template< typename LatticeModel_T >
typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T *
MyBoundaryHandling<LatticeModel_T>::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   using PdfField_T = typename Types< LatticeModel_T >::PdfField_T;

   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagField_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfField_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                                           NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                                           UBB_T( "velocity bounce back", UBB_Flag, pdfField, velocity_, real_c(0), real_c(0) ) );

   CellInterval domainBB = storage->getDomainCellBB();
   storage->transformGlobalToBlockLocalCellInterval( domainBB, *block );

   // no slip WEST
   CellInterval west( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, west );

   // no slip EAST
   CellInterval east( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, east );

   // no slip SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, south );

   // no slip NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( NoSlip_Flag, north );

   // no slip BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   handling->forceBoundary( NoSlip_Flag, bottom );

   // velocity bounce back TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( UBB_Flag, top );

   handling->fillWithDomain( domainBB );

   return handling;
}






/////////
// VTK //
/////////

template< typename LatticeModel_T >
class MyVTKOutput {

public:

   MyVTKOutput( const ConstBlockDataID & pdfField, const ConstBlockDataID & flagField,
                const vtk::VTKOutput::BeforeFunction& pdfGhostLayerSync ) :
      pdfField_( pdfField ), flagField_( flagField ), pdfGhostLayerSync_( pdfGhostLayerSync ) {}

   void operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                    std::map< std::string, vtk::VTKOutput::CellFilter > &          filters,
                    std::map< std::string, vtk::VTKOutput::BeforeFunction > &      beforeFunctions );

private:

   const ConstBlockDataID pdfField_;
   const ConstBlockDataID flagField_;

   vtk::VTKOutput::BeforeFunction pdfGhostLayerSync_;

}; // class MyVTKOutput

template< typename LatticeModel_T >
void MyVTKOutput<LatticeModel_T>::operator()( std::vector< shared_ptr<vtk::BlockCellDataWriterInterface> > & writers,
                                              std::map< std::string, vtk::VTKOutput::CellFilter > &          filters,
                                              std::map< std::string, vtk::VTKOutput::BeforeFunction > &      beforeFunctions )
{
   // block data writers

   writers.push_back( make_shared< lbm::VelocityVTKWriter<LatticeModel_T> >( pdfField_, "VelocityFromPDF" ) );
   writers.push_back( make_shared< lbm::DensityVTKWriter<LatticeModel_T> >( pdfField_, "DensityFromPDF" ) );

   writers.push_back( make_shared< field::VTKWriter< FlagField_T > >( flagField_, "FlagField" ) );

   // cell filters

   field::FlagFieldCellFilter<FlagField_T> fluidFilter( flagField_ );
   fluidFilter.addFlag( Fluid_Flag );
   filters[ "FluidFilter" ] = fluidFilter;

   field::FlagFieldCellFilter<FlagField_T> obstacleFilter( flagField_ );
   obstacleFilter.addFlag( NoSlip_Flag );
   obstacleFilter.addFlag(    UBB_Flag );
   filters[ "ObstacleFilter" ] = obstacleFilter;

   // before functions

   beforeFunctions[ "PDFGhostLayerSync" ] = pdfGhostLayerSync_;
}






////////////////////
// THE SIMULATION //
////////////////////



template< typename LatticeModel_T, class Enable = void >
struct AddLB
{
   using PdfField = typename Types< LatticeModel_T >::PdfField_T;
   using CommunicationStencil = typename Types< LatticeModel_T >::CommunicationStencil_T;

   static void add( shared_ptr< blockforest::StructuredBlockForest > & blocks, SweepTimeloop & timeloop,
                    const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const BlockDataID & boundaryHandlingId,
                    const bool split, const bool pure, const bool fullComm, const bool fused, const bool directComm )
   {
      // setup of the LB communication for synchronizing the pdf field between neighboring blocks

      std::function< void () > commFunction;
      if( directComm )
      {
         if( fullComm )
         {
            blockforest::communication::UniformDirectScheme< stencil::D3Q27 > comm( blocks );
            auto mpiDatatypeInfo = make_shared<field::communication::UniformMPIDatatypeInfo< PdfField > >( pdfFieldId );
            comm.addDataToCommunicate( mpiDatatypeInfo );
            commFunction = comm;
         }
         else
         {
            blockforest::communication::UniformDirectScheme< CommunicationStencil > comm( blocks );
            auto mpiDatatypeInfo = make_shared<lbm::communication::PdfFieldMPIDatatypeInfo< PdfField > >( pdfFieldId );
            comm.addDataToCommunicate( mpiDatatypeInfo );
            commFunction = comm;
         }
      }
      else
      {
         if( fullComm )
         {
            blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > scheme( blocks );
            scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField > >( pdfFieldId ) );
            commFunction = scheme;
         }
         else
         {
            blockforest::communication::UniformBufferedScheme< CommunicationStencil > scheme( blocks );
            scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );
            commFunction = scheme;
         }
      }

      if( fused )
      {
         timeloop.add() << BeforeFunction( commFunction, "LB communication" )
                        << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "LB boundary sweep" );

         if( split )
         {
            if( pure )
               timeloop.add() << Sweep( lbm::SplitPureSweep< LatticeModel_T >( pdfFieldId ), "split pure LB sweep (stream & collide)" );
            else
               timeloop.add() << Sweep( lbm::SplitSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag ), "split LB sweep (stream & collide)" );
         }
         else
            timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag ) ), "cell-wise LB sweep (stream & collide)" );
      }
      else
      {
         if( split )
         {
            if( pure )
            {
               using Sweep_T = lbm::SplitPureSweep< LatticeModel_T >;
               auto sweep = make_shared< Sweep_T >( pdfFieldId );

               timeloop.add() << Sweep( lbm::CollideSweep< Sweep_T >( sweep ), "split pure LB sweep (collide)" );
               timeloop.add() << BeforeFunction( commFunction, "LB communication" )
                              << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "LB boundary sweep" );
               timeloop.add() << Sweep( lbm::StreamSweep< Sweep_T >( sweep ), "split pure LB sweep (stream)" );
            }
            else
            {
               typedef lbm::SplitSweep< LatticeModel_T, FlagField_T > Sweep_T;
               auto sweep = make_shared< Sweep_T >( pdfFieldId, flagFieldId, Fluid_Flag );

               timeloop.add() << Sweep( lbm::CollideSweep< Sweep_T >( sweep ), "split LB sweep (collide)" );
               timeloop.add() << BeforeFunction( commFunction, "LB communication" )
                              << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "LB boundary sweep" );
               timeloop.add() << Sweep( lbm::StreamSweep< Sweep_T >( sweep ), "split LB sweep (stream)" );
            }
         }
         else
         {
            auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );

            timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );
            timeloop.add() << BeforeFunction( commFunction, "LB communication" )
                           << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "LB boundary sweep" );
            timeloop.add() << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
         }
      }
   }
};

template< typename LatticeModel_T  >
struct AddLB< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag, lbm::collision_model::MRT_tag >::value ||
                                                       std::is_same< typename LatticeModel_T::CollisionModel::tag, lbm::collision_model::Cumulant_tag >::value
                                                       >::type >
{
   using PdfField = typename Types< LatticeModel_T >::PdfField_T;
   using CommunicationStencil = typename Types< LatticeModel_T >::CommunicationStencil_T;

   static void add( shared_ptr< blockforest::StructuredBlockForest > & blocks, SweepTimeloop & timeloop,
                    const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const BlockDataID & boundaryHandlingId,
                    const bool /*split*/, const bool /*pure*/, const bool fullComm, const bool fused, const bool directComm )
   {
      // setup of the LB communication for synchronizing the pdf field between neighboring blocks

      std::function< void () > commFunction;
      if( directComm )
      {
         
         if( fullComm )
         {
            blockforest::communication::UniformDirectScheme< stencil::D3Q27 > comm( blocks );
            auto mpiDatatypeInfo = make_shared<field::communication::UniformMPIDatatypeInfo< PdfField > >( pdfFieldId );
            comm.addDataToCommunicate( mpiDatatypeInfo );
            commFunction = comm;
         }
         else
         {
            blockforest::communication::UniformDirectScheme< CommunicationStencil > comm( blocks );
            auto mpiDatatypeInfo = make_shared<lbm::communication::PdfFieldMPIDatatypeInfo< PdfField > >( pdfFieldId );
            comm.addDataToCommunicate( mpiDatatypeInfo );
            commFunction = comm;
         }
      }
      else
      {
         if( fullComm )
         {
             blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > scheme( blocks );
             scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField > >( pdfFieldId ) );
             commFunction = scheme;
         }
         else
         {
            blockforest::communication::UniformBufferedScheme< CommunicationStencil > scheme( blocks );
            scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId ) );
            commFunction = scheme;
         }
      }

      if( fused )
      {
         timeloop.add() << BeforeFunction( commFunction, "LB communication" )
                        << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "LB boundary sweep" );
         timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag ) ), "cell-wise LB sweep (stream & collide)" );
      }
      else
      {
         auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldId, flagFieldId, Fluid_Flag );

         timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );
         timeloop.add() << BeforeFunction( commFunction, "LB communication" )
                        << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), "LB boundary sweep" );
         timeloop.add() << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
      }
   }
};



template< typename LatticeModel_T >
void run( const shared_ptr< Config > & config, const LatticeModel_T & latticeModel,
          const bool split, const bool pure, const bool fzyx, const bool fullComm, const bool fused, const bool directComm )
{
   using PdfField = typename Types<LatticeModel_T>::PdfField_T;

   Config::BlockHandle configBlock = config->getBlock( "UniformGrid" );

   // creating the block structure

   auto blocks = createStructuredBlockForest( configBlock );

   // add pdf field to blocks

   BlockDataID pdfFieldId = fzyx ? lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                              Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::fzyx ) :
                                   lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel,
                                                              Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_t(1),
                                                              FieldGhostLayers, field::zyxf );

   // add flag field to blocks

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // add LB boundary handling to blocks

   const real_t velocity = configBlock.getParameter< real_t >( "velocity", real_t(0.05) );

   BlockDataID boundaryHandlingId = blocks->template addStructuredBlockData< typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T >(
            MyBoundaryHandling< LatticeModel_T >( flagFieldId, pdfFieldId, velocity ), "boundary handling" );

   // creating the time loop

   const uint_t outerTimeSteps = configBlock.getParameter< uint_t >( "outerTimeSteps", uint_c(1 ) );
   const uint_t innerTimeSteps = configBlock.getParameter< uint_t >( "innerTimeSteps", uint_c(10) );

   SweepTimeloop timeloop( blocks->getBlockStorage(), outerTimeSteps * innerTimeSteps );

   // VTK

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
   pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo< PdfField > >( pdfFieldId ) );

   MyVTKOutput< LatticeModel_T > myVTKOutput( pdfFieldId, flagFieldId, pdfGhostLayerSync );

   std::map< std::string, vtk::SelectableOutputFunction > vtkOutputFunctions;
   vtk::initializeVTKOutput( vtkOutputFunctions, myVTKOutput, blocks, config );

   for( auto output = vtkOutputFunctions.begin(); output != vtkOutputFunctions.end(); ++output )
      timeloop.addFuncBeforeTimeStep( output->second.outputFunction, std::string("VTK: ") + output->first,
                                      output->second.requiredGlobalStates, output->second.incompatibleGlobalStates );

   // add LB kernel, boundary handling, and communication to time loop

   AddLB< LatticeModel_T >::add( blocks, timeloop, pdfFieldId, flagFieldId, boundaryHandlingId, split, pure, fullComm, fused, directComm );

   // logging right before the benchmark starts

   lbm::BlockForestEvaluation< FlagField_T > blockForest( blocks, flagFieldId, Fluid_Flag );
   blockForest.logInfoOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark parameters:"
                              "\n- collision model:                 " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n- stencil:                         " << StencilString< LatticeModel_T >::str() <<
                              "\n- compressible:                    " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n- fused (stream & collide) kernel: " << ( fused ? "yes" : "no" ) <<
                              "\n- split (collision) kernel:        " << ( split ? "yes" : "no" ) <<
                              "\n- pure kernel:                     " << ( pure ? "yes (collision is also performed within obstacle cells)" : "no" ) <<
                              "\n- data layout:                     " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n- communication:                   " << ( fullComm ? "full synchronization" : "direction-aware optimizations" ) <<
                              "\n- direct communication:            " << ( directComm ? "enabled" : "disabled" ) );

   // run the benchmark

   lbm::PerformanceEvaluation< FlagField_T > performance( blocks, flagFieldId, Fluid_Flag );

   for( uint_t outerRun = 0; outerRun < outerTimeSteps; ++outerRun )
   {
      WcTimingPool timeloopTiming;

      WALBERLA_MPI_WORLD_BARRIER();
      WcTimer timer;
      timer.start();

      for( uint_t innerRun = 0; innerRun < innerTimeSteps; ++innerRun )
         timeloop.singleStep( timeloopTiming );

      timer.end();

      double time = timer.max();
      mpi::reduceInplace( time, mpi::MAX );

      const auto reducedTimeloopTiming = timeloopTiming.getReduced();

      WALBERLA_LOG_RESULT_ON_ROOT( "Time loop timing:\n" << *reducedTimeloopTiming );

      performance.logResultOnRoot( innerTimeSteps, time );

      WALBERLA_ROOT_SECTION()
      {
         // logging in SQL database

         if( configBlock.getParameter< bool >( "logToSqlDB", true ) )
         {
            const std::string sqlFile = configBlock.getParameter< std::string >( "sqlFile", "performance.sqlite" );

            std::map< std::string, int >         integerProperties;
            std::map< std::string, double >      realProperties;
            std::map< std::string, std::string > stringProperties;

            performance.getResultsForSQLOnRoot( integerProperties, realProperties, stringProperties, innerTimeSteps, time );
            blockForest.getResultsForSQLOnRoot( integerProperties, realProperties, stringProperties );

            stringProperties[ "collisionModel" ]    = CollisionModelString< LatticeModel_T >::str();
            stringProperties[ "stencil" ]           = StencilString< LatticeModel_T >::str();
            stringProperties[ "compressible" ]      = ( LatticeModel_T::compressible ? "yes" : "no" );
            stringProperties[ "fusedKernel" ]       = ( fused ? "yes" : "no" );
            stringProperties[ "splitKernel" ]       = ( split ? "yes" : "no" );
            stringProperties[ "pureKernel" ]        = ( pure ? "yes" : "no" );
            stringProperties[ "dataLayout" ]        = ( fzyx ? "fzyx" : "zyxf" );
            stringProperties[ "fullCommunication" ] = ( fullComm ? "yes" : "no" );
            stringProperties[ "directComm"]         = ( directComm ? "yes" : "no" );

            auto runId = sqlite::storeRunInSqliteDB( sqlFile, integerProperties, stringProperties, realProperties );
            sqlite::storeTimingPoolInSqliteDB( sqlFile, runId, *reducedTimeloopTiming, "Timeloop" );
         }
      }
   }

   // logging once again at the end of the simulation, identical to logging at the beginning :-)

   blockForest.logInfoOnRoot();

   WALBERLA_LOG_INFO_ON_ROOT( "Benchmark parameters:"
                              "\n- collision model:                 " << CollisionModelString< LatticeModel_T >::str() <<
                              "\n- stencil:                         " << StencilString< LatticeModel_T >::str() <<
                              "\n- compressible:                    " << ( LatticeModel_T::compressible ? "yes" : "no" ) <<
                              "\n- fused (stream & collide) kernel: " << ( fused ? "yes" : "no" ) <<
                              "\n- split (collision) kernel:        " << ( split ? "yes" : "no" ) <<
                              "\n- pure kernel:                     " << ( pure ? "yes (collision is also performed within obstacle cells)" : "no" ) <<
                              "\n- data layout:                     " << ( fzyx ? "fzyx (structure of arrays [SoA])" : "zyxf (array of structures [AoS])" ) <<
                              "\n- communication:                   " << ( fullComm ? "full synchronization" : "direction-aware optimizations" ) <<
                              "\n- direct communication:            " << ( directComm ? "enabled" : "disabled" ) );

}


//////////
// MAIN //
//////////

enum CM { CMSRT, CMTRT, CMMRT, CMCUM };

int main( int argc, char **argv )
{
   mpi::Environment env( argc, argv );

   if( argc < 2 )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::cout << "Usage: " << argv[0] << " path-to-configuration-file [--trt | --mrt] [--comp] [--not-split] [--not-pure] [--zyxf] [--full-comm] [--not-fused] [--direct-comm]\n"
                      "\n"
                      "By default, SRT is selected as collision model, a communication with direction-aware optimizations is chosen, and an\n"
                      "incompressible, split, pure LB kernel is executed on a PDF field with layout 'fzyx' (= structure of arrays [SoA]).\n"
                      "\n"
                      "Optional arguments:\n"
                      " --trt:         collision model = TRT\n"
                      " --mrt:         collision model = MRT\n"
                      " --cumulant     collision model = cumulant\n"
                      " --comp:        LB kernel is switched from incompressible to compressible\n"
                      " --not-split:   LB kernel NOT split by PDF direction but executed cell by cell\n"
                      " --not-pure:    LB kernel is only executed in fluid cells, not in obstacle/boundary cells.\n"
                      "                Will be automatically selected for non-split LB kernels.\n"
                      " --zyxf:        data layout switched to 'zyxf' (array of structures [AoS])\n"
                      "                Probably the best layout for non-split kernels.\n"
                      " --full-comm:   A full synchronization of neighboring blocks is performed instead of using a communication\n"
                      "                that uses direction-aware optimizations.\n"
                      " --not-fused:   Selects separate LB kernels for collision and streaming.\n"
                      "                By default, a 'fused' stream & collide kernel is used.\n"
                      " --direct-comm: Enables bufferless direct communication\n"
                      "\n"
                      "Please note: For small/very small blocks (i.e., blocks with only few cells), the best performance may be achieved with\n"
                      "             basic (non-split), incompressible LB kernels combined with an array of structures ('zyxf') data layout!" << std::endl;
      }
      return EXIT_SUCCESS;
   }

   logging::Logging::printHeaderOnStream();
   //WALBERLA_ROOT_SECTION() { logging::Logging::instance()->setLogLevel( logging::Logging::PROGRESS ); }

#ifdef _OPENMP
   if( std::getenv( "OMP_NUM_THREADS" ) == nullptr )
      WALBERLA_ABORT( "If you are using a version of the benchmark that was compiled with OpenMP you have to "
                      "specify the environment variable \'OMP_NUM_THREADS\' accordingly!" );
#endif

   // open configuration file

   shared_ptr< Config > config = make_shared< Config >();
   config->readParameterFile( argv[1] );

   Config::BlockHandle configBlock = config->getBlock( "UniformGrid" );

   if( !configBlock )
      WALBERLA_ABORT( "You have to specify a \"UniformGrid\" block in the configuration file!" );

   // In case 'sbffile' and 'processes' are specified in the configuration file:
   // -> just create the block structure and save it to file for later use

   if( configBlock.isDefined( "sbffile" ) && configBlock.isDefined( "processes" ) )
   {
      std::string  sbffile           = configBlock.getParameter< std::string >( "sbffile" );
      const uint_t numberOfProcesses = configBlock.getParameter< uint_t >( "processes" );
      const uint_t blocksPerProcess  = configBlock.getParameter< uint_t >( "blocksPerProcess", uint_t( 1 ) );

      std::ostringstream infoString;
      infoString << "You have selected the option of just creating the block structure (= domain decomposition) and saving the result to file\n"
                    "by specifying the output file name \'" << sbffile << "\' AND providing a targeted number of processes ("
                 << numberOfProcesses << ").\n";

      if( MPIManager::instance()->numProcesses() > 1 )
         WALBERLA_ABORT( infoString.str() << "In this mode you need to start " << argv[0] << " with just one process!" );

      WALBERLA_LOG_INFO_ON_ROOT( infoString.str() << "Creating the block structure ..." );

      blockforest::SetupBlockForest sforest;
      createSetupBlockForest( sforest, configBlock, numberOfProcesses, blocksPerProcess );

      sforest.saveToFile( sbffile.c_str() );

      logging::Logging::printFooterOnStream();
      return EXIT_SUCCESS;
   }

   // reading optional parameters from passed arguments

   CM collisionModel = CMSRT;
   bool compressible = false;
   bool split        = true;
   bool pure         = true;
   bool fzyx         = true;
   bool fullComm     = false;
   bool fused        = true;
   bool directComm   = false;

   for( int i = 2; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--trt" )         == 0 ) collisionModel = CMTRT;
      if( std::strcmp( argv[i], "--mrt" )         == 0 ) collisionModel = CMMRT;
      if( std::strcmp( argv[i], "--cumulant" )    == 0 ) collisionModel = CMCUM;
      if( std::strcmp( argv[i], "--comp" )        == 0 ) compressible   = true;
      if( std::strcmp( argv[i], "--not-split" )   == 0 ) split          = false;
      if( std::strcmp( argv[i], "--not-pure" )    == 0 ) pure           = false;
      if( std::strcmp( argv[i], "--zyxf" )        == 0 ) fzyx           = false;
      if( std::strcmp( argv[i], "--full-comm" )   == 0 ) fullComm       = true;
      if( std::strcmp( argv[i], "--not-fused" )   == 0 ) fused          = false;
      if( std::strcmp( argv[i], "--direct-comm" ) == 0 ) directComm     = true;
   }

   if( pure && !split )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "You called the benchmark with \"--not-split\" but without \"--not-pure\".\n"
                                    "\"Pure\" kernels are only available for \"split\" kernels! Setting \"pure\" to false ..." );
      pure = pure && split; // pure only works in combination with split
   }

   if( collisionModel == CMMRT && ( compressible || split || pure ) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Option \"--comp\" is not available for MRT! Also, MRT requires options \"--not-split\" and \"--not-pure\"!\n"
                                    "Setting \"compressible\", \"split\", and \"pure\" to false ..." );
      compressible = false;
      split        = false;
      pure         = false;
   }
   if( collisionModel == CMCUM && ( compressible || split || pure ) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Option \"--comp\" has to be set for Cumulant! Also, Cumulant requires options \"--not-split\" and \"--not-pure\"!\n"
                                            "Setting \"compressible\", to true, \"split\", and \"pure\" to false ..." );
      compressible = true;
      split        = false;
      pure         = false;
   }

   WALBERLA_NON_MPI_SECTION()
   {
      if( directComm )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "Direct (bufferless) communication is not available when building with MPI disabled.\n"
                                       "Switching to buffered comm..." );
      }
      directComm = false;
   }

   const real_t omega = configBlock.getParameter< real_t >( "omega", real_t(1.4) ); // on the coarsest grid!

   // executing benchmark

   if( collisionModel == CMSRT ) // SRT
   {
      if( compressible )
      {
         D3Q19_SRT_COMP latticeModel = D3Q19_SRT_COMP( lbm::collision_model::SRT( omega ) );
         run( config, latticeModel, split, pure, fzyx, fullComm, fused, directComm );
      }
      else
      {
         D3Q19_SRT_INCOMP latticeModel = D3Q19_SRT_INCOMP( lbm::collision_model::SRT( omega ) );
         run( config, latticeModel, split, pure, fzyx, fullComm, fused, directComm );
      }
   }
   else if( collisionModel == CMTRT ) // TRT
   {
      if( compressible )
      {
         D3Q19_TRT_COMP latticeModel = D3Q19_TRT_COMP( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
         run( config, latticeModel, split, pure, fzyx, fullComm, fused, directComm );
      }
      else
      {
         D3Q19_TRT_INCOMP latticeModel = D3Q19_TRT_INCOMP( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
         run( config, latticeModel, split, pure, fzyx, fullComm, fused, directComm );
      }
   }
   else if( collisionModel == CMMRT ) // MRT
   {
      D3Q19_MRT_INCOMP latticeModel = D3Q19_MRT_INCOMP( lbm::collision_model::D3Q19MRT::constructTRTWithMagicNumber( omega ) );
      run( config, latticeModel, split, pure, fzyx, fullComm, fused, directComm );
   }
   else  // Cumulant
   {
      D3Q27_CUMULANT_COMP latticeModel = D3Q27_CUMULANT_COMP( lbm::collision_model::D3Q27Cumulant(omega) );
      run( config, latticeModel, split, pure, fzyx, fullComm, fused, directComm );
   }

   logging::Logging::printFooterOnStream();
   return EXIT_SUCCESS;
}


} // namespace uniform_grid

int main( int argc, char ** argv )
{
	return uniform_grid::main( argc, argv );
}
