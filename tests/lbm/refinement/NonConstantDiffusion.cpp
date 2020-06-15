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
//! \file NonConstantDiffusion.cpp
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//! \brief Calculates temporal evolution of a temperature profile between two plates
//!
//! Configuration:
//!  - CHANNEL_LENGTH  distance between both plates (x)
//!  - CHANNEL_WIDTH   y and z size
//!  - OMEGA           relaxation parameter -> depending on thermal diffusivity k: omega = 1/(3.k+0.5)
//!  - DELTA_SCALAR    temperature difference between plates
//!  - NUM_TIMESTEPS   amount of timesteps to simulate
//! Setup:
//!  - periodic in y and z coordinate
//!  - plates in the x-z planes
//!  - no hydrodynamic computations
//!  - comparison with analytical solution in each timestep
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/SharedFunctor.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/config/Config.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/IntegerFactorization.h"
#include "core/math/Limits.h"
#include "core/mpi/Environment.h"
#include "core/stringToNum.h"

#include "domain_decomposition/SharedSweep.h"

#include "lbm/boundary/factories/DefaultDiffusionBoundaryHandling.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/refinement/PdfFieldSyncPackInfo.h"
#include "lbm/refinement/TimeStep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/GhostLayerField.h"
#include "field/communication/PackInfo.h"
#include "field/iterators/FieldPointer.h"

#include "stencil/D3Q19.h"

#include "timeloop/SweepTimeloop.h"

#include <stdexcept>
#include <functional>
#include <string>

#include "gather/GnuPlotGraphWriter.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"
#include "vtk/VTKOutput.h"



namespace walberla {

typedef GhostLayerField< real_t, 1 >          ScalarField;
typedef GhostLayerField< Vector3<real_t>, 1 > VectorField;

typedef lbm::D3Q19< lbm::collision_model::SRTField<ScalarField>, true, lbm::force_model::None, 1 >  LM;

using Stencil = LM::Stencil;
using MyPdfField = lbm::PdfField<LM>;

using flag_t = uint8_t;
using MyFlagField = FlagField<flag_t>;

typedef lbm::DefaultDiffusionBoundaryHandlingFactory< LM, MyFlagField > MyBoundaryHandling;

/////////////////////////////////////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, const uint_t levels )
{
   const AABB & domain = forest.getDomain();

   const real_t domainxMax = domain.xMax() / real_c( pow( real_t(2), int_c( levels - uint_t(1u) ) ) );

   AABB left( domain.xMin(), domain.yMin(), domain.zMin(),
              domainxMax,    domain.yMax(), domain.zMax() );

   for( auto block = forest.begin(); block != forest.end(); ++block )
      if( block->getAABB().intersects( left ) )
         if( block->getLevel() < ( levels - uint_t(1u) ) )
            block->setMarker( true );

}

static void workloadAndMemoryAssignment( SetupBlockForest& forest ) {

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
      block->setMemory( numeric_cast< memory_t >(1) );
   }
}


shared_ptr< StructuredBlockForest > makeStructuredBlockStorage( uint_t length, uint_t width, uint_t refinement )
{
   // initialize SetupBlockForest = determine domain decomposition
    SetupBlockForest sforest;

    uint_t cells[]  = { length, width, width  };
    uint_t blocks[] = { uint_t(1u), uint_t(1u), uint_t(1u) };
    sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, refinement ) );
    sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

    sforest.init(
       AABB( real_t(0),        real_t(0),        real_t(0),             // blocks/processes in x/y/z direction
             real_c(cells[0]), real_c(cells[1]), real_c(cells[2]) ),    // cells per block in x/y/z direction
             blocks[0]  , blocks[1]  , blocks[2],                       // one block per process
             false      , true       , true);                           // periodicity

    // calculate process distribution
    const memory_t memoryLimit = math::Limits< memory_t >::inf();

    sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

    MPIManager::instance()->useWorldComm();

    // create StructuredBlockForest (encapsulates a newly created BlockForest)
    shared_ptr< StructuredBlockForest > sbf =
          make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, false ),
                   cells[0]/blocks[0], cells[1]/blocks[1], cells[2]/blocks[2] );
    sbf->createCellBoundingBoxes();
    return sbf;
}

/////////////////////////////////////////////////////
// Boundary        //
/////////////////////

const FlagUID& getFluidFlag()           { static FlagUID uid( "Fluid" );           return uid; }

void setFlags( shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & boundaryHandlingId, uint_t ghostLayers, bool closed, real_t lval, real_t rval )
{
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      auto * boundaryHandling = block->getData< MyBoundaryHandling::BoundaryHandling_T >( boundaryHandlingId );

      const uint_t level = blocks->getLevel(*block);

      CellInterval domainBB = blocks->getDomainCellBB( level );
      blocks->transformGlobalToBlockLocalCellInterval( domainBB, *block );

      //const cell_idx_t width = cell_idx_c( uint_t(1) << level ) - cell_idx_t(1);

      domainBB.xMin() -= cell_idx_c( 1 );
      domainBB.xMax() += cell_idx_c( 1 );
      domainBB.yMin() -= cell_idx_c( ghostLayers );
      domainBB.yMax() += cell_idx_c( ghostLayers );
      domainBB.zMin() -= cell_idx_c( ghostLayers );
      domainBB.zMax() += cell_idx_c( ghostLayers );

      MyBoundaryHandling::SimpleDiffusionDirichlet_T::ScalarConfiguration scl( lval );
      MyBoundaryHandling::SimpleDiffusionDirichlet_T::ScalarConfiguration scr( rval );

      // LEFT
      CellInterval left( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(),
                         domainBB.xMin(), domainBB.yMax(), domainBB.zMax());
      boundaryHandling->forceBoundary( MyBoundaryHandling::getSimpleDiffusionDirichletFlagUID1(), left, scl );

      // RIGHT
      CellInterval right( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(),
                          domainBB.xMax(), domainBB.yMax(), domainBB.zMax());
      boundaryHandling->forceBoundary( MyBoundaryHandling::getSimpleDiffusionDirichletFlagUID2(), right, scr );


      if( closed)
      {
         domainBB.xMin() -= cell_idx_c( ghostLayers );
         domainBB.xMax() += cell_idx_c( ghostLayers );
         domainBB.yMin() -= cell_idx_c( ghostLayers );
         domainBB.yMax() += cell_idx_c( ghostLayers );
         domainBB.zMin() -= cell_idx_c( 1 );
         domainBB.zMax() += cell_idx_c( 1 );

         // TOP
         CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
         boundaryHandling->forceBoundary(MyBoundaryHandling::getFreeDiffusionFlagUID(), top );

         // BOTTOM
         CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
         boundaryHandling->forceBoundary( MyBoundaryHandling::getFreeDiffusionFlagUID(), bottom );
      }

       boundaryHandling->fillWithDomain( ghostLayers );
   }
}

class OmegaSweep
{
public:
   OmegaSweep( BlockDataID advDiffFieldID, real_t minOmega, real_t maxOmega, real_t minDensity, real_t maxDensity ) : 
      advDiffFieldID_(advDiffFieldID), minOmega_(minOmega), maxOmega_(maxOmega), minDensity_(minDensity), maxDensity_(maxDensity) {}

   void operator()( IBlock * block )
   {
      auto advDiffField = block->getData< MyPdfField >( advDiffFieldID_ );

      auto cm = advDiffField->latticeModel().collisionModel();
      for( auto it = advDiffField->begin(); it != advDiffField->end(); ++it ){
         real_t density = advDiffField->getDensity( it.cell() );
         real_t omega = ( maxOmega_ ) * ( density - minDensity_ )
                      - ( minOmega_ ) * ( density - maxDensity_ );
         cm.reset( it.x(), it.y(), it.z(), omega );
      }
   }

private:
   const BlockDataID advDiffFieldID_;
   const real_t minOmega_;
   const real_t maxOmega_;
   const real_t minDensity_;
   const real_t maxDensity_;
};

} // namespace walberla

int main( int argc, char **argv )
{
   using namespace walberla;
   debug::enterTestMode();
   mpi::Environment env( argc, argv );
   logging::Logging::printHeaderOnStream();

   uint_t length = uint_t( 16u  );
   uint_t width  = uint_t( 16u  );
   uint_t time   = uint_t( 10u  );
   real_t dv     = real_t(  1   );
   real_t v      = real_t(  1   );
   real_t domega = real_t(  0.8 );
   real_t omega  = real_t(  1.1 );
   bool   closed = false;
   uint_t levels = uint_t(  1u );

   bool useVTK = false;

   if( argc > 1 ) {
      std::vector<std::string> args( argv, argv + argc );
      for( uint_t i = 1; i < uint_c(argc); ++i ) {
              if( std::string(argv[i]) == "-l"    )   length  = stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-w"    )   width   = stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-t"    )   time    = stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "-dv"   )   dv      = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-v"    )   v       = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-do"   )   domega  = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-o"    )   omega   = stringToNum<real_t>( args[++i] );
         else if( std::string(argv[i]) == "-c"    )   closed  = stringToNum<int>( args[++i] ) != 0;
         else if( std::string(argv[i]) == "-r"    )   levels += stringToNum<uint_t>( args[++i] );
         else if( std::string(argv[i]) == "--vtk" )   useVTK  = true;
         else if( argv[i][0] != '-' ){
            std::cerr << "Usage: -option value" << std::endl; return EXIT_FAILURE;
         } else
            std::cerr << "Ignore unknown option " << argv[i++] << std::endl;
      }
   }

   uint_t ghostLayers = uint_t(4u);
   if(closed || levels == uint_t(1u))
   {
      levels      = uint_t(1u);
      ghostLayers = uint_t(1u);
   }

   auto blockStorage = makeStructuredBlockStorage( length, width, levels);

   BlockDataID omegaFieldID = field::addToStorage<ScalarField>( blockStorage, "Flag field", omega+real_c(0.5)*domega );

   LM          lm              = LM( lbm::collision_model::SRTField<ScalarField>( omegaFieldID ) );
   BlockDataID advDiffFieldID  = lbm::addPdfFieldToStorage( blockStorage, "PDF field", lm, Vector3<real_t>(), v+real_c(0.5)*dv, ghostLayers );

   BlockDataID velFieldID      = field::addToStorage<VectorField>( blockStorage, "Velocity field", Vector3<real_t>() );
   BlockDataID flagFieldID     = field::addFlagFieldToStorage<MyFlagField>( blockStorage, "Flag field", ghostLayers );

   BlockDataID boundaryHandling = MyBoundaryHandling::addDefaultDiffusionBoundaryHandlingToStorage( blockStorage, "DiffusionboundaryHandling", flagFieldID, getFluidFlag(), advDiffFieldID );

   setFlags( blockStorage, boundaryHandling, ghostLayers, closed, v+dv, v );

   SweepTimeloop timeloop( blockStorage->getBlockStorage(), time );

   if( levels == uint_t(1u) )
   {
      blockforest::communication::UniformBufferedScheme< Stencil > scheme( blockStorage );
      scheme.addPackInfo( make_shared< field::communication::PackInfo< MyPdfField > >( advDiffFieldID ) );
      timeloop.addFuncBeforeTimeStep( scheme, "Communication" );

      timeloop.add() << Sweep( MyBoundaryHandling::BoundaryHandling_T::getBlockSweep(boundaryHandling), "Boundary" );
      timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseAdvectionDiffusionSweep< LM, VectorField, MyFlagField >(
                                                   advDiffFieldID, velFieldID, flagFieldID, getFluidFlag() ) ), "LBM_SRT" );
   }
   else
   {
      auto mySweep = lbm::makeCellwiseAdvectionDiffusionSweep< LM, VectorField, MyFlagField >( advDiffFieldID, velFieldID, flagFieldID, getFluidFlag() );

      timeloop.addFuncBeforeTimeStep( makeSharedFunctor( lbm::refinement::makeTimeStep< LM, MyBoundaryHandling::BoundaryHandling_T >(
                                                            blockStorage, mySweep, advDiffFieldID, boundaryHandling  ) ), "LBM refinement time step" );
   }

   timeloop.add() << Sweep( OmegaSweep( advDiffFieldID, omega, omega+domega, v, v+dv ), "OmegaSweep" );

   WcTimingPool timeloopTiming;
   timeloop.run(timeloopTiming );
   timeloopTiming.logResultOnRoot();

   if( useVTK )
   {
      auto vtkOut = vtk::createVTKOutput_BlockData( *blockStorage, "block_data", uint_t(1u), uint_t(0u), false, "vtk_out/NonConstantDiffusion" );
      auto densityWriter  = make_shared< lbm::DensityVTKWriter<LM> >( advDiffFieldID, "E" );
      auto omegaWriter    = make_shared< field::VTKWriter<ScalarField>              >( omegaFieldID,   "omega" );
      vtkOut->addCellDataWriter( densityWriter );
      vtkOut->addCellDataWriter( omegaWriter   );
      vtkOut->write();
      vtk::writeDomainDecomposition( blockStorage, "domain_decomposition", "vtk_out/NonConstantDiffusion" );
      field::createVTKOutput<MyFlagField>( flagFieldID, *blockStorage, "flag_field", uint_t(1u), uint_t(1u), false, "vtk_out/NonConstantDiffusion" )();
   }

   logging::Logging::printFooterOnStream();

   return 0;
}

