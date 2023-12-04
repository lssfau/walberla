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
//! \file Initialization.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "BlockNeighborhoodSection.h"
#include "Initialization.h"
#include "SetupBlockForest.h"
#include "loadbalancing/Cartesian.h"
#include "loadbalancing/StaticCurve.h"

#include "core/Abort.h"
#include "core/cell/CellInterval.h"
#include "core/math/IntegerFactorization.h"
#include "core/mpi/MPIManager.h"

#include "stencil/D3Q19.h"

#include <functional>
#include <memory>

namespace walberla {
namespace blockforest {


//**********************************************************************************************************************
/*!
* Parses config block called 'DomainSetup' and creates a StructuredBlockForest
*
* For more information see function below.
*/
//**********************************************************************************************************************
shared_ptr< StructuredBlockForest > createUniformBlockGridFromConfig( const shared_ptr< Config > & config,
                                                                      CellInterval * requestedDomainSize,
                                                                      const bool keepGlobalBlockInformation )
{
   if( config != nullptr )
   {
      auto block = config->getGlobalBlock();
      if( block ) {
         auto subBlock = block.getBlock( "DomainSetup" );
         if ( ! subBlock ) {
            WALBERLA_ABORT_NO_DEBUG_INFO( "Unable to create uniform block grid from configuration file."
                                          "\nDid you forget to specify a \"DomainSetup\" block in the configuration file?" );
         }
         return createUniformBlockGridFromConfig( subBlock, requestedDomainSize, keepGlobalBlockInformation );
      }
   }
   WALBERLA_ABORT_NO_DEBUG_INFO( "No Configuration specified" );
   return shared_ptr<StructuredBlockForest>();
}



//**********************************************************************************************************************
/*!
* \brief Parses config block and creates a StructuredBlockForest
*
* Two possibilities:
*    1)  Using the cells per block and number of blocks for each direction
\verbatim
          {
             cellsPerBlock < 10, 20, 30 > ;  // required
             blocks        < 1,   2,  3 > ;  // required
             periodic      < 0, 0, 1 >;      // not required, defaults to no periodicity

             dx  0.01;  // defaults to 1.0
          }
\endverbatim
*          An optional config parameter 'cartesianSetup' is available. Its default, true, causes one block to be
*          assigned to each process. Setting it to false allows multiple blocks to be assigned to each process.
*          If the number of blocks is not divisible by the number of processes, the loadbalancer tries to assign
*          the blocks to processes as evenly as possible.
*    2) Using the number of global cells, \#blocks = \#processes, if this does not fit, extend the domain
\verbatim
          {
             cells <    10,40,90>;    // required
             periodic   < 0, 0, 1 >;  // not required, defaults to no periodicity
             dx  0.01;                // defaults to 1.0
          }
\endverbatim
*          An optional config parameter 'oneBlockPerProcess' is available. Setting it to false forces all
*          blocks to be assigned to a single process, which may be useful for debugging purposes. Otherwise,
*          one block is assigned to each process.
*          Example:  cells < 31,31,31> started using 8 processors <BR>
*                    calculated processor distribution <2,2,2>    <BR>
*                    real domain is then extended to <32,32,32> and every processor gets a block of <16,16,16>
*
*          When this setup is used and requestedDomainSize is not the null pointer, it is set to <BR>
*          the requested domain size ( in the example <31,31,31> )
*
*/
//**********************************************************************************************************************
shared_ptr< StructuredBlockForest > createUniformBlockGridFromConfig( const Config::BlockHandle & configBlock,
                                                                      CellInterval * requestedDomainSize,
                                                                      const bool keepGlobalBlockInformation )
{
   const Vector3<bool> periodic = configBlock.getParameter<Vector3<bool> >( "periodic",  Vector3<bool> (false) );
   const real_t        dx       = configBlock.getParameter<real_t        >( "dx",        real_t(1)             );

   Vector3<uint_t> cellsPerBlock;
   Vector3<uint_t> blocks;

   if ( configBlock.isDefined("cells") )
   {
      if ( configBlock.isDefined("cellsPerBlock") ||  configBlock.isDefined("blocks")  )
         WALBERLA_ABORT_NO_DEBUG_INFO("Config Error:  Use either ('cellsPerBlock' and 'blocks') or 'cells', not both!");

      Vector3<uint_t> cells = configBlock.getParameter<Vector3<uint_t>  >( "cells" );

      if ( requestedDomainSize )
         *requestedDomainSize = CellInterval( 0,0,0,
                                              cell_idx_c(cells[0])-1,
                                              cell_idx_c(cells[1])-1,
                                              cell_idx_c(cells[2])-1 );

      const uint_t nrOfProcesses = uint_c( MPIManager::instance()->numProcesses() );

      calculateCellDistribution( cells, nrOfProcesses, blocks, cellsPerBlock );
   }
   else
   {
      cellsPerBlock = configBlock.getParameter<Vector3<uint_t>  >( "cellsPerBlock" );
      blocks        = configBlock.getParameter<Vector3<uint_t>  >( "blocks"        );

      if ( requestedDomainSize )
         *requestedDomainSize = CellInterval( 0, 0, 0, cell_idx_c( cellsPerBlock[0] * blocks[0] ),
                                                       cell_idx_c( cellsPerBlock[1] * blocks[1] ),
                                                       cell_idx_c( cellsPerBlock[2] * blocks[2] ) );
   }

   const bool oneBlockPerProcess = configBlock.getParameter<bool> ( "oneBlockPerProcess", true );
   const bool cartesian = configBlock.getParameter<bool> ( "cartesianSetup", true );

   if ( !cartesian )
   {
      if ( configBlock.isDefined("oneBlockPerProcess") )
         WALBERLA_ABORT_NO_DEBUG_INFO("Config Error:  Set either 'oneBlockPerProcess' or set 'cartesianSetup' to false, not both!");

      return createUniformBlockGrid(
               blocks[0],        blocks[1],        blocks[2],         // blocks in x/y/z direction
               cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],  // cells per block in x/y/z direction
               dx,                                                    // cell size
               uint_t(0),                                             // maximum number of blocks per process
               true, false,                                           // include but don't force Metis
               periodic[0],      periodic[1],      periodic[2],       // periodicity
               keepGlobalBlockInformation                             // keep global block information
               );
   }

   return createUniformBlockGrid(
            blocks[0],        blocks[1],        blocks[2],         // blocks/processes in x/y/z direction
            cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],  // cells per block in x/y/z direction
            dx,                                                    // cell size
            oneBlockPerProcess,                                    // one block per process
            periodic[0], periodic[1], periodic[2],                 // periodicity
            keepGlobalBlockInformation                             // keep global block information
            );
}



//**********************************************************************************************************************
/*!
*   \brief Function for creating a block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks, each block has the same size.
*   The distribution of blocks to processes also follows a Cartesian decomposition.
*
*   \param domainAABB                 An axis-aligned bounding box that spans the entire simulation space/domain
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXProcesses         Number of processes the blocks are distributed to in x direction
*   \param numberOfYProcesses         Number of processes the blocks are distributed to in y direction
*   \param numberOfZProcesses         Number of processes the blocks are distributed to in z direction
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< BlockForest >
   createBlockForest(      const AABB& domainAABB,
                     const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,         const uint_t numberOfZBlocks,
                     const uint_t numberOfXProcesses,      const uint_t numberOfYProcesses,      const uint_t numberOfZProcesses,
                     const bool   xPeriodic /* = false */, const bool   yPeriodic /* = false */, const bool   zPeriodic /* = false */,
                     const bool keepGlobalBlockInformation /* = false */ ) {

   const uint_t numberOfProcesses = numberOfXProcesses * numberOfYProcesses * numberOfZProcesses;

   if( numeric_cast< int >( numberOfProcesses ) != MPIManager::instance()->numProcesses() )
      WALBERLA_ABORT( "The number of requested processes (" << numberOfProcesses << ") doesn't match the number "
                                                                                   "of active MPI processes (" << MPIManager::instance()->numProcesses() << ")!" );

   // initialize SetupBlockForest = determine domain decomposition

   SetupBlockForest sforest;

   sforest.addWorkloadMemorySUIDAssignmentFunction( uniformWorkloadAndMemoryAssignment );

   sforest.init( domainAABB, numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic );

   // if possible, create Cartesian MPI communicator

   std::vector< uint_t > processIdMap(0);

   WALBERLA_MPI_SECTION()
   {
      auto mpiManager = MPIManager::instance();
      if (!mpiManager->hasWorldCommSetup())
      {
         //create cartesian communicator only if not yet a cartesian communicator (or other communicator was created)
         if ( ! mpiManager->rankValid() )
         {
            mpiManager->createCartesianComm(numberOfXProcesses, numberOfYProcesses, numberOfZProcesses, xPeriodic,
                                            yPeriodic, zPeriodic);
         }

         processIdMap.resize( numberOfProcesses );

         for( uint_t z = 0; z != numberOfZProcesses; ++z ) {
            for( uint_t y = 0; y != numberOfYProcesses; ++y ) {
               for( uint_t x = 0; x != numberOfXProcesses; ++x ) {

                  processIdMap[ z * numberOfXProcesses * numberOfYProcesses + y * numberOfXProcesses + x ] =
                     uint_c( MPIManager::instance()->cartesianRank(x,y,z) );
               }
            }
         }
      }

   }

   // calculate process distribution

   sforest.balanceLoad( blockforest::CartesianDistribution( numberOfXProcesses, numberOfYProcesses, numberOfZProcesses, &processIdMap ),
                       numberOfXProcesses * numberOfYProcesses * numberOfZProcesses );

   // create StructuredBlockForest (encapsulates a newly created BlockForest)

   return std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation );
}




class FixedRefinementLevelSelector
{
public:
   FixedRefinementLevelSelector( const uint_t level ) : level_(level) {}
   void operator()( SetupBlockForest& forest )
   {
      for( auto block = forest.begin(); block != forest.end(); ++block )
      {
         if( block->getLevel() < level_ )
            block->setMarker( true );
      }
   }
private:
   uint_t level_;
};

std::unique_ptr<SetupBlockForest> createSetupBlockForest(const math::AABB& simulationDomain,
                                                         Vector3<uint_t> blocks,
                                                         const Vector3<bool>& isPeriodic,
                                                         const uint_t numberOfProcesses,
                                                         const uint_t initialRefinementLevel)
{
   WALBERLA_MPI_SECTION()
   {
      if (!MPIManager::instance()->rankValid())
         MPIManager::instance()->useWorldComm();
   }

   auto sforest = std::make_unique<SetupBlockForest>( );
   sforest->addWorkloadMemorySUIDAssignmentFunction( blockforest::uniformWorkloadAndMemoryAssignment );
   sforest->addRefinementSelectionFunction( FixedRefinementLevelSelector(initialRefinementLevel) );
   sforest->init( simulationDomain, blocks[0], blocks[1], blocks[2], isPeriodic[0], isPeriodic[1], isPeriodic[2] );

   WALBERLA_LOG_INFO_ON_ROOT( "Balancing " << sforest->getNumberOfBlocks() << " blocks for " << numberOfProcesses << " processes...");

   sforest->balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), numberOfProcesses, real_t(0), memory_t(0), false, true );
   return sforest;
}

shared_ptr<BlockForest> createBlockForest(const math::AABB& simulationDomain,
                                          const Vector3<uint_t>& blocks,
                                          const Vector3<bool>& isPeriodic,
                                          const uint_t numberOfProcesses,
                                          const uint_t initialRefinementLevel,
                                          const bool keepGlobalBlockInformation)
{
   WALBERLA_MPI_SECTION()
   {
      if (!MPIManager::instance()->rankValid())
         MPIManager::instance()->useWorldComm();
   }

   std::unique_ptr<SetupBlockForest> sforest( createSetupBlockForest( simulationDomain, blocks, isPeriodic, numberOfProcesses, initialRefinementLevel ));
   return std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), *sforest, keepGlobalBlockInformation );
}

shared_ptr<BlockForest> createBlockForest(const math::AABB& simulationDomain,
                                          const Vector3<uint_t>& blocks,
                                          const Vector3<bool>& isPeriodic,
                                          const bool setupRun,
                                          const std::string& sbffile,
                                          const uint_t numberOfProcesses,
                                          const uint_t initialRefinementLevel,
                                          const bool keepGlobalBlockInformation )
{
   if (setupRun)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Setup run. For production run specify 'setupRun = false'");

      if( MPIManager::instance()->numProcesses() > 1 )
         WALBERLA_LOG_WARNING_ON_ROOT( "Setup run with more than one process! Only root is doing work! I hope you know what you are doing!" );

      WALBERLA_ROOT_SECTION()
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure ..." );

         std::unique_ptr<SetupBlockForest> sforest( createSetupBlockForest(simulationDomain, blocks, isPeriodic, numberOfProcesses, initialRefinementLevel) );
         sforest->saveToFile( sbffile.c_str() );

         WALBERLA_LOG_INFO_ON_ROOT( "SetupBlockForest successfully saved to file!" );
      }

      return shared_ptr<BlockForest>();
   }

   WALBERLA_MPI_SECTION()
   {
      if (!MPIManager::instance()->rankValid())
         MPIManager::instance()->useWorldComm();
   }

   WALBERLA_LOG_INFO_ON_ROOT( "Production Run!" );
   WALBERLA_LOG_INFO_ON_ROOT( "Creating the block structure: loading from file \'" << sbffile << "\' ..." );
   return std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sbffile.c_str(), true, keepGlobalBlockInformation );
}


shared_ptr<BlockForest> createBlockForestFromConfig(const Config::BlockHandle& mainConf,
                                                    const bool keepGlobalBlockInformation)
{
   bool setupRun                    = mainConf.getParameter< bool >( "setupRun", false );
   Vector3<real_t> simulationCorner = mainConf.getParameter<Vector3<real_t>>("simulationCorner", Vector3<real_t>(0, 0, 0));
   Vector3<real_t> simulationSize   = mainConf.getParameter<Vector3<real_t>>("simulationDomain", Vector3<real_t>(10, 10, 10));
   math::AABB simulationDomain      = math::AABB( simulationCorner, simulationCorner + simulationSize );
   Vector3<uint_t> blocks           = mainConf.getParameter<Vector3<uint_t>>("blocks", Vector3<uint_t>(3, 3, 3));
   Vector3<bool> isPeriodic         = mainConf.getParameter<Vector3<bool>>("isPeriodic", Vector3<bool>(true, true, true));
   uint_t numberOfProcesses         = mainConf.getParameter<uint_t>( "numberOfProcesses",
                                                                  setupRun ? blocks[0] * blocks[1] * blocks[2] : uint_c(mpi::MPIManager::instance()->numProcesses()) );
   uint_t initialRefinementLevel = mainConf.getParameter<uint_t>( "initialRefinementLevel", uint_t(0) );

   if( !mainConf.isDefined( "sbfFile" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "No setup file specified: Creation without setup file!" );
      return createBlockForest(simulationDomain, blocks, isPeriodic, numberOfProcesses, initialRefinementLevel, keepGlobalBlockInformation);
   }

   // sbf file given -> try to load or save domain decomposition
   std::string sbffile = mainConf.getParameter< std::string >( "sbfFile" );
   WALBERLA_LOG_INFO_ON_ROOT( "Setup file specified: Using " << sbffile );
   return createBlockForest(simulationDomain, blocks, isPeriodic, setupRun, sbffile, numberOfProcesses, initialRefinementLevel, keepGlobalBlockInformation);
}



//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   The distribution of blocks to processes also follows a Cartesian decomposition.
*
*   \param domainAABB                 An axis-aligned bounding box that spans the entire simulation space/domain
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param numberOfXProcesses         Number of processes the blocks are distributed to in x direction
*   \param numberOfYProcesses         Number of processes the blocks are distributed to in y direction
*   \param numberOfZProcesses         Number of processes the blocks are distributed to in z direction
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const AABB& domainAABB,
                        const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,         const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock,  const uint_t numberOfZCellsPerBlock,
                        const uint_t numberOfXProcesses,      const uint_t numberOfYProcesses,      const uint_t numberOfZProcesses,
                        const bool   xPeriodic /* = false */, const bool   yPeriodic /* = false */, const bool   zPeriodic /* = false */,
                        const bool keepGlobalBlockInformation /* = false */ )
{
   auto bf = createBlockForest(
            domainAABB,
            numberOfXBlocks,
            numberOfYBlocks,
            numberOfZBlocks,
            numberOfXProcesses,
            numberOfYProcesses,
            numberOfZProcesses,
            xPeriodic,
            yPeriodic,
            zPeriodic,
            keepGlobalBlockInformation);

   auto sbf = std::make_shared< StructuredBlockForest >( bf, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();

   return sbf;
}



//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   The distribution of blocks to processes also follows a Cartesian decomposition.
*
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param dx                         Edge length of each cell (cells are assumed to be cubes)
*   \param numberOfXProcesses         Number of processes the blocks are distributed to in x direction
*   \param numberOfYProcesses         Number of processes the blocks are distributed to in y direction
*   \param numberOfZProcesses         Number of processes the blocks are distributed to in z direction
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,         const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock,  const uint_t numberOfZCellsPerBlock,
                        const real_t dx,
                        const uint_t numberOfXProcesses,      const uint_t numberOfYProcesses,      const uint_t numberOfZProcesses,
                        const bool   xPeriodic /* = false */, const bool   yPeriodic /* = false */, const bool   zPeriodic /* = false */,
                        const bool keepGlobalBlockInformation /* = false */ ) {

   return createUniformBlockGrid( AABB( real_c(0), real_c(0), real_c(0), dx * real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                                         dx * real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                                         dx * real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                                  numberOfXBlocks, numberOfYBlocks, numberOfZBlocks,
                                  numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock,
                                  numberOfXProcesses, numberOfYProcesses, numberOfZProcesses,
                                  xPeriodic, yPeriodic, zPeriodic, keepGlobalBlockInformation );
}



//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   Either all blocks are assigned to the same process (useful for non-parallel simulations) or each blocks is assigned
*   to a different process (useful if only one block shall be assigned to each process).
*
*   \param domainAABB                 An axis-aligned bounding box that spans the entire simulation space/domain
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param oneBlockPerProcess         If true, each block is assigned to a different process. If false, all blocks are
*                                     assigned to the same process (process 0).
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const AABB& domainAABB,
                        const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,         const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock,  const uint_t numberOfZCellsPerBlock,
                        const bool oneBlockPerProcess,
                        const bool   xPeriodic /* = false */, const bool   yPeriodic /* = false */, const bool   zPeriodic /* = false */,
                        const bool keepGlobalBlockInformation /* = false */ ) {

   if( oneBlockPerProcess )
      return createUniformBlockGrid( domainAABB, numberOfXBlocks, numberOfYBlocks, numberOfZBlocks,
                                     numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock,
                                     numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic, keepGlobalBlockInformation );

   // all blocks on the same process
   return createUniformBlockGrid( domainAABB, numberOfXBlocks, numberOfYBlocks, numberOfZBlocks,
                                  numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock, uint_c(1), uint_c(1), uint_c(1),
                                  xPeriodic, yPeriodic, zPeriodic, keepGlobalBlockInformation );
}

//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   Either all blocks are assigned to the same process (useful for non-parallel simulations) or each blocks is assigned
*   to a different process (useful if only one block shall be assigned to each process).
*
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param dx                         Edge length of each cell (cells are assumed to be cubes)
*   \param oneBlockPerProcess         If true, each block is assigned to a different process. If false, all blocks are
*                                     assigned to the same process (process 0).
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,         const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock,  const uint_t numberOfZCellsPerBlock,
                        const real_t dx,
                        const bool oneBlockPerProcess,
                        const bool   xPeriodic /* = false */, const bool   yPeriodic /* = false */, const bool   zPeriodic /* = false */,
                        const bool keepGlobalBlockInformation /* = false */ ) {

   return createUniformBlockGrid( AABB( real_c(0), real_c(0), real_c(0), dx * real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                                         dx * real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                                         dx * real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                                  numberOfXBlocks, numberOfYBlocks, numberOfZBlocks,
                                  numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock,
                                  oneBlockPerProcess, xPeriodic, yPeriodic, zPeriodic, keepGlobalBlockInformation );
}

//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   The number of active MPI processes is used in order to determine the process distribution = in order to perform the
*   initial, static load balancing. Each block is assumed to generate the same amount of work and to require the same
*   amount of memory.
*
*   \param domainAABB                 An axis-aligned bounding box that spans the entire simulation space/domain
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param maxBlocksPerProcess        Maximum number of blocks that are allowed to be assigned to one process. If a
*                                     value of '0' is provided, any number of blocks are allowed to be located on one
*                                     process - meaning static load balancing doesn't try to obey any memory limit. ['0' by default]
*   \param includeMetis               If true (and if available!), METIS is also used during load balancing. [true by default]
*   \param forceMetis                 If true, METIS is always preferred over space filling curves [false by default]
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const AABB& domainAABB,
                        const uint_t numberOfXBlocks,             const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,      const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const uint_t maxBlocksPerProcess /*= 0*/, const bool includeMetis /*= true*/,  const bool forceMetis /*= false*/,
                        const bool   xPeriodic /*= false*/,       const bool   yPeriodic /*= false*/,  const bool   zPeriodic /*= false*/,
                        const bool keepGlobalBlockInformation /*= false*/ ) {

   // initialize SetupBlockForest = determine domain decomposition

   SetupBlockForest sforest;

   sforest.addWorkloadMemorySUIDAssignmentFunction( uniformWorkloadAndMemoryAssignment );

   sforest.init( domainAABB, numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic );

   // calculate process distribution

   const memory_t memoryLimit = ( maxBlocksPerProcess == 0 ) ? numeric_cast< memory_t >( sforest.getNumberOfBlocks() ) :
                                                               numeric_cast< memory_t >( maxBlocksPerProcess );

   GlobalLoadBalancing::MetisConfiguration< SetupBlock > metisConfig( includeMetis, forceMetis,
                                                                      std::bind( cellWeightedCommunicationCost, std::placeholders::_1, std::placeholders::_2,
                                                                                   numberOfXCellsPerBlock,
                                                                                   numberOfYCellsPerBlock,
                                                                                   numberOfZCellsPerBlock ) );

   sforest.calculateProcessDistribution_Default( uint_c( MPIManager::instance()->numProcesses() ), memoryLimit, "hilbert", 10, false, metisConfig );

   if( !MPIManager::instance()->rankValid() )
      MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)

   auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation );

   auto sbf = std::make_shared< StructuredBlockForest >( bf, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();

   return sbf;
}



//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   The number of active MPI processes is used in order to determine the process distribution = in order to perform the
*   initial, static load balancing. Each block is assumed to generate the same amount of work and to require the same
*   amount of memory.
*
*   \param numberOfXBlocks            Number of blocks in x direction
*   \param numberOfYBlocks            Number of blocks in y direction
*   \param numberOfZBlocks            Number of blocks in z direction
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param dx                         Edge length of each cell (cells are assumed to be cubes)
*   \param maxBlocksPerProcess        Maximum number of blocks that are allowed to be assigned to one process. If a
*                                     value of '0' is provided, any number of blocks are allowed to be located on one
*                                     process - meaning static load balancing doesn't try to obey any memory limit. ['0' by default]
*   \param includeMetis               If true (and if available!), METIS is also used during load balancing. [true by default]
*   \param forceMetis                 If true, METIS is always preferred over space filling curves [false by default]
*   \param xPeriodic                  If true, the block structure is periodic in x direction [false by default]
*   \param yPeriodic                  If true, the block structure is periodic in y direction [false by default]
*   \param zPeriodic                  If true, the block structure is periodic in z direction [false by default]
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const uint_t numberOfXBlocks,             const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,      const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const real_t dx,
                        const uint_t maxBlocksPerProcess /*= 0*/, const bool includeMetis /*= true*/,  const bool forceMetis /*= false*/,
                        const bool   xPeriodic /*= false*/,       const bool   yPeriodic /*= false*/,  const bool   zPeriodic /*= false*/,
                        const bool keepGlobalBlockInformation /*= false*/ ) {

   return createUniformBlockGrid( AABB( real_c(0), real_c(0), real_c(0), dx * real_c( numberOfXBlocks * numberOfXCellsPerBlock ),
                                                                         dx * real_c( numberOfYBlocks * numberOfYCellsPerBlock ),
                                                                         dx * real_c( numberOfZBlocks * numberOfZCellsPerBlock ) ),
                                  numberOfXBlocks, numberOfYBlocks, numberOfZBlocks,
                                  numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock,
                                  maxBlocksPerProcess, includeMetis, forceMetis, xPeriodic, yPeriodic, zPeriodic, keepGlobalBlockInformation );
}



//**********************************************************************************************************************
/*!
*   \brief Function for creating a structured block forest that represents a uniform block grid.
*
*   Uniform block grid: Cartesian domain decomposition into blocks of cells, each block has the same size and contains
*                       the same number of cells.
*   The entire block structure and its corresponding process distribution are loaded from file.
*
*   \param filename                   A file that stores a block structure and its corresponding process distribution
*   \param numberOfXCellsPerBlock     Number of cells of each block in x direction
*   \param numberOfYCellsPerBlock     Number of cells of each block in y direction
*   \param numberOfZCellsPerBlock     Number of cells of each block in z direction
*   \param keepGlobalBlockInformation If true, each process keeps information about remote blocks (blocks that reside
*                                     on other processes). This information includes the process rank, the state, and
*                                     the axis-aligned bounding box of any block (local or remote). [false by default]
*/
//**********************************************************************************************************************

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const std::string& filename,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const bool keepGlobalBlockInformation /*= false*/ )
{

   if( !MPIManager::instance()->rankValid() )
      MPIManager::instance()->useWorldComm();

   auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), filename.c_str(), true, keepGlobalBlockInformation );

   if( !bf->storesUniformBlockGrid() )
      WALBERLA_ABORT( "The block forest loaded from file \'" << filename << "\' does not contain a uniform block grid!" );

   auto sbf = std::make_shared< StructuredBlockForest >( bf, numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );
   sbf->createCellBoundingBoxes();

   return sbf;
}



///////////////////////////////////////////
// HELPER FUNCTIONS                      //
///////////////////////////////////////////



//*******************************************************************************************************************
/*! Tries to distribute a given amount of total cells to a given amount of blocks
*
* It may happen that divisibility of the nr of cells requested prevents a distribution
* in this case the number of cells is chosen bigger than requested
*
*
* \param cells              total number of cells requested
* \param nrOfBlocks         total number of blocks to distribute the cells to
* \param[out] blocksOut     calculated number of blocks in x/y/z
* \param[out] cellsPerBlock how many cells to put on each block
*                           it may happen that divisibility of the number of cells requested prevents a distribution
*                           in this case the number of cells is chosen (slightly) bigger than requested
*
* Example: in:  cells = (10,15,16)
*          in:  blocks = 8
*          out: blocks = (2,2,2)
*          out: cellsPerBlock = (5,8,8)
*          out: newCells = (10,16,16)
*/
//*******************************************************************************************************************
void calculateCellDistribution( const Vector3<uint_t> & cells, uint_t nrOfBlocks,
                                Vector3<uint_t> & blocksOut, Vector3<uint_t> & cellsPerBlock)
{
   std::vector< real_t > weighting;
   weighting.push_back( real_c( cells[0]) );
   weighting.push_back( real_c( cells[1]) );
   weighting.push_back( real_c( cells[2]) );
   std::vector<uint_t> blocks = math::getFactors( nrOfBlocks, 3, weighting );

   for( uint_t i = 0; i < 3; ++i )
   {
      if ( uint_c( cells[i] ) % blocks[i] == 0 )
         cellsPerBlock[i] = cells[i] /  blocks[i];
      else // extend the domain if processesCount does not divide the cell count in this direction
         cellsPerBlock[i] = ( cells[i] +   blocks[i] ) /  blocks[i];
   }
   for( uint_t i = 0; i < 3; ++i )
      blocksOut[i] = blocks[i];
}


void uniformWorkloadAndMemoryAssignment( SetupBlockForest& forest ) {

   std::vector< SetupBlock* > blocks;
   forest.getBlocks( blocks );

   for( uint_t i = 0; i != blocks.size(); ++i ) {
      blocks[i]->setWorkload( numeric_cast< workload_t >(1) );
      blocks[i]->setMemory( numeric_cast< memory_t >(1) );
   }
}


memory_t cellWeightedCommunicationCost( const SetupBlock* const a, const SetupBlock* const b,
                                        uint_t xCellsPerBlock, uint_t yCellsPerBlock, uint_t zCellsPerBlock )
{
   for ( auto dIter = stencil::D3Q19::beginNoCenter(); dIter != stencil::D3Q19::end(); ++dIter )
   {
      auto neighborHoodIdxA = getBlockNeighborhoodSectionIndex( *dIter );
      for ( uint_t j = 0; j != a->getNeighborhoodSectionSize( neighborHoodIdxA ); ++j )
      {
         if( a->getNeighbor(neighborHoodIdxA,j) == b )
         {
            switch ( *dIter )
            {
               //faces
               case stencil::W: return memory_c( yCellsPerBlock * zCellsPerBlock );
               case stencil::E: return memory_c( yCellsPerBlock * zCellsPerBlock );
               case stencil::N: return memory_c( xCellsPerBlock * zCellsPerBlock );
               case stencil::S: return memory_c( xCellsPerBlock * zCellsPerBlock );
               case stencil::T: return memory_c( xCellsPerBlock * yCellsPerBlock );
               case stencil::B: return memory_c( xCellsPerBlock * yCellsPerBlock );
               //edges
               case stencil::NW: return memory_c( zCellsPerBlock );
               case stencil::NE: return memory_c( zCellsPerBlock );
               case stencil::SW: return memory_c( zCellsPerBlock );
               case stencil::SE: return memory_c( zCellsPerBlock );
               case stencil::TN: return memory_c( xCellsPerBlock );
               case stencil::TS: return memory_c( xCellsPerBlock );
               case stencil::TW: return memory_c( yCellsPerBlock );
               case stencil::TE: return memory_c( yCellsPerBlock );
               case stencil::BN: return memory_c( xCellsPerBlock );
               case stencil::BS: return memory_c( xCellsPerBlock );
               case stencil::BW: return memory_c( yCellsPerBlock );
               case stencil::BE: return memory_c( yCellsPerBlock );
               default:
                  WALBERLA_ABORT( "Unknown direction. Should not happen!" )
            }
         }
      }
   }

   // Return 1 for corners
   return numeric_cast< memory_t >(1);

}


memory_t uniformFacesDominantCommunication( const SetupBlock* const a, const SetupBlock* const b ) {

   uint_t faces[] = { 4, 10, 12, 13, 15, 21 };

   for( uint_t i = 0; i != 6; ++i ) {
      for( uint_t j = 0; j != a->getNeighborhoodSectionSize(faces[i]); ++j )
         if( a->getNeighbor(faces[i],j) == b )
            return numeric_cast< memory_t >(1000);
   }

   return numeric_cast< memory_t >(1);
}



} // namespace blockforest
} // namespace walberla
