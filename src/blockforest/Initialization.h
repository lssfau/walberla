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
//! \file Initialization.h
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StructuredBlockForest.h"

#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"



namespace walberla {
namespace blockforest {



shared_ptr< StructuredBlockForest >  createUniformBlockGridFromConfig( const shared_ptr< Config > & config,
                                                                       CellInterval * requestedDomainSize = NULL,
                                                                       const bool keepGlobalBlockInformation = false );

shared_ptr< StructuredBlockForest >  createUniformBlockGridFromConfig( const Config::BlockHandle & configBlock,
                                                                       CellInterval * requestedDomainSize = NULL,
                                                                       const bool keepGlobalBlockInformation = false );



shared_ptr< BlockForest >
createBlockForest(      const AABB& domainAABB,
                        const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXProcesses,     const uint_t numberOfYProcesses,     const uint_t numberOfZProcesses,
                        const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );
shared_ptr<BlockForest>
createBlockForest(const math::AABB& simulationDomain,
                  const Vector3<uint_t>& blocks,
                  const Vector3<bool>& isPeriodic,
                  const uint_t numberOfProcesses = uint_c(mpi::MPIManager::instance()->numProcesses()),
                  const uint_t initialRefinementLevel = uint_t(0),
                  const bool keepGlobalBlockInformation = false);
shared_ptr<BlockForest>
createBlockForest(const math::AABB& simulationDomain,
                  const Vector3<uint_t>& blocks,
                  const Vector3<bool>& isPeriodic,
                  const bool setupRun,
                  const std::string& sbffile,
                  const uint_t numberOfProcesses = uint_c(mpi::MPIManager::instance()->numProcesses()),
                  const uint_t initialRefinementLevel = uint_t(0),
                  const bool keepGlobalBlockInformation = false);
shared_ptr<BlockForest>
createBlockForestFromConfig(const Config::BlockHandle& mainConf,
                            const bool keepGlobalBlockInformation = false);

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const AABB& domainAABB,
                        const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const uint_t numberOfXProcesses,     const uint_t numberOfYProcesses,     const uint_t numberOfZProcesses,
                        const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const real_t dx,
                        const uint_t numberOfXProcesses,     const uint_t numberOfYProcesses,     const uint_t numberOfZProcesses,
                        const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );



shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const AABB& domainAABB,
                        const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const bool oneBlockPerProcess,
                        const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const uint_t numberOfXBlocks,        const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const real_t dx,
                        const bool oneBlockPerProcess,
                        const bool   xPeriodic = false,      const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );



shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const AABB& domainAABB,
                        const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const uint_t maxBlocksPerProcess = 0, const bool includeMetis = true,      const bool forceMetis = false,
                        const bool   xPeriodic = false,       const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );

shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const uint_t numberOfXBlocks,         const uint_t numberOfYBlocks,        const uint_t numberOfZBlocks,
                        const uint_t numberOfXCellsPerBlock,  const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const real_t dx,
                        const uint_t maxBlocksPerProcess = 0, const bool includeMetis = true,      const bool forceMetis = false,
                        const bool   xPeriodic = false,       const bool   yPeriodic = false,      const bool   zPeriodic = false,
                        const bool keepGlobalBlockInformation = false );



shared_ptr< StructuredBlockForest >
createUniformBlockGrid( const std::string& filename,
                        const uint_t numberOfXCellsPerBlock, const uint_t numberOfYCellsPerBlock, const uint_t numberOfZCellsPerBlock,
                        const bool keepGlobalBlockInformation = false );





void calculateCellDistribution( const Vector3<uint_t> & cells, uint_t nrOfBlocks,
                                Vector3<uint_t> & blocksOut, Vector3<uint_t> & cellsPerBlock);

void     uniformWorkloadAndMemoryAssignment( SetupBlockForest& forest );
memory_t uniformFacesDominantCommunication ( const SetupBlock* const a, const SetupBlock* const b );
memory_t cellWeightedCommunicationCost( const SetupBlock* const a, const SetupBlock* const b,
                                        uint_t xCellsPerBlock, uint_t yCellsPerBlock, uint_t zCellsPerBlock );


} // namespace blockforest
} // namespace walberla


