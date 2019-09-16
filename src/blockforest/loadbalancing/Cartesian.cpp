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
//! \file Cartesian.cpp
//! \ingroup blockforest
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "Cartesian.h"

#include "core/debug/CheckFunctions.h"



namespace walberla {
namespace blockforest {



uint_t CartesianDistribution::operator()( SetupBlockForest & forest, const uint_t numberOfProcesses, const memory_t /*perProcessMemoryLimit*/ )
{
   if( numberOfProcesses != ( numberOfXProcesses_ * numberOfYProcesses_ * numberOfZProcesses_ ) )
      WALBERLA_ABORT( "Load balancing failed: The total number of processes must be identical to the product "
                      "of the \'number of processes in x-, y-, and z-direction\'." );

   if( numberOfXProcesses_ > forest.getXSize() )
      WALBERLA_ABORT( "Load balancing failed: \'Number of processes in x-direction\' must be in (0," << forest.getXSize() << "]. "
                      "You specified \'" << numberOfXProcesses_ << "\'." );

   if( numberOfYProcesses_ > forest.getYSize() )
      WALBERLA_ABORT( "Load balancing failed: \'Number of processes in y-direction\' must be in (0," << forest.getYSize() << "]. "
                      "You specified \'" << numberOfYProcesses_ << "\'." );

   if( numberOfZProcesses_ > forest.getZSize() )
      WALBERLA_ABORT( "Load balancing failed: \'Number of processes in z-direction\' must be in (0," << forest.getZSize() << "]. "
                      "You specified \'" << numberOfZProcesses_ << "\'." );

   if( !processIdMap_->empty() )
      WALBERLA_CHECK_EQUAL( processIdMap_->size(), numberOfProcesses );

   uint_t partitions[3];
   partitions[0] = numberOfXProcesses_;
   partitions[1] = numberOfYProcesses_;
   partitions[2] = numberOfZProcesses_;

   std::vector< uint_t > indices[3];

   for( uint_t i = 0; i != 3; ++i )
   {
      const uint_t div = forest.getSize(i) / partitions[i];
      const uint_t mod = forest.getSize(i) % partitions[i];

      indices[i].resize( partitions[i] + 1, div );
      indices[i][0] = 0;

      for( uint_t j = 0; j != mod; ++j )
         ++indices[i][j+1];
      for( uint_t j = 1; j != indices[i].size(); ++j )
         indices[i][j] += indices[i][j-1];
   }

   for( uint_t z = 0; z != partitions[2]; ++z ) {
      for( uint_t y = 0; y != partitions[1]; ++y ) {
         for( uint_t x = 0; x != partitions[0]; ++x )
         {
            std::vector< SetupBlock * > partitionBlocks;

            forest.getBlocks( partitionBlocks, indices[0][x], indices[1][y], indices[2][z], indices[0][x+1], indices[1][y+1], indices[2][z+1] );

            for( auto block = partitionBlocks.begin(); block != partitionBlocks.end(); ++block )
            {
               const uint_t index = z * partitions[0] * partitions[1] + y * partitions[0] + x;

               (*block)->assignTargetProcess( ( !processIdMap_->empty() ) ? (*processIdMap_)[ index ] : index );
            }
         }
      }
   }

   return numberOfProcesses;
}



} // namespace blockforest
} // namespace walberla
