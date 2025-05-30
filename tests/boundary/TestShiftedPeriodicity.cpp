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
//! \file TestShiftedPeriodicity.cpp
//! \ingroup boundary
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "blockforest/StructuredBlockForest.h"

#include "core/cell/Cell.h"
#include "core/cell/CellInterval.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Initialization.h"
#include "core/math/Vector3.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Operation.h"
#include "core/mpi/Reduce.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"

#include "field/Layout.h"
//#include "field/vtk/VTKWriter.h"

#include "python_coupling/CreateConfig.h"

#include "stencil/D3Q27.h"
#include "stencil/Directions.h"

#include <blockforest/Initialization.h>
#include <boundary/ShiftedPeriodicity.h>
#include <core/DataTypes.h>
#include <core/debug/Debug.h>
#include <core/debug/TestSubsystem.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <memory>
#include <vector>

namespace walberla {

using Stencil_T = stencil::D3Q27;

using ValueType_T = real_t;
using Field_T = GhostLayerField< ValueType_T, 3 >;

constexpr cell_idx_t fieldGhostLayers = 1;

//////////
// MAIN //
//////////

template< typename FieldType_T >
class FieldInitialiser {

 public:
   FieldInitialiser( const std::weak_ptr< StructuredBlockForest > & sbf, const BlockDataID fieldId )
   : sbf_(sbf), fieldId_(fieldId)
   {}

   void operator()() {


      const auto blocks = sbf_.lock();
      WALBERLA_ASSERT_NOT_NULLPTR(blocks)

      for ( auto & block : *blocks )
      {
         // get field
         auto * field = block.template getData<FieldType_T>(fieldId_);
         WALBERLA_ASSERT_NOT_NULLPTR(field)

         // get cell interval
         auto ci = field->xyzSizeWithGhostLayer();

         for (const auto& cell : ci) {

            // get global coordinates
            Cell globalCell;
            blocks->transformBlockLocalToGlobalCell(globalCell, block, cell);

            for (uint_t d = 0; d < FieldType_T::F_SIZE; ++d)
            {
               field->get(cell, d) = real_c(  (globalCell.x() + 2 * fieldGhostLayers)
                                            * (globalCell.y() + 2 * fieldGhostLayers)
                                            * (globalCell.z() + 2 * fieldGhostLayers)
                                            * int_c(d + 1));
            }
         }
      }

   }

 private:
   std::weak_ptr< StructuredBlockForest > sbf_{};
   const BlockDataID fieldId_{};

};


int main( int argc, char **argv ) {

   const mpi::Environment env(argc, argv);

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      auto config = *cfg;
      logging::configureLogging(config);

      // create the domain, flag field and vector field (non-uniform initialisation)
      auto configBlock = config->getBlock("DomainSetup");

      const Vector3<bool> periodicity = configBlock.getParameter<Vector3<bool>  >( "periodic" );
      const Vector3<uint_t> cellsPerBlock = configBlock.getParameter<Vector3<uint_t>  >( "cellsPerBlock" );
      const Vector3<uint_t> nBlocks        = configBlock.getParameter<Vector3<uint_t>  >( "blocks" );
      const bool zeroCenteredDomain = configBlock.getParameter<bool>( "zeroCenteredDomain" );

      Vector3<real_t> domainSize{};
      Vector3<real_t> domainCorner{};
      for(uint_t d = 0; d < 3; ++d) {
         domainSize[d] = real_c(nBlocks[d] * cellsPerBlock[d]);
         domainCorner[d] = zeroCenteredDomain ? real_t(0) : - domainSize[d] / real_t(2);
      }

      auto blocks = blockforest::createUniformBlockGrid(
         AABB(domainCorner, domainCorner + domainSize),         // AABB spanning the entire domain
         nBlocks[0], nBlocks[1], nBlocks[2],                    // number of blocks in x/y/z direction
         cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2],  // cells per block in x/y/z direction
         uint_t(0),                                             // maximum number of blocks per process
         true, false,                                           // include but don't force Metis
         periodicity[0], periodicity[1], periodicity[2],        // periodicity
         true                                                   // keep global block information
      );

      const auto fieldID = field::addToStorage< Field_T >(blocks, "test field", real_t(0), field::Layout::fzyx, fieldGhostLayers);
      FieldInitialiser< Field_T > initialiser(blocks, fieldID);

      // re-initialise fields
      initialiser();

      // create periodic shift boundary condition

      const auto spConfig = config->getBlock("Boundaries").getBlock("ShiftedPeriodicity");
      const uint_t shiftDir = spConfig.getParameter<uint_t>("shiftDir");
      const int shiftValue = spConfig.getParameter<int>("shiftValue");
      const uint_t normalDir = spConfig.getParameter<uint_t>("normalDir");

      boundary::ShiftedPeriodicity<Field_T> shiftedPeriodicity(
         blocks, fieldID, fieldGhostLayers, normalDir, shiftDir, shiftValue
      );

//      auto vtkOutput = field::createVTKOutput<Field_T>(fieldID, *blocks, "test_field", 1, fieldGhostLayers,
//                                                         false, "vtk_out", "simulation_step", false, false);
//      vtkOutput();

      // apply boundary condition and standard communication
      shiftedPeriodicity();

//      vtkOutput();

      /// compare values
      // create local domain slices to compare values before and after shift

      const uint_t remDir = 3 - shiftDir - normalDir;
      const auto shift = shiftedPeriodicity.shift();

      const uint_t shiftSize = blocks->getNumberOfCells(shiftDir);
      const uint_t remSize   = blocks->getNumberOfCells(remDir);
      const uint_t sizeSlice = shiftSize * remSize * Field_T::F_SIZE;

      std::vector<ValueType_T> innerMin(sizeSlice, ValueType_T(0));
      std::vector<ValueType_T> innerMax(sizeSlice, ValueType_T(0));
      std::vector<ValueType_T> glMin(sizeSlice, ValueType_T(0));
      std::vector<ValueType_T> glMax(sizeSlice, ValueType_T(0));

      auto getIdx = [&remSize](const cell_idx_t shiftIdx, const cell_idx_t remIdx){
         return shiftIdx * cell_idx_c(Field_T::F_SIZE * remSize)
                + remIdx * cell_idx_c(Field_T::F_SIZE);
      };

      // fill slices for comparing values
      for(auto & block : *blocks) {
         const bool atMin = blocks->atDomainMinBorder(normalDir, block);
         const bool atMax = blocks->atDomainMaxBorder(normalDir, block);
         if(!atMin && !atMax)
            continue;

         auto * field = block.getData<Field_T>(fieldID);
         WALBERLA_ASSERT_NOT_NULLPTR(field)

         // fill innerMin and glMin
         if(atMin) {
            Vector3<int> dirVector{};
            dirVector[normalDir] = -1;
            const auto dir = stencil::vectorToDirection(dirVector);

            CellInterval innerMinCI;
            CellInterval glMinCI;
            field->getSliceBeforeGhostLayer(dir, innerMinCI, fieldGhostLayers, false);
            field->getGhostRegion(dir, glMinCI, fieldGhostLayers, false);

            // fill inner min
            for(const auto & cell : innerMinCI){
               Cell globalCell;
               blocks->transformBlockLocalToGlobalCell(globalCell, block, cell);

               cell_idx_t idxShiftDir = globalCell[shiftDir] - shift[shiftDir];
               if(idxShiftDir >= cell_idx_c(shiftSize)) idxShiftDir -= cell_idx_c(shiftSize);
               if(idxShiftDir <= - int_c(fieldGhostLayers)) idxShiftDir += cell_idx_c(shiftSize);

               const cell_idx_t idx = getIdx(idxShiftDir, cell_idx_c(globalCell[remDir]));

               for(cell_idx_t f = 0; f < cell_idx_c(Field_T::F_SIZE); ++f) {
                  WALBERLA_ASSERT(field->coordinatesValid(cell.x(), cell.y(), cell.z(), f))
                  WALBERLA_ASSERT_LESS(idx + f, innerMin.size())

                  innerMin[uint_c(idx + f)] = field->get(cell, f);
               }
            }

            // fill gl min
            for(const auto & cell : glMinCI){
               Cell globalCell;
               blocks->transformBlockLocalToGlobalCell(globalCell, block, cell);

               const cell_idx_t idx = getIdx(cell_idx_c(globalCell[shiftDir]), globalCell[remDir]);

               for(cell_idx_t f = 0; f < cell_idx_c(Field_T::F_SIZE); ++f) {
                  WALBERLA_ASSERT(field->coordinatesValid(cell.x(), cell.y(), cell.z(), f))
                  WALBERLA_ASSERT_LESS(idx + f, glMin.size())

                  glMin[uint_c(idx + f)] = field->get(cell, f);
               }
            }
         }

         // fill innerMax and glMax
         if(atMax) {
            Vector3<int> dirVector{};
            dirVector[normalDir] = 1;
            const auto dir = stencil::vectorToDirection(dirVector);

            CellInterval innerMaxCI;
            CellInterval glMaxCI;
            field->getSliceBeforeGhostLayer(dir, innerMaxCI, fieldGhostLayers, false);
            field->getGhostRegion(dir, glMaxCI, fieldGhostLayers, false);

            // fill inner max
            for(const auto & cell : innerMaxCI){
               Cell globalCell;
               blocks->transformBlockLocalToGlobalCell(globalCell, block, cell);

               cell_idx_t idxShiftDir = globalCell[shiftDir] + shift[shiftDir];
               if(idxShiftDir >= cell_idx_c(shiftSize)) idxShiftDir -= cell_idx_c(shiftSize);
               if(idxShiftDir <= - int_c(fieldGhostLayers)) idxShiftDir += cell_idx_c(shiftSize);

               const cell_idx_t idx = getIdx(idxShiftDir, cell_idx_c(globalCell[remDir]));

               for(cell_idx_t f = 0; f < cell_idx_c(Field_T::F_SIZE); ++f) {
                  WALBERLA_ASSERT(field->coordinatesValid(cell.x(), cell.y(), cell.z(), f))
                  WALBERLA_ASSERT_LESS(idx + f, innerMax.size())

                  innerMax[uint_c(idx + f)] = field->get(cell, f);
               }
            }

            // fill gl max
            for(const auto & cell : glMaxCI){
               Cell globalCell;
               blocks->transformBlockLocalToGlobalCell(globalCell, block, cell);

               const cell_idx_t idx = getIdx(cell_idx_c(globalCell[shiftDir]), cell_idx_c(globalCell[remDir]));

               for(cell_idx_t f = 0; f < cell_idx_c(Field_T::F_SIZE); ++f) {
                  WALBERLA_ASSERT(field->coordinatesValid(cell.x(), cell.y(), cell.z(), f))
                  WALBERLA_ASSERT_LESS(idx + f, glMax.size())

                  glMax[uint_c(idx + f)] = field->get(cell, f);
               }
            }
         }


      }

      WALBERLA_MPI_SECTION() {

         mpi::reduceInplace(innerMin, mpi::SUM);
         mpi::reduceInplace(innerMax, mpi::SUM);
         mpi::reduceInplace(glMin, mpi::SUM);
         mpi::reduceInplace(glMax, mpi::SUM);

      }

      WALBERLA_ROOT_SECTION() {
         for(uint_t i = 0; i < sizeSlice; ++i) {
            WALBERLA_CHECK_FLOAT_EQUAL(innerMin[i], glMax[i])
            WALBERLA_CHECK_FLOAT_EQUAL(innerMax[i], glMin[i])
         }
      }

   }

   return 0;
}

} // namespace walberla



int main(int argc, char **argv) {

   walberla::debug::enterTestMode();

   return walberla::main(argc,argv);
}
