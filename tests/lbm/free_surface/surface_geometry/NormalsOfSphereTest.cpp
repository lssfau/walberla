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
//! \file NormalsOfSphereTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Initialize spherical bubble and compare calculated normals to analytical (radial) normals.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/DistributedSample.h"

#include "field/AddToStorage.h"

#include "geometry/bodies/Sphere.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>

namespace walberla
{
namespace free_surface
{
namespace NormalsOfSphereTest
{
// define types
using flag_t = uint32_t;

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

// compute and evaluate the absolute error of the angle of the interface normal in each interface cell;
// reference: vector that points from the gas bubble's, i.e., sphere's center to the current interface cell's center
template< typename Stencil_T >
class EvaluateNormalError
{
 public:
   EvaluateNormalError(const std::weak_ptr< const StructuredBlockForest >& blockForest, const geometry::Sphere& sphere,
                       const ConstBlockDataID& normalFieldID, const ConstBlockDataID& flagFieldID,
                       const field::FlagUID& interfaceFlagID, bool computeNormalsInInterfaceNeighbors,
                       bool smoothFillLevel)
      : blockForest_(blockForest), sphere_(sphere), normalFieldID_(normalFieldID), flagFieldID_(flagFieldID),
        interfaceFlagID_(interfaceFlagID), computeNormalsInInterfaceNeighbors_(computeNormalsInInterfaceNeighbors),
        smoothFillLevel_(smoothFillLevel)
   {}

   void operator()(IBlock* const block)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // get fields
      const FlagField_T* const flagField     = block->getData< const FlagField_T >(flagFieldID_);
      const VectorField_T* const normalField = block->getData< const VectorField_T >(normalFieldID_);

      const flag_t interfaceFlag = flagField->getFlag(interfaceFlagID_);

      math::DistributedSample errorSample;

      WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, normalFieldIt, normalField, omp critical, {
         // evaluate normal only in interface cell (and possibly in an interface cell's neighborhood)
         if (isFlagSet(flagFieldIt, interfaceFlag) ||
             (computeNormalsInInterfaceNeighbors_ && isFlagInNeighborhood< Stencil_T >(flagFieldIt, interfaceFlag)))
         {
            // get a vector pointing from sphere's center to current cell's center
            Vector3< real_t > r;

            // get the global coordinate of the current cell's center
            blockForest->getBlockLocalCellCenter(*block, flagFieldIt.cell(), r);
            r -= sphere_.midpoint();

            // invert normal direction: in this FSLBM implementation, the normal vector is defined to point from liquid
            // to gas
            r *= real_c(-1);
            normalize(r);

            // calculate the error in the angle of the interface normal
            // the domain of arccosine is [-1,1]; due to numerical inaccuracies, the dot product "*normalFieldIt*r"
            // might be slightly out of this range; these cases are ignored here
            const real_t dotProduct = *normalFieldIt * r;
            if (dotProduct >= real_c(-1) && dotProduct <= real_c(1))
            {
               const real_t angleDiff = std::acos(dotProduct);
               errorSample.insert(angleDiff);
            }
         }
      }) // WALBERLA_FOR_ALL_CELLS_OMP

      errorSample.mpiAllGather();

      const real_t maxError  = errorSample.max();
      const real_t meanError = errorSample.mean();

      WALBERLA_LOG_RESULT("Mean absolute error in angle of normal = " << meanError);
      WALBERLA_LOG_RESULT("Maximum absolute error in angle of normal = " << maxError);
      WALBERLA_LOG_RESULT("Minimum absolute error in angle of normal = " << errorSample.min());

      // the following reference errors have been obtained with a version of the code that is believed to be correct
      if (computeNormalsInInterfaceNeighbors_ && smoothFillLevel_)
      {
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q27 >)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.036));
            WALBERLA_CHECK_LESS(maxError, real_c(0.2));
         }
         else
         {
            if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >)
            {
               WALBERLA_CHECK_LESS(meanError, real_c(0.036));
               WALBERLA_CHECK_LESS(maxError, real_c(0.19));
            }
         }
      }

      if (!computeNormalsInInterfaceNeighbors_ && smoothFillLevel_)
      {
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q27 >)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.016));
            WALBERLA_CHECK_LESS(maxError, real_c(0.045));
         }
         else
         {
            if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >)
            {
               WALBERLA_CHECK_LESS(meanError, real_c(0.0147));
               WALBERLA_CHECK_LESS(maxError, real_c(0.0457));
            }
         }
      }

      if (!computeNormalsInInterfaceNeighbors_ && !smoothFillLevel_)
      {
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q27 >)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.021));
            WALBERLA_CHECK_LESS(maxError, real_c(0.08));
         }
         else
         {
            if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >)
            {
               WALBERLA_CHECK_LESS(meanError, real_c(0.037));
               WALBERLA_CHECK_LESS(maxError, real_c(0.096));
            }
         }
      }

      if (computeNormalsInInterfaceNeighbors_ && !smoothFillLevel_)
      {
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q27 >)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.139));
            WALBERLA_CHECK_LESS(maxError, real_c(1.58));
         }
         else
         {
            if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >)
            {
               WALBERLA_CHECK_LESS(meanError, real_c(0.121));
               WALBERLA_CHECK_LESS(maxError, real_c(1.58));
            }
         }
      }
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   geometry::Sphere sphere_;

   ConstBlockDataID normalFieldID_;
   ConstBlockDataID flagFieldID_;
   field::FlagUID interfaceFlagID_;

   bool computeNormalsInInterfaceNeighbors_;
   bool smoothFillLevel_;
}; // class EvaluateNormalError

template< typename LatticeModel_T >
void test(uint_t numCells, Vector3< real_t > offset, bool computeNormalsInInterfaceNeighbors, bool smoothFillLevel)
{
   // define types
   using Stencil_T                     = typename LatticeModel_T::Stencil;
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;
   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::CommunicationStencil >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(numCells + uint_c(6));
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          false, false, false);                                 // periodicity

   // create lattice model with omega=1
   LatticeModel_T latticeModel(real_c(1.0));

   // add fields
   BlockDataID pdfFieldID =
      lbm::addPdfFieldToStorage< LatticeModel_T >(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(1.0), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID smoothFillFieldID;

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();

   // add gas sphere, i.e., bubble
   Vector3< real_t > midpoint(real_c(domainSize[0]) * real_c(0.5), real_c(domainSize[1]) * real_c(0.5),
                              real_c(domainSize[2]) * real_c(0.5));
   geometry::Sphere sphere(midpoint + offset, real_c(numCells) * real_c(0.5));
   freeSurfaceBoundaryHandling->addFreeSurfaceObject(sphere);
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   BlockDataID* relevantFillFieldID = &fillFieldID;

   if (smoothFillLevel)
   {
      smoothFillFieldID =
         field::addToStorage< ScalarField_T >(blockForest, "Smooth fill levels", real_c(1.0), field::fzyx, uint_c(1));

      // add sweep for smoothing the fill level field
      SmoothingSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > smoothingSweep(
         smoothFillFieldID, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs,
         freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), false);
      timeloop.add() << Sweep(smoothingSweep, "Smoothing sweep")
                     << AfterFunction(Communication_T(blockForest, smoothFillFieldID),
                                      "Communication after smoothing sweep");

      relevantFillFieldID = &smoothFillFieldID;
   }

   // add sweep for computing interface normals
   NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalsSweep(
      normalFieldID, *relevantFillFieldID, flagFieldID, flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs,
      freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), computeNormalsInInterfaceNeighbors, false, false,
      false);
   timeloop.add() << Sweep(normalsSweep, "Normals sweep")
                  << AfterFunction(Communication_T(blockForest, normalFieldID), "Communication after normal sweep");

   // add sweep for evaluating the error of the interface normals (by comparing with analytical normal)
   EvaluateNormalError< Stencil_T > errorEvaluationSweep(blockForest, sphere, normalFieldID, flagFieldID,
                                                         flagIDs::interfaceFlagID, computeNormalsInInterfaceNeighbors,
                                                         smoothFillLevel);
   timeloop.add() << Sweep(errorEvaluationSweep, "Error evaluation sweep");

   // perform a single time step
   timeloop.singleStep();

   MPIManager::instance()->resetMPI();
}

template< typename LatticeModel_T >
void runAllTests()
{
   // test with various bubble diameters
   for (uint_t i = uint_c(10); i <= uint_c(30); i += uint_c(10))
   {
      WALBERLA_LOG_RESULT("Bubble diameter " << i << " cells; normals computed only in interface cells;");
      test< LatticeModel_T >(i, Vector3< real_t >(real_c(0.5), real_c(0), real_c(0)), false, false);

      WALBERLA_LOG_RESULT("Bubble diameter "
                          << i << " cells; normals computed only in interface cells; fill level field smoothed;");
      test< LatticeModel_T >(i, Vector3< real_t >(real_c(0.5), real_c(0), real_c(0)), false, true);

      WALBERLA_LOG_RESULT("Bubble diameter " << i
                                             << " cells; normals computed in D3Q27 neighborhood of interface cells;");
      test< LatticeModel_T >(i, Vector3< real_t >(real_c(0.5), real_c(0), real_c(0)), true, false);

      WALBERLA_LOG_RESULT(
         "Bubble diameter "
         << i << " cells; normals computed in D3Q27 neighborhood of interface cells; fill level field smoothed;");
      test< LatticeModel_T >(i, Vector3< real_t >(real_c(0.5), real_c(0), real_c(0)), true, true);
   }

   // test with various offsets of sphere's center
   for (real_t off = real_c(0.0); off < real_c(1.0); off += real_c(0.2))
   {
      WALBERLA_LOG_RESULT("Bubble offset " << off << " cells; normals computed only in interface cells;");
      test< LatticeModel_T >(uint_c(10), Vector3< real_t >(off, real_c(0), real_c(0)), false, false);

      WALBERLA_LOG_RESULT("Bubble offset "
                          << off << " cells; normals computed only in interface cells; fill level field smoothed;");
      test< LatticeModel_T >(uint_c(10), Vector3< real_t >(off, real_c(0), real_c(0)), false, true);

      WALBERLA_LOG_RESULT("Bubble offset " << off
                                           << " cells; normals computed in D3Q27 neighborhood of interface cells;");
      test< LatticeModel_T >(uint_c(10), Vector3< real_t >(off, real_c(0), real_c(0)), true, false);

      WALBERLA_LOG_RESULT(
         "Bubble offset "
         << off << " cells; normals computed in D3Q27 neighborhood of interface cells; fill level field smoothed;");
      test< LatticeModel_T >(uint_c(10), Vector3< real_t >(off, real_c(0), real_c(0)), true, true);
   }
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_RESULT("###################################");
   WALBERLA_LOG_RESULT("### Testing with D3Q19 stencil. ###");
   WALBERLA_LOG_RESULT("###################################");
   runAllTests< lbm::D3Q19< lbm::collision_model::SRT > >();

   WALBERLA_LOG_RESULT("###################################");
   WALBERLA_LOG_RESULT("### Testing with D3Q27 stencil. ###");
   WALBERLA_LOG_RESULT("###################################");
   runAllTests< lbm::D3Q27< lbm::collision_model::SRT > >();

   return EXIT_SUCCESS;
}
} // namespace NormalsOfSphereTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::NormalsOfSphereTest::main(argc, argv); }