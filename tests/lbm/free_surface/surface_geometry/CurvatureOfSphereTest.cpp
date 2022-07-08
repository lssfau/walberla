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
//! \file CurvatureOfSphereTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Compare computed mean curvature of spherical gas bubbles with analytical curvature.
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
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/ContactAngle.h"
#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"
#include "lbm/free_surface/surface_geometry/ExtrapolateNormalsSweep.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

namespace walberla
{
namespace free_surface
{
namespace CurvatureOfSphereTest
{
// define types
using flag_t = uint32_t;

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

template< typename Stencil_T >
class ComputeAnalyticalNormal
{
 public:
   ComputeAnalyticalNormal(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                           const BlockDataID& normalField, const ConstBlockDataID& flagField,
                           const field::FlagUID& interfaceFlag, const Vector3< real_t > sphereMidpoint)
      : blockForest_(blockForest), normalFieldID_(normalField), flagFieldID_(flagField),
        interfaceFlagID_(interfaceFlag), sphereMidpoint_(sphereMidpoint)
   {}

   void operator()(IBlock* const block)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // get fields
      const FlagField_T* const flagField = block->getData< const FlagField_T >(flagFieldID_);
      VectorField_T* const normalField   = block->getData< VectorField_T >(normalFieldID_);

      flag_t interfaceFlag = flagField->getFlag(interfaceFlagID_);

      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, normalFieldIt, normalField, {
         // evaluate normal only in interface cell (and possibly in an interface cell's neighborhood)
         if (isFlagSet(flagFieldIt, interfaceFlag) ||
             field::isFlagInNeighborhood< Stencil_T >(flagFieldIt, interfaceFlag))
         {
            // get a vector pointing from sphere's center to current cell's center
            Vector3< real_t > r;

            // get the global coordinate of the current cell's center
            blockForest->getBlockLocalCellCenter(*block, flagFieldIt.cell(), r);
            r -= sphereMidpoint_;

            // invert normal direction: in this FSLBM implementation, the normal vector is defined to point from liquid
            // to gas
            r *= real_c(-1);
            normalize(r);

            // store analytical normal
            *normalFieldIt = r;
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;

   BlockDataID normalFieldID_;
   ConstBlockDataID flagFieldID_;
   field::FlagUID interfaceFlagID_;

   Vector3< real_t > sphereMidpoint_;
}; // class ComputeAnalyticalNormal

class ComputeCurvatureError
{
 public:
   ComputeCurvatureError(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                         const ConstBlockDataID& curvatureFieldID, const ConstBlockDataID& flagFieldID,
                         const field::FlagUID& interfaceFlagID, real_t sphereDiameter,
                         const Vector3< uint_t >& domainSize, const std::shared_ptr< real_t >& l2Error)
      : blockForest_(blockForest), curvatureFieldID_(curvatureFieldID), flagFieldID_(flagFieldID),
        interfaceFlagID_(interfaceFlagID), sphereDiameter_(sphereDiameter), domainSize_(domainSize), l2Error_(l2Error)
   {}

   void operator()(IBlock* const block)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // get fields
      const ScalarField_T* const curvatureField = block->getData< const ScalarField_T >(curvatureFieldID_);
      const FlagField_T* const flagField        = block->getData< const FlagField_T >(flagFieldID_);

      const flag_t interfaceFlag = flagField->getFlag(interfaceFlagID_);

      const real_t analyticalCurvature = real_c(-1) / (real_c(0.5) * real_c(sphereDiameter_));

      real_t curvDiffSum2 = real_c(0);
      real_t analCurvSum2 = real_c(0);

      WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, curvatureFieldIt, curvatureField,
                                 omp parallel for schedule(static) reduction(+:curvDiffSum2) reduction(+:analCurvSum2),
      {
         // skip non-interface cells
         if (isFlagSet(flagFieldIt, interfaceFlag))
         {
            const real_t curvDiff = *curvatureFieldIt - analyticalCurvature;
            curvDiffSum2 += curvDiff * curvDiff;
            analCurvSum2 += analyticalCurvature * analyticalCurvature;
         }
      }) // WALBERLA_FOR_ALL_CELLS_OMP

      mpi::allReduceInplace< real_t >(curvDiffSum2, mpi::SUM);
      mpi::allReduceInplace< real_t >(analCurvSum2, mpi::SUM);

      *l2Error_ = std::pow(curvDiffSum2 / analCurvSum2, real_c(0.5));

      WALBERLA_LOG_RESULT("Relative error in curvature according to L2 norm = " << *l2Error_);
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;

   ConstBlockDataID curvatureFieldID_;
   ConstBlockDataID flagFieldID_;
   field::FlagUID interfaceFlagID_;

   real_t sphereDiameter_;
   Vector3< uint_t > domainSize_;

   std::shared_ptr< real_t > l2Error_;
}; // class ComputeCurvatureError

template< typename LatticeModel_T >
real_t test(uint_t sphereDiameter, Vector3< real_t > offset, bool useTriangulation, bool useAnalyticalNormal)
{
   // define types
   using Stencil_T                     = typename LatticeModel_T::Stencil;
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;
   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::CommunicationStencil >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   const Vector3< uint_t > domainSize(sphereDiameter + uint_c(6));
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
   BlockDataID curvatureFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Curvature", real_c(0), field::fzyx, uint_c(1));
   BlockDataID smoothFillFieldID;

   // obstacle normal field is only a dummy field here (wetting effects are not considered in this test)
   BlockDataID dummyObstNormalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Dummy obstacle normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();

   // add gas sphere, i.e., bubble
   Vector3< real_t > midpoint(real_c(domainSize[0]) * real_c(0.5), real_c(domainSize[1]) * real_c(0.5),
                              real_c(domainSize[2]) * real_c(0.5));
   geometry::Sphere sphere(midpoint + offset, real_c(sphereDiameter) * real_c(0.5));
   freeSurfaceBoundaryHandling->addFreeSurfaceObject(sphere);
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // contact angle is only a dummy object here (wetting effects are not considered in this test)
   ContactAngle dummyContactAngle(real_c(0));

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(1));

   BlockDataID* relevantFillFieldID = &fillFieldID;

   if (useAnalyticalNormal)
   {
      // add sweep for getting analytical interface normal
      ComputeAnalyticalNormal< Stencil_T > analyticalNormalSweep(blockForest, normalFieldID, flagFieldID,
                                                                 flagIDs::interfaceFlagID, midpoint);
      timeloop.add() << Sweep(analyticalNormalSweep, "Analytical normal sweep")
                     << AfterFunction(Communication_T(blockForest, normalFieldID),
                                      "Communication after analytical normal sweep");
   }
   else
   {
      if (!useTriangulation)
      {
         smoothFillFieldID = field::addToStorage< ScalarField_T >(blockForest, "Smooth fill levels", real_c(1.0),
                                                                  field::fzyx, uint_c(1));

         // add sweep for smoothing the fill level field when not using local triangulation
         SmoothingSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > smoothingSweep(
            smoothFillFieldID, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs,
            freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), false);
         timeloop.add() << Sweep(smoothingSweep, "Smoothing sweep")
                        << AfterFunction(Communication_T(blockForest, smoothFillFieldID),
                                         "Communication after smoothing sweep");

         relevantFillFieldID = &smoothFillFieldID;
      }

      NormalSweep< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > normalsSweep(
         normalFieldID, *relevantFillFieldID, flagFieldID, flagIDs::interfaceFlagID, flagIDs::liquidInterfaceGasFlagIDs,
         freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), true, false, false, false);
      timeloop.add() << Sweep(normalsSweep, "Normal sweep")
                     << AfterFunction(Communication_T(blockForest, normalFieldID), "Communication after normal sweep");
   }

   if (useTriangulation) // use local triangulation for curvature computation
   {
      CurvatureSweepLocalTriangulation< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvatureSweep(
         blockForest, curvatureFieldID, normalFieldID, fillFieldID, flagFieldID, dummyObstNormalFieldID,
         flagIDs::interfaceFlagID, freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), false,
         dummyContactAngle);
      timeloop.add() << Sweep(curvatureSweep, "Curvature sweep");
   }
   else // use finite differences for curvature computation
   {
      CurvatureSweepFiniteDifferences< Stencil_T, FlagField_T, ScalarField_T, VectorField_T > curvatureSweep(
         curvatureFieldID, normalFieldID, dummyObstNormalFieldID, flagFieldID, flagIDs::interfaceFlagID,
         flagIDs::liquidInterfaceGasFlagIDs, freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), false,
         dummyContactAngle);
      timeloop.add() << Sweep(curvatureSweep, "Curvature sweep");
   }

   const std::shared_ptr< real_t > l2Error = std::make_shared< real_t >(real_c(0));
   ComputeCurvatureError errorEvaluationSweep(blockForest, curvatureFieldID, flagFieldID, flagIDs::interfaceFlagID,
                                              real_c(sphereDiameter), domainSize, l2Error);
   timeloop.add() << Sweep(errorEvaluationSweep, "Error evaluation sweep");

   // perform a single time step
   timeloop.singleStep();

   MPIManager::instance()->resetMPI();

   return *l2Error;
}

template< typename LatticeModel_T >
void runAllTests()
{
   real_t l2Error;

   // test with various bubble diameters
   for (uint_t diameter = uint_c(10); diameter <= uint_c(50); diameter += uint_c(10))
   {
      WALBERLA_LOG_RESULT("Bubble diameter " << diameter << " cells; curvature with finite difference method;")
      l2Error = test< LatticeModel_T >(diameter, Vector3< real_t >(real_c(0), real_c(0), real_c(0)), false, false);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.2));

      WALBERLA_LOG_RESULT("Bubble diameter "
                          << diameter
                          << " cells; curvature with finite difference method; normal from analytical solution;")
      l2Error = test< LatticeModel_T >(diameter, Vector3< real_t >(real_c(0), real_c(0), real_c(0)), false, true);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.19));

      WALBERLA_LOG_RESULT("Bubble diameter " << diameter << " cells; curvature with local triangulation;")
      l2Error = test< LatticeModel_T >(diameter, Vector3< real_t >(real_c(0), real_c(0), real_c(0)), true, false);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.29));

      WALBERLA_LOG_RESULT("Bubble diameter "
                          << diameter << " cells; curvature with local triangulation; normal from analytical solution;")
      l2Error = test< LatticeModel_T >(diameter, Vector3< real_t >(real_c(0), real_c(0), real_c(0)), true, true);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.29));
   }

   // test with various offsets of sphere's center
   for (real_t off = real_c(0.1); off < real_c(1.0); off += real_c(0.2))
   {
      WALBERLA_LOG_RESULT("Sphere center offset " << off << " cells; curvature with finite difference method;")
      l2Error = test< LatticeModel_T >(uint_c(20), Vector3< real_t >(off, off, real_c(0)), false, false);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.82));

      WALBERLA_LOG_RESULT("Sphere center offset "
                          << off << " cells; curvature with finite difference method; normal from analytical solution;")
      l2Error = test< LatticeModel_T >(uint_c(20), Vector3< real_t >(off, off, real_c(0)), false, true);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.83));

      WALBERLA_LOG_RESULT("Sphere center offset " << off << " cells; curvature with local triangulation;")
      l2Error = test< LatticeModel_T >(uint_c(20), Vector3< real_t >(off, off, real_c(0)), true, false);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.13));

      WALBERLA_LOG_RESULT("Sphere center offset "
                          << off << " cells; curvature with local triangulation; normal from analytical solution;")
      l2Error = test< LatticeModel_T >(uint_c(20), Vector3< real_t >(off, off, real_c(0)), true, true);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.16));
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
} // namespace CurvatureOfSphereTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::CurvatureOfSphereTest::main(argc, argv); }