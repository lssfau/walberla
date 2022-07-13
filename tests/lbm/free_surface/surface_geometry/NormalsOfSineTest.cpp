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
//! \file NormalsOfSineTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Initialize sine profile and compare calculated normals to analytical normals.
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/math/DistributedSample.h"

#include "field/AddToStorage.h"

#include "lbm/blockforest/communication/SimpleCommunication.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/NormalSweep.h"
#include "lbm/free_surface/surface_geometry/SmoothingSweep.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/D2Q9.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q19.h"
#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <cmath>

namespace walberla
{
namespace free_surface
{
namespace NormalsOfSineTest
{
// define types
using flag_t = uint32_t;

// define fields
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;
using FlagField_T   = FlagField< flag_t >;

// function describing the global sine profile
inline real_t function(real_t x, real_t amplitude, real_t offset, uint_t domainWidth)
{
   return amplitude * std::sin(x / real_c(domainWidth) * real_c(2) * math::pi) + offset;
}

// derivative of the function that is describing the sine profile
inline real_t derivative(real_t x, real_t amplitude, uint_t domainWidth)
{
   const real_t domainWidthInv = real_c(1) / real_c(domainWidth);
   return amplitude * std::cos(x * domainWidthInv * real_c(2) * math::pi) * real_c(2) * math::pi * domainWidthInv;
}

// compute and evaluate the absolute error of the angle of the interface normal in each interface cell;
// reference: derivative of sine function, i.e., analytically computed normal
template< typename Stencil_T >
class EvaluateNormalError
{
 public:
   EvaluateNormalError(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                       const ConstBlockDataID& normalFieldID, const ConstBlockDataID& flagFieldID,
                       const field::FlagUID& interfaceFlagID, const Vector3< uint_t >& domainSize, real_t amplitude,
                       bool computeNormalsInInterfaceNeighbors, bool smoothFillLevel)
      : blockForest_(blockForest), normalFieldID_(normalFieldID), flagFieldID_(flagFieldID),
        interfaceFlagID_(interfaceFlagID), domainSize_(domainSize), amplitude_(amplitude),
        computeNormalsInInterfaceNeighbors_(computeNormalsInInterfaceNeighbors), smoothFillLevel_(smoothFillLevel)
   {}

   void operator()(IBlock* const block)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // compute analytical normals and compare their directions with the computed normals
      const FlagField_T* const flagField     = block->getData< const FlagField_T >(flagFieldID_);
      const VectorField_T* const normalField = block->getData< const VectorField_T >(normalFieldID_);

      const flag_t interfaceFlag = flagField->getFlag(interfaceFlagID_);

      math::DistributedSample errorSample;

      // avoid OpenMP parallelization as DistributedSample can not be reduced by OpenMP
      WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, normalFieldIt, normalField, omp critical, {
         if (isFlagSet(flagFieldIt, interfaceFlag) ||
             (computeNormalsInInterfaceNeighbors_ && isFlagInNeighborhood< Stencil_T >(flagFieldIt, interfaceFlag)))
         {
            Cell globalCell = flagFieldIt.cell();
            blockForest->transformBlockLocalToGlobalCell(globalCell, *block, flagFieldIt.cell());

            // global x-location of this cell's center
            const real_t globalXCenter = real_c(globalCell[0]) + real_c(0.5);

            // get normal vector (slope of the negative inverse derivative gives direction of the normal)
            Vector3< real_t > analyticalNormal = Vector3< real_t >(
               real_c(1), -real_c(1) / derivative(globalXCenter, amplitude_, domainSize_[0]), real_c(0));
            analyticalNormal = analyticalNormal.getNormalized();

            // mirror vectors that are pointing upwards, i.e., from gas to liquid; in this FSLBM implementation, the
            // normal vector is defined to point from liquid to gas
            if (analyticalNormal[1] < real_c(0)) { analyticalNormal *= -real_c(1); }

            // calculate the error in the angle of the interface normal
            // the domain of arccosine is [-1,1]; due to numerical inaccuracies, the dot product
            // "*normalFieldIt*analyticalNormal" might be slightly out of this range; these cases are ignored here
            const real_t dotProduct = *normalFieldIt * analyticalNormal;
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
      if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 > || std::is_same_v< Stencil_T, stencil::D3Q27 >)
      {
         if (computeNormalsInInterfaceNeighbors_ && smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.036));
            WALBERLA_CHECK_LESS(maxError, real_c(0.135));
         }

         if (!computeNormalsInInterfaceNeighbors_ && smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.019));
            WALBERLA_CHECK_LESS(maxError, real_c(0.044));
         }

         if (!computeNormalsInInterfaceNeighbors_ && !smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.027));
            WALBERLA_CHECK_LESS(maxError, real_c(0.07));
         }

         if (computeNormalsInInterfaceNeighbors_ && !smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.158));
            WALBERLA_CHECK_LESS(maxError, real_c(1.58));
         }
      }

      // normal computation is less accurate if corner directions are not included (as opposed to D2Q9/D3Q27)
      if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >)
      {
         if (computeNormalsInInterfaceNeighbors_ && smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.036));
            WALBERLA_CHECK_LESS(maxError, real_c(0.135));
         }

         if (!computeNormalsInInterfaceNeighbors_ && smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.023));
            WALBERLA_CHECK_LESS(maxError, real_c(0.046));
         }

         if (!computeNormalsInInterfaceNeighbors_ && !smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.048));
            WALBERLA_CHECK_LESS(maxError, real_c(0.101));
         }

         if (computeNormalsInInterfaceNeighbors_ && !smoothFillLevel_)
         {
            WALBERLA_CHECK_LESS(meanError, real_c(0.158));
            WALBERLA_CHECK_LESS(maxError, real_c(1.58));
         }
      }
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   ConstBlockDataID normalFieldID_;
   ConstBlockDataID flagFieldID_;
   field::FlagUID interfaceFlagID_;

   Vector3< uint_t > domainSize_;
   real_t amplitude_;

   bool computeNormalsInInterfaceNeighbors_;
   bool smoothFillLevel_;
}; // class EvaluateNormalError

template< typename LatticeModel_T >
void test(uint_t domainWidth, real_t amplitude, real_t offset, uint_t fillLevelInitSamples,
          bool computeNormalsInInterfaceNeighbors, bool smoothFillLevel)
{
   // define types
   using Stencil_T       = typename LatticeModel_T::Stencil;
   using Communication_T = blockforest::SimpleCommunication< typename LatticeModel_T::CommunicationStencil >;
   using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

   // define the domain size
   const Vector3< uint_t > numBlocks(uint_c(1));
   Vector3< uint_t > domainSize(domainWidth);
   domainSize[2] = uint_c(1);
   const Vector3< uint_t > cellsPerBlock(domainSize[0] / numBlocks[0], domainSize[1] / numBlocks[1],
                                         domainSize[2] / numBlocks[2]);

   // create block grid
   const std::shared_ptr< StructuredBlockForest > blockForest =
      blockforest::createUniformBlockGrid(numBlocks[0], numBlocks[1], numBlocks[2],             // blocks
                                          cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2], // cells
                                          real_c(1.0),                                          // dx
                                          true,                                                 // one block per process
                                          true, false, true);                                   // periodicity

   // create (dummy) lattice model with omega=1
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

   // initialize sine profile such that there is exactly one period; every length is normalized with domainSize[0]
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      uint_t numTotalPoints = fillLevelInitSamples * fillLevelInitSamples;
      const real_t stepsize = real_c(1) / real_c(fillLevelInitSamples);

      ScalarField_T* const fillField = blockIt->getData< ScalarField_T >(fillFieldID);

      WALBERLA_FOR_ALL_CELLS(fillFieldIt, fillField, {
         // cell in block-local coordinates
         const Cell localCell = fillFieldIt.cell();

         // get cell in global coordinates
         Cell globalCell = fillFieldIt.cell();
         blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

         // Monte-Carlo like estimation of the fill level:
         // create uniformly-distributed sample points in each cell and count the number of points below the sine
         // profile; this fraction of points is used as the fill level to initialize the profile
         uint_t numPointsBelow = uint_c(0);

         for (uint_t xSample = uint_c(0); xSample < fillLevelInitSamples; ++xSample)
         {
            // value of the sine-function
            const real_t functionValue =
               function(real_c(globalCell[0]) + real_c(xSample) * stepsize, amplitude, offset, domainSize[0]);

            for (uint_t ySample = uint_c(0); ySample < fillLevelInitSamples; ++ySample)
            {
               const real_t yPoint = real_c(globalCell[1]) + real_c(ySample) * stepsize;
               if (yPoint < functionValue) { ++numPointsBelow; }
            }
         }

         // fill level is fraction of points below sine profile
         fillField->get(localCell) = real_c(numPointsBelow) / real_c(numTotalPoints);
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // communicate fill level field to have meaningful values in ghost layer cells in periodic directions
   Communication_T(blockForest, fillFieldID)();

   // initialize fill level field
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // communicate initialized flag field
   Communication_T(blockForest, flagFieldID)();

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(50));

   BlockDataID* relevantFillFieldID = &fillFieldID;

   if (smoothFillLevel)
   {
      smoothFillFieldID =
         field::addToStorage< ScalarField_T >(blockForest, "Smooth fill levels", real_c(0.0), field::fzyx, uint_c(1));

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
   EvaluateNormalError< Stencil_T > errorEvaluationSweep(blockForest, normalFieldID, flagFieldID,
                                                         flagIDs::interfaceFlagID, domainSize, amplitude,
                                                         computeNormalsInInterfaceNeighbors, smoothFillLevel);
   timeloop.add() << Sweep(errorEvaluationSweep, "Error evaluation sweep");

   // perform a single time step
   timeloop.singleStep();

   MPIManager::instance()->resetMPI();
}

template< typename LatticeModel_T >
void runAllTests()
{
   // used for initializing the fill level in a Monte-Carlo-like fashion; each cell is sampled by the specified value in
   // x- and y-direction
   const uint_t fillLevelInitSamples = uint_c(100);

   // test with various domain sizes, i.e., resolutions
   for (uint_t i = uint_c(50); i <= uint_c(200); i += uint_c(50))
   {
      const real_t amplitude = real_c(0.1) * real_c(i); // amplitude of the sine profile
      const real_t offset    = real_c(0.5) * real_c(i); // offset the sine profile in y-direction

      WALBERLA_LOG_RESULT("Domain size " << i << " cells; normals computed only in interface cells;");
      test< LatticeModel_T >(i, amplitude, offset, fillLevelInitSamples, false, false);

      WALBERLA_LOG_RESULT(
         "Domain size " << i << " cells; normals computed only in interface cells; fill level field smoothed;");
      test< LatticeModel_T >(i, amplitude, offset, fillLevelInitSamples, false, true);

      WALBERLA_LOG_RESULT("Domain size " << i
                                         << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells;");
      test< LatticeModel_T >(i, amplitude, offset, fillLevelInitSamples, true, false);

      WALBERLA_LOG_RESULT(
         "Domain size "
         << i << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells; fill level field smoothed;");
      test< LatticeModel_T >(i, amplitude, offset, fillLevelInitSamples, true, true);
   }

   // default values
   const uint_t domainWidth = uint_c(100);                       // size of the domain in x- and y-direction
   const real_t amplitude   = real_c(0.1) * real_c(domainWidth); // amplitude of the sine profile
   const real_t offset      = real_c(0.5) * real_c(domainWidth); // offset the sine profile in y-direction

   // test with various amplitudes
   for (real_t i = real_c(0.01); i < real_c(0.3); i += real_c(0.03))
   {
      WALBERLA_LOG_RESULT("Amplitude " << i << " cells; normals computed only in interface cells;");
      test< LatticeModel_T >(domainWidth, i * real_c(domainWidth), offset, fillLevelInitSamples, false, false);

      WALBERLA_LOG_RESULT(
         "Amplitude " << i << " cells; normals computed only in interface cells; fill level field smoothed;");
      test< LatticeModel_T >(domainWidth, i * real_c(domainWidth), offset, fillLevelInitSamples, false, true);

      WALBERLA_LOG_RESULT("Amplitude " << i
                                       << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells;");
      test< LatticeModel_T >(domainWidth, i * real_c(domainWidth), offset, fillLevelInitSamples, true, false);

      WALBERLA_LOG_RESULT(
         "Amplitude "
         << i << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells; fill level field smoothed;");
      test< LatticeModel_T >(domainWidth, i * real_c(domainWidth), offset, fillLevelInitSamples, true, true);
   }

   // test with various offsets, i.e., position of the sine-profile's zero-line
   for (real_t i = real_c(0.4); i < real_c(0.6); i += real_c(0.03))
   {
      WALBERLA_LOG_RESULT("Offset " << i << " cells; normals computed only in interface cells;");
      test< LatticeModel_T >(domainWidth, amplitude, i * real_c(domainWidth), fillLevelInitSamples, false, false);

      WALBERLA_LOG_RESULT("Offset " << i
                                    << " cells; normals computed only in interface cells; fill level field smoothed;");
      test< LatticeModel_T >(domainWidth, amplitude, i * real_c(domainWidth), fillLevelInitSamples, false, true);

      WALBERLA_LOG_RESULT("Offset " << i << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells;");
      test< LatticeModel_T >(domainWidth, amplitude, i * real_c(domainWidth), fillLevelInitSamples, true, false);

      WALBERLA_LOG_RESULT(
         "Offset "
         << i << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells; fill level field smoothed;");
      test< LatticeModel_T >(domainWidth, amplitude, i * real_c(domainWidth), fillLevelInitSamples, true, true);
   }

   // test with different fill level initialization samples (more samples => initialization of fill levels is
   // closer to the real sine-profile)
   for (uint_t i = uint_c(50); i <= uint_c(200); i += uint_c(50))
   {
      WALBERLA_LOG_RESULT("Fill level initialization sample " << i
                                                              << " cells; normals computed only in interface cells;");
      test< LatticeModel_T >(domainWidth, amplitude, offset, i, false, false);

      WALBERLA_LOG_RESULT("Fill level initialization sample "
                          << i << " cells; normals computed only in interface cells; fill level field smoothed;");
      test< LatticeModel_T >(domainWidth, amplitude, offset, i, false, true);

      WALBERLA_LOG_RESULT("Fill level initialization sample "
                          << i << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells;");
      test< LatticeModel_T >(domainWidth, amplitude, offset, i, true, false);

      WALBERLA_LOG_RESULT(
         "Fill level initialization sample "
         << i << " cells; normals computed in D2Q9/D3Q27 neighborhood of interface cells; fill level field smoothed;");
      test< LatticeModel_T >(domainWidth, amplitude, offset, i, true, true);
   }
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   WALBERLA_LOG_RESULT("##################################");
   WALBERLA_LOG_RESULT("### Testing with D2Q9 stencil. ###");
   WALBERLA_LOG_RESULT("##################################");
   runAllTests< lbm::D2Q9< lbm::collision_model::SRT > >();

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
} // namespace NormalsOfSineTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::NormalsOfSineTest::main(argc, argv); }