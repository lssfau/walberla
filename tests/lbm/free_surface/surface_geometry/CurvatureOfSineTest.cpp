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
//! \file CurvatureOfSineTest.cpp
//! \ingroup lbm/free_surface/surface_geometry
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Initialize sine profile and compare computed curvature with analytical curvature.
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
#include "lbm/free_surface/boundary/FreeSurfaceBoundaryHandling.h"
#include "lbm/free_surface/surface_geometry/CurvatureSweep.h"
#include "lbm/free_surface/surface_geometry/ExtrapolateNormalsSweep.h"
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
namespace CurvatureOfSineTest
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

// second derivative of the function that is describing the sine profile
inline real_t secondDerivative(real_t x, real_t amplitude, uint_t domainWidth)
{
   const real_t domainWidthInv = real_c(1) / real_c(domainWidth);

   // the minus has been neglected on purpose: due to the definition of the normal (pointing from liquid to gas), the
   // curvature also has a different sign
   return amplitude * std::sin(x * domainWidthInv * real_c(2) * math::pi) * real_c(4) * math::pi * math::pi *
          domainWidthInv * domainWidthInv;
}

// compute the analytical normal vector and evaluate the error of the computed normal (absolute error in angle)
template< typename Stencil_T >
class ComputeAnalyticalNormal
{
 public:
   ComputeAnalyticalNormal(const std::weak_ptr< const StructuredBlockForest >& blockForest, BlockDataID& normalFieldID,
                           const ConstBlockDataID& flagFieldID, const field::FlagUID& interfaceFlagID,
                           const Vector3< uint_t >& domainSize, real_t amplitude)
      : blockForest_(blockForest), normalFieldID_(normalFieldID), flagFieldID_(flagFieldID),
        interfaceFlagID_(interfaceFlagID), domainSize_(domainSize), amplitude_(amplitude)
   {}

   void operator()(IBlock* const block)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // compute analytical normals
      const FlagField_T* const flagField = block->getData< const FlagField_T >(flagFieldID_);
      VectorField_T* const normalField   = block->getData< VectorField_T >(normalFieldID_);

      const flag_t interfaceFlag = flagField->getFlag(interfaceFlagID_);

      // explicitly use either D2Q9 or D3Q27 here, as the geometry operations require (or are most accurate with) the
      // full neighborhood;
      using NeighborhoodStencil_T =
         typename std::conditional< Stencil_T::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;

      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, normalFieldIt, normalField, {
         // only treat interface cells
         if (!isFlagSet(flagFieldIt, interfaceFlag) &&
             !field::isFlagInNeighborhood< NeighborhoodStencil_T >(flagFieldIt, interfaceFlag))
         {
            continue;
         }

         Cell globalCell = flagFieldIt.cell();
         blockForest->transformBlockLocalToGlobalCell(globalCell, *block, flagFieldIt.cell());

         // global x-location of this cell's center
         const real_t globalXCenter = real_c(globalCell[0]) + real_c(0.5);

         // get normal vector (slope of the negative inverse derivative gives direction of the normal)
         Vector3< real_t > normalVec =
            Vector3< real_t >(real_c(1), -real_c(1) / derivative(globalXCenter, amplitude_, domainSize_[0]), real_c(0));
         normalVec = normalVec.getNormalized();

         // mirror vectors that are pointing downwards, as normal vector is defined to point from fluid to gas in FSLBM
         if (normalVec[1] < real_c(0)) { normalVec *= -real_c(1); }

         *normalFieldIt = normalVec;
      }) // WALBERLA_FOR_ALL_CELLS
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   BlockDataID normalFieldID_;
   ConstBlockDataID flagFieldID_;

   field::FlagUID interfaceFlagID_;

   Vector3< uint_t > domainSize_;
   real_t amplitude_;
}; // class ComputeAnalyticalNormal

// compute the analytical normal vector and evaluate the error of the computed normal (absolute error in angle)
template< typename Stencil_T >
class ComputeCurvatureError
{
 public:
   ComputeCurvatureError(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                         const ConstBlockDataID& curvatureFieldID, const ConstBlockDataID& flagFieldID,
                         const field::FlagUID& interfaceFlagID, const Vector3< uint_t >& domainSize, real_t amplitude,
                         const std::shared_ptr< real_t >& l2Error)
      : blockForest_(blockForest), curvatureFieldID_(curvatureFieldID), flagFieldID_(flagFieldID),
        interfaceFlagID_(interfaceFlagID), domainSize_(domainSize), amplitude_(amplitude), l2Error_(l2Error)
   {}

   void operator()(IBlock* const block)
   {
      const auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      // compute analytical normals and compare their directions with the computed normals
      const FlagField_T* const flagField        = block->getData< const FlagField_T >(flagFieldID_);
      const ScalarField_T* const curvatureField = block->getData< const ScalarField_T >(curvatureFieldID_);

      const flag_t interfaceFlag = flagField->getFlag(interfaceFlagID_);

      real_t curvDiffSum2 = real_c(0);
      real_t analCurvSum2 = real_c(0);

      WALBERLA_FOR_ALL_CELLS_OMP(flagFieldIt, flagField, curvatureFieldIt, curvatureField,
                                 omp parallel for schedule(static) reduction(+:curvDiffSum2) reduction(+:analCurvSum2),
      {
         // only treat interface cells
         if (isFlagSet(flagFieldIt, interfaceFlag))
         {
            Cell globalCell = flagFieldIt.cell();
            blockForest->transformBlockLocalToGlobalCell(globalCell, *block, flagFieldIt.cell());

            // normalized global x-location of this cell's center
            const real_t globalXCenter = real_c(globalCell[0]) + real_c(0.5);

            // get analytical curvature
            const real_t analyticalCurvature = secondDerivative(globalXCenter, amplitude_, domainSize_[0]);

            // calculate the relative error in curvature (dx=1/domainSize[0] for converting from LBM to physical units)
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

   Vector3< uint_t > domainSize_;
   real_t amplitude_;

   std::shared_ptr< real_t > l2Error_;
}; // class ComputeCurvatureError

template< typename LatticeModel_T >
real_t test(uint_t domainWidth, real_t amplitude, real_t offset, uint_t fillLevelInitSamples, bool useTriangulation,
            bool useAnalyticalNormal)
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

   // create lattice model with omega=1
   LatticeModel_T latticeModel(real_c(1.0));

   // add fields
   BlockDataID pdfFieldID =
      lbm::addPdfFieldToStorage< LatticeModel_T >(blockForest, "PDF field", latticeModel, field::fzyx);
   BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill levels", real_c(1), field::fzyx, uint_c(1));
   BlockDataID normalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Normals", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));
   BlockDataID curvatureFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Curvature", real_c(0), field::fzyx, uint_c(1));
   BlockDataID smoothFillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Smooth fill levels", real_c(1.0), field::fzyx, uint_c(1));

   // obstacle normal field is only a dummy field here (wetting effects are not considered in this test)
   BlockDataID dummyObstNormalFieldID = field::addToStorage< VectorField_T >(
      blockForest, "Dummy obstacle normal field", Vector3< real_t >(real_c(0)), field::fzyx, uint_c(1));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();

   // initialize sine profile such that there is exactly one period; every length is normalized with domainSize[0]
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const uint_t numTotalPoints = fillLevelInitSamples * fillLevelInitSamples;
      const real_t stepsize       = real_c(1) / real_c(fillLevelInitSamples);

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

   // contact angle is only a dummy object here (wetting effects are not considered in this test)
   ContactAngle dummyContactAngle(real_c(0));

   // create timeloop
   SweepTimeloop timeloop(blockForest, uint_c(50));

   BlockDataID* relevantFillFieldID = &fillFieldID;

   if (useAnalyticalNormal)
   {
      // add sweep for getting analytical interface normal
      ComputeAnalyticalNormal< Stencil_T > analyticalNormalSweep(blockForest, normalFieldID, flagFieldID,
                                                                 flagIDs::interfaceFlagID, domainSize, amplitude);
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
         freeSurfaceBoundaryHandling->getFlagInfo().getObstacleIDSet(), !useTriangulation, false, false, false);
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

   // add sweep for computing the error (with respect to analytical solution) of the interface curvature
   std::shared_ptr< real_t > l2Error = std::make_shared< real_t >(real_c(0));
   ComputeCurvatureError< Stencil_T > errorEvaluationSweep(blockForest, curvatureFieldID, flagFieldID,
                                                           flagIDs::interfaceFlagID, domainSize, amplitude, l2Error);
   timeloop.add() << Sweep(errorEvaluationSweep, "Error evaluation sweep");

   // perform a single time step
   timeloop.singleStep();

   MPIManager::instance()->resetMPI();

   return *l2Error;
}

template< typename LatticeModel_T >
void runAllTests()
{
   using Stencil_T = typename LatticeModel_T::Stencil;

   real_t l2Error;

   // used for initializing the fill level in a Monte-Carlo-like fashion; each cell is sampled by the specified value in
   // x- and y-direction
   const uint_t fillLevelInitSamples = uint_c(100);

   // test with various domain sizes, i.e., resolutions
   for (uint_t domainWidth = uint_c(50); domainWidth <= uint_c(200); domainWidth += uint_c(50))
   {
      const real_t amplitude = real_c(0.1) * real_c(domainWidth); // amplitude of the sine profile
      const real_t offset    = real_c(0.5) * real_c(domainWidth); // offset the sine profile in y-direction

      WALBERLA_LOG_RESULT("Domain width " << domainWidth << " cells; curvature with finite difference method;");
      l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, false, false);
      if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 > || std::is_same_v< Stencil_T, stencil::D3Q27 >)
      {
         WALBERLA_CHECK_LESS(l2Error, real_c(0.63));
      }
      // computation is less accurate if corner directions are not included (as opposed to D2Q9/D3Q27)
      if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >) { WALBERLA_CHECK_LESS(l2Error, real_c(0.756)); }

      WALBERLA_LOG_RESULT("Domain width "
                          << domainWidth
                          << " cells; curvature with finite difference method; normal from analytical solution;");
      l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, false, true);
      if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 > || std::is_same_v< Stencil_T, stencil::D3Q27 >)
      {
         WALBERLA_CHECK_LESS(l2Error, real_c(0.52));
      }
      // computation is less accurate if corner directions are not included (as opposed to D2Q9/D3Q27)
      if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >) { WALBERLA_CHECK_LESS(l2Error, real_c(0.58)); }

      if constexpr (!std::is_same_v< Stencil_T, stencil::D2Q9 >)
      {
         WALBERLA_LOG_RESULT("Domain width " << domainWidth << " cells; curvature with local triangulation;");
         l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, true, false);
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q27 >) { WALBERLA_CHECK_LESS(l2Error, real_c(0.79)); }
         // computation is less accurate if corner directions are not included (as opposed to D2Q9/D3Q27)
         if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >) { WALBERLA_CHECK_LESS(l2Error, real_c(0.845)); }

         WALBERLA_LOG_RESULT("Domain width "
                             << domainWidth
                             << " cells; curvature with local triangulation; normal from analytical solution;");
         l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, true, true);
         WALBERLA_CHECK_LESS(l2Error, real_c(0.71));
      }
   }

   // test with various amplitudes
   for (real_t i = real_c(0.1); i <= real_c(0.3); i += real_c(0.05))
   {
      const uint_t domainWidth = uint_c(100);                       // size of the domain in x- and y-direction
      const real_t amplitude   = i * real_c(domainWidth);           // amplitude of the sine profile
      const real_t offset      = real_c(0.5) * real_c(domainWidth); // offset the sine profile in y-direction

      WALBERLA_LOG_RESULT("Amplitude " << amplitude << " cells; curvature with finite difference method;");
      l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, false, false);
      if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 > || std::is_same_v< Stencil_T, stencil::D3Q27 >)
      {
         WALBERLA_CHECK_LESS(l2Error, real_c(0.76));
      }
      // computation is less accurate if corner directions are not included (as opposed to D2Q9/D3Q27)
      if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >) { WALBERLA_CHECK_LESS(l2Error, real_c(0.79)); }

      WALBERLA_LOG_RESULT(
         "Amplitude " << amplitude
                      << " cells; curvature with finite difference method; normal from analytical solution;");
      l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, false, true);
      if constexpr (std::is_same_v< Stencil_T, stencil::D2Q9 > || std::is_same_v< Stencil_T, stencil::D3Q27 >)
      {
         WALBERLA_CHECK_LESS(l2Error, real_c(0.77));
      }
      // computation is less accurate if corner directions are not included (as opposed to D2Q9/D3Q27)
      if constexpr (std::is_same_v< Stencil_T, stencil::D3Q19 >) { WALBERLA_CHECK_LESS(l2Error, real_c(0.79)); }

      if constexpr (!std::is_same_v< Stencil_T, stencil::D2Q9 >)
      {
         WALBERLA_LOG_RESULT("Amplitude " << amplitude << " cells; curvature with local triangulation;");
         l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, true, false);
         WALBERLA_CHECK_LESS(l2Error, real_c(0.8));

         WALBERLA_LOG_RESULT("Amplitude "
                             << amplitude
                             << " cells; curvature with local triangulation; normal from analytical solution;");
         l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, true, true);
         WALBERLA_CHECK_LESS(l2Error, real_c(0.8));
      }
   }

   // test with various offsets
   for (real_t i = real_c(0.4); i <= real_c(0.6); i += real_c(0.04))
   {
      const uint_t domainWidth = uint_c(100);                       // size of the domain in x- and y-direction
      const real_t amplitude   = real_c(0.2) * real_c(domainWidth); // amplitude of the sine profile
      const real_t offset      = i * real_c(domainWidth);           // offset the sine profile in y-direction

      WALBERLA_LOG_RESULT("Offset " << offset << " cells; curvature with finite difference method;");
      l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, false, false);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.72));

      WALBERLA_LOG_RESULT(
         "Offset " << offset << " cells; curvature with finite difference method; normal from analytical solution;");
      l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, false, true);
      WALBERLA_CHECK_LESS(l2Error, real_c(0.73));

      if constexpr (!std::is_same_v< Stencil_T, stencil::D2Q9 >)
      {
         WALBERLA_LOG_RESULT("Offset " << offset << " cells; curvature with local triangulation;");
         l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, true, false);
         WALBERLA_CHECK_LESS(l2Error, real_c(0.699));

         WALBERLA_LOG_RESULT(
            "Offset " << offset << " cells; curvature with local triangulation; normal from analytical solution;");
         l2Error = test< LatticeModel_T >(domainWidth, amplitude, offset, fillLevelInitSamples, true, true);
         WALBERLA_CHECK_LESS(l2Error, real_c(0.706));
      }
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
} // namespace CurvatureOfSineTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::CurvatureOfSineTest::main(argc, argv); }