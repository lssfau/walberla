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
//! \file TaylorBubble.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief This showcase simulates a rising gas bubble in a cylindrical tube.
//
// This showcase simulates an elongated gas bubble in a cylindrical tube. Reference experiments are available from
// Bugg, Saad (2002), doi:10.1016/S0301-9322(02)00002-2
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/StringUtility.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/LoadBalancing.h"
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/TotalMassComputer.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/free_surface/surface_geometry/Utility.h"
#include "lbm/lattice_model/D3Q19.h"

namespace walberla
{
namespace free_surface
{
namespace TaylorBubble
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using CollisionModel_T      = lbm::collision_model::SRT;
using ForceModel_T          = lbm::force_model::GuoField< VectorField_T >;
using LatticeModel_T        = lbm::D3Q19< CollisionModel_T, true, ForceModel_T, 2 >;
using LatticeModelStencil_T = LatticeModel_T::Stencil;
using PdfField_T            = lbm::PdfField< LatticeModel_T >;
using PdfCommunication_T    = blockforest::SimpleCommunication< LatticeModelStencil_T >;

// the geometry computations in SurfaceGeometryHandler require meaningful values in the ghost layers in corner
// directions (flag field and fill level field); this holds, even if the lattice model uses a D3Q19 stencil
using CommunicationStencil_T =
   std::conditional_t< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >;
using Communication_T = blockforest::SimpleCommunication< CommunicationStencil_T >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

// write each entry in "vector" to line in a file; the first column contains the current timestep; all values are
// separated by tabs
template< typename T >
void writeVectorToFile(const std::vector< T >& vector, uint_t timestep, const std::string& filename);

// IMPORTANT REMARK: this does not work when multiple bubbles are in the system (e.g. after splitting)
template< typename FreeSurfaceBoundaryHandling_T >
class CenterOfMassComputer
{
 public:
   CenterOfMassComputer(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                        const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                        uint_t frequency, const std::shared_ptr< Vector3< real_t > >& centerOfMass)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling),
        centerOfMass_(centerOfMass), frequency_(frequency), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % frequency_ != uint_c(0)) { return; }

      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      *centerOfMass_ = Vector3< real_t >(real_c(0));

      // use real_t to avoid integer overflow in sum below when using high grid resolutions
      Vector3< real_t > cellSumVector = Vector3< real_t >(real_c(0));
      real_t cellCount                = real_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID);

         WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, {
            if (flagInfo.isGas(flagFieldIt))
            {
               Cell globalCell = flagFieldIt.cell();
               // transform local cell to global coordinates
               blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt);
               cellSumVector += Vector3< real_t >(real_c(globalCell[0]), real_c(globalCell[1]), real_c(globalCell[2]));
               cellCount += real_c(1);
            }
         }) // WALBERLA_FOR_ALL_CELLS
      }

      mpi::allReduceInplace< real_t >(cellCount, mpi::SUM);
      mpi::allReduceInplace< real_t >(cellSumVector, mpi::SUM);

      (*centerOfMass_)[0] = cellSumVector[0] / cellCount;
      (*centerOfMass_)[1] = cellSumVector[1] / cellCount;
      (*centerOfMass_)[2] = cellSumVector[2] / cellCount;
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;
   std::shared_ptr< Vector3< real_t > > centerOfMass_;

   uint_t frequency_;
   uint_t executionCounter_;
}; // class CenterOfMassComputer

int main(int argc, char** argv)
{
   Environment walberlaEnv(argc, argv);

   if (argc < 2) { WALBERLA_ABORT("Please specify a parameter file as input argument.") }

   // print content of parameter file
   WALBERLA_LOG_INFO_ON_ROOT(*walberlaEnv.config());

   // get block forest parameters from parameter file
   auto blockForestParameters              = walberlaEnv.config()->getOneBlock("BlockForestParameters");
   const Vector3< uint_t > cellsPerBlock   = blockForestParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");
   const Vector3< bool > periodicity       = blockForestParameters.getParameter< Vector3< bool > >("periodicity");
   const uint_t loadBalancingFrequency     = blockForestParameters.getParameter< uint_t >("loadBalancingFrequency");
   const bool printLoadBalancingStatistics = blockForestParameters.getParameter< bool >("printLoadBalancingStatistics");

   // read domain parameters from parameter file
   const auto domainParameters              = walberlaEnv.config()->getOneBlock("DomainParameters");
   const real_t tubeDiameter                = domainParameters.getParameter< real_t >("tubeDiameter");
   const real_t bubbleDiameter              = domainParameters.getParameter< real_t >("bubbleDiameter");
   const real_t bubbleHeight                = domainParameters.getParameter< real_t >("bubbleHeight");
   const Vector3< real_t > bubbleBottomEnd  = domainParameters.getParameter< Vector3< real_t > >("bubbleBottomEnd");
   const Vector3< real_t > domainSizeFactor = domainParameters.getParameter< Vector3< real_t > >("domainSizeFactor");

   // define domain size
   Vector3< uint_t > domainSize;
   domainSize[0] = uint_c(domainSizeFactor[0] * tubeDiameter);
   domainSize[1] = uint_c(domainSizeFactor[1] * tubeDiameter);
   domainSize[2] = uint_c(domainSizeFactor[2] * tubeDiameter);

   // compute number of blocks as defined by domainSize and cellsPerBlock
   Vector3< uint_t > numBlocks;
   numBlocks[0] = uint_c(std::ceil(real_c(domainSize[0]) / real_c(cellsPerBlock[0])));
   numBlocks[1] = uint_c(std::ceil(real_c(domainSize[1]) / real_c(cellsPerBlock[1])));
   numBlocks[2] = uint_c(std::ceil(real_c(domainSize[2]) / real_c(cellsPerBlock[2])));

   // get number of (MPI) processes
   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   WALBERLA_CHECK_LESS_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
                             "The number of MPI processes is greater than the number of blocks as defined by "
                             "\"domainSize/cellsPerBlock\". This would result in unused MPI processes. Either decrease "
                             "the number of MPI processes or increase \"cellsPerBlock\".")

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numProcesses);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellsPerBlock);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numBlocks);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(tubeDiameter);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(bubbleDiameter);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(bubbleHeight);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(bubbleBottomEnd);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSizeFactor);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(periodicity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSize);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellsPerBlock);

   Vector3< uint_t > realDomainSize;
   realDomainSize[0] = cellsPerBlock[0] * numBlocks[0];
   realDomainSize[1] = cellsPerBlock[1] * numBlocks[1];
   realDomainSize[2] = cellsPerBlock[2] * numBlocks[2];

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(realDomainSize);

   if (domainSize[0] != realDomainSize[0] && periodicity[0])
   {
      WALBERLA_ABORT(
         "The specified domain size in x-direction can not be obtained with the numer of blocks you specified.")
   }
   if (domainSize[1] != realDomainSize[1] && periodicity[1])
   {
      WALBERLA_ABORT(
         "The specified domain size in y-direction can not be obtained with the numer of blocks you specified.")
   }
   // the z-direction must be no slip in this setup and is simply extended (see below) to obtain the specified size
   int boundaryThicknessZ = int_c(realDomainSize[2]) - int_c(domainSize[2]);
   if (boundaryThicknessZ < 0)
   {
      WALBERLA_ABORT("Something went wrong: the resulting domain size in z-direction is less than specified.")
   }

   // read physics parameters from parameter file
   const auto physicsParameters = walberlaEnv.config()->getOneBlock("PhysicsParameters");
   const real_t bondNumber      = physicsParameters.getParameter< real_t >("bondNumber");
   const real_t mortonNumber    = physicsParameters.getParameter< real_t >("mortonNumber");
   const real_t relaxationRate  = physicsParameters.getParameter< real_t >("relaxationRate");
   const bool enableWetting     = physicsParameters.getParameter< bool >("enableWetting");
   const real_t contactAngle    = physicsParameters.getParameter< real_t >("contactAngle");
   const uint_t timesteps       = physicsParameters.getParameter< uint_t >("timesteps");

   const CollisionModel_T collisionModel = CollisionModel_T(relaxationRate);
   const real_t viscosity                = collisionModel.viscosity();

   const real_t surfaceTension =
      real_c(std::pow(bondNumber * std::pow(viscosity, 4) / mortonNumber / real_c(tubeDiameter * tubeDiameter), 0.5));

   const real_t gravitationalAccelerationZ =
      mortonNumber * real_c(std::pow(surfaceTension, 3)) / real_c(std::pow(viscosity, 4));

   const Vector3< real_t > force = Vector3< real_t >(real_c(0), real_c(0), -gravitationalAccelerationZ);

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(mortonNumber);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(bondNumber);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableWetting);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(contactAngle);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(timesteps);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(force);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(surfaceTension);

   // read model parameters from parameter file
   const auto modelParameters               = walberlaEnv.config()->getOneBlock("ModelParameters");
   const std::string pdfReconstructionModel = modelParameters.getParameter< std::string >("pdfReconstructionModel");
   const std::string pdfRefillingModel      = modelParameters.getParameter< std::string >("pdfRefillingModel");
   const std::string excessMassDistributionModel =
      modelParameters.getParameter< std::string >("excessMassDistributionModel");
   const std::string curvatureModel          = modelParameters.getParameter< std::string >("curvatureModel");
   const bool useSimpleMassExchange          = modelParameters.getParameter< bool >("useSimpleMassExchange");
   const real_t cellConversionThreshold      = modelParameters.getParameter< real_t >("cellConversionThreshold");
   const real_t cellConversionForceThreshold = modelParameters.getParameter< real_t >("cellConversionForceThreshold");
   const bool enableBubbleModel              = modelParameters.getParameter< bool >("enableBubbleModel");
   const bool enableBubbleSplits             = modelParameters.getParameter< bool >("enableBubbleSplits");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfReconstructionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfRefillingModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(excessMassDistributionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(curvatureModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(useSimpleMassExchange);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionForceThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleSplits);

   // read evaluation parameters from parameter file
   const auto evaluationParameters      = walberlaEnv.config()->getOneBlock("EvaluationParameters");
   const uint_t performanceLogFrequency = evaluationParameters.getParameter< uint_t >("performanceLogFrequency");
   const uint_t evaluationFrequency     = evaluationParameters.getParameter< uint_t >("evaluationFrequency");
   const std::string filename           = evaluationParameters.getParameter< std::string >("filename");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(performanceLogFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evaluationFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(filename);

   // create non-uniform block forest (non-uniformity required for load balancing)
   const std::shared_ptr< StructuredBlockForest > blockForest =
      createNonUniformBlockForest(domainSize, cellsPerBlock, numBlocks, periodicity);

   // add force field
   const BlockDataID forceDensityFieldID =
      field::addToStorage< VectorField_T >(blockForest, "Force field", force, field::fzyx, uint_c(1));

   // create lattice model
   const LatticeModel_T latticeModel = LatticeModel_T(collisionModel, ForceModel_T(forceDensityFieldID));

   // add various fields
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(1.0), field::fzyx, uint_c(2));

   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   const Vector3< real_t > bubbleTopEnd =
      Vector3< real_t >(bubbleBottomEnd[0], bubbleBottomEnd[1], bubbleBottomEnd[2] + bubbleHeight);
   const geometry::Cylinder cylinderBubble(bubbleBottomEnd * tubeDiameter, bubbleTopEnd * tubeDiameter,
                                           bubbleDiameter * tubeDiameter * real_c(0.5));

   bubble_model::addBodyToFillLevelField< geometry::Cylinder >(*blockForest, fillFieldID, cylinderBubble, true);

   // add boundary conditions from config file
   const auto boundaryParameters = walberlaEnv.config()->getOneBlock("BoundaryParameters");
   freeSurfaceBoundaryHandling->initFromConfig(boundaryParameters);
   for (int i = -1; i != boundaryThicknessZ; ++i)
   {
      freeSurfaceBoundaryHandling->setNoSlipAtBorder(stencil::T, cell_idx_c(i));
   }

   // initialize cylindrical domain walls
   const Vector3< real_t > domainCylinderBottomEnd =
      Vector3< real_t >(real_c(domainSize[0]) * real_c(0.5), real_c(domainSize[1]) * real_c(0.5), real_c(0));
   const Vector3< real_t > domainCylinderTopEnd = Vector3< real_t >(
      real_c(domainSize[0]) * real_c(0.5), real_c(domainSize[1]) * real_c(0.5), real_c(domainSize[2]));
   const geometry::Cylinder cylinderTube(domainCylinderBottomEnd, domainCylinderTopEnd, tubeDiameter * real_c(0.5));

   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS(flagFieldIt, flagField, {
         Cell globalCell = flagFieldIt.cell();
         blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt);
         const Vector3< real_t > globalPoint = blockForest->getCellCenter(globalCell);

         if (!geometry::contains(cylinderTube, globalPoint))
         {
            freeSurfaceBoundaryHandling->setNoSlipInCell(globalCell);
         }
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // IMPORTANT REMARK: this must be only called after every solid flag has been set; otherwise, the boundary handling
   // might not detect solid flags correctly
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // communication after initialization
   Communication_T communication(blockForest, flagFieldID, fillFieldID, forceDensityFieldID);
   communication();

   PdfCommunication_T pdfCommunication(blockForest, pdfFieldID);
   pdfCommunication();

   // add bubble model
   std::shared_ptr< bubble_model::BubbleModelBase > bubbleModel = nullptr;
   if (enableBubbleModel)
   {
      const std::shared_ptr< bubble_model::BubbleModel< LatticeModelStencil_T > > bubbleModelDerived =
         std::make_shared< bubble_model::BubbleModel< LatticeModelStencil_T > >(blockForest, enableBubbleSplits);
      bubbleModelDerived->initFromFillLevelField(fillFieldID);

      bubbleModel = std::static_pointer_cast< bubble_model::BubbleModelBase >(bubbleModelDerived);
   }
   else { bubbleModel = std::make_shared< bubble_model::BubbleModelConstantPressure >(real_c(1)); }

   // initialize hydrostatic pressure
   initHydrostaticPressure< PdfField_T >(blockForest, pdfFieldID, force, real_c(domainSize[2]) * real_c(0.5));

   // set density in non-liquid or non-interface cells to 1 (after initializing with hydrostatic pressure)
   setDensityInNonFluidCellsToOne< FlagField_T, PdfField_T >(blockForest, flagInfo, flagFieldID, pdfFieldID);

   // create timeloop
   SweepTimeloop timeloop(blockForest, timesteps);

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with 0 surface tension
   bool computeCurvature = false;
   if (!realIsEqual(surfaceTension, real_c(0), real_c(1e-14))) { computeCurvature = true; }

   // add surface geometry handler
   const SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, curvatureModel, computeCurvature, enableWetting,
      contactAngle);

   geometryHandler.addSweeps(timeloop);

   const ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
   const ConstBlockDataID normalFieldID    = geometryHandler.getConstNormalFieldID();

   // add boundary handling for standard boundaries and free surface boundaries
   const SurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID, normalFieldID, curvatureFieldID,
      freeSurfaceBoundaryHandling, bubbleModel, pdfReconstructionModel, pdfRefillingModel, excessMassDistributionModel,
      relaxationRate, force, surfaceTension, useSimpleMassExchange, cellConversionThreshold,
      cellConversionForceThreshold);

   dynamicsHandler.addSweeps(timeloop);

   // add load balancing
   const LoadBalancer< FlagField_T, CommunicationStencil_T, LatticeModelStencil_T > loadBalancer(
      blockForest, communication, pdfCommunication, bubbleModel, uint_c(50), uint_c(10), uint_c(5),
      loadBalancingFrequency, printLoadBalancingStatistics);
   timeloop.addFuncAfterTimeStep(loadBalancer, "Sweep: load balancing");

   const std::shared_ptr< Vector3< real_t > > centerOfMass =
      std::make_shared< Vector3< real_t > >(Vector3< real_t >(real_c(0)));
   const CenterOfMassComputer< FreeSurfaceBoundaryHandling_T > centerOfMassComputer(
      blockForest, freeSurfaceBoundaryHandling, evaluationFrequency, centerOfMass);
   timeloop.addFuncAfterTimeStep(centerOfMassComputer, "Evaluator: center of mass");

   // add evaluator for total and excessive mass (mass that is currently undistributed)
   const std::shared_ptr< real_t > totalMass  = std::make_shared< real_t >(real_c(0));
   const std::shared_ptr< real_t > excessMass = std::make_shared< real_t >(real_c(0));
   const TotalMassComputer< FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T > totalMassComputer(
      blockForest, freeSurfaceBoundaryHandling, pdfFieldID, fillFieldID, dynamicsHandler.getConstExcessMassFieldID(),
      evaluationFrequency, totalMass, excessMass);
   timeloop.addFuncAfterTimeStep(totalMassComputer, "Evaluator: total mass");

   // add VTK output
   addVTKOutput< LatticeModel_T, FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T, VectorField_T >(
      blockForest, timeloop, walberlaEnv.config(), flagInfo, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID,
      geometryHandler.getCurvatureFieldID(), geometryHandler.getNormalFieldID(),
      geometryHandler.getObstNormalFieldID());

   // add triangle mesh output of free surface
   SurfaceMeshWriter< ScalarField_T, FlagField_T > surfaceMeshWriter(
      blockForest, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs, real_c(1), walberlaEnv.config());
   surfaceMeshWriter(); // write initial mesh
   timeloop.addFuncAfterTimeStep(surfaceMeshWriter, "Writer: surface mesh");

   // add performance logging
   const lbm::PerformanceLogger< FlagField_T > perfLogger(blockForest, flagFieldID, flagIDs::liquidInterfaceFlagIDs,
                                                          performanceLogFrequency);
   timeloop.addFuncAfterTimeStep(perfLogger, "Evaluator: performance logging");

   WcTimingPool timingPool;

   const real_t stoppingHeight = real_c(domainSize[2]) - bubbleHeight * tubeDiameter;

   uint_t timestepOld                = uint_c(0);
   Vector3< real_t > centerOfMassOld = Vector3< real_t >(real_c(0));

   for (uint_t t = uint_c(0); t != timesteps; ++t)
   {
      timeloop.singleStep(timingPool, true);

      if (t % evaluationFrequency == uint_c(0) && t > uint_c(0))
      {
         WALBERLA_ROOT_SECTION()
         {
            const real_t riseVelocity = ((*centerOfMass)[2] - centerOfMassOld[2]) / real_c(t - timestepOld);
            const real_t dragForce    = real_c(4) / real_c(3) * gravitationalAccelerationZ * real_c(bubbleDiameter) /
                                     (riseVelocity * riseVelocity);

            WALBERLA_LOG_DEVEL("time step = " << t << "\n\t\tcenterOfMass = " << *centerOfMass
                                              << "\n\t\triseVelocity = " << riseVelocity
                                              << "\n\t\tdragForce = " << dragForce << "\n\t\ttotalMass = " << *totalMass
                                              << "\n\t\texcessMass = " << *excessMass);

            const std::vector< real_t > resultVector{ (*centerOfMass)[2], riseVelocity, dragForce };

            writeVectorToFile(resultVector, t, filename);

            timestepOld     = t;
            centerOfMassOld = *centerOfMass;
         }
      }

      // stop simulation before bubble hits the top wall
      if ((*centerOfMass)[2] > stoppingHeight) { break; }

      if (t % performanceLogFrequency == uint_c(0) && t > uint_c(0)) { timingPool.logResultOnRoot(); }
   }
   return EXIT_SUCCESS;
}

template< typename T >
void writeVectorToFile(const std::vector< T >& vector, uint_t timestep, const std::string& filename)
{
   std::fstream file;
   file.open(filename, std::fstream::app);

   file << timestep;
   for (const auto i : vector)
   {
      file << "\t" << i;
   }

   file << "\n";
   file.close();
}

} // namespace TaylorBubble
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::TaylorBubble::main(argc, argv); }