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
//! \file ZalesakDisk.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// This showcase simulates a slotted disk of gas in a constant rotating velocity field. This benchmark is commonly
// referred to as Zalesak's rotating disk (see doi: 10.1016/0021-9991(79)90051-2).
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"

#include "field/Gather.h"

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
#include "lbm/lattice_model/D2Q9.h"

namespace walberla
{
namespace free_surface
{
namespace ZalesakDisk
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using CollisionModel_T      = lbm::collision_model::SRT;
using ForceModel_T          = lbm::force_model::GuoField< VectorField_T >;
using LatticeModel_T        = lbm::D2Q9< CollisionModel_T, true, ForceModel_T, 2 >;
using LatticeModelStencil_T = LatticeModel_T::Stencil;
using PdfField_T            = lbm::PdfField< LatticeModel_T >;
using PdfCommunication_T    = blockforest::SimpleCommunication< LatticeModelStencil_T >;

// the geometry computations in SurfaceGeometryHandler require meaningful values in the ghost layers in corner
// directions (flag field and fill level field); this holds, even if the lattice model uses a D3Q19 stencil
using CommunicationStencil_T =
   typename std::conditional< LatticeModel_T::Stencil::D == uint_t(2), stencil::D2Q9, stencil::D3Q27 >::type;
using Communication_T = blockforest::SimpleCommunication< CommunicationStencil_T >;

using flag_t                        = uint32_t;
using FlagField_T                   = FlagField< flag_t >;
using FreeSurfaceBoundaryHandling_T = FreeSurfaceBoundaryHandling< LatticeModel_T, FlagField_T, ScalarField_T >;

// function describing the initialization velocity profile (in global cell coordinates)
inline Vector3< real_t > velocityProfile(real_t angularVelocity, Cell globalCell, Vector3< real_t > domainCenter)
{
   // add 0.5 to get Cell's center
   const real_t velocityX = -angularVelocity * ((real_c(globalCell.y()) + real_c(0.5)) - domainCenter[0]);
   const real_t velocityY = angularVelocity * ((real_c(globalCell.x()) + real_c(0.5)) - domainCenter[1]);

   return Vector3< real_t >(velocityX, velocityY, real_c(0));
}

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

   // get domain parameters from parameter file
   auto domainParameters    = walberlaEnv.config()->getOneBlock("DomainParameters");
   const uint_t domainWidth = domainParameters.getParameter< uint_t >("domainWidth");

   const real_t diskRadius = real_c(domainWidth) * real_c(0.15);
   const Vector3< real_t > diskCenter =
      Vector3< real_t >(real_c(domainWidth) * real_c(0.5), real_c(domainWidth) * real_c(0.75), real_c(0.5));
   const real_t diskSlotLength = real_c(0.25) * real_c(domainWidth);
   const real_t diskSlotWidth  = real_c(0.05) * real_c(domainWidth);

   // define domain size
   Vector3< uint_t > domainSize;
   domainSize[0] = domainWidth;
   domainSize[1] = domainWidth;
   domainSize[2] = uint_c(1);
   const Vector3< real_t > domainCenter =
      real_c(0.5) * Vector3< real_t >(real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]));

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
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSize);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numBlocks);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainWidth);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(periodicity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(loadBalancingFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(printLoadBalancingStatistics);

   // get physics parameters from parameter file
   auto physicsParameters = walberlaEnv.config()->getOneBlock("PhysicsParameters");
   const uint_t timesteps = physicsParameters.getParameter< uint_t >("timesteps");

   const real_t relaxationRate           = physicsParameters.getParameter< real_t >("relaxationRate");
   const CollisionModel_T collisionModel = CollisionModel_T(relaxationRate);
   const real_t viscosity                = collisionModel.viscosity();

   const real_t reynoldsNumber = physicsParameters.getParameter< real_t >("reynoldsNumber");
   const real_t angularVelocity =
      reynoldsNumber * viscosity / (real_c(0.5) * real_c(domainWidth) * real_c(domainWidth));
   const Vector3< real_t > force(real_c(0), real_c(0), real_c(0));

   const bool enableWetting  = physicsParameters.getParameter< bool >("enableWetting");
   const real_t contactAngle = physicsParameters.getParameter< real_t >("contactAngle");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(reynoldsNumber);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableWetting);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(contactAngle);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(timesteps);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(viscosity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(force);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(angularVelocity);
   WALBERLA_LOG_DEVEL_ON_ROOT("Timesteps for full rotation " << real_c(2) * math::pi / angularVelocity);

   // read model parameters from parameter file
   const auto modelParameters               = walberlaEnv.config()->getOneBlock("ModelParameters");
   const std::string pdfReconstructionModel = modelParameters.getParameter< std::string >("pdfReconstructionModel");
   const std::string pdfRefillingModel      = modelParameters.getParameter< std::string >("pdfRefillingModel");
   const std::string excessMassDistributionModel =
      modelParameters.getParameter< std::string >("excessMassDistributionModel");
   const std::string curvatureModel          = modelParameters.getParameter< std::string >("curvatureModel");
   const bool enableForceWeighting           = modelParameters.getParameter< bool >("enableForceWeighting");
   const bool useSimpleMassExchange          = modelParameters.getParameter< bool >("useSimpleMassExchange");
   const real_t cellConversionThreshold      = modelParameters.getParameter< real_t >("cellConversionThreshold");
   const real_t cellConversionForceThreshold = modelParameters.getParameter< real_t >("cellConversionForceThreshold");
   const bool enableBubbleModel              = modelParameters.getParameter< bool >("enableBubbleModel");
   const bool enableBubbleSplits             = modelParameters.getParameter< bool >("enableBubbleSplits");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfReconstructionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfRefillingModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(excessMassDistributionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(curvatureModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableForceWeighting);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(useSimpleMassExchange);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionForceThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleSplits);

   // read evaluation parameters from parameter file
   const auto evaluationParameters      = walberlaEnv.config()->getOneBlock("EvaluationParameters");
   const uint_t performanceLogFrequency = evaluationParameters.getParameter< uint_t >("performanceLogFrequency");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(performanceLogFrequency);

   // create non-uniform block forest (non-uniformity required for load balancing)
   const std::shared_ptr< StructuredBlockForest > blockForest =
      createNonUniformBlockForest(domainSize, cellsPerBlock, numBlocks, periodicity);

   // add force field
   const BlockDataID forceFieldID =
      field::addToStorage< VectorField_T >(blockForest, "Force field", force, field::fzyx, uint_c(1));

   // create lattice model
   const LatticeModel_T latticeModel = LatticeModel_T(collisionModel, ForceModel_T(forceFieldID));

   // add pdf field
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

   // add fill level field (initialized with 0, i.e., gas everywhere)
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(1.0), field::fzyx, uint_c(2));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   // initialize the velocity profile
   for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
   {
      PdfField_T* const pdfField   = blockIt->getData< PdfField_T >(pdfFieldID);
      FlagField_T* const flagField = blockIt->getData< FlagField_T >(flagFieldID);

      WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
         //  cell in block-local coordinates
         const Cell localCell = pdfFieldIt.cell();

         // get cell in global coordinates
         Cell globalCell = pdfFieldIt.cell();
         blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

         // set velocity profile
         const Vector3< real_t > initialVelocity = velocityProfile(angularVelocity, globalCell, domainCenter);
         pdfField->setDensityAndVelocity(localCell, initialVelocity, real_c(1));
      }) // WALBERLA_FOR_ALL_CELLS
   }

   // create the disk
   const Vector3< real_t > diskBottomEnd = diskCenter - Vector3< real_t >(real_c(0), real_c(0), -real_c(domainSize[2]));
   const Vector3< real_t > diskTopEnd    = diskCenter - Vector3< real_t >(real_c(0), real_c(0), real_c(domainSize[2]));
   const geometry::Cylinder disk(diskBottomEnd, diskTopEnd, diskRadius);
   bubble_model::addBodyToFillLevelField< geometry::Cylinder >(*blockForest, fillFieldID, disk, true);

   // create the disk's slot
   const Vector3< real_t > slotMinCorner = Vector3< real_t >(diskCenter[0] - diskSlotWidth * real_c(0.5),
                                                             diskCenter[1] - diskRadius, -real_c(domainSize[2]));
   const Vector3< real_t > slotMaxCorner = Vector3< real_t >(
      diskCenter[0] + diskSlotWidth * real_c(0.5), diskCenter[1] - diskRadius + diskSlotLength, real_c(domainSize[2]));
   AABB slotAABB(slotMinCorner, slotMaxCorner);
   bubble_model::addBodyToFillLevelField< AABB >(*blockForest, fillFieldID, slotAABB, false);

   // initialize domain boundary conditions from config file
   const auto boundaryParameters = walberlaEnv.config()->getOneBlock("BoundaryParameters");
   freeSurfaceBoundaryHandling->initFromConfig(boundaryParameters);

   // IMPORTANT REMARK: this must be called only after every solid flag has been set; otherwise, the boundary handling
   // might not detect solid flags correctly
   freeSurfaceBoundaryHandling->initFlagsFromFillLevel();

   // communication after initialization
   Communication_T communication(blockForest, flagFieldID, fillFieldID, forceFieldID);
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

   // set density in non-liquid or non-interface cells to one (after initializing with hydrostatic pressure)
   setDensityInNonFluidCellsToOne< FlagField_T, PdfField_T >(blockForest, flagInfo, flagFieldID, pdfFieldID);

   // create timeloop
   SweepTimeloop timeloop(blockForest, timesteps);

   const real_t surfaceTension = real_c(0);

   // Laplace pressure = 2 * surface tension * curvature; curvature computation is not necessary with zero surface
   // tension
   bool computeCurvature = false;
   if (!realIsEqual(surfaceTension, real_c(0), real_c(1e-14))) { computeCurvature = true; }

   // add surface geometry handler
   const SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, curvatureModel, computeCurvature, enableWetting,
      contactAngle);

   geometryHandler.addSweeps(timeloop);

   // get fields created by surface geometry handler
   const ConstBlockDataID curvatureFieldID = geometryHandler.getConstCurvatureFieldID();
   const ConstBlockDataID normalFieldID    = geometryHandler.getConstNormalFieldID();

   // add boundary handling for standard boundaries and free surface boundaries
   const SurfaceDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, forceFieldID, normalFieldID, curvatureFieldID,
      freeSurfaceBoundaryHandling, bubbleModel, pdfReconstructionModel, pdfRefillingModel, excessMassDistributionModel,
      relaxationRate, force, surfaceTension, enableForceWeighting, useSimpleMassExchange, cellConversionThreshold,
      cellConversionForceThreshold);

   dynamicsHandler.addSweeps(timeloop);

   // add load balancing
   const LoadBalancer< FlagField_T, CommunicationStencil_T, LatticeModelStencil_T > loadBalancer(
      blockForest, communication, pdfCommunication, bubbleModel, uint_c(50), uint_c(10), uint_c(5),
      loadBalancingFrequency, printLoadBalancingStatistics);
   timeloop.addFuncAfterTimeStep(loadBalancer, "Sweep: load balancing");

   // add VTK output
   addVTKOutput< LatticeModel_T, FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T, VectorField_T >(
      blockForest, timeloop, walberlaEnv.config(), flagInfo, pdfFieldID, flagFieldID, fillFieldID, forceFieldID,
      geometryHandler.getCurvatureFieldID(), geometryHandler.getNormalFieldID(),
      geometryHandler.getObstNormalFieldID());

   // add triangle mesh output of free surface
   SurfaceMeshWriter< ScalarField_T, FlagField_T > surfaceMeshWriter(
      blockForest, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs, real_c(0), walberlaEnv.config());
   surfaceMeshWriter(); // write initial mesh
   timeloop.addFuncAfterTimeStep(surfaceMeshWriter, "Writer: surface mesh");

   // add logging for computational performance
   const lbm::PerformanceLogger< FlagField_T > performanceLogger(
      blockForest, flagFieldID, flagIDs::liquidInterfaceFlagIDs, performanceLogFrequency);
   timeloop.addFuncAfterTimeStep(performanceLogger, "Evaluator: performance logging");

   WcTimingPool timingPool;

   for (uint_t t = uint_c(0); t != timesteps; ++t)
   {
      timeloop.singleStep(timingPool, true);

      // set the constant velocity profile
      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         PdfField_T* const pdfField   = blockIt->getData< PdfField_T >(pdfFieldID);
         FlagField_T* const flagField = blockIt->getData< FlagField_T >(flagFieldID);

         WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
            if (flagInfo.isInterface(flagFieldIt) || flagInfo.isLiquid(flagFieldIt))
            {
               //  cell in block-local coordinates
               const Cell localCell = pdfFieldIt.cell();

               // get cell in global coordinates
               Cell globalCell = pdfFieldIt.cell();
               blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

               // set velocity profile
               const Vector3< real_t > initialVelocity = velocityProfile(angularVelocity, globalCell, domainCenter);
               pdfField->setDensityAndVelocity(localCell, initialVelocity, real_c(1));
            }
         }) // WALBERLA_FOR_ALL_CELLS
      }

      if (t % performanceLogFrequency == uint_c(0) && t > uint_c(0)) { timingPool.logResultOnRoot(); }
   }

   return EXIT_SUCCESS;
}

} // namespace ZalesakDisk
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::ZalesakDisk::main(argc, argv); }