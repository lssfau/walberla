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
// This benchmark simulates the advection of a slotted disk of gas in a constant rotating velocity field. The disk
// returns to its initial position, where it should take its initial form. The relative geometrical error of the
// bubble's shape after one rotation is evaluated. There is no LBM flow simulation performed here, it is a test case for
// the FSLBM's mass advection. This benchmark is commonly referred to as Zalesak's rotating disk (see
// doi: 10.1016/0021-9991(79)90051-2). The setup chosen here is identical to the one used by Janssen (see
// doi: 10.1016/j.camwa.2009.08.064).
//======================================================================================================================

#include "core/Environment.h"

#include "field/Gather.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/LoadBalancing.h"
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/TotalMassComputer.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/free_surface/surface_geometry/Utility.h"
#include "lbm/lattice_model/D2Q9.h"

#include "functionality/AdvectionDynamicsHandler.h"
#include "functionality/GeometricalErrorEvaluator.h"

namespace walberla
{
namespace free_surface
{
namespace ZalesakDisk
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

// Lattice model is only created for dummy purposes; no LBM simulation is performed
using CollisionModel_T      = lbm::collision_model::SRT;
using LatticeModel_T        = lbm::D2Q9< CollisionModel_T, true >;
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
inline Vector3< real_t > velocityProfile(real_t angularVelocity, Cell globalCell, const Vector3< real_t >& domainCenter)
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
   auto blockForestParameters            = walberlaEnv.config()->getOneBlock("BlockForestParameters");
   const Vector3< uint_t > cellsPerBlock = blockForestParameters.getParameter< Vector3< uint_t > >("cellsPerBlock");
   const Vector3< bool > periodicity     = blockForestParameters.getParameter< Vector3< bool > >("periodicity");

   // get domain parameters from parameter file
   auto domainParameters    = walberlaEnv.config()->getOneBlock("DomainParameters");
   const uint_t domainWidth = domainParameters.getParameter< uint_t >("domainWidth");

   const real_t diskRadius = real_c(domainWidth) * real_c(0.125);
   const Vector3< real_t > diskCenter =
      Vector3< real_t >(real_c(domainWidth) * real_c(0.5), real_c(domainWidth) * real_c(0.75), real_c(0.5));
   const real_t diskSlotLength = real_c(2) * diskRadius - real_c(0.1) * real_c(domainWidth);
   const real_t diskSlotWidth  = real_c(0.06) * real_c(domainWidth);

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
   const uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
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

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(diskRadius);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(diskCenter);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(diskSlotLength);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(diskSlotWidth);

   // get physics parameters from parameter file
   auto physicsParameters             = walberlaEnv.config()->getOneBlock("PhysicsParameters");
   const uint_t timesteps             = physicsParameters.getParameter< uint_t >("timesteps");
   const uint_t timestepsFullRotation = physicsParameters.getParameter< uint_t >("timestepsFullRotation");

   // compute CFL number
   const real_t dx_SI = real_c(4) / real_c(domainWidth);
   const real_t dt_SI = real_c(12.59652) / real_c(timestepsFullRotation);
   const real_t CFL   = dt_SI / dx_SI; // with velocity_SI = 1
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(CFL);

   // dummy collision model (LBM not simulated in this benchmark)
   const CollisionModel_T collisionModel = CollisionModel_T(real_c(1));

   const real_t angularVelocity = real_c(2) * math::pi / real_c(timestepsFullRotation);

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(timestepsFullRotation);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(timesteps);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(angularVelocity);

   // read model parameters from parameter file
   const auto modelParameters               = walberlaEnv.config()->getOneBlock("ModelParameters");
   const std::string pdfReconstructionModel = modelParameters.getParameter< std::string >("pdfReconstructionModel");
   const std::string excessMassDistributionModel =
      modelParameters.getParameter< std::string >("excessMassDistributionModel");
   const std::string curvatureModel          = modelParameters.getParameter< std::string >("curvatureModel");
   const bool useSimpleMassExchange          = modelParameters.getParameter< bool >("useSimpleMassExchange");
   const real_t cellConversionThreshold      = modelParameters.getParameter< real_t >("cellConversionThreshold");
   const real_t cellConversionForceThreshold = modelParameters.getParameter< real_t >("cellConversionForceThreshold");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfReconstructionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(excessMassDistributionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(curvatureModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(useSimpleMassExchange);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionForceThreshold);

   // read evaluation parameters from parameter file
   const auto evaluationParameters      = walberlaEnv.config()->getOneBlock("EvaluationParameters");
   const uint_t performanceLogFrequency = evaluationParameters.getParameter< uint_t >("performanceLogFrequency");
   const uint_t evaluationFrequency     = evaluationParameters.getParameter< uint_t >("evaluationFrequency");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(performanceLogFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evaluationFrequency);

   // create non-uniform block forest (non-uniformity required for load balancing)
   const std::shared_ptr< StructuredBlockForest > blockForest =
      createNonUniformBlockForest(domainSize, cellsPerBlock, numBlocks, periodicity);

   // create lattice model
   const LatticeModel_T latticeModel = LatticeModel_T(collisionModel);

   // add pdf field
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

   // add fill level field (initialized with 1, i.e., liquid everywhere)
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
      PdfField_T* const pdfField = blockIt->getData< PdfField_T >(pdfFieldID);

      WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, {
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
   Communication_T communication(blockForest, flagFieldID, fillFieldID);
   communication();

   PdfCommunication_T pdfCommunication(blockForest, pdfFieldID);
   pdfCommunication();

   const ConstBlockDataID initialFillFieldID =
      field::addCloneToStorage< ScalarField_T >(blockForest, fillFieldID, "Initial fill level field");

   // add bubble model
   const std::shared_ptr< bubble_model::BubbleModelBase > bubbleModel =
      std::make_shared< bubble_model::BubbleModelConstantPressure >(real_c(1));

   // create timeloop
   SweepTimeloop timeloop(blockForest, timesteps);

   // add surface geometry handler
   const SurfaceGeometryHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > geometryHandler(
      blockForest, freeSurfaceBoundaryHandling, fillFieldID, curvatureModel, false, false, real_c(0));

   geometryHandler.addSweeps(timeloop);

   // get fields created by surface geometry handler
   const ConstBlockDataID normalFieldID = geometryHandler.getConstNormalFieldID();

   // add boundary handling for standard boundaries and free surface boundaries
   const AdvectionDynamicsHandler< LatticeModel_T, FlagField_T, ScalarField_T, VectorField_T > dynamicsHandler(
      blockForest, pdfFieldID, flagFieldID, fillFieldID, normalFieldID, freeSurfaceBoundaryHandling, bubbleModel,
      pdfReconstructionModel, excessMassDistributionModel, useSimpleMassExchange, cellConversionThreshold,
      cellConversionForceThreshold);

   dynamicsHandler.addSweeps(timeloop);

   // add evaluator for geometrical
   const std::shared_ptr< real_t > geometricalError = std::make_shared< real_t >(real_c(0));
   const GeometricalErrorEvaluator< FreeSurfaceBoundaryHandling_T, FlagField_T, ScalarField_T >
      geometricalErrorEvaluator(blockForest, freeSurfaceBoundaryHandling, initialFillFieldID, fillFieldID,
                                evaluationFrequency, geometricalError);
   timeloop.addFuncAfterTimeStep(geometricalErrorEvaluator, "Evaluator: geometrical errors");

   // add evaluator for total and excessive mass (mass that is currently undistributed)
   const std::shared_ptr< real_t > totalMass  = std::make_shared< real_t >(real_c(0));
   const std::shared_ptr< real_t > excessMass = std::make_shared< real_t >(real_c(0));
   const TotalMassComputer< FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T > totalMassComputer(
      blockForest, freeSurfaceBoundaryHandling, pdfFieldID, fillFieldID, dynamicsHandler.getConstExcessMassFieldID(),
      evaluationFrequency, totalMass, excessMass);
   timeloop.addFuncAfterTimeStep(totalMassComputer, "Evaluator: total mass");

   // add VTK output
   addVTKOutput< LatticeModel_T, FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T, VectorField_T >(
      blockForest, timeloop, walberlaEnv.config(), flagInfo, pdfFieldID, flagFieldID, fillFieldID, BlockDataID(),
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

      if (t % evaluationFrequency == uint_c(0))
      {
         WALBERLA_LOG_DEVEL_ON_ROOT("time step = " << t << "\n\t\ttotal mass = " << *totalMass << "\n\t\texcess mass = "
                                                   << *excessMass << "\n\t\tgeometrical error = " << *geometricalError);
      }

      // set the constant velocity profile
      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         PdfField_T* const pdfField   = blockIt->getData< PdfField_T >(pdfFieldID);
         FlagField_T* const flagField = blockIt->getData< FlagField_T >(flagFieldID);

         WALBERLA_FOR_ALL_CELLS(pdfFieldIt, pdfField, flagFieldIt, flagField, {
            const Cell localCell = pdfFieldIt.cell();

            // get cell in global coordinates
            Cell globalCell = pdfFieldIt.cell();
            blockForest->transformBlockLocalToGlobalCell(globalCell, *blockIt, localCell);

            // set velocity profile
            const Vector3< real_t > initialVelocity = velocityProfile(angularVelocity, globalCell, domainCenter);
            pdfField->setDensityAndVelocity(localCell, initialVelocity, real_c(1));
         }) // WALBERLA_FOR_ALL_CELLS
      }

      pdfCommunication();

      if (t % performanceLogFrequency == uint_c(0) && t > uint_c(0)) { timingPool.logResultOnRoot(); }
   }

   return EXIT_SUCCESS;
}

} // namespace ZalesakDisk
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::ZalesakDisk::main(argc, argv); }