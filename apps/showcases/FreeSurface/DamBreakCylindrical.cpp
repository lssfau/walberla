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
//! \file DamBreakCylindrical.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
// This showcase simulates the collapse of a cylindrical liquid column in 3D. Reference experiments are available from
// Martin, Moyce (1952), doi:10.1098/rsta.1952.0006
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "core/Environment.h"
#include "core/math/DistributedSample.h"
#include "core/math/Sample.h"

#include "field/Gather.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/free_surface/LoadBalancing.h"
#include "lbm/free_surface/SurfaceMeshWriter.h"
#include "lbm/free_surface/VtkWriter.h"
#include "lbm/free_surface/bubble_model/Geometry.h"
#include "lbm/free_surface/dynamics/SurfaceDynamicsHandler.h"
#include "lbm/free_surface/surface_geometry/SurfaceGeometryHandler.h"
#include "lbm/free_surface/surface_geometry/Utility.h"
#include "lbm/lattice_model/D3Q19.h"

#include "stencil/D3Q27.h"

namespace walberla
{
namespace free_surface
{
namespace DamBreakCylindrical
{
using ScalarField_T = GhostLayerField< real_t, 1 >;
using VectorField_T = GhostLayerField< Vector3< real_t >, 1 >;

using CollisionModel_T      = lbm::collision_model::SRTField< ScalarField_T >;
using ForceModel_T          = lbm::force_model::GuoField< VectorField_T >;
using LatticeModel_T        = lbm::D3Q19< CollisionModel_T, true, ForceModel_T, 2 >;
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

template< typename T >
void writeNumberVector(const std::vector< T >& numberVector, const uint_t& timestep, const std::string& filename)
{
   std::fstream file;
   file.open(filename, std::fstream::app);

   file << timestep;
   for (const auto number : numberVector)
   {
      file << "\t" << number;
   }

   file << "\n";
   file.close();
}

// get statistical Sample of column width, i.e., radius of the liquid front at the bottom layer (y=0); the center point
// of the cylindrical column is assumed to be constant throughout the simulations
// IMPORTANT REMARK: This implementation is not very efficient, as it gathers a slice of the whole flag field on a
// single process to perform the evaluation.
template< typename FreeSurfaceBoundaryHandling_T >
class columnRadiusEvaluator
{
 public:
   columnRadiusEvaluator(const std::weak_ptr< StructuredBlockForest >& blockForest,
                         const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                         const Vector3< uint_t >& domainSize, const Vector3< real_t >& initialOrigin, uint_t interval,
                         const std::shared_ptr< math::Sample >& columnRadiusSample)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), domainSize_(domainSize),
        initialOrigin_(initialOrigin), columnRadiusSample_(columnRadiusSample), interval_(interval),
        executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % interval_ != uint_c(0)) { return; }

      // gather a slice of the flag field on rank 0 (WARNING: simple, but inefficient)
      std::shared_ptr< FlagField_T > flagFieldGathered = nullptr;
      WALBERLA_ROOT_SECTION()
      {
         flagFieldGathered = std::make_shared< FlagField_T >(domainSize_[0], uint_c(1), domainSize_[2], uint_c(0));
      }
      field::gatherSlice< FlagField_T, FlagField_T >(
         *flagFieldGathered, blockForest, freeSurfaceBoundaryHandling->getFlagFieldID(), 1, cell_idx_c(0), 0);

      WALBERLA_ROOT_SECTION()
      {
         columnRadiusSample_->clear();

         const Cell startCell(cell_idx_c(0), cell_idx_c(0), cell_idx_c(0));

         const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo =
            freeSurfaceBoundaryHandling->getFlagInfo();
         WALBERLA_CHECK(flagInfo.isGas(flagFieldGathered->get(startCell)),
                        "The \"startCell\" in columnRadiusEvaluator's flood fill algorithm must be a gas cell.");

         getRadiusFromFloodFill(flagFieldGathered, startCell, *columnRadiusSample_);
      }
   }

   void getRadiusFromFloodFill(const std::shared_ptr< FlagField_T >& flagFieldGathered, const Cell& startCell,
                               math::Sample& columnRadiusSample)
   {
      WALBERLA_CHECK_EQUAL(startCell.y(), cell_idx_c(0),
                           "columnRadiusEvaluator is meant to be search at the domain's bottom layer only (at y=0).");

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      std::queue< Cell > cellQueue;
      cellQueue.push(startCell);

      while (!cellQueue.empty())
      {
         Cell cell = cellQueue.front();
         cellQueue.pop();

         if (flagInfo.isGas(flagFieldGathered->get(cell)))
         {
            // remove flag such that cell is not detected as gas when searching from neighboring cells
            flagFieldGathered->get(cell) = flag_t(0);

            for (auto dir = stencil::D3Q27::beginNoCenter(); dir != stencil::D3Q27::end(); ++dir)
            {
               // only search at y=0
               if (dir.cy() != 0) { continue; }

               const Cell neighborCell = Cell(cell.x() + dir.cx(), cell.y() + dir.cy(), cell.z() + dir.cz());

               // make sure that the algorithm stays within the field dimensions
               if (neighborCell.x() >= cell_idx_c(0) && neighborCell.x() < cell_idx_c(flagFieldGathered->xSize()) &&
                   neighborCell.y() >= cell_idx_c(0) && neighborCell.y() < cell_idx_c(flagFieldGathered->ySize()) &&
                   neighborCell.z() >= cell_idx_c(0) && neighborCell.z() < cell_idx_c(flagFieldGathered->zSize()))
               {
                  // add neighboring cell to queue
                  cellQueue.push(neighborCell);
               }
            }
         }
         else
         {
            if (flagInfo.isInterface(flagFieldGathered->get(cell)))
            {
               // get center point of this cell; directly use y=0 here
               const Vector3< real_t > cellCenter(real_c(cell.x()) + real_c(0.5), real_c(0),
                                                  real_c(cell.z()) + real_c(0.5));

               const real_t radius = (initialOrigin_ - cellCenter).length();

               columnRadiusSample.castToRealAndInsert(radius);
            } // else: cell is neither gas, nor interface => do nothing
         }
      }
   }

 private:
   std::weak_ptr< StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;
   Vector3< uint_t > domainSize_;
   Vector3< real_t > initialOrigin_;
   std::shared_ptr< math::Sample > columnRadiusSample_;

   uint_t interval_;
   uint_t executionCounter_;
}; // class columnRadiusEvaluator

// get height of residual liquid column, i.e., height of liquid at initial center
template< typename FreeSurfaceBoundaryHandling_T >
class ColumnHeightEvaluator
{
 public:
   ColumnHeightEvaluator(const std::weak_ptr< const StructuredBlockForest >& blockForest,
                         const std::weak_ptr< const FreeSurfaceBoundaryHandling_T >& freeSurfaceBoundaryHandling,
                         const Vector3< uint_t >& domainSize, const Vector3< real_t >& initialOrigin, uint_t interval,
                         const std::shared_ptr< cell_idx_t >& currentColumnHeight)
      : blockForest_(blockForest), freeSurfaceBoundaryHandling_(freeSurfaceBoundaryHandling), domainSize_(domainSize),
        initialOrigin_(initialOrigin), currentColumnHeight_(currentColumnHeight), interval_(interval),
        executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      auto blockForest = blockForest_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blockForest);

      auto freeSurfaceBoundaryHandling = freeSurfaceBoundaryHandling_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(freeSurfaceBoundaryHandling);

      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % interval_ != uint_c(0)) { return; }

      const BlockDataID flagFieldID = freeSurfaceBoundaryHandling->getFlagFieldID();
      const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

      *currentColumnHeight_ = cell_idx_c(0);

      for (auto blockIt = blockForest->begin(); blockIt != blockForest->end(); ++blockIt)
      {
         cell_idx_t maxColumnHeight = cell_idx_c(0);
         bool isInterfaceFound      = false;

         const CellInterval globalSearchInterval(cell_idx_c(initialOrigin_[0]), cell_idx_c(0),
                                                 cell_idx_c(initialOrigin_[2]), cell_idx_c(initialOrigin_[0]),
                                                 cell_idx_c(domainSize_[1]), cell_idx_c(initialOrigin_[2]));

         if (blockForest->getBlockCellBB(*blockIt).overlaps(globalSearchInterval))
         {
            CellInterval localSearchInterval = globalSearchInterval;

            // get intersection of globalSearchInterval and this block's bounding box (both in global coordinates)
            localSearchInterval.intersect(blockForest->getBlockCellBB(*blockIt));

            blockForest->transformGlobalToBlockLocalCellInterval(localSearchInterval, *blockIt);

            const FlagField_T* const flagField = blockIt->template getData< const FlagField_T >(flagFieldID);

            for (auto c = localSearchInterval.begin(); c != localSearchInterval.end(); ++c)
            {
               if (flagInfo.isInterface(flagField->get(*c)))
               {
                  if (c->y() >= maxColumnHeight)
                  {
                     maxColumnHeight  = c->y();
                     isInterfaceFound = true;
                  }
               }
            }

            if (isInterfaceFound)
            {
               // transform local y-coordinate to global coordinate
               Cell localResultCell = Cell(cell_idx_c(0), maxColumnHeight, cell_idx_c(0));
               blockForest->transformBlockLocalToGlobalCell(localResultCell, *blockIt);
               if (localResultCell[1] > *currentColumnHeight_) { *currentColumnHeight_ = localResultCell[1]; }
            }
         }
      }
      mpi::allReduceInplace< cell_idx_t >(*currentColumnHeight_, mpi::MAX);
   }

 private:
   std::weak_ptr< const StructuredBlockForest > blockForest_;
   std::weak_ptr< const FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling_;
   Vector3< uint_t > domainSize_;
   Vector3< real_t > initialOrigin_;
   std::shared_ptr< cell_idx_t > currentColumnHeight_;

   uint_t interval_;
   uint_t executionCounter_;
}; // class ColumnHeightEvaluator

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
   auto domainParameters     = walberlaEnv.config()->getOneBlock("DomainParameters");
   const real_t columnRadius = domainParameters.getParameter< real_t >("columnRadius");
   const real_t columnRatio  = domainParameters.getParameter< real_t >("columnRatio");
   const real_t columnHeight = columnRadius * columnRatio;

   // define domain size
   Vector3< uint_t > domainSize;
   domainSize[0] = uint_c(real_c(12) * columnRadius);
   domainSize[1] = uint_c(real_c(2) * columnHeight);
   domainSize[2] = domainSize[0];

   // compute number of blocks as defined by domainSize and cellsPerBlock
   Vector3< uint_t > numBlocks;
   numBlocks[0] = uint_c(std::ceil(real_c(domainSize[0]) / real_c(cellsPerBlock[0])));
   numBlocks[1] = uint_c(std::ceil(real_c(domainSize[1]) / real_c(cellsPerBlock[1])));
   numBlocks[2] = uint_c(std::ceil(real_c(domainSize[2]) / real_c(cellsPerBlock[2])));

   uint_t numProcesses = uint_c(MPIManager::instance()->numProcesses());
   WALBERLA_CHECK_LESS_EQUAL(numProcesses, numBlocks[0] * numBlocks[1] * numBlocks[2],
                             "The number of MPI processes is greater than the number of blocks as defined by "
                             "\"domainSize/cellsPerBlock\". This would result in unused MPI processes. Either decrease "
                             "the number of MPI processes or increase \"cellsPerBlock\".")

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numProcesses);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellsPerBlock);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSize);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numBlocks);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(columnRadius);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(columnHeight);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(columnRatio);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(periodicity);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(loadBalancingFrequency);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(printLoadBalancingStatistics);

   // read physics parameters from parameter file
   auto physicsParameters = walberlaEnv.config()->getOneBlock("PhysicsParameters");
   const uint_t timesteps = physicsParameters.getParameter< uint_t >("timesteps");

   // Galilei number: Ga = density * force * L^3 / kinematicViscosity^2
   const real_t galileiNumber = physicsParameters.getParameter< real_t >("galileiNumber");

   // Bond (Eötvös) number: Bo = (density_liquid - density_gas) * force * L^2 / surfaceTension
   const real_t bondNumber = physicsParameters.getParameter< real_t >("bondNumber");

   const real_t relaxationRate = physicsParameters.getParameter< real_t >("relaxationRate");
   const real_t viscosity      = (real_c(1) / relaxationRate - real_c(0.5)) / real_c(3);
   const real_t forceY         = galileiNumber * viscosity * viscosity / real_c(std::pow(columnRadius, real_c(3)));
   const Vector3< real_t > force(real_c(0), -forceY, real_c(0));
   const real_t surfaceTension = std::abs(forceY) * columnRadius * columnRadius / bondNumber;

   const bool enableWetting  = physicsParameters.getParameter< bool >("enableWetting");
   const real_t contactAngle = physicsParameters.getParameter< real_t >("contactAngle");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(galileiNumber);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(bondNumber);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(relaxationRate);
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
   const real_t smagorinskyConstant          = modelParameters.getParameter< real_t >("smagorinskyConstant");

   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfReconstructionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfRefillingModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(excessMassDistributionModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(curvatureModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(useSimpleMassExchange);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellConversionForceThreshold);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleModel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(enableBubbleSplits);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(smagorinskyConstant);

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

   // add relaxationRate field (initialized with relaxationRate from parameter file)
   const BlockDataID relaxationRateFieldID = field::addToStorage< ScalarField_T >(
      blockForest, "Relaxation rate field", relaxationRate, field::fzyx, uint_c(1));

   const CollisionModel_T collisionModel = lbm::collision_model::SRTField< ScalarField_T >(relaxationRateFieldID);

   // add force field
   const BlockDataID forceDensityFieldID =
      field::addToStorage< VectorField_T >(blockForest, "Force field", force, field::fzyx, uint_c(1));

   // create lattice model
   const LatticeModel_T latticeModel =
      LatticeModel_T(collisionModel, lbm::force_model::GuoField< VectorField_T >(forceDensityFieldID));

   // add pdf field
   const BlockDataID pdfFieldID = lbm::addPdfFieldToStorage(blockForest, "PDF field", latticeModel, field::fzyx);

   // add fill level field (initialized with 0, i.e., gas everywhere)
   const BlockDataID fillFieldID =
      field::addToStorage< ScalarField_T >(blockForest, "Fill level field", real_c(0.0), field::fzyx, uint_c(2));

   // add boundary handling
   const std::shared_ptr< FreeSurfaceBoundaryHandling_T > freeSurfaceBoundaryHandling =
      std::make_shared< FreeSurfaceBoundaryHandling_T >(blockForest, pdfFieldID, fillFieldID);
   const BlockDataID flagFieldID                                      = freeSurfaceBoundaryHandling->getFlagFieldID();
   const typename FreeSurfaceBoundaryHandling_T::FlagInfo_T& flagInfo = freeSurfaceBoundaryHandling->getFlagInfo();

   const geometry::Cylinder cylinderColumn(
      Vector3< real_t >(real_c(0.5) * real_c(domainSize[0]), real_c(-1), real_c(0.5) * real_c(domainSize[2])),
      Vector3< real_t >(real_c(0.5) * real_c(domainSize[0]), columnHeight, real_c(0.5) * real_c(domainSize[2])),
      columnRadius);
   bubble_model::addBodyToFillLevelField< geometry::Cylinder >(*blockForest, fillFieldID, cylinderColumn, false);

   // initialize boundary conditions from config file
   const auto boundaryParameters = walberlaEnv.config()->getOneBlock("BoundaryParameters");
   freeSurfaceBoundaryHandling->initFromConfig(boundaryParameters);

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
      bubbleModelDerived->setAtmosphere(
         Cell(domainSize[0] - uint_c(1), domainSize[1] - uint_c(1), domainSize[2] - uint_c(1)), real_c(1));

      bubbleModel = std::static_pointer_cast< bubble_model::BubbleModelBase >(bubbleModelDerived);
   }
   else { bubbleModel = std::make_shared< bubble_model::BubbleModelConstantPressure >(real_c(1)); }

   // initialize hydrostatic pressure
   initHydrostaticPressure< PdfField_T >(blockForest, pdfFieldID, force, columnHeight);

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
      cellConversionForceThreshold, relaxationRateFieldID, smagorinskyConstant);

   dynamicsHandler.addSweeps(timeloop);

   // add load balancing
   LoadBalancer< FlagField_T, CommunicationStencil_T, LatticeModelStencil_T > loadBalancer(
      blockForest, communication, pdfCommunication, bubbleModel, uint_c(50), uint_c(10), uint_c(5),
      loadBalancingFrequency, printLoadBalancingStatistics);
   timeloop.addFuncAfterTimeStep(loadBalancer, "Sweep: load balancing");

   // add sweep for evaluating the column height at the origin
   const std::shared_ptr< cell_idx_t > currentColumnHeight = std::make_shared< cell_idx_t >(columnHeight);
   const ColumnHeightEvaluator< FreeSurfaceBoundaryHandling_T > heightEvaluator(
      blockForest, freeSurfaceBoundaryHandling, domainSize,
      Vector3< real_t >(real_c(0.5) * real_c(domainSize[0]), real_c(0), real_c(0.5) * real_c(domainSize[2])),
      evaluationFrequency, currentColumnHeight);
   timeloop.addFuncAfterTimeStep(heightEvaluator, "Evaluator: column height");

   // add sweep for evaluating the column width (distance of front to origin)
   const std::shared_ptr< math::Sample > columnRadiusSample = std::make_shared< math::Sample >();
   const columnRadiusEvaluator< FreeSurfaceBoundaryHandling_T > columnRadiusEvaluator(
      blockForest, freeSurfaceBoundaryHandling, domainSize,
      Vector3< real_t >(real_c(0.5) * real_c(domainSize[0]), real_c(0), real_c(0.5) * real_c(domainSize[2])),
      evaluationFrequency, columnRadiusSample);
   timeloop.addFuncAfterTimeStep(columnRadiusEvaluator, "Evaluator: radius");

   // add VTK output
   addVTKOutput< LatticeModel_T, FreeSurfaceBoundaryHandling_T, PdfField_T, FlagField_T, ScalarField_T, VectorField_T >(
      blockForest, timeloop, walberlaEnv.config(), flagInfo, pdfFieldID, flagFieldID, fillFieldID, forceDensityFieldID,
      geometryHandler.getCurvatureFieldID(), geometryHandler.getNormalFieldID(),
      geometryHandler.getObstNormalFieldID());

   // add triangle mesh output of free surface
   SurfaceMeshWriter< ScalarField_T, FlagField_T > surfaceMeshWriter(
      blockForest, fillFieldID, flagFieldID, flagIDs::liquidInterfaceGasFlagIDs, real_c(0), walberlaEnv.config());
   surfaceMeshWriter(); // write initial mesh
   timeloop.addFuncAfterTimeStep(surfaceMeshWriter, "Writer: surface mesh");

   // add logging for computational performance
   const lbm::PerformanceLogger< FlagField_T > perfLogger(blockForest, flagFieldID, flagIDs::liquidInterfaceFlagIDs,
                                                          performanceLogFrequency);
   timeloop.addFuncAfterTimeStep(perfLogger, "Evaluator: performance logging");

   WcTimingPool timingPool;

   for (uint_t t = uint_c(0); t != timesteps; ++t)
   {
      timeloop.singleStep(timingPool, true);

      if (t % evaluationFrequency == uint_c(0))
      {
         // set initial values
         real_t T              = real_c(0);
         real_t Z_mean         = real_c(1);
         real_t Z_max          = real_c(1);
         real_t Z_min          = real_c(1);
         real_t Z_stdDeviation = real_c(0);
         real_t H              = real_c(1);

         if (!columnRadiusSample->empty())
         {
            // compute dimensionless quantities as defined in paper from Martin and Moyce (1952)
            T              = real_c(t) * std::sqrt(columnRatio * std::abs(force[1]) / columnRadius);
            Z_mean         = columnRadiusSample->mean() / columnRadius;
            Z_max          = columnRadiusSample->max() / columnRadius;
            Z_min          = columnRadiusSample->min() / columnRadius;
            Z_stdDeviation = columnRadiusSample->stdDeviation() / columnRadius;
            H              = real_c(*currentColumnHeight) / columnHeight;
         }

         WALBERLA_LOG_DEVEL_ON_ROOT("time step =" << t);
         WALBERLA_LOG_DEVEL_ON_ROOT("\t\tT = " << T << "\n\t\tZ_mean = " << Z_mean << "\n\t\tZ_max = " << Z_max
                                               << "\n\t\tZ_min = " << Z_min
                                               << "\n\t\tZ_stdDeviation = " << Z_stdDeviation << "\n\t\tH = " << H);

         WALBERLA_ROOT_SECTION()
         {
            const std::vector< real_t > resultVector{ T, Z_mean, Z_max, Z_min, Z_stdDeviation, H };
            writeNumberVector(resultVector, t, filename);
         }

         mpi::broadcastObject(Z_max, 0);

         // simulation is considered converged
         if (Z_max >= real_c(domainSize[0]) * real_c(0.5) / columnRadius - real_c(0.5))
         {
            WALBERLA_LOG_DEVEL_ON_ROOT("Liquid has reached domain borders.");
            break;
         }
      }

      if (t % performanceLogFrequency == uint_c(0) && t > uint_c(0)) { timingPool.logResultOnRoot(); }
   }

   return EXIT_SUCCESS;
}
} // namespace DamBreakCylindrical
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::DamBreakCylindrical::main(argc, argv); }