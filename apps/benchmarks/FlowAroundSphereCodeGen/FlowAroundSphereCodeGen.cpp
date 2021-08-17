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
//! \file FlowAroundSphereCodeGen.cpp
//! \author Frederik Hennig <frederik.hennig@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
#include "blockforest/all.h"

#include "core/all.h"

#include "domain_decomposition/all.h"

#include "field/all.h"

#include "geometry/all.h"

#include "lbm/inplace_streaming/TimestepTracker.h"
#include "lbm/vtk/QCriterion.h"

#include "python_coupling/CreateConfig.h"
#include "python_coupling/PythonCallback.h"

#include "timeloop/all.h"

#if defined(WALBERLA_BUILD_WITH_CUDA)
#   include "cuda/AddGPUFieldToStorage.h"
#   include "cuda/DeviceSelectMPI.h"
#   include "cuda/HostFieldAllocator.h"
#   include "cuda/NVTX.h"
#   include "cuda/ParallelStreams.h"
#   include "cuda/communication/GPUPackInfo.h"
#   include "cuda/communication/UniformGPUScheme.h"
#endif

// CodeGen includes
#include "FlowAroundSphereCodeGen_InfoHeader.h"

namespace walberla
{
typedef lbm::FlowAroundSphereCodeGen_PackInfoEven PackInfoEven_T;
typedef lbm::FlowAroundSphereCodeGen_PackInfoOdd PackInfoOdd_T;

typedef walberla::uint8_t flag_t;
typedef FlagField< flag_t > FlagField_T;

#if defined(WALBERLA_BUILD_WITH_CUDA)
typedef cuda::GPUField< real_t > GPUField;
#endif

using namespace std::placeholders;

auto pdfFieldAdder = [](IBlock* const block, StructuredBlockStorage* const storage) {
   return new PdfField_T(storage->getNumberOfXCells(*block), storage->getNumberOfYCells(*block),
                         storage->getNumberOfZCells(*block), uint_t(1), field::fzyx,
                         make_shared< field::AllocateAligned< real_t, 64 > >());
};

auto VelocityCallback = [](const Cell& pos, const shared_ptr< StructuredBlockForest >& SbF, IBlock& block,
                           real_t inflow_velocity, const bool constant_inflow = true) {
   if (constant_inflow)
   {
      Vector3< real_t > result(inflow_velocity, 0.0, 0.0);
      return result;
   }
   else
   {
      Cell globalCell;
      CellInterval domain = SbF->getDomainCellBB();
      real_t h_y          = real_c(domain.ySize());
      real_t h_z          = real_c(domain.zSize());
      SbF->transformBlockLocalToGlobalCell(globalCell, block, pos);

      real_t y1 = real_c(globalCell[1] - (h_y / 2.0 - 0.5));
      real_t z1 = real_c(globalCell[2] - (h_z / 2.0 - 0.5));

      real_t u = (inflow_velocity * real_c(16)) / (h_y * h_y * h_z * h_z) * (h_y / real_c(2.0) - y1) *
                 (h_y / real_c(2.0) + y1) * (h_z / real_c(2.0) - z1) * (h_z / real_c(2.0) + z1);

      Vector3< real_t > result(u, 0.0, 0.0);
      return result;
   }
};

class AlternatingBeforeFunction
{
 public:
   typedef std::function< void() > BeforeFunction;

   AlternatingBeforeFunction(BeforeFunction evenFunc, BeforeFunction oddFunc,
                             std::shared_ptr< lbm::TimestepTracker >& tracker)
      : tracker_(tracker), funcs_{ evenFunc, oddFunc } {};

   void operator()() { funcs_[tracker_->getCounter()](); }

 private:
   std::shared_ptr< lbm::TimestepTracker > tracker_;
   std::vector< BeforeFunction > funcs_;
};

class Filter
{
 public:
   explicit Filter(Vector3< uint_t > numberOfCells) : numberOfCells_(numberOfCells) {}

   void operator()(const IBlock& /*block*/) {}

   bool operator()(const cell_idx_t x, const cell_idx_t y, const cell_idx_t z) const
   {
      return x >= -1 && x <= cell_idx_t(numberOfCells_[0]) && y >= -1 && y <= cell_idx_t(numberOfCells_[1]) &&
             z >= -1 && z <= cell_idx_t(numberOfCells_[2]);
   }

 private:
   Vector3< uint_t > numberOfCells_;
};

using FluidFilter_T = Filter;

int main(int argc, char** argv)
{
   walberla::Environment walberlaEnv(argc, argv);
#if defined(WALBERLA_BUILD_WITH_CUDA)
   cuda::selectDeviceBasedOnMpiRank();
#endif

   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      WALBERLA_MPI_WORLD_BARRIER()

      auto config = *cfg;
      logging::configureLogging(config);
      auto blocks = blockforest::createUniformBlockGridFromConfig(config);

      // read parameters
      Vector3< uint_t > cellsPerBlock =
         config->getBlock("DomainSetup").getParameter< Vector3< uint_t > >("cellsPerBlock");
      auto parameters = config->getOneBlock("Parameters");

      const uint_t timesteps       = parameters.getParameter< uint_t >("timesteps", uint_c(10));
      const real_t omega           = parameters.getParameter< real_t >("omega", real_t(1.9));
      const real_t u_max           = parameters.getParameter< real_t >("u_max", real_t(0.05));
      const real_t reynolds_number = parameters.getParameter< real_t >("reynolds_number", real_t(1000));
      const uint_t diameter_sphere = parameters.getParameter< uint_t >("diameter_sphere", uint_t(5));
      const bool constant_inflow = parameters.getParameter< bool >("constant_inflow", true);

      const double remainingTimeLoggerFrequency =
         parameters.getParameter< double >("remainingTimeLoggerFrequency", 3.0); // in seconds

      // create fields
      BlockDataID pdfFieldID     = blocks->addStructuredBlockData< PdfField_T >(pdfFieldAdder, "PDFs");
      BlockDataID velFieldID     = field::addToStorage< VelocityField_T >(blocks, "velocity", real_t(0), field::fzyx);
      BlockDataID densityFieldID = field::addToStorage< ScalarField_T >(blocks, "density", real_t(0), field::fzyx);

#if defined(WALBERLA_BUILD_WITH_CUDA)
      BlockDataID pdfFieldIDGPU = cuda::addGPUFieldToStorage< PdfField_T >(blocks, pdfFieldID, "PDFs on GPU", true);
      BlockDataID velFieldIDGPU =
         cuda::addGPUFieldToStorage< VelocityField_T >(blocks, velFieldID, "velocity on GPU", true);
      BlockDataID densityFieldIDGPU =
         cuda::addGPUFieldToStorage< ScalarField_T >(blocks, densityFieldID, "density on GPU", true);
#endif

      BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

      // initialise all PDFs
#if defined(WALBERLA_BUILD_WITH_CUDA)
      pystencils::FlowAroundSphereCodeGen_MacroSetter setterSweep(pdfFieldIDGPU, velFieldIDGPU);
      for (auto& block : *blocks)
         setterSweep(&block);
      cuda::fieldCpy< PdfField_T, GPUField >(blocks, pdfFieldID, pdfFieldIDGPU);
#else
      pystencils::FlowAroundSphereCodeGen_MacroSetter setterSweep(pdfFieldID, velFieldID);
      for (auto& block : *blocks)
         setterSweep(&block);
#endif
      // Create communication

#if defined(WALBERLA_BUILD_WITH_CUDA)
      // This way of using alternating pack infos is temporary and will soon be replaced
      // by something more straight-forward

      cuda::communication::UniformGPUScheme< Stencil_T > comEven(blocks, false);
      comEven.addPackInfo(make_shared< PackInfoEven_T >(pdfFieldIDGPU));
      auto evenComm = std::function< void() >([&]() { comEven.communicate(nullptr); });

      cuda::communication::UniformGPUScheme< Stencil_T > comODD(blocks, false);
      comODD.addPackInfo(make_shared< PackInfoOdd_T >(pdfFieldIDGPU));
      auto oddComm = std::function< void() >([&]() { comODD.communicate(nullptr); });
#else
      blockforest::communication::UniformBufferedScheme< Stencil_T > evenComm(blocks);
      evenComm.addPackInfo(make_shared< PackInfoEven_T >(pdfFieldID));

      blockforest::communication::UniformBufferedScheme< Stencil_T > oddComm(blocks);
      oddComm.addPackInfo(make_shared< PackInfoOdd_T >(pdfFieldID));
#endif

      // create and initialize boundary handling
      const FlagUID fluidFlagUID("Fluid");

      auto boundariesConfig = config->getOneBlock("Boundaries");

      std::function< Vector3< real_t >(const Cell&, const shared_ptr< StructuredBlockForest >&, IBlock&) >
         velocity_initialisation = std::bind(VelocityCallback, _1, _2, _3, u_max, constant_inflow);

#if defined(WALBERLA_BUILD_WITH_CUDA)
      lbm::FlowAroundSphereCodeGen_UBB ubb(blocks, pdfFieldIDGPU, velocity_initialisation);
      lbm::FlowAroundSphereCodeGen_NoSlip noSlip(blocks, pdfFieldIDGPU);
      lbm::FlowAroundSphereCodeGen_Outflow outflow(blocks, pdfFieldIDGPU, pdfFieldID);

      lbm::FlowAroundSphereCodeGen_LbSweep lbSweep(densityFieldIDGPU, pdfFieldIDGPU, velFieldIDGPU, omega);
#else
      lbm::FlowAroundSphereCodeGen_UBB ubb(blocks, pdfFieldID, velocity_initialisation);
      lbm::FlowAroundSphereCodeGen_NoSlip noSlip(blocks, pdfFieldID);
      lbm::FlowAroundSphereCodeGen_Outflow outflow(blocks, pdfFieldID);

      lbm::FlowAroundSphereCodeGen_LbSweep lbSweep(densityFieldID, pdfFieldID, velFieldID, omega);
#endif

      geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, boundariesConfig);
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

      ubb.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("UBB"), fluidFlagUID);
      noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("NoSlip"), fluidFlagUID);
      outflow.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("Outflow"), fluidFlagUID);

      // create time loop
      SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

      // Timestep Tracking and Sweeps
      auto tracker = make_shared< lbm::TimestepTracker >(0);

      AlternatingBeforeFunction communication(evenComm, oddComm, tracker);

      // add LBM sweep and communication to time loop
      timeloop.add() << BeforeFunction(communication, "communication") << Sweep(ubb.getSweep(tracker), "ubb boundary");
      timeloop.add() << Sweep(outflow.getSweep(tracker), "outflow boundary");
      timeloop.add() << Sweep(noSlip.getSweep(tracker), "noSlip boundary");
      timeloop.add() << BeforeFunction(tracker->getAdvancementFunction(), "Timestep Advancement")
                     << Sweep(lbSweep.getSweep(tracker), "LB update rule");

      // LBM stability check
      timeloop.addFuncAfterTimeStep(makeSharedFunctor(field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                       config, blocks, pdfFieldID, flagFieldId, fluidFlagUID)),
                                    "LBM stability check");

      // log remaining time
      timeloop.addFuncAfterTimeStep(
         timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), remainingTimeLoggerFrequency),
         "remaining time logger");

      // add VTK output to time loop
      uint_t vtkWriteFrequency = parameters.getParameter< uint_t >("vtkWriteFrequency", 0);
      if (vtkWriteFrequency > 0)
      {
         auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "vtk", vtkWriteFrequency, 0, false, "vtk_out",
                                                         "simulation_step", false, true, true, false, 0);

#if defined(WALBERLA_BUILD_WITH_CUDA)
         vtkOutput->addBeforeFunction([&]() {
            cuda::fieldCpy< VelocityField_T, GPUField >(blocks, velFieldID, velFieldIDGPU);
            cuda::fieldCpy< ScalarField_T, GPUField >(blocks, densityFieldID, densityFieldIDGPU);
         });
#endif
         auto velWriter     = make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "velocity");
         auto densityWriter = make_shared< field::VTKWriter< ScalarField_T > >(densityFieldID, "density");
         FluidFilter_T filter(cellsPerBlock);

         auto QCriterionWriter = make_shared< lbm::QCriterionVTKWriter< VelocityField_T, FluidFilter_T > >(
            blocks, filter, velFieldID, "Q-Criterion");

         vtkOutput->addCellDataWriter(velWriter);
         vtkOutput->addCellDataWriter(densityWriter);
         vtkOutput->addCellDataWriter(QCriterionWriter);

         timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
      }

      WcTimer simTimer;

      WALBERLA_LOG_INFO_ON_ROOT("Simulating flow around sphere:"
                                "\n timesteps:               "
                                << timesteps << "\n reynolds number:         " << reynolds_number
                                << "\n relaxation rate:         " << omega << "\n maximum inflow velocity: " << u_max
                                << "\n diameter_sphere:         " << diameter_sphere)

      simTimer.start();
      timeloop.run();
      simTimer.end();
      WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")
      auto time            = simTimer.last();
      auto nrOfCells       = real_c(cellsPerBlock[0] * cellsPerBlock[1] * cellsPerBlock[2]);
      auto mlupsPerProcess = nrOfCells * real_c(timesteps) / time * 1e-6;
      WALBERLA_LOG_RESULT_ON_ROOT("MLUPS per process " << mlupsPerProcess)
      WALBERLA_LOG_RESULT_ON_ROOT("Time per time step " << time / real_c(timesteps))
   }

   return EXIT_SUCCESS;
}

} // namespace walberla

int main(int argc, char** argv) { walberla::main(argc, argv); }
