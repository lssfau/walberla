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
//! \file TurbulentChannel.cpp
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include <memory>
#include <cmath>
#include <string>
#include <iostream>

#include <blockforest/all.h>
#include <core/all.h>
#include <domain_decomposition/all.h>
#include <field/all.h>
#include <field/vtk/VTKWriter.h>
#include <geometry/all.h>
#include <timeloop/all.h>
#include <lbm/all.h>

//    Codegen Includes
#include "CodegenIncludes.h"

namespace walberla {
   ///////////////////////
   /// Typedef Aliases ///
   ///////////////////////

   using Stencil_T = codegen::Stencil_T;

   // Communication Pack Info
   using PackInfo_T = pystencils::TurbulentChannel_PackInfo;

   // PDF field type
   using PdfField_T = field::GhostLayerField< real_t, Stencil_T::Size >;

   // Field Types
   using ScalarField_T = field::GhostLayerField< real_t, 1 >;
   using VectorField_T = field::GhostLayerField< real_t, Stencil_T::D >;
   using TensorField_T = field::GhostLayerField< real_t, Stencil_T::D*Stencil_T::D >;

   using Setter_T = pystencils::TurbulentChannel_Setter;

   using StreamCollideSweep_T = pystencils::TurbulentChannel_Sweep;
   using WelfordSweep_T = pystencils::TurbulentChannel_Welford;
   using TKEWelfordSweep_T = pystencils::TurbulentChannel_Welford_TKE_SGS;

   using TkeSgsWriter_T = pystencils::TurbulentChannel_TKE_SGS_Writer;

   // Boundary Handling
   using flag_t = uint8_t;
   using FlagField_T = FlagField< flag_t >;
   using NoSlip_T = lbm::TurbulentChannel_NoSlip;
   using FreeSlip_top_T = lbm::TurbulentChannel_FreeSlip_top;
   using WFB_bottom_T = lbm::TurbulentChannel_WFB_bottom;
   using WFB_top_T = lbm::TurbulentChannel_WFB_top;

   /// DEAN CORRELATIONS

   namespace dean_correlation {

      real_t calculateFrictionReynoldsNumber(const real_t reynoldsBulk) {
         return std::pow(0.073_r / 8_r, 1_r / 2_r) * std::pow(reynoldsBulk, 7_r / 8_r);
      }

      real_t calculateBulkReynoldsNumber(const real_t reynoldsFriction) {
         return std::pow(8_r / 0.073_r, 4_r / 7_r) * std::pow(reynoldsFriction, 8_r / 7_r);
      }
   } // namespace dean_correlation


   /// VELOCITY FIELD INITIALISATION

   /*
    * Initialises the velocity field with a logarithmic profile and sinusoidal perturbations to trigger turbulence.
    * This initialisation is provided by Henrik Asmuth.
    */
   template<typename VelocityField_T>
   void setVelocityFieldsAsmuth( const std::weak_ptr<StructuredBlockStorage>& forest,
                                 const BlockDataID & velocityFieldId, const BlockDataID & meanVelocityFieldId,
                                 const real_t frictionVelocity, const uint_t channel_half_width,
                                 const real_t B, const real_t kappa, const real_t viscosity,
                                 const uint_t wallAxis, const uint_t flowAxis ) {

      auto blocks = forest.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks)

      const auto domainSize = blocks->getDomain().max();
      const auto delta = real_c(channel_half_width);
      const auto remAxis = 3 - wallAxis - flowAxis;

      for( auto block = blocks->begin(); block != blocks->end(); ++block ) {

         auto * velocityField = block->template getData<VelocityField_T>(velocityFieldId);
         WALBERLA_CHECK_NOT_NULLPTR(velocityField)

         auto * meanVelocityField = block->template getData<VelocityField_T>(meanVelocityFieldId);
         WALBERLA_CHECK_NOT_NULLPTR(meanVelocityField)

         const auto ci = velocityField->xyzSizeWithGhostLayer();
         for(auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt) {

            Cell globalCell(*cellIt);
            blocks->transformBlockLocalToGlobalCell(globalCell, *block);
            Vector3<real_t> cellCenter;
            blocks->getCellCenter(cellCenter, globalCell);

            const auto y = cellCenter[wallAxis];
            const auto rel_x = cellCenter[flowAxis] / domainSize[flowAxis];
            const auto rel_z = cellCenter[remAxis] / domainSize[remAxis];

            const real_t pos = std::max(delta - std::abs(y - delta - 1_r), 0.05_r);
            const auto rel_y = pos / delta;

            auto initialVel = frictionVelocity * (std::log(frictionVelocity * pos / viscosity) / kappa + B);

            Vector3<real_t> vel;
            vel[flowAxis] = initialVel;

            vel[remAxis] = 2_r * frictionVelocity / kappa * std::sin(math::pi * 16_r * rel_x) *
                           std::sin(math::pi * 8_r * rel_y) / (std::pow(rel_y, 2_r) + 1_r);

            vel[wallAxis] = 8_r * frictionVelocity / kappa *
                            (std::sin(math::pi * 8_r * rel_z) * std::sin(math::pi * 8_r * rel_y) +
                             std::sin(math::pi * 8_r * rel_x)) / (std::pow(0.5_r * delta - pos, 2_r) + 1_r);

            for(uint_t d = 0; d < 3; ++d) {
               velocityField->get(*cellIt, d) = vel[d];
               meanVelocityField->get(*cellIt, d) = vel[d];
            }

         }
      }

   } // function setVelocityFieldsHenrik

   /// SIMULATION PARAMETERS

   struct SimulationParameters {

      SimulationParameters(const Config::BlockHandle & config)
      {
         channelHalfWidth = config.getParameter<uint_t>("channel_half_width");
         fullChannel = config.getParameter<bool>("full_channel", false);

         /// TARGET QUANTITIES

         targetFrictionReynolds = config.getParameter<real_t>("target_friction_reynolds");
         targetBulkVelocity = config.getParameter<real_t>("target_bulk_velocity", 0.05_r);

         targetBulkReynolds = dean_correlation::calculateBulkReynoldsNumber(targetFrictionReynolds);
         viscosity = 2_r * real_c(channelHalfWidth) * targetBulkVelocity / targetBulkReynolds;
         targetFrictionVelocity = targetFrictionReynolds * viscosity / real_c(channelHalfWidth);

         /// TIMESTEPS

         timesteps = config.getParameter<uint_t>("timesteps", 0);
         const uint_t turnoverPeriods = config.getParameter<uint_t>("turnover_periods", 0);

         WALBERLA_ASSERT((timesteps != 0) != (turnoverPeriods != 0),
                         "Either timesteps OR turnover periods must be given.")

         if(turnoverPeriods != 0) {
            // turnover period defined by T = delta / u_tau
            timesteps = turnoverPeriods * uint_c((real_c(channelHalfWidth) / targetFrictionVelocity));
         }

         /// DOMAIN DEFINITIONS

         // obtained from codegen script -> adapt there
         wallAxis = codegen::wallAxis;
         flowAxis = codegen::flowAxis;

         uint_t sizeFlowAxis      = config.getParameter<uint_t>("size_flow_axis", 0);
         uint_t sizeRemainingAxis = config.getParameter<uint_t>("size_remaining_axis", 0);

         WALBERLA_ASSERT_NOT_IDENTICAL(wallAxis, flowAxis, "Wall and flow axis must be different.")

         const auto sizeFactor = channelHalfWidth / uint_t(10);
         if( !sizeFlowAxis) sizeFlowAxis = sizeFactor * 64;
         if( !sizeRemainingAxis) sizeRemainingAxis = sizeFactor * 32;

         domainSize[wallAxis] = fullChannel ? 2 * channelHalfWidth : channelHalfWidth;
         domainSize[flowAxis] = sizeFlowAxis;
         domainSize[3- wallAxis - flowAxis] = sizeRemainingAxis;

         periodicity[wallAxis] = false;

         boundaryCondition = config.getParameter<std::string>("wall_boundary_condition", "WFB");

         /// OUTPUT

         auto tsPerPeriod = uint_c((real_c(channelHalfWidth) / targetFrictionVelocity));

         vtkFrequency = config.getParameter<uint_t>("vtk_frequency", 0);
         vtkStart = config.getParameter<uint_t>("vtk_start", 0);
         plotFrequency = config.getParameter<uint_t>("plot_frequency", 0);
         plotStart = config.getParameter<uint_t>("plot_start", 0);

         // vtk start
         vtkStart = config.getParameter<uint_t>("vtk_start_timesteps", 0);
         const uint_t vtkStartPeriods = config.getParameter<uint_t>("vtk_start_periods", 0);

         if(vtkStart || vtkStartPeriods) {
            WALBERLA_ASSERT((vtkStart != 0) != (vtkStartPeriods != 0),
                            "VTK start must be given in timesteps OR periods, not both.")
         }

         if(vtkStartPeriods != 0) {
            // turnover period defined by T = delta / u_tau
            vtkStart = vtkStartPeriods * tsPerPeriod;
         }

         // plot start
         plotStart = config.getParameter<uint_t>("plot_start_timesteps", 0);
         const uint_t plotStartPeriods = config.getParameter<uint_t>("plot_start_periods", 0);

         if(plotStart || plotStartPeriods) {
            WALBERLA_ASSERT((plotStart != 0) != (plotStartPeriods != 0),
                            "Plotting start must be given in timesteps OR periods, not both.")
         }

         if(plotStartPeriods != 0) {
            // turnover period defined by T = delta / u_tau
            plotStart = plotStartPeriods * tsPerPeriod;
         }

         // frequencies
         if(plotFrequency) {
            timesteps = uint_c(std::ceil(real_c(timesteps) / real_c(plotFrequency))) * plotFrequency;
         }

         // sampling start & interval
         samplingStart = config.getParameter<uint_t>("sampling_start_timesteps", 0);
         const uint_t samplingStartPeriods = config.getParameter<uint_t>("sampling_start_periods", 0);

         if(samplingStart || samplingStartPeriods) {
            WALBERLA_ASSERT((samplingStart != 0) != (samplingStartPeriods != 0),
                            "Sampling start must be given in timesteps OR periods, not both.")
         }

         if(samplingStartPeriods != 0) {
            // turnover period defined by T = delta / u_tau
            samplingStart = samplingStartPeriods * tsPerPeriod;
         }

         samplingInterval = config.getParameter<uint_t>("sampling_interval_timesteps", 0);
         const uint_t samplingIntervalPeriods = config.getParameter<uint_t>("sampling_interval_periods", 0);

         if(samplingInterval || samplingIntervalPeriods) {
            WALBERLA_ASSERT((samplingInterval != 0) != (samplingIntervalPeriods != 0),
                            "Sampling start must be given in timesteps OR periods, not both.")
         }

         if(samplingStartPeriods != 0) {
            // turnover period defined by T = delta / u_tau
            samplingInterval = samplingIntervalPeriods * tsPerPeriod;
         }

         timesteps += 1;

      }

      uint_t channelHalfWidth{};
      bool fullChannel{};
      Vector3<uint_t> domainSize{};
      Vector3<uint_t> periodicity{true};

      real_t targetFrictionReynolds{};
      real_t targetBulkReynolds{};
      real_t targetFrictionVelocity{};
      real_t targetBulkVelocity{};

      real_t viscosity{};
      const real_t density{1.0};

      uint_t timesteps{};

      std::string boundaryCondition{};

      uint_t wallAxis{};
      uint_t flowAxis{};

      /// output
      uint_t vtkFrequency{};
      uint_t vtkStart{};
      uint_t plotFrequency{};
      uint_t plotStart{};
      uint_t samplingStart{};
      uint_t samplingInterval{};

   };


   namespace boundaries {
      void createBoundaryConfig(const SimulationParameters & parameters, Config::Block & boundaryBlock) {

         auto & bottomWall = boundaryBlock.createBlock("Border");
         bottomWall.addParameter("direction", stencil::dirToString[stencil::directionFromAxis(parameters.wallAxis, true)]);
         bottomWall.addParameter("walldistance", "-1");
         if(parameters.boundaryCondition == "NoSlip") {
            bottomWall.addParameter("flag", "NoSlip");
         } else if(parameters.boundaryCondition == "WFB") {
            bottomWall.addParameter("flag", "WFB_bottom");
         }

         auto & topWall = boundaryBlock.createBlock("Border");
         topWall.addParameter("direction", stencil::dirToString[stencil::directionFromAxis(parameters.wallAxis, false)]);
         topWall.addParameter("walldistance", "-1");
         if(parameters.fullChannel) {
            if (parameters.boundaryCondition == "NoSlip") {
               topWall.addParameter("flag", "NoSlip");
            } else if (parameters.boundaryCondition == "WFB") {
               topWall.addParameter("flag", "WFB_top");
            }
         } else {
            topWall.addParameter("flag", "FreeSlip");
         }

      }
   }

   /// BULK VELOCITY CALCULATION
   template< typename VelocityField_T >
   class ForceCalculator {

    public:
      ForceCalculator(const std::weak_ptr<StructuredBlockStorage> & blocks, const BlockDataID meanVelocityId,
                      const SimulationParameters & parameter)
         : blocks_(blocks), meanVelocityId_(meanVelocityId), channelHalfWidth_(real_c(parameter.channelHalfWidth)),
           targetBulkVelocity_(parameter.targetBulkVelocity), targetFrictionVelocity_(parameter.targetFrictionVelocity)
      {
         const auto & domainSize = parameter.domainSize;

         Cell maxCell;
         maxCell[parameter.wallAxis] = int_c(parameter.channelHalfWidth) - 1;
         maxCell[flowDirection_] = int_c(domainSize[flowDirection_]) - 1;
         const auto remainingIdx = 3 - parameter.wallAxis - flowDirection_;
         maxCell[remainingIdx] = int_c(domainSize[remainingIdx]) - 1;
         ci_ = CellInterval(Cell{}, maxCell);

         numCells_ = real_c(parameter.channelHalfWidth * domainSize[flowDirection_] * domainSize[remainingIdx]);
      }

      real_t bulkVelocity() const { return bulkVelocity_; }
      void setBulkVelocity(const real_t bulkVelocity) { bulkVelocity_ = bulkVelocity; }

      void calculateBulkVelocity() {

         // reset bulk velocity
         bulkVelocity_ = 0_r;

         auto blocks = blocks_.lock();
         WALBERLA_CHECK_NOT_NULLPTR(blocks)

         for( auto block = blocks->begin(); block != blocks->end(); ++block) {

            auto * meanVelocityField = block->template getData<VelocityField_T>(meanVelocityId_);
            WALBERLA_CHECK_NOT_NULLPTR(meanVelocityField)

            auto fieldSize = meanVelocityField->xyzSize();
            CellInterval localCi;
            blocks->transformGlobalToBlockLocalCellInterval(localCi, *block, ci_);
            fieldSize.intersect(localCi);

            auto * slicedField = meanVelocityField->getSlicedField(fieldSize);
            WALBERLA_CHECK_NOT_NULLPTR(meanVelocityField)

            for(auto fieldIt = slicedField->beginXYZ(); fieldIt != slicedField->end(); ++fieldIt) {
               const auto localMean = fieldIt[flowDirection_];
               bulkVelocity_ += localMean;
            }

         }

         mpi::allReduceInplace< real_t >(bulkVelocity_, mpi::SUM);
         bulkVelocity_ /= numCells_;

      }

      real_t calculateDrivingForce() const {

         // forcing term as in Malaspinas (2014) "Wall model for large-eddy simulation based on the lattice Boltzmann method"
         const auto force = targetFrictionVelocity_ * targetFrictionVelocity_ / channelHalfWidth_
                            + (targetBulkVelocity_ - bulkVelocity_) * targetBulkVelocity_ / channelHalfWidth_;

         return force;
      }

    private:
      const std::weak_ptr<StructuredBlockStorage> blocks_{};
      const BlockDataID meanVelocityId_{};

      const uint_t flowDirection_{};
      const real_t channelHalfWidth_{};
      const real_t targetBulkVelocity_{};
      const real_t targetFrictionVelocity_{};

      CellInterval ci_{};

      real_t numCells_{};
      real_t bulkVelocity_{};
   };

   template< typename Welford_T >
   class TurbulentChannelPlotter {

    public:
      TurbulentChannelPlotter(SimulationParameters const * const parameters, Timeloop * const timeloop,
                              ForceCalculator<VectorField_T> const * const forceCalculator,
                              const std::weak_ptr<StructuredBlockStorage> & blocks,
                              const BlockDataID velocityFieldId, const BlockDataID meanVelocityFieldId,
                              const BlockDataID meanTkeSGSFieldId, Welford_T * velocityWelford,
                              const bool separateFile = false)
         : parameters_(parameters), forceCalculator_(forceCalculator), timeloop_(timeloop), blocks_(blocks),
           velocityWelford_(velocityWelford), velocityFieldId_(velocityFieldId), meanVelocityFieldId_(meanVelocityFieldId),
           meanTkeSGSFieldId_(meanTkeSGSFieldId), plotFrequency_(parameters->plotFrequency), plotStart_(parameters->plotStart),
           separateFiles_(separateFile)
      {
         if(!plotFrequency_)
            return;

         // prepare output folder
         const filesystem::path path(baseFolder_);
         std::string fileSuffix = parameters->boundaryCondition + "_";

         if(parameters->fullChannel)
            fileSuffix += "full_D";
         else
            fileSuffix += "half_D";

         fileSuffix += std::to_string(parameters->channelHalfWidth) + "_Re" +
                       std::to_string(int(parameters->targetFrictionReynolds)) ;

         velocityProfilesFilePath_ = path / ("velocity_profiles_" + fileSuffix);
         forcingDataFilePath_ = path / ("forcing_data_" + fileSuffix + "_t" +
                                        std::to_string(parameters->timesteps-1) + ".txt");

         WALBERLA_ROOT_SECTION() {
            // create directory if not existent; empty if existent
            if( !filesystem::exists(path) ) {
               filesystem::create_directories(path);
            } else {
               for (const auto& entry : filesystem::directory_iterator(path))
                  std::filesystem::remove_all(entry.path());
            }
         }

         // write force header
         std::ofstream os (forcingDataFilePath_, std::ios::out | std::ios::trunc);
         if(os.is_open()) {
            os << "# timestep\t bulk_velocity\t driving_force\n";
            os.close();
         } else {
            WALBERLA_ABORT("Could not open forcing data file.")
         }

      }

      void operator()() {

         const auto ts = timeloop_->getCurrentTimeStep();
         if(ts < plotStart_)
            return;

         if(!plotFrequency_ || (ts % plotFrequency_))
            return;

         const auto channelHalfWidth = real_c(parameters_->channelHalfWidth);
         const auto bulkVelocity = forceCalculator_->bulkVelocity();

         /// write force data

         WALBERLA_ROOT_SECTION() {
            std::ofstream forceOS(forcingDataFilePath_, std::ios::out | std::ios::app);
            if (forceOS.is_open())
            {
               forceOS << ts << "\t" << bulkVelocity << "\t" << forceCalculator_->calculateDrivingForce() << "\n";
               forceOS.close();
            }
         }

         /// write velocity data

         // gather velocity data
         std::vector<real_t> instantaneousVelocityVector(parameters_->channelHalfWidth, 0_r);
         std::vector<real_t> meanVelocityVector(parameters_->channelHalfWidth, 0_r);
         std::vector<real_t> tkeSGSVector(parameters_->channelHalfWidth, 0_r);
         std::vector<real_t> tkeResolvedVector(parameters_->channelHalfWidth, 0_r);
         std::vector<real_t> reynoldsStressVector(parameters_->channelHalfWidth * TensorField_T::F_SIZE, 0_r);

         const auto idxFlow = int_c(parameters_->domainSize[parameters_->flowAxis] / uint_t(2));
         const auto idxRem = int_c(parameters_->domainSize[3 - parameters_->flowAxis - parameters_->wallAxis] / uint_t(2));

         Cell point;
         point[parameters_->flowAxis] = idxFlow;
         point[3 - parameters_->flowAxis - parameters_->wallAxis] = idxRem;

         const auto flowAxis = int_c(parameters_->flowAxis);

         auto blocks = blocks_.lock();
         WALBERLA_CHECK_NOT_NULLPTR(blocks)

         for(auto block = blocks->begin(); block != blocks->end(); ++block) {

            const auto * const velocity = block->template getData<VectorField_T>(velocityFieldId_);
            WALBERLA_CHECK_NOT_NULLPTR(velocity)

            const auto * const meanVelocity = block->template getData<VectorField_T>(meanVelocityFieldId_);
            WALBERLA_CHECK_NOT_NULLPTR(meanVelocity)

            const auto * const tkeSGS = block->template getData<ScalarField_T>(meanTkeSGSFieldId_);
            WALBERLA_CHECK_NOT_NULLPTR(tkeSGS)

            const auto * const sop = block->template getData<TensorField_T>(velocityWelford_->sum_of_productsID);
            WALBERLA_CHECK_NOT_NULLPTR(sop)

            for(uint_t idx = 0; idx < parameters_->channelHalfWidth; ++idx) {

               point[parameters_->wallAxis] = int_c(idx);

               Cell localCell;
               blocks->transformGlobalToBlockLocalCell(localCell, *block, point);

               if(velocity->xyzSize().contains(localCell)){
                  instantaneousVelocityVector[idx] = velocity->get(localCell, flowAxis);
                  meanVelocityVector[idx] = meanVelocity->get(localCell, flowAxis);
                  tkeSGSVector[idx] = tkeSGS->get(localCell);
                  for(uint_t i = 0; i < TensorField_T::F_SIZE; ++i) {
                     reynoldsStressVector[idx*TensorField_T::F_SIZE+i] = sop->get(localCell,i) / velocityWelford_->counter_;
                  }
                  tkeResolvedVector[idx] = real_c(0.5) * (
                     reynoldsStressVector[idx*TensorField_T::F_SIZE+0] +
                     reynoldsStressVector[idx*TensorField_T::F_SIZE+4] +
                     reynoldsStressVector[idx*TensorField_T::F_SIZE+8]
                  );
               }
            }
         }

         // MPI exchange information
         mpi::reduceInplace(instantaneousVelocityVector, mpi::SUM);
         mpi::reduceInplace(meanVelocityVector, mpi::SUM);
         mpi::reduceInplace(tkeSGSVector, mpi::SUM);
         mpi::reduceInplace(tkeResolvedVector, mpi::SUM);
         mpi::reduceInplace(reynoldsStressVector, mpi::SUM);

         WALBERLA_ROOT_SECTION()
         {
            std::ofstream velocityOS;
            filesystem::path path = velocityProfilesFilePath_;
            if (separateFiles_) {
               path.concat("_t" + std::to_string(timeloop_->getCurrentTimeStep()) + ".txt");
               velocityOS.open(path, std::ios::out | std::ios::trunc);
            } else {
               path.concat("_t" + std::to_string(parameters_->timesteps-1) + ".txt");
               velocityOS.open(path, std::ios::out | std::ios::trunc);
            }

            if (velocityOS.is_open()) {
               if (!separateFiles_) velocityOS << "# t = " << ts << "\n";
               velocityOS << "# y/delta\t y+\t u+\t u_mean\t u_instantaneous\t TKE_SGS\t TKE_resolved\t uu_rms\t uv_rms\t uw_rms\t vu_rms\t vv_rms\t vw_rms\t wu_rms\t wv_rms\t ww_rms\n";

               const auto & viscosity = parameters_->viscosity;
               const auto bulkReynolds = 2_r * channelHalfWidth * bulkVelocity / viscosity;
               const auto frictionReynolds = dean_correlation::calculateFrictionReynoldsNumber(bulkReynolds);
               const auto frictionVelocity = frictionReynolds * viscosity / channelHalfWidth;

               for(uint_t idx = 0; idx < parameters_->channelHalfWidth; ++idx) {
                  // relative position
                  velocityOS << (real_c(idx)+0.5_r) / channelHalfWidth << "\t";
                  // y+
                  velocityOS << (real_c(idx)+0.5_r) * frictionVelocity / viscosity << "\t";
                  // u+
                  velocityOS << meanVelocityVector[idx] / frictionVelocity << "\t";
                  // mean velocity
                  velocityOS << meanVelocityVector[idx] << "\t";
                  // instantaneous velocity
                  velocityOS << instantaneousVelocityVector[idx] << "\t";
                  // subgrid-scale TKE
                  velocityOS << tkeSGSVector[idx] << "\t";
                  // resolved TKE
                  velocityOS << tkeResolvedVector[idx] << "\t";
                  // Reynolds stresses
                  for(uint_t i = 0; i < TensorField_T::F_SIZE; ++i) {
                     velocityOS << reynoldsStressVector[idx*TensorField_T::F_SIZE+i] << "\t";
                  }

                  velocityOS << "\n";
               }
               velocityOS.close();
            } else{
               WALBERLA_ABORT("Could not open velocity plot file " << path.generic_string())
            }
         }
      }

    private:


      SimulationParameters const * const parameters_{};
      ForceCalculator<VectorField_T> const * const forceCalculator_{};

      Timeloop * const timeloop_{};
      const std::weak_ptr<StructuredBlockStorage> blocks_;
      Welford_T * const velocityWelford_{};

      const BlockDataID velocityFieldId_{};
      const BlockDataID meanVelocityFieldId_{};
      const BlockDataID meanTkeSGSFieldId_{};

      const uint_t plotFrequency_{};
      const uint_t plotStart_{};

      const bool separateFiles_{false};
      const std::string baseFolder_{"output"};
      filesystem::path velocityProfilesFilePath_;
      filesystem::path forcingDataFilePath_;
   };

   /////////////////////
   /// Main Function ///
   /////////////////////

   int main(int argc, char** argv) {

      Environment walberlaEnv(argc, argv);

      if (!walberlaEnv.config()) { WALBERLA_ABORT("No configuration file specified!") }

      ///////////////////////////////////////////////////////
      /// Block Storage Creation and Simulation Parameter ///
      ///////////////////////////////////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Creating block forest...")

      const auto channelParameter = walberlaEnv.config()->getOneBlock("TurbulentChannel");
      const SimulationParameters simulationParameters(channelParameter);

      // domain creation
      std::shared_ptr<StructuredBlockForest> blocks;
      {
         Vector3< uint_t > numBlocks;
         Vector3< uint_t > cellsPerBlock;
         blockforest::calculateCellDistribution(simulationParameters.domainSize,
                                                uint_c(mpi::MPIManager::instance()->numProcesses()),
                                                numBlocks, cellsPerBlock);

         const auto & periodicity = simulationParameters.periodicity;
         const auto & domainSize = simulationParameters.domainSize;

         SetupBlockForest sforest;

         sforest.addWorkloadMemorySUIDAssignmentFunction( blockforest::uniformWorkloadAndMemoryAssignment );

         sforest.init( AABB(0_r, 0_r, 0_r, real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2])),
                       numBlocks[0], numBlocks[1], numBlocks[2], periodicity[0], periodicity[1], periodicity[2] );

         // calculate process distribution

         const memory_t memoryLimit = numeric_cast< memory_t >( sforest.getNumberOfBlocks() );

         const blockforest::GlobalLoadBalancing::MetisConfiguration< SetupBlock > metisConfig(
            true, false, std::bind( blockforest::cellWeightedCommunicationCost, std::placeholders::_1, std::placeholders::_2,
                                    cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] ) );

         sforest.calculateProcessDistribution_Default( uint_c( MPIManager::instance()->numProcesses() ), memoryLimit,
                                                       "hilbert", 10, false, metisConfig );

         if( !MPIManager::instance()->rankValid() )
            MPIManager::instance()->useWorldComm();

         // create StructuredBlockForest (encapsulates a newly created BlockForest)

         WALBERLA_LOG_INFO_ON_ROOT("SetupBlockForest created successfully:\n" << sforest)

         sforest.writeVTKOutput("domain_decomposition");

         auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, false );

         blocks = std::make_shared< StructuredBlockForest >( bf, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] );
         blocks->createCellBoundingBoxes();

      }

      ////////////////////////////////////
      /// PDF Field and Velocity Setup ///
      ////////////////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Creating fields...")

      // Common Fields
      const BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0.0), codegen::layout);
      const BlockDataID meanVelocityFieldId = field::addToStorage< VectorField_T >(blocks, "mean velocity", real_c(0.0), codegen::layout);
      const BlockDataID sopFieldId = field::addToStorage< TensorField_T >(blocks, "sum of products", real_c(0.0), codegen::layout);

      const BlockDataID tkeSgsFieldId = field::addToStorage< ScalarField_T >(blocks, "tke_SGS", real_c(0.0), codegen::layout);
      const BlockDataID meanTkeSgsFieldId = field::addToStorage< ScalarField_T >(blocks, "mean_tke_SGS", real_c(0.0), codegen::layout);

      const BlockDataID omegaFieldId = field::addToStorage< ScalarField_T >(blocks, "omega_out", real_c(0.0), codegen::layout);

      const BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

      // CPU Field for PDFs
      const BlockDataID pdfFieldId = field::addToStorage< PdfField_T >(blocks, "pdf field", real_c(0.0), codegen::layout);

      ///////////////////////////////////////////
      /// Force and bulk velocity calculation ///
      ///////////////////////////////////////////

      ForceCalculator<VectorField_T> forceCalculator(blocks, velocityFieldId, simulationParameters);

      //////////////
      /// Setter ///
      //////////////

      WALBERLA_LOG_INFO_ON_ROOT("Setting up fields...")

      // Velocity field setup
      setVelocityFieldsAsmuth<VectorField_T>(
         blocks, velocityFieldId, meanVelocityFieldId,
         simulationParameters.targetFrictionVelocity, simulationParameters.channelHalfWidth,
         5.5_r, 0.41_r, simulationParameters.viscosity,
         simulationParameters.wallAxis, simulationParameters.flowAxis );

      forceCalculator.setBulkVelocity(simulationParameters.targetBulkVelocity);
      const auto initialForce = forceCalculator.calculateDrivingForce();

      // pdfs setup
      Setter_T pdfSetter(pdfFieldId, velocityFieldId, initialForce, simulationParameters.density);

      for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
         pdfSetter(blockIt.get());

      /////////////
      /// Sweep ///
      /////////////

      WALBERLA_LOG_INFO_ON_ROOT("Creating sweeps...")

      const auto omega = lbm::collision_model::omegaFromViscosity(simulationParameters.viscosity);
      StreamCollideSweep_T streamCollideSweep(omegaFieldId, pdfFieldId, velocityFieldId, initialForce, omega);

      WelfordSweep_T welfordSweep(meanVelocityFieldId, sopFieldId, velocityFieldId, 0_r);
      TKEWelfordSweep_T welfordTKESweep(meanTkeSgsFieldId, tkeSgsFieldId, 0_r);

      TkeSgsWriter_T tkeSgsWriter(omegaFieldId, pdfFieldId, tkeSgsFieldId, initialForce, omega);

      /////////////////////////
      /// Boundary Handling ///
      /////////////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Creating boundary handling...")

      const FlagUID fluidFlagUID("Fluid");

      Config::Block boundaryBlock;
      boundaries::createBoundaryConfig(simulationParameters, boundaryBlock);

      std::unique_ptr<WFB_bottom_T> wfb_bottom_ptr = std::make_unique<WFB_bottom_T>(blocks, meanVelocityFieldId, pdfFieldId, omega, simulationParameters.targetFrictionVelocity);
      std::unique_ptr<WFB_top_T > wfb_top_ptr = std::make_unique<WFB_top_T>(blocks, meanVelocityFieldId, pdfFieldId, omega, simulationParameters.targetFrictionVelocity);

      NoSlip_T noSlip(blocks, pdfFieldId);
      FreeSlip_top_T freeSlip_top(blocks, pdfFieldId);

      geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, Config::BlockHandle(&boundaryBlock));
      geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, fluidFlagUID);

      noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("NoSlip"), fluidFlagUID);
      freeSlip_top.fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("FreeSlip"), fluidFlagUID);
      wfb_bottom_ptr->fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("WFB_bottom"), fluidFlagUID);
      wfb_top_ptr->fillFromFlagField< FlagField_T >(blocks, flagFieldId, FlagUID("WFB_top"), fluidFlagUID);

      //////////////
      /// Output ///
      //////////////

      WALBERLA_LOG_INFO_ON_ROOT("Creating field output...")

      // vtk output
      auto vtkWriter = vtk::createVTKOutput_BlockData(
         blocks, "field_writer", simulationParameters.vtkFrequency, 0, false, "vtk_out", "simulation_step",
         false, false, true, false
      );
      vtkWriter->setInitialWriteCallsToSkip(simulationParameters.vtkStart);

      // velocity field writer
      auto velocityWriter = std::make_shared<field::VTKWriter<VectorField_T>>(velocityFieldId, "instantaneous velocity");
      vtkWriter->addCellDataWriter(velocityWriter);

      auto meanVelocityFieldWriter = std::make_shared<field::VTKWriter<VectorField_T>>(meanVelocityFieldId, "mean velocity");
      vtkWriter->addCellDataWriter(meanVelocityFieldWriter);

      // vtk writer
      {
         auto flagOutput = vtk::createVTKOutput_BlockData(
            blocks, "flag_writer", 1, 1, false, "vtk_out", "simulation_step",
            false, true, true, false
         );
         auto flagWriter = std::make_shared<field::VTKWriter<FlagField_T>>(flagFieldId, "flag field");
         flagOutput->addCellDataWriter(flagWriter);
         flagOutput->write();
      }


      /////////////////
      /// Time Loop ///
      /////////////////

      WALBERLA_LOG_INFO_ON_ROOT("Creating timeloop...")

      SweepTimeloop timeloop(blocks->getBlockStorage(), simulationParameters.timesteps);

      // Communication
      blockforest::communication::UniformBufferedScheme< Stencil_T > communication(blocks);
      communication.addPackInfo(make_shared< PackInfo_T >(pdfFieldId));

      auto setNewForce = [&](const real_t newForce) {
         streamCollideSweep.F_x_ = newForce;
         tkeSgsWriter.F_x_ = newForce;
         tkeSgsWriter.F_x_ = newForce;
      };

      // plotting
      const bool outputSeparateFiles = channelParameter.getParameter<bool>("separate_files", false);
      const TurbulentChannelPlotter<WelfordSweep_T > plotter(&simulationParameters, &timeloop, &forceCalculator, blocks,
                                                             velocityFieldId, meanVelocityFieldId,
                                                             meanTkeSgsFieldId, &welfordSweep,
                                                             outputSeparateFiles);

      //NOTE must convert sweeps that are altered to lambdas, otherwise copy and counter will stay 0
      auto welfordLambda = [&welfordSweep, &welfordTKESweep](IBlock * block) {
         welfordSweep(block);
         welfordTKESweep(block);
      };

      auto wfbLambda = [&wfb_bottom_ptr, &wfb_top_ptr](IBlock * block) {
         wfb_bottom_ptr->operator()(block);
         wfb_top_ptr->operator()(block);
      };

      auto streamCollideLambda = [&streamCollideSweep](IBlock * block) {
         streamCollideSweep(block);
      };

      // Timeloop
      timeloop.add() << BeforeFunction(communication, "communication")
                     << BeforeFunction([&](){forceCalculator.calculateBulkVelocity();}, "bulk velocity calculation")
                     << BeforeFunction([&](){
                           const auto newForce = forceCalculator.calculateDrivingForce();
                           setNewForce(newForce);
                        }, "new force setter")
                     << Sweep([](IBlock *){}, "new force setter");
      timeloop.add() << Sweep(freeSlip_top, "freeSlip");
      timeloop.add() << Sweep(noSlip, "noSlip");
      timeloop.add() << Sweep(wfbLambda, "wall function bounce");
      timeloop.add() << Sweep(streamCollideLambda, "stream and collide");
      timeloop.add() << BeforeFunction([&](){
                           const uint_t velCtr = uint_c(welfordSweep.counter_);
                           if((timeloop.getCurrentTimeStep() == simulationParameters.samplingStart) ||
                              (timeloop.getCurrentTimeStep() > simulationParameters.samplingStart && simulationParameters.samplingInterval && (velCtr % simulationParameters.samplingInterval == 0))) {
                              welfordSweep.counter_ = real_t(0);
                              welfordTKESweep.counter_ = real_t(0);
                               for(auto & block : *blocks) {
                                  auto * sopField = block.template getData<TensorField_T >(sopFieldId);
                                  sopField->setWithGhostLayer(0.0);

                                  auto * tkeField = block.template getData<ScalarField_T>(tkeSgsFieldId);
                                  tkeField->setWithGhostLayer(0.0);
                               }
                           }

                           welfordSweep.counter_ = welfordSweep.counter_ + real_c(1);
                           welfordTKESweep.counter_ = welfordTKESweep.counter_ + real_c(1);
                        }, "welford sweep")
                     << Sweep(welfordLambda, "welford sweep");
      timeloop.add() << Sweep(tkeSgsWriter, "TKE_SGS writer");

      timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkWriter), "VTK field output");
      timeloop.addFuncAfterTimeStep(plotter, "Turbulent quantity plotting");

      // LBM stability check
      timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                       walberlaEnv.config(), blocks, pdfFieldId, flagFieldId, fluidFlagUID ) ),
                                    "LBM stability check" );

      // Time logger
      timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), 5_r),
                                    "remaining time logger");


      WALBERLA_LOG_INFO_ON_ROOT("Running timeloop with " << timeloop.getNrOfTimeSteps() - 1 << " timesteps...")

      WcTimingPool timing;

      WcTimer timer;
      timer.start();

      timeloop.run(timing);

      timer.end();

      double time = timer.max();
      walberla::mpi::reduceInplace(time, walberla::mpi::MAX);

      const auto timeloopTiming = timing.getReduced();
      WALBERLA_LOG_INFO_ON_ROOT("Timeloop timing:\n" << *timeloopTiming)

      const walberla::lbm::PerformanceEvaluation<FlagField_T> performance(blocks, flagFieldId, fluidFlagUID);
      performance.logResultOnRoot(simulationParameters.timesteps, time);

      return EXIT_SUCCESS;
   }

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
