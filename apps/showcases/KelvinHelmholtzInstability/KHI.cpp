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
//! \file KHI.cpp
//! \author Brendan Waters <brendan.waters@sydney.edu.au>
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

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   #include "gpu/AddGPUFieldToStorage.h"
   #include "gpu/DeviceSelectMPI.h"
   #include "gpu/FieldCopy.h"
   #include "gpu/GPUWrapper.h"
   #include "gpu/HostFieldAllocator.h"
   #include "gpu/ParallelStreams.h"
   #include "gpu/communication/UniformGPUScheme.h"
   #include "gpu/timing/DeviceTimingPool.h"
   #include "gpu/timeloop/DeviceSweepTimeloop.h"
#endif

#include <lbm/all.h>

#include "lbm_generated/field/PdfField.h"
#include "lbm_generated/field/AddToStorage.h"
#include "lbm_generated/evaluation/PerformanceEvaluation.h"
#include "lbm_generated/communication/UniformGeneratedPdfPackInfo.h"

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   #include "lbm_generated/gpu/UniformGeneratedGPUPdfPackInfo.h"
   #include "lbm_generated/gpu/GPUPdfField.h"
   #include "lbm_generated/gpu/AddToStorage.h"
#endif

#include "python_coupling/CreateConfig.h"
#include "python_coupling/DictWrapper.h"

#include <timeloop/all.h>

//    Codegen Includes
#include "InfoHeader.h"

namespace walberla {
   // See InfoHeader.h for typedef aliases 
   constexpr uint_t FieldGhostLayer{1};

   /*------------------------------------------------------------------------------------------------------*/
   struct SimulationParameters {

      SimulationParameters(const Config::BlockHandle & config)
      {
         reNumber   = config.getParameter<real_t>( "reNumber");
         machNumber = config.getParameter<real_t>( "machNumber");

         const real_t viscosity_SI = config.getParameter<real_t>("viscosity_SI");
         Href_SI      = config.getParameter<real_t>("Href_SI");

         const auto velocity_SI  = (reNumber * viscosity_SI) / Href_SI ;

         cellsPerLength = config.getParameter<uint_t>("cellsPerLength");

         const auto speedOfSound = 1.0_r / std::sqrt( 3.0_r );
         velocity_LB  = machNumber * speedOfSound;

         dx_SI = Href_SI / real_c(cellsPerLength);
         const auto dt_SI = (velocity_LB / velocity_SI) * dx_SI;

         const auto C_l   = dx_SI;
         const auto dx_LB = dx_SI/C_l; 
         
         const auto Href_LB = Href_SI/C_l;

         const auto C_t   = dt_SI;
         const auto dt_LB = dt_SI/C_t;

         const auto viscosity_LB = viscosity_SI *  C_t/(C_l * C_l);    // Physical units are m^2/s

         omega = lbm::collision_model::omegaFromViscosity(viscosity_LB);

         const real_t theta_SI = config.getParameter<real_t>( "momentumThickness_SI"); 
         theta_LB = theta_SI/C_l; 

         flowAxis = codegen::flowAxis;
         wallAxis = codegen::wallAxis;

         remAxis = 3 - flowAxis - wallAxis;

         periodicity[wallAxis] = false;

         domainSize = Vector3<uint_t>(cellsPerLength);
         domainSize[remAxis] = 1;

         timestepsPerPeriod = Href_LB/velocity_LB ;
         const real_t simulationPeriods = config.getParameter<real_t>("simulationPeriods");
         timesteps = uint_c( std::ceil( simulationPeriods * timestepsPerPeriod ) ) + 1;

         const real_t vtkPeriods = config.getParameter<real_t>("vtk_periods",  real_c(0.0));
         vtkFrequency = uint_c( std::ceil( vtkPeriods * timestepsPerPeriod ) );
         vtkStart     = config.getParameter<uint_t>("vtk_start", 0);

         const real_t plotPeriods = config.getParameter<real_t>( "plot_periods",  real_c(0.05) );
         plotStart = config.getParameter<uint_t>("plot_start", 0);
         plotFrequency = uint_c( std::ceil( plotPeriods * timestepsPerPeriod ) );

         outputBaseFolder = "PostProcessing/SimData/Re" + std::to_string(int(reNumber)) + "/" + \
                            codegen::configID + "/nCells_" + std::to_string(int(cellsPerLength));

         WALBERLA_ROOT_SECTION() {
            if( !filesystem::exists(outputBaseFolder) ) {
               filesystem::create_directories(outputBaseFolder);
            }
         }

         std::ostringstream oss;
         oss <<"- Time Parameters:"
             << "\n   + No. Course Timesteps: " << timesteps;
            
         oss <<"\n- Simulation Parameters:"  
             << "\n   + dx (-):              " << dx_LB
             << "\n   + dt (-):              " << dt_LB
             << "\n   + velocity (-):        " << velocity_LB
             << "\n   + viscosity (-):       " << viscosity_LB
             << "\n   + L x W x H (-):       " << Href_LB
             << "\n   + Interface Width(-):  " << theta_LB
             << "\n   + Reynolds Number:     " << reNumber  << " (LB: "<< velocity_LB *Href_LB/viscosity_LB << " )"
             << "\n   + Mach Number:         " << machNumber
             << "\n   + Relaxation Rate:     " << omega;
            
         oss <<"\n- Physical Parameters:"  
             << "\n   + dx (m):              " << dx_SI
             << "\n   + dt (s):              " << dt_SI
             << "\n   + velocity (m/s):      " << velocity_SI
             << "\n   + viscosity (m^2/s):   " << viscosity_SI
             << "\n   + L x W x H (m):       " << Href_SI
             << "\n   + Interface Width:     " << theta_SI;
             
         std::string str = oss.str();

         WALBERLA_LOG_INFO_ON_ROOT( "Parameters:\n" << str);
      }

      real_t reNumber{};
      real_t machNumber{};
      uint_t cellsPerLength{};
      
      real_t dx_SI {};
      real_t Href_SI {};

      Vector3<uint_t> domainSize{};

      uint_t wallAxis{};
      uint_t flowAxis{};
      uint_t remAxis{};

      Vector3<uint_t> periodicity{true};
      const real_t density{1.0_r};

      real_t velocity_LB {};
      real_t omega{};
      real_t theta_LB {};

      real_t timestepsPerPeriod{};
      uint_t timesteps{};

      uint_t vtkFrequency{};
      uint_t vtkStart{};
      uint_t plotFrequency{};
      uint_t plotStart{};

      std::string outputBaseFolder{};
   };

   /*------------------------------------------------------------------------------------------------------*/
   /*
    * Initialises the velocity field as a quasi-2D shear-layer with stream function perterbation for single mode KHI
    * Uses 6th order Gaussian Quadrature
    */
   template<typename VelocityField_T>
   class VelocityFieldSetter {
   public:
      VelocityFieldSetter( const std::weak_ptr<StructuredBlockStorage>& forest,
                           const BlockDataID & velocityFieldId,  const SimulationParameters & parameters)
         : forest_(forest), velocityFieldId_(velocityFieldId), parameters_(parameters),
           wallAxis_(parameters_.wallAxis), flowAxis_(parameters_.flowAxis), remAxis_( 3 - wallAxis_ - flowAxis_),
           rV1_( 0.025_r * parameters_.velocity_LB ),rV2_( 0.05_r * parameters_.velocity_LB ) {}

      void operator()() {
         auto blocks = forest_.lock();
         WALBERLA_ASSERT_NOT_NULLPTR(blocks)

         const Vector3<real_t> dxdydz {blocks->dx(), blocks->dy(), blocks->dz()};

         const auto domainSize = blocks->getDomain().max();

         rLength_ = domainSize[flowAxis_];
         interfacePosition_ = 0.5_r * domainSize[wallAxis_];

         rK1_ = ( 2.0_r * math::pi ) / domainSize[flowAxis_];
         rK2_ = ( 4.0_r * math::pi ) / domainSize[flowAxis_];

         for( auto block = blocks->begin(); block != blocks->end(); ++block ) {

            auto * velocityField = block->template getData<VelocityField_T>(velocityFieldId_);
            WALBERLA_ASSERT_NOT_NULLPTR(velocityField)

            const auto ci = velocityField->xyzSize();
            for(auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt) {

               Cell globalCell(*cellIt);
               blocks->transformBlockLocalToGlobalCell(globalCell, *block);
               Vector3<real_t> cellCenter;
               blocks->getCellCenter(cellCenter, globalCell);

               Vector3<real_t> vel; 
               vel[flowAxis_] = quadrature( cellCenter[flowAxis_], cellCenter[wallAxis_], dxdydz[flowAxis_], dxdydz[wallAxis_], flowAxis_ );
               vel[wallAxis_] = quadrature( cellCenter[flowAxis_], cellCenter[wallAxis_], dxdydz[flowAxis_], dxdydz[wallAxis_], wallAxis_ );
               vel[remAxis_]  = 0.0_r;
               
               for (uint_t idx = 0; idx < VelocityField_T::F_SIZE; ++idx)  {
                  velocityField->get(*cellIt, idx) = vel[idx];
               }
            }
         }
      }

   private:
      const std::weak_ptr<StructuredBlockStorage> forest_;
      const BlockDataID velocityFieldId_;
      const SimulationParameters parameters_;
      const uint_t wallAxis_, flowAxis_, remAxis_;

      const real_t rV1_, rV2_;

      real_t interfacePosition_;
      real_t rK1_, rK2_;
      real_t rLength_;
      
      real_t velocityProfile(const real_t & x, const real_t & y, const uint_t & axis)
      {  
         const auto yShifted = y - interfacePosition_;

         const auto rA1 = (1.0_r - std::exp(- 2.0_r * rK1_ * (0.5_r * rLength_ - std::abs(yShifted))))/(1.0_r - std::exp(-rK1_ * rLength_));
         const auto rA2 = (1.0_r - std::exp(- 2.0_r * rK2_ * (0.5_r * rLength_ - std::abs(yShifted))))/(1.0_r - std::exp(-rK2_ * rLength_));
         
         if (axis == flowAxis_){
            return -0.5_r * parameters_.velocity_LB * std::tanh( yShifted /( 2.0_r * (2.0_r * interfacePosition_) * parameters_.theta_LB )) + \
                           std::copysign(1.0_r, yShifted) * real_c(  rV1_ * rA1 * std::exp( -rK1_ * std::abs( yShifted ) ) * std::cos( rK1_ * x ) + \
                                                                     rV2_ * rA2 * std::exp( -rK2_ * std::abs( yShifted ) ) * std::cos( rK2_ * x ) + \
                                                                     2.0_r * rV1_ * std::exp( -rK1_ * std::abs( yShifted ) ) * \
                                                                                    std::exp( -2.0_r *  rK1_* ( rLength_/2.0_r - std::abs( yShifted ) )) * \
                                                                                    std::cos( rK1_ * x ) / (1.0_r - std::exp( -rLength_ * rK1_ )) + \
                                                                     2.0_r * rV2_ * std::exp( -rK2_ * std::abs( yShifted ) ) * \
                                                                                    std::exp( -2.0_r * rK2_ * ( rLength_/2.0_r - std::abs( yShifted ) )) * \
                                                                                    std::cos( rK2_ * x ) / (1.0_r - std::exp( -rLength_ * rK2_ )));
         } else if (axis == wallAxis_){
            return - real_c(  rA1 * rV1_ * std::sin( rK1_ * x ) * std::exp( - rK1_ * std::abs(yShifted) )+
                              rA2 * rV2_ * std::sin( rK2_ * x ) * std::exp( - rK2_ * std::abs(yShifted) ));
            
         } else { WALBERLA_ABORT("Incompatible Initial Condition!!") }

      }

      real_t quadrature( const real_t & x_cc, const real_t & y_cc, 
                         const real_t & dx, const real_t & dy, const uint_t & axis )
      {
         // Sixth order Gaussian rule
         const auto s  =  std::sqrt(3.0_r) / std::sqrt(5.0_r);
         const auto w1 =  5.0_r / 9.0_r;
         const auto w2 =  8.0_r / 9.0_r;

         const auto x  = std::array{ x_cc - dx * s / 2.0_r, x_cc + dx * s / 2.0_r, x_cc };
         const auto y  = std::array{ y_cc - dy * s / 2.0_r,   y_cc + dy * s / 2.0_r, y_cc };
         const auto wx = std::array{ w1, w1, w2 };
         const auto wy = std::array{ w1, w1, w2 };
         
         auto result { 0.0_r };
         for (uint_t i = 0; i < 3; ++i) {
            for (uint_t j = 0; j < 3; ++j){
               result += wx[i] * wy[j] * velocityProfile(x[i], y[j], axis);
            }
         }

         return (0.25_r * result);
      }
   };

   /*
    * Initialises the scalar field as a contact surface Φ ∈ [1,2]
    */
   template<typename ScalarField_T>
   void setScalarField( const std::weak_ptr<StructuredBlockStorage>& forest,
                         const BlockDataID & scalarFieldId,  const SimulationParameters & parameters ) {

      auto blocks = forest.lock();
      WALBERLA_ASSERT_NOT_NULLPTR(blocks)

      const auto domainSize = blocks->getDomain().max();

      const auto wallAxis = parameters.wallAxis;
      
      for( auto block = blocks->begin(); block != blocks->end(); ++block ) {

         auto * scalarField = block->template getData<ScalarField_T>(scalarFieldId);
         WALBERLA_ASSERT_NOT_NULLPTR(scalarField)

         const auto ci = scalarField->xyzSizeWithGhostLayer();
         for(auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt) {

            Cell globalCell(*cellIt);
            blocks->transformBlockLocalToGlobalCell(globalCell, *block);
            Vector3<real_t> cellCenter;
            blocks->getCellCenter(cellCenter, globalCell);

            if (cellCenter[wallAxis] <= 0.5_r * domainSize[wallAxis]){
               scalarField->get(*cellIt, 0) = 2.0_r;
            } else {
               scalarField->get(*cellIt, 0) = 1.0_r;
            }
         }
      }

   }

   /*------------------------------------------------------------------------------------------------------*/
   namespace boundaries {
      void createBoundaryConfig(const SimulationParameters & parameters,  Config::Block & boundaryBlock, const bool & scalar_config = false) {
         if ( !scalar_config ){
            auto & bottomWall = boundaryBlock.createBlock("Border");
            bottomWall.addParameter("direction", stencil::dirToString[stencil::directionFromAxis(parameters.wallAxis, true)]);
            bottomWall.addParameter("walldistance", "-1");
            bottomWall.addParameter("flag", "FreeSlip");

            auto & topWall = boundaryBlock.createBlock("Border");
            topWall.addParameter("direction", stencil::dirToString[stencil::directionFromAxis(parameters.wallAxis, false)]);
            topWall.addParameter("walldistance", "-1");
            topWall.addParameter("flag", "FreeSlip");
         } else {
            auto & bottomWall = boundaryBlock.createBlock("Border");
            bottomWall.addParameter("direction", stencil::dirToString[stencil::directionFromAxis(parameters.wallAxis, true)]);
            bottomWall.addParameter("walldistance", "-1");
            bottomWall.addParameter("flag", "Neumann");

            auto & topWall = boundaryBlock.createBlock("Border");
            topWall.addParameter("direction", stencil::dirToString[stencil::directionFromAxis(parameters.wallAxis, false)]);
            topWall.addParameter("walldistance", "-1");
            topWall.addParameter("flag", "Neumann");
         }
      }
   }

   /*------------------------------------------------------------------------------------------------------*/
   template< typename Timeloop_T >
   class MomentumThicknessPlotter {
    public:
      MomentumThicknessPlotter(  const std::weak_ptr<StructuredBlockStorage> & blocks, SimulationParameters const * const parameters, 
                                 Timeloop_T * const timeloop, const BlockDataID velocityFieldId)
                              : blocks_(blocks), parameters_(parameters), timeloop_(timeloop), velocityFieldId_(velocityFieldId),
                              plotFrequency_(parameters->plotFrequency), plotStart_(parameters->plotStart), baseFolder_(parameters->outputBaseFolder)
      {
         if(!plotFrequency_)
            return;

         const filesystem::path path(baseFolder_);
         momentumThicknessFilePath_ = path / ("momentum_thickness.txt");
         tkeFilePath_ = path / ("tke.txt");

         WALBERLA_ROOT_SECTION() {
            filesystem::remove(momentumThicknessFilePath_);
            filesystem::remove(tkeFilePath_);
         }

         WALBERLA_MPI_WORLD_BARRIER();

         std::ofstream mtOs(momentumThicknessFilePath_, std::ios::out | std::ios::trunc);
         if(mtOs.is_open()) {
            mtOs << "# KHI timescale (t*U)/L \t momentum thickness (theta/L)\n";
            mtOs.close();
         } else {
            WALBERLA_ABORT("Could not open momentum thickness data file.")
         }

         std::ofstream tkeOs(tkeFilePath_, std::ios::out | std::ios::trunc);
         if(tkeOs.is_open()) {
            tkeOs << "# y/L \t tke \n";
            tkeOs.close();
         } else {
            WALBERLA_ABORT("Could not open tke data file.")
         }

      }

      using BeforeFunction = std::function<void ()>;
      void addBeforeFunction( BeforeFunction f ) { beforeFunctions_.push_back( f ); }

      real_t getMaxMinVelocity( const real_t cell_center, const real_t domain_max ){

         const auto yShifted = cell_center - 0.5_r * domain_max; // Interface at 0.5 * ymax
         const auto rel_y = yShifted / domain_max;

         return - 0.5_r * parameters_->velocity_LB * std::tanh( rel_y /( 2.0_r * parameters_->theta_LB ));
      }

      void operator()() {

         const auto ts = timeloop_->getCurrentTimeStep();
         if(ts < plotStart_)
            return;

         if (!plotFrequency_ || (ts - plotStart_) % plotFrequency_ != 0)
            return;

         for( auto func = beforeFunctions_.begin(); func != beforeFunctions_.end(); ++func )
            ( *func )( );

         auto blocks = blocks_.lock();
         WALBERLA_ASSERT_NOT_NULLPTR(blocks)

         const auto domainSize = blocks->getDomain().max();

         const auto yMax = uint_c(domainSize[parameters_->wallAxis]);
         
         const auto rU1 =  getMaxMinVelocity( 0.5_r, real_c(yMax) );
         const auto rU2 =  getMaxMinVelocity( real_c(yMax) - 0.5_r, real_c(yMax) );

         std::vector<real_t> planeAverageVelocityVectorX(yMax, 0.0_r);
         std::vector<real_t> planeAverageVelocityVectorY(yMax, 0.0_r);

         std::vector<real_t> planeAverageTKEVector(yMax, 0.0_r);

         const auto idxRem  = int_c(domainSize[parameters_->remAxis] / uint_t(2));

         Cell point;
         point[parameters_->remAxis] = idxRem;

         // Loop for mean velocity
         for(auto block = blocks->begin(); block != blocks->end(); ++block) {

            const auto * const velocity = block->template getData<VectorField_T>(velocityFieldId_);
            WALBERLA_ASSERT_NOT_NULLPTR(velocity)

            for(uint_t idx_y = 0; idx_y < uint_c(domainSize[parameters_->wallAxis]) ; ++idx_y) {
               point[parameters_->wallAxis] = int_c(idx_y);

               auto planeAveVelx {0.0_r};
               auto planeAveVely {0.0_r};

               for(uint_t idx_x = 0; idx_x < uint_c(domainSize[parameters_->flowAxis]) ; ++idx_x) {
                  point[parameters_->flowAxis] = int_c(idx_x);

                  Cell localCell;
                  blocks->transformGlobalToBlockLocalCell(localCell, *block, point);

                  if(velocity->xyzSize().contains(localCell)){
                     planeAveVelx += velocity->get(localCell, parameters_->flowAxis);
                     planeAveVely += velocity->get(localCell, parameters_->wallAxis);
                  }
               }

               planeAverageVelocityVectorX[idx_y] = planeAveVelx / real_c( domainSize[parameters_->flowAxis] );
               planeAverageVelocityVectorY[idx_y] = planeAveVely / real_c( domainSize[parameters_->wallAxis] );

            }
         }

         mpi::reduceInplace(planeAverageVelocityVectorX, mpi::SUM);
         mpi::reduceInplace(planeAverageVelocityVectorY, mpi::SUM);

         // Loop for tke
         for(auto block = blocks->begin(); block != blocks->end(); ++block) {

            const auto * const velocity = block->template getData<VectorField_T>(velocityFieldId_);
            WALBERLA_ASSERT_NOT_NULLPTR(velocity)

            for(uint_t idx_y = 0; idx_y < uint_c(domainSize[parameters_->wallAxis]) ; ++idx_y) {
               point[parameters_->wallAxis] = int_c(idx_y);

               const auto rUmean = planeAverageVelocityVectorX[idx_y];
               const auto rVmean = planeAverageVelocityVectorY[idx_y];

               auto planeAveTKE {0.0_r};

               for(uint_t idx_x = 0; idx_x < uint_c(domainSize[parameters_->flowAxis]) ; ++idx_x) {
                  point[parameters_->flowAxis] = int_c(idx_x);

                  Cell localCell;
                  blocks->transformGlobalToBlockLocalCell(localCell, *block, point);

                  if(velocity->xyzSize().contains(localCell)){
                     const auto u_prime = velocity->get(localCell, parameters_->flowAxis) - rUmean;
                     const auto v_prime = velocity->get(localCell, parameters_->wallAxis) - rVmean;
                     planeAveTKE += 0.5_r * (u_prime*u_prime + v_prime*v_prime);
                  }
               }

               planeAverageTKEVector[idx_y] = planeAveTKE / real_c( domainSize[parameters_->flowAxis] );

            }
         }

         mpi::reduceInplace(planeAverageTKEVector, mpi::SUM);

         WALBERLA_ROOT_SECTION() {
            std::ofstream mtOS(momentumThicknessFilePath_, std::ios::out | std::ios::app);
            if (mtOS.is_open())
            {
               auto theta { 0.0_r };

               for(uint_t idx = 0; idx < yMax ; ++idx) {
                  const auto rUmean = planeAverageVelocityVectorX[idx];

                  // Approximate integration for cell center values
                  theta += ((rU1 - rUmean)*(rUmean - rU2))/std::pow((rU1 - rU2), 2.0_r);
               }

               mtOS << real_c(ts) / (parameters_ -> timestepsPerPeriod) << "\t" << theta / real_c( domainSize[parameters_->wallAxis] ) << "\n"; 
               mtOS.close();
            }

            std::ofstream tkeOS(tkeFilePath_, std::ios::out | std::ios::app);
            if (tkeOS.is_open())
            {
               for(uint_t idx = 0; idx < yMax ; ++idx) {
                  const auto rTKEmean = planeAverageTKEVector[idx];

                  tkeOS << ( real_c(idx) + 0.5_r )  / ( real_c(yMax) - 0.5_r ) << "\t" << rTKEmean / (parameters_->velocity_LB * parameters_->velocity_LB)<< "\n";  
               }
               tkeOS.close();
            }
         }

      }

    private:
      const std::weak_ptr<StructuredBlockStorage> blocks_;
      
      SimulationParameters const * const parameters_{};

      Timeloop_T * const timeloop_{};

      const BlockDataID velocityFieldId_{};

      const uint_t plotFrequency_{};
      const uint_t plotStart_{};

      const std::string baseFolder_{};
      filesystem::path momentumThicknessFilePath_;
      filesystem::path tkeFilePath_;

      std::vector< BeforeFunction > beforeFunctions_;
   };

   /*------------------------------------------------------------------------------------------------------*
    * Main Function
    *------------------------------------------------------------------------------------------------------*/
   int main(int argc, char** argv) {

      Environment walberlaEnv(argc, argv);
      
      const std::string input_filename(argv[1]);
      const bool inputIsPython = string_ends_with(input_filename, ".py");

      if (!inputIsPython) { WALBERLA_ABORT("Configuration file must be a python script!") }

      mpi::MPIManager::instance()->useWorldComm();

      #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         gpu::selectDeviceBasedOnMpiRank();
         WALBERLA_GPU_CHECK(gpuPeekAtLastError())
      #endif

      for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
      {
         WALBERLA_MPI_WORLD_BARRIER()

         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            WALBERLA_GPU_CHECK(gpuPeekAtLastError())
         #endif

         auto config = *cfg;
         logging::configureLogging(config);

         /*------------------------------------------------------------------------------------------------------*
         *-----------------------------------CONFIGURE FROM INPUT PARAMETERS------------------------------------*
         *------------------------------------------------------------------------------------------------------*/
         const auto simulation_parameters = config->getOneBlock("Parameters");
         SimulationParameters simulationParameters(simulation_parameters);

         /*------------------------------------------------------------------------------------------------------
         *-------------------------------------------Block Storage Creation--------------------------------------
         *------------------------------------------------------------------------------------------------------*/
         WALBERLA_LOG_INFO_ON_ROOT("Creating block forest...")

         std::shared_ptr<StructuredBlockForest> blocks;
         {
            Vector3< uint_t > numBlocks;
            Vector3< uint_t > cellsPerBlock;
            blockforest::calculateCellDistribution(simulationParameters.domainSize,
                                                   uint_c(mpi::MPIManager::instance()->numProcesses()),
                                                   numBlocks, cellsPerBlock);

            const auto & periodicity = simulationParameters.periodicity;
            auto & domainSize = simulationParameters.domainSize;
            const Vector3<uint_t> newDomainSize(numBlocks[0] * cellsPerBlock[0], numBlocks[1] * cellsPerBlock[1], numBlocks[2] * cellsPerBlock[2]);

            if(domainSize != newDomainSize) {
               WALBERLA_LOG_WARNING_ON_ROOT("\nWARNING: Domain size has changed due to the chosen domain decomposition. \n Some Physical -> LBM scaling laws might not be valid\n")
            }

            SetupBlockForest sforest;

            sforest.addWorkloadMemorySUIDAssignmentFunction( blockforest::uniformWorkloadAndMemoryAssignment );

            sforest.init( AABB(0_r, 0_r, 0_r, real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2])),
                        numBlocks[0], numBlocks[1], numBlocks[2], periodicity[0], periodicity[1], periodicity[2] );

            const memory_t memoryLimit = numeric_cast< memory_t >( sforest.getNumberOfBlocks() );

            const blockforest::GlobalLoadBalancing::MetisConfiguration< SetupBlock > metisConfig(
               true, false, std::bind( blockforest::cellWeightedCommunicationCost, std::placeholders::_1, std::placeholders::_2,
                                       cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] ) );

            sforest.calculateProcessDistribution_Default( uint_c( MPIManager::instance()->numProcesses() ), memoryLimit,
                                                         "hilbert", 10, false, metisConfig );

            if( !MPIManager::instance()->rankValid() )
               MPIManager::instance()->useWorldComm();

            WALBERLA_LOG_INFO_ON_ROOT("SetupBlockForest created successfully:\n" << sforest)

            sforest.writeVTKOutput("domain_decomposition");

            auto bf = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, true );

            blocks = std::make_shared< StructuredBlockForest >( bf, cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2] );
            blocks->createCellBoundingBoxes();
         }

         /*------------------------------------------------------------------------------------------------------
         *------------------------------------------- DATA FIELDS ----------------------------------------------
         *------------------------------------------------------------------------------------------------------*/
         WALBERLA_LOG_INFO_ON_ROOT("Creating CPU fields...")

         const StorageSpecification_T StorageSpec = StorageSpecification_T();
         const BlockDataID pdfFieldId  = lbm_generated::addPdfFieldToStorage(blocks, "pdf field", StorageSpec, FieldGhostLayer, codegen::layout);

         const BlockDataID densityFieldId  = field::addToStorage< ScalarField_T >(blocks, "density", real_c(simulationParameters.density), codegen::layout);
         const BlockDataID velocityFieldId = field::addToStorage< VectorField_T >(blocks, "velocity", real_c(0.0), codegen::layout, FieldGhostLayer);
         
         const BlockDataID flagFieldId     = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

         const ScalarStorageSpecification_T ScalarStorageSpec = ScalarStorageSpecification_T();
         const BlockDataID scalarPdfFieldId  = lbm_generated::addPdfFieldToStorage(blocks, "scalar pdf field", ScalarStorageSpec, FieldGhostLayer, codegen::layout);
         const BlockDataID scalarFieldId     = field::addToStorage< ScalarField_T >(blocks, "scalar", real_c(0.0), codegen::layout);
         const BlockDataID scalarFlagFieldId = field::addFlagFieldToStorage< FlagField_T >(blocks, "scalar flag field");
         
         //------------------------------------------------------------------------------------------------------//
         // Initialise fields before transferring to GPU (if required)
         //------------------------------------------------------------------------------------------------------//
         VelocityFieldSetter<VectorField_T> setter(blocks, velocityFieldId, simulationParameters);
         setter();

         setScalarField<ScalarField_T>( blocks,  scalarFieldId, simulationParameters );

         //------------------------------------------------------------------------------------------------------//
         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            WALBERLA_LOG_INFO_ON_ROOT("Creating GPU fields...")
            const bool usedPitchedMem = simulation_parameters.getParameter< bool >("usedPitchedMem", false);

            const BlockDataID pdfFieldGpuId           = lbm_generated::addGPUPdfFieldToStorage< PdfField_T >(blocks, pdfFieldId, StorageSpec, "pdfs GPU", usedPitchedMem);
            const BlockDataID densityFieldGpuId       = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, densityFieldId, "density GPU", usedPitchedMem);
            const BlockDataID velocityFieldGpuId      = gpu::addGPUFieldToStorage< VectorField_T >(blocks, velocityFieldId,  "velocity GPU" , usedPitchedMem);

            const BlockDataID scalarPdfFieldGpuId     = lbm_generated::addGPUPdfFieldToStorage< ScalarPdfField_T >(blocks, scalarPdfFieldId, ScalarStorageSpec, "scalar pdfs GPU", usedPitchedMem);
            const BlockDataID scalarFieldGpuId        = gpu::addGPUFieldToStorage< ScalarField_T >(blocks, scalarFieldId, "scalar GPU", usedPitchedMem);
         #endif
         
         //------------------------------------------------------------------------------------------------------//

         WALBERLA_LOG_INFO_ON_ROOT("Creating sweeps...")
         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            const Vector3< int32_t >  gpuBlockSize = simulation_parameters.getParameter< Vector3< int32_t > >( "gpuBlockSize", Vector3< int32_t >(64, 1, 1) );

            SweepCollection_T streamCollideSweep( blocks, densityFieldGpuId, pdfFieldGpuId, velocityFieldGpuId,
                                                gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], simulationParameters.omega);

            ScalarSweepCollection_T scalarStreamCollideSweep(  blocks, scalarPdfFieldGpuId, scalarFieldGpuId, velocityFieldGpuId,
                                                               gpuBlockSize[0], gpuBlockSize[1], gpuBlockSize[2], simulationParameters.omega );

            int streamHighPriority = 0;
            int streamLowPriority  = 0;

            WALBERLA_GPU_CHECK(gpuDeviceGetStreamPriorityRange(&streamLowPriority, &streamHighPriority))
            
            streamCollideSweep.setOuterPriority(streamHighPriority);
            scalarStreamCollideSweep.setOuterPriority(streamHighPriority);
         #else 
            SweepCollection_T streamCollideSweep( blocks, densityFieldId, pdfFieldId, velocityFieldId, simulationParameters.omega);
            ScalarSweepCollection_T scalarStreamCollideSweep(  blocks, scalarPdfFieldId, scalarFieldId, velocityFieldId, simulationParameters.omega );
         #endif
         
         for (auto& block : *blocks)
         {
            streamCollideSweep.initialise(&block);
            scalarStreamCollideSweep.initialise(&block);
         }

         /*------------------------------------------------------------------------------------------------------
         *------------------------------------------- BOUNDARY HANDLING -----------------------------------------
         *------------------------------------------------------------------------------------------------------*/
         WALBERLA_LOG_INFO_ON_ROOT("Creating boundary handling...")

         Config::Block boundaryBlock;
         boundaries::createBoundaryConfig(simulationParameters, boundaryBlock);

         geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldId, Config::BlockHandle(&boundaryBlock));
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldId, FluidFlagUID);

         Config::Block scalarBoundaryBlock;
         boundaries::createBoundaryConfig(simulationParameters, scalarBoundaryBlock, true);

         geometry::initBoundaryHandling< FlagField_T >(*blocks, scalarFlagFieldId, Config::BlockHandle(&scalarBoundaryBlock));
         geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, scalarFlagFieldId, FluidFlagUID);

         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            BoundaryCollection_T boundaryCollection( blocks, flagFieldId, pdfFieldGpuId, FluidFlagUID );
            ScalarBoundaryCollection_T scalarBoundaryCollection( blocks, scalarFlagFieldId, scalarPdfFieldGpuId, FluidFlagUID );
         #else 
            BoundaryCollection_T boundaryCollection( blocks, flagFieldId, pdfFieldId, FluidFlagUID ); 
            ScalarBoundaryCollection_T scalarBoundaryCollection( blocks, scalarFlagFieldId, scalarPdfFieldId, FluidFlagUID );
         #endif

         /*------------------------------------------------------------------------------------------------------
         *--------------------------------------------- COMMUNICATION -------------------------------------------
         *------------------------------------------------------------------------------------------------------*/
         WALBERLA_LOG_INFO_ON_ROOT("Creating communicators...")
         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            const bool gpuEnabledMPI = simulation_parameters.getParameter< bool >("gpuEnabledMPI", false);

            UniformGPUScheme< LBMCommunicationStencil_T > communication(blocks, gpuEnabledMPI, false);
            communication.addPackInfo( std::make_shared< PackInfo_T >( pdfFieldGpuId ) );

            UniformGPUScheme< ScalarLBMCommunicationStencil_T > scalarCommunication(blocks, gpuEnabledMPI, false);
            scalarCommunication.addPackInfo( std::make_shared< ScalarPackInfo_T >( scalarPdfFieldGpuId ) );
         #else
            UniformBufferedScheme< LBMCommunicationStencil_T > communication(blocks);
            communication.addPackInfo( std::make_shared< PackInfo_T >( pdfFieldId ) );

            UniformBufferedScheme< ScalarLBMCommunicationStencil_T > scalarCommunication(blocks);
            scalarCommunication.addPackInfo( std::make_shared< ScalarPackInfo_T >( scalarPdfFieldId ) );
         #endif

         /*------------------------------------------------------------------------------------------------------
         *------------------------------------------------- TIMELOOP --------------------------------------------
         *------------------------------------------------------------------------------------------------------*/
         WALBERLA_LOG_INFO_ON_ROOT("Creating timeloop...")

         Timeloop_T timeloop(blocks->getBlockStorage(), simulationParameters.timesteps);

         MomentumThicknessPlotter< Timeloop_T > plotter( blocks, &simulationParameters, &timeloop,  velocityFieldId);
         plotter.addBeforeFunction([&]() {
            #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
               gpu::fieldCpy< VectorField_T, GPUField_T >(blocks, velocityFieldId, velocityFieldGpuId);
            #endif  
         });

         timeloop.addFuncBeforeTimeStep(plotter, "Momentum thickness plotter");

         timeloop.add() << BeforeFunction(communication.getCommunicateFunctor(),                         "lbm communication")
                        << Sweep(boundaryCollection.getSweep(BoundaryCollection_T::ALL),                 "lbm boundary conditions");
         timeloop.add() << Sweep(streamCollideSweep.streamCollide(SweepCollection_T::ALL),               "lbm stream and collide");
         timeloop.add() << BeforeFunction(scalarCommunication.getCommunicateFunctor(),                   "scalar communication")
                        << Sweep(scalarBoundaryCollection.getSweep(ScalarBoundaryCollection_T::ALL),     "scalar boundary conditions");
         timeloop.add() << Sweep(scalarStreamCollideSweep.streamCollide(ScalarSweepCollection_T::ALL),   "scalar stream and collide");

         //------------------------------------------------------------------------------------------------------//
         // vtk output
         auto vtkWriter = vtk::createVTKOutput_BlockData(
            blocks, "field_writer", simulationParameters.vtkFrequency, 0, false, simulationParameters.outputBaseFolder + "/vtk_out", "simulation_step",
            false, false, true, false
         );
         vtkWriter->setInitialWriteCallsToSkip(simulationParameters.vtkStart);

         // velocity field writer
         auto velocityWriter = std::make_shared<field::VTKWriter<VectorField_T>>(velocityFieldId, "velocity");
         vtkWriter->addCellDataWriter(velocityWriter);

         auto densityFieldWriter = std::make_shared<field::VTKWriter<ScalarField_T>>(densityFieldId, "density");
         vtkWriter->addCellDataWriter(densityFieldWriter);

         auto scalarFieldWriter = std::make_shared<field::VTKWriter<ScalarField_T>>(scalarFieldId, "scalar");
         vtkWriter->addCellDataWriter(scalarFieldWriter);

         vtkWriter->addBeforeFunction([&]() {
            #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
               gpu::fieldCpy< ScalarField_T, GPUField_T >(blocks, densityFieldId, densityFieldGpuId);
               gpu::fieldCpy< VectorField_T, GPUField_T >(blocks, velocityFieldId, velocityFieldGpuId);
               gpu::fieldCpy< ScalarField_T, GPUField_T >(blocks, scalarFieldId, scalarFieldGpuId);
            #endif  
         });

         timeloop.addFuncAfterTimeStep(vtk::writeFiles(vtkWriter), "VTK field output");

         {
         auto flagOutput = vtk::createVTKOutput_BlockData(
            blocks, "flag_writer", 1, 1, false, simulationParameters.outputBaseFolder + "/vtk_out", "simulation_step",
            false, true, true, false
         );
         auto flagWriter = std::make_shared<field::VTKWriter<FlagField_T>>(flagFieldId, "flag field");
         flagOutput->addCellDataWriter(flagWriter);
         flagOutput->write();
         }

         //------------------------------------------------------------------------------------------------------//
         // LBM stability check
         timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< PdfField_T, FlagField_T >(
                                       config, blocks, pdfFieldId, flagFieldId, FluidFlagUID ) ),
                                       "LBM stability check" );

         // Time logger
         timeloop.addFuncAfterTimeStep(timing::RemainingTimeLogger(timeloop.getNrOfTimeSteps(), 30_r),
                                       "remaining time logger");

         /*------------------------------------------------------------------------------------------------------
         *------------------------------------------------- RUN SIM --------------------------------------------
         *------------------------------------------------------------------------------------------------------*/
         WALBERLA_LOG_INFO_ON_ROOT("Running timeloop with " << timeloop.getNrOfTimeSteps() - 1 << " timesteps...")

         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            WALBERLA_GPU_CHECK(gpuPeekAtLastError())

            DeviceSynchronizeTimingPool timing;

            WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
            WALBERLA_GPU_CHECK( gpuPeekAtLastError() )
         #else
            WcTimingPool timing;
         #endif

         WcTimer timer;

         timer.start();
         timeloop.run(timing);
         timer.end();

         #ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
            WALBERLA_GPU_CHECK( gpuDeviceSynchronize() )
         #endif

         double time = timer.max();
         walberla::mpi::reduceInplace(time, walberla::mpi::MAX);

         WALBERLA_LOG_INFO_ON_ROOT("Simulation finished")

         const walberla::lbm::PerformanceEvaluation<FlagField_T> performance(blocks, flagFieldId, FluidFlagUID);
         performance.logResultOnRoot( simulationParameters.timesteps, time );

         timing.unifyRegisteredTimersAcrossProcesses();
         timing.logResultOnRoot( timing::REDUCE_TOTAL, true );
      }
            
      return EXIT_SUCCESS;
   }

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
