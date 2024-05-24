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
//! \file DragForceSpherePSMCPUGPU.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Modification of partially_saturated_cells_method/DragForceSphere.cpp
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/vtk/VTKWriter.h"

#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/ParticleSelector.h"

#include <iostream>

// codegen
#include "InitializeDomainForPSM.h"
#include "PSMPackInfo.h"
#include "PSMSweep.h"
#include "PSM_InfoHeader.h"
#include "PSM_MacroGetter.h"

namespace drag_force_sphere_psm
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;
using namespace lbm_mesapd_coupling::psm::gpu;

typedef pystencils::PSMPackInfo PackInfo_T;

////////////////
// PARAMETERS //
////////////////

struct Setup
{
   uint_t checkFrequency;
   real_t visc;
   real_t tau;
   real_t radius;
   uint_t length;
   real_t chi;
   real_t extForce;
   real_t analyticalDrag;
};

template< typename ParticleAccessor_T >
class DragForceEvaluator
{
 public:
   DragForceEvaluator(SweepTimeloop* timeloop, Setup* setup, const shared_ptr< StructuredBlockStorage >& blocks,
                      const BlockDataID& velocityFieldID, const shared_ptr< ParticleAccessor_T >& ac,
                      walberla::id_t sphereID)
      : timeloop_(timeloop), setup_(setup), blocks_(blocks), velocityFieldID_(velocityFieldID), ac_(ac),
        sphereID_(sphereID), normalizedDragOld_(0.0), normalizedDragNew_(0.0)
   {
      // calculate the analytical drag force value based on the series expansion of chi
      // see also Sangani - Slow flow through a periodic array of spheres, IJMF 1982. Eq. 60 and Table 1
      real_t analyticalDrag = real_c(0);
      real_t tempChiPowS    = real_c(1);

      // coefficients to calculate the drag in a series expansion
      real_t dragCoefficients[31] = { real_c(1.000000),  real_c(1.418649),  real_c(2.012564),   real_c(2.331523),
                                      real_c(2.564809),  real_c(2.584787),  real_c(2.873609),   real_c(3.340163),
                                      real_c(3.536763),  real_c(3.504092),  real_c(3.253622),   real_c(2.689757),
                                      real_c(2.037769),  real_c(1.809341),  real_c(1.877347),   real_c(1.534685),
                                      real_c(0.9034708), real_c(0.2857896), real_c(-0.5512626), real_c(-1.278724),
                                      real_c(1.013350),  real_c(5.492491),  real_c(4.615388),   real_c(-0.5736023),
                                      real_c(-2.865924), real_c(-4.709215), real_c(-6.870076),  real_c(0.1455304),
                                      real_c(12.51891),  real_c(9.742811),  real_c(-5.566269) };

      for (uint_t s = 0; s <= uint_t(30); ++s)
      {
         analyticalDrag += dragCoefficients[s] * tempChiPowS;
         tempChiPowS *= setup->chi;
      }
      setup_->analyticalDrag = analyticalDrag;
   }

   // evaluate the acting drag force
   void operator()()
   {
      const uint_t timestep(timeloop_->getCurrentTimeStep() + 1);

      if (timestep % setup_->checkFrequency != 0) return;

      // get force in x-direction acting on the sphere
      real_t forceX = computeDragForce();
      // get average volumetric flowrate in the domain
      real_t uBar = computeAverageVel();

      // f_total = f_drag + f_buoyancy
      real_t totalForce =
         forceX + real_c(4.0 / 3.0) * math::pi * setup_->radius * setup_->radius * setup_->radius * setup_->extForce;

      real_t normalizedDragForce = totalForce / real_c(6.0 * math::pi * setup_->visc * setup_->radius * uBar);

      // update drag force values
      normalizedDragOld_ = normalizedDragNew_;
      normalizedDragNew_ = normalizedDragForce;
   }

   // return the relative temporal change in the normalized drag
   real_t getDragForceDiff() const { return std::fabs((normalizedDragNew_ - normalizedDragOld_) / normalizedDragNew_); }

   // return the drag force
   real_t getDragForce() const { return normalizedDragNew_; }

   void logResultToFile(const std::string& filename) const
   {
      // write to file if desired
      // format: length tau viscosity simulatedDrag analyticalDrag\n
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open(filename.c_str(), std::ofstream::app);
         file.precision(8);
         file << setup_->length << " " << setup_->tau << " " << setup_->visc << " " << normalizedDragNew_ << " "
              << setup_->analyticalDrag << "\n";
         file.close();
      }
   }

 private:
   // obtain the drag force acting on the sphere by summing up all the process local parts of fX
   real_t computeDragForce()
   {
      size_t idx = ac_->uidToIdx(sphereID_);
      WALBERLA_ASSERT_UNEQUAL(idx, ac_->getInvalidIdx(), "Index of particle is invalid!");
      real_t force = real_t(0);
      if (idx != ac_->getInvalidIdx()) { force = ac_->getHydrodynamicForce(idx)[0]; }

      WALBERLA_MPI_SECTION() { mpi::allReduceInplace(force, mpi::SUM); }

      return force;
   }

   // calculate the average velocity in forcing direction (here: x) inside the domain (assuming dx=1)
   real_t computeAverageVel()
   {
      auto velocity_sum = real_t(0);
      // iterate all blocks stored locally on this process
      for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt)
      {
         // retrieve the pdf field and the flag field from the block
         VelocityField_T* velocityField = blockIt->getData< VelocityField_T >(velocityFieldID_);

         // get the flag that marks a cell as being fluid

         auto xyzField = velocityField->xyzSize();
         for (auto cell : xyzField)
         {
            // TODO: fix velocity computation by using getPSMMacroscopicVelocity
            velocity_sum += velocityField->get(cell, 0);
         }
      }

      WALBERLA_MPI_SECTION() { mpi::allReduceInplace(velocity_sum, mpi::SUM); }

      return velocity_sum / real_c(setup_->length * setup_->length * setup_->length);
   }

   SweepTimeloop* timeloop_;

   Setup* setup_;

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID velocityFieldID_;

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereID_;

   // drag coefficient
   real_t normalizedDragOld_;
   real_t normalizedDragNew_;
};

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Testcase that checks the drag force acting on a fixed sphere in the center of a cubic domain in Stokes flow
 *
 * The drag force for this problem (often denoted as Simple Cubic setup) is given by a semi-analytical series expansion.
 * The cubic domain is periodic in all directions, making it a physically infinite periodic array of spheres.
   \verbatim
           _______________
        ->|               |->
        ->|      ___      |->
      W ->|     /   \     |-> E
      E ->|    |  x  |    |-> A
      S ->|     \___/     |-> S
      T ->|               |-> T
        ->|_______________|->

   \endverbatim
 *
 * The collision model used for the LBM is TRT with a relaxation parameter tau=1.5 and the magic parameter 3/16.
 * The Stokes approximation of the equilibrium PDFs is used.
 * The flow is driven by a constant particle force of 1e-5.
 * The domain is 32x32x32, and the sphere has a diameter of 16 cells ( chi * domainlength )
 * The simulation is run until the relative change in the dragforce between 100 time steps is less than 1e-5.
 * The RPD is not used since the sphere is kept fixed and the force is explicitly reset after each time step.
 * To avoid periodicity constrain problems, the sphere is set as global.
 *
 */
//*******************************************************************************************************************

int main(int argc, char** argv)
{
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   auto processes = MPIManager::instance()->numProcesses();

   if (processes != 1 && processes != 2 && processes != 4 && processes != 8)
   {
      std::cerr << "Number of processes must be equal to either 1, 2, 4, or 8!" << std::endl;
      return EXIT_FAILURE;
   }

   ///////////////////
   // Customization //
   ///////////////////

   bool shortrun       = false;
   bool funcTest       = false;
   bool logging        = false;
   uint_t vtkFrequency = uint_c(0);
   real_t tau          = real_c(1.5);
   uint_t length       = uint_c(32);

   for (int i = 1; i < argc; ++i)
   {
      if (std::strcmp(argv[i], "--shortrun") == 0)
      {
         shortrun = true;
         continue;
      }
      if (std::strcmp(argv[i], "--funcTest") == 0)
      {
         funcTest = true;
         continue;
      }
      if (std::strcmp(argv[i], "--logging") == 0)
      {
         logging = true;
         continue;
      }
      if (std::strcmp(argv[i], "--tau") == 0)
      {
         tau = real_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--length") == 0)
      {
         length = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--vtkFrequency") == 0)
      {
         vtkFrequency = uint_c(std::atof(argv[++i]));
         continue;
      }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   Setup setup;

   setup.length                  = length;       // length of the cubic domain in lattice cells
   setup.chi                     = real_c(0.5);  // porosity parameter: diameter / length
   setup.tau                     = tau;          // relaxation time
   setup.extForce                = real_c(1e-7); // constant particle force in lattice units
   setup.checkFrequency          = uint_t(100);  // evaluate the drag force only every checkFrequency time steps
   setup.radius                  = real_c(0.5) * setup.chi * real_c(setup.length); // sphere radius
   setup.visc                    = (setup.tau - real_c(0.5)) / real_c(3);          // viscosity in lattice units
   const real_t omega            = real_c(1) / setup.tau;                          // relaxation rate
   const real_t dx               = real_c(1);                                      // lattice dx
   const real_t convergenceLimit = real_c(1e-7); // tolerance for relative change in drag force
   const uint_t timesteps =
      funcTest ? 1 : (shortrun ? uint_c(150) : uint_c(50000)); // maximum number of time steps for the whole simulation

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = (processes >= 2) ? uint_t(2) : uint_t(1);
   const uint_t YBlocks = (processes >= 4) ? uint_t(2) : uint_t(1);
   const uint_t ZBlocks = (processes == 8) ? uint_t(2) : uint_t(1);
   const uint_t XCells  = setup.length / XBlocks;
   const uint_t YCells  = setup.length / YBlocks;
   const uint_t ZCells  = setup.length / ZBlocks;

   // create fully periodic domain
   auto blocks = blockforest::createUniformBlockGrid(XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells, dx, true, true,
                                                     true, true);

   /////////
   // RPD //
   /////////

   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());

   // init data structures
   auto ps                  = std::make_shared< mesa_pd::data::ParticleStorage >(1);
   auto ss                  = std::make_shared< mesa_pd::data::ShapeStorage >();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor            = make_shared< ParticleAccessor_T >(ps, ss);
   auto sphereShape         = ss->create< mesa_pd::data::Sphere >(setup.radius);

   //////////////////
   // RPD COUPLING //
   //////////////////

   // connect to pe
   const real_t overlap = real_t(1.5) * dx;

   if (setup.radius > real_c(setup.length) * real_t(0.5) - overlap)
   {
      std::cerr << "Periodic sphere is too large and would lead to incorrect mapping!" << std::endl;
      // solution: create the periodic copies explicitly
      return EXIT_FAILURE;
   }

   // create the sphere in the middle of the domain
   // it is global and thus present on all processes
   Vector3< real_t > position(real_c(setup.length) * real_c(0.5));
   walberla::id_t sphereID;
   {
      mesa_pd::data::Particle&& p = *ps->create(true);
      p.setPosition(position);
      p.setInteractionRadius(setup.radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereID = p.getUid();
   }

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add fields
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   BlockDataID pdfFieldID =
      field::addToStorage< PdfField_T >(blocks, "pdf field (fzyx)", real_c(std::nan("")), field::fzyx);
   BlockDataID BFieldID         = field::addToStorage< BField_T >(blocks, "B", real_t(0), field::fzyx);
   BlockDataID pdfFieldCPUGPUID = gpu::addGPUFieldToStorage< PdfField_T >(blocks, pdfFieldID, "pdf field GPU");
#else
   BlockDataID pdfFieldCPUGPUID =
      field::addToStorage< PdfField_T >(blocks, "pdf field CPU", real_c(std::nan("")), field::fzyx);
#endif

   BlockDataID densityFieldID = field::addToStorage< DensityField_T >(blocks, "Density", real_t(0), field::fzyx);
   BlockDataID velFieldID     = field::addToStorage< VelocityField_T >(blocks, "Velocity", real_t(0), field::fzyx);

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   gpu::communication::UniformGPUScheme< Stencil_T > com(blocks, 0, false);
#else
   walberla::blockforest::communication::UniformBufferedScheme< Stencil_T > com(blocks);
#endif
   com.addPackInfo(make_shared< PackInfo_T >(pdfFieldCPUGPUID));
   auto communication = std::function< void() >([&]() { com.communicate(); });

   // add particle and volume fraction data structures
   ParticleAndVolumeFractionSoA_T< 1 > particleAndVolumeFractionSoA(blocks, omega);

   // map particles and calculate solid volume fraction initially
   PSMSweepCollection psmSweepCollection(blocks, accessor, lbm_mesapd_coupling::GlobalParticlesSelector(),
                                         particleAndVolumeFractionSoA, Vector3(8));
   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      psmSweepCollection.particleMappingSweep(&(*blockIt));
   }

   pystencils::InitializeDomainForPSM pdfSetter(
      particleAndVolumeFractionSoA.BsFieldID, particleAndVolumeFractionSoA.BFieldID,
      particleAndVolumeFractionSoA.particleVelocitiesFieldID, pdfFieldCPUGPUID, real_t(0), real_t(0), real_t(0),
      real_t(1.0), real_t(0), real_t(0), real_t(0));

   for (auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      // pdfSetter requires particle velocities at cell centers
      psmSweepCollection.setParticleVelocitiesSweep(&(*blockIt));
      pdfSetter(&(*blockIt));
   }

   // since external forcing is applied, the evaluation of the velocity has to be carried out directly after the
   // streaming step however, the default sweep is a  stream - collide step, i.e. after the sweep, the velocity
   // evaluation is not correct solution: split the sweep explicitly into collide and stream
   pystencils::PSMSweep PSMSweep(particleAndVolumeFractionSoA.BsFieldID, particleAndVolumeFractionSoA.BFieldID,
                                 particleAndVolumeFractionSoA.particleForcesFieldID,
                                 particleAndVolumeFractionSoA.particleVelocitiesFieldID, pdfFieldCPUGPUID,
                                 setup.extForce, real_t(0.0), real_t(0.0), omega);

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils::PSM_MacroGetter getterSweep(BFieldID, densityFieldID, pdfFieldID, velFieldID, setup.extForce,
                                           real_t(0.0), real_t(0.0));
#else
   pystencils::PSM_MacroGetter getterSweep(particleAndVolumeFractionSoA.BFieldID, densityFieldID, pdfFieldCPUGPUID,
                                           velFieldID, setup.extForce, real_t(0.0), real_t(0.0));
#endif

   // add LBM communication function and streaming & force evaluation
   using DragForceEval_T = DragForceEvaluator< ParticleAccessor_T >;
   auto forceEval        = make_shared< DragForceEval_T >(&timeloop, &setup, blocks, velFieldID, accessor, sphereID);
   timeloop.add() << BeforeFunction(communication, "LBM Communication")
                  << Sweep(deviceSyncWrapper(psmSweepCollection.setParticleVelocitiesSweep), "Set particle velocities");
   timeloop.add() << Sweep(deviceSyncWrapper(PSMSweep), "cell-wise PSM sweep");
   timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.reduceParticleForcesSweep), "Reduce particle forces");
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   timeloop.add() << Sweep(gpu::fieldCpyFunctor< PdfField_T, gpu::GPUField< real_t > >(pdfFieldID, pdfFieldCPUGPUID),
                           "copy pdf from GPU to CPU");
   timeloop.add() << Sweep(
      gpu::fieldCpyFunctor< BField_T, gpu::GPUField< real_t > >(BFieldID, particleAndVolumeFractionSoA.BFieldID),
      "copy B field from GPU to CPU");
#endif
   timeloop.add() << Sweep(getterSweep, "compute velocity")
                  << AfterFunction(SharedFunctor< DragForceEval_T >(forceEval), "drag force evaluation");

   // resetting force
   timeloop.addFuncAfterTimeStep(
      [ps, accessor]() {
         ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor,
                             lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel(), *accessor);
      },
      "reset force on sphere");

   timeloop.addFuncAfterTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   if (vtkFrequency > 0)
   {
      const std::string path = "vtk_out/dragForceSphereCPUGPU";
      auto vtkOutput = vtk::createVTKOutput_BlockData(*blocks, "psm_velocity_fieldCPUGPU", vtkFrequency, 0, false, path,
                                                      "simulation_step", false, true, true, false, 0);

      vtkOutput->addBeforeFunction([&]() {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         gpu::fieldCpy< PdfField_T, gpu::GPUField< real_t > >(blocks, pdfFieldID, pdfFieldCPUGPUID);
         gpu::fieldCpy< BField_T, gpu::GPUField< real_t > >(blocks, BFieldID, particleAndVolumeFractionSoA.BFieldID);
#endif
         for (auto& block : *blocks)
            getterSweep(&block);
      });
      vtkOutput->addCellDataWriter(make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "Velocity"));

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(vtkOutput), "VTK Output");
   }

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   // time loop
   for (uint_t i = 0; i < timesteps; ++i)
   {
      // perform a single simulation step
      timeloop.singleStep(timeloopTiming);

      // check if the relative change in the normalized drag force is below the specified convergence criterion
      if (i > setup.checkFrequency && forceEval->getDragForceDiff() < convergenceLimit)
      {
         // if simulation has converged, terminate simulation
         break;
      }
   }

   timeloopTiming.logResultOnRoot();

   if (!funcTest && !shortrun)
   {
      // check the result
      real_t relErr = std::fabs((setup.analyticalDrag - forceEval->getDragForce()) / setup.analyticalDrag);
      if (logging)
      {
         WALBERLA_ROOT_SECTION()
         {
            std::cout << "Analytical drag: " << setup.analyticalDrag << "\n"
                      << "Simulated drag: " << forceEval->getDragForce() << "\n"
                      << "Relative error: " << relErr << "\n";
         }
         forceEval->logResultToFile("log_DragForceSphere.txt");
      }
      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS(relErr, real_c(0.1));
   }

   return 0;
}

} // namespace drag_force_sphere_psm

int main(int argc, char** argv) { drag_force_sphere_psm::main(argc, argv); }
