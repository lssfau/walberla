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
//! \file TorqueSpherePSMGPU.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Modification of pe_coupling/partially_saturated_cells_method/TorqueSpherePSM.cpp
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

#include "field/AddToStorage.h"

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
#include "mesa_pd/mpi/SyncNextNeighbors.h"

#include <iostream>
// codegen
#include "InitializeDomainForPSM.h"
#include "PSMPackInfo.h"
#include "PSMSweep.h"
#include "PSM_InfoHeader.h"

namespace torque_sphere_psm
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;
using namespace lbm_mesapd_coupling::psm::gpu;

using flag_t      = walberla::uint8_t;
using FlagField_T = FlagField< flag_t >;

using PackInfo_T = pystencils::PSMPackInfo;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag("fluid");

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
   real_t phi;
   real_t angularVel;
   real_t analyticalTorque;
};

//*******************************************************************************************************************
/*!\brief Evaluating the torque on a sphere, rotating with a constant angular velocity
 */
//*******************************************************************************************************************
template< typename ParticleAccessor_T >
class TorqueEval
{
 public:
   TorqueEval(SweepTimeloop* timeloop, Setup* setup, const shared_ptr< StructuredBlockStorage >& blocks,
              const shared_ptr< ParticleAccessor_T >& ac, bool fileIO)
      : timeloop_(timeloop), setup_(setup), blocks_(blocks), ac_(ac), fileIO_(fileIO)
   {
      // calculate the (semi)analytical torque value
      // see also Hofmann et al. - Hydrodynamic interactions in colloidal crystals:(II). Application to dense cubic and
      // tetragonal arrays (1999), Eqs. 5.1 and 5.5
      const real_t S = real_c(1.95708);
      setup_->analyticalTorque =
         -setup_->visc *
         (real_c(6) * setup_->phi / (real_c(1) - setup_->phi - S * std::pow(setup_->phi, real_c(10. / 3.)))) *
         setup_->angularVel * real_c(setup_->length * setup_->length * setup_->length);

      if (fileIO_)
      {
         std::ofstream file;
         filename_ = "TorqueSpherePSMGPU.txt";
         WALBERLA_ROOT_SECTION()
         {
            file.open(filename_.c_str());
            file << "#\t torqueSim\t torqueAnaly\n";
            file.close();
         }
      }
   }

   // evaluate the acting torque
   void operator()()
   {
      const uint_t timestep(timeloop_->getCurrentTimeStep() + 1);

      if (timestep % setup_->checkFrequency != 0) return;

      // update torque values
      torqueOld_ = torqueNew_;
      torqueNew_ = calculateTorque();

      // write to file if desired
      WALBERLA_ROOT_SECTION()
      {
         if (fileIO_)
         {
            std::ofstream file;
            file.open(filename_.c_str(), std::ofstream::app);
            file.setf(std::ios::unitbuf);
            file.precision(15);
            file << timestep << " " << torqueNew_ << " " << setup_->analyticalTorque << "\n";
            file.close();
         }
      }
   }

   // obtain the torque acting on the sphere by summing up all the process local parts
   real_t calculateTorque()
   {
      real_t torque = real_c(0);
      for (auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt)
      {
         for (size_t idx = 0; idx < ac_->size(); ++idx)
         {
            torque += ac_->getHydrodynamicTorque(idx)[1];
         }
      }

      WALBERLA_MPI_SECTION() { mpi::allReduceInplace(torque, mpi::SUM); }
      return torque;
   }

   // return the relative temporal change in the torque
   real_t getTorqueDiff() const { return std::fabs((torqueNew_ - torqueOld_) / torqueNew_); }

   // return the torque
   real_t getTorque() const { return torqueNew_; }

 private:
   SweepTimeloop* timeloop_;

   Setup* setup_;

   shared_ptr< StructuredBlockStorage > blocks_;
   shared_ptr< ParticleAccessor_T > ac_;

   bool fileIO_;
   std::string filename_;

   real_t torqueOld_{ 0.0 };
   real_t torqueNew_{ 0.0 };
};

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Testcase that checks the torque acting on a constantly rotating sphere in the center of a cubic domain
 *
 * The torque for this problem (often denoted as Simple Cubic setup) is given by a semi-analytical formula.
 * The cubic domain is periodic in all directions, making it a physically infinite periodic array of spheres.
   \verbatim
         _______________
        |       <-      |
        |      ___      |
        |     /   \     |
        |    |  x  |    |
        |     \___/     |
        |      ->       |
        |_______________|

   \endverbatim
 *
 * The collision model used for the LBM is TRT with a relaxation parameter tau=1.5 and the magic parameter 3/16.
 * The Stokes approximation of the equilibrium PDFs is used.
 * The sphere rotates with a angular velocity of 1e-5.
 * The domain is 32x32x32, and the sphere has a diameter of around 27 cells ( chi * domainlength )
 * The simulation is run until the relative change in the torque between 100 time steps is less than 1e-5.
 * The pe is not used since the angular velocity is kept constant and the force is explicitly reset after each time
 step.
 *
 */
//*******************************************************************************************************************

int main(int argc, char** argv)
{
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   logging::Logging::instance()->setLogLevel(logging::Logging::INFO);

   auto processes = MPIManager::instance()->numProcesses();

   if (processes != 1 && processes != 2 && processes != 4 && processes != 8)
   {
      std::cerr << "Number of processes must be equal to either 1, 2, 4, or 8!" << std::endl;
      return EXIT_FAILURE;
   }

   ///////////////////
   // Customization //
   ///////////////////

   bool shortrun = false;
   bool funcTest = false;
   bool fileIO   = false;

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
      if (std::strcmp(argv[i], "--fileIO") == 0)
      {
         fileIO = true;
         continue;
      }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   Setup setup;

   setup.length         = uint_t(32);   // length of the cubic domain in lattice cells
   const real_t chi     = real_c(0.85); // porosity parameter: diameter / length
   setup.tau            = real_c(1.5);  // relaxation time
   setup.angularVel     = real_c(1e-5); // angular velocity of the sphere
   setup.checkFrequency = uint_t(100);  // evaluate the torque only every checkFrequency time steps
   setup.radius         = real_c(0.5) * chi * real_c(setup.length); // sphere radius
   setup.visc           = (setup.tau - real_c(0.5)) / real_c(3);    // viscosity in lattice units
   setup.phi            = real_c(4.0 / 3.0) * math::pi * setup.radius * setup.radius * setup.radius /
               (real_c(setup.length * setup.length * setup.length)); // solid volume fraction
   const real_t omega            = real_c(1) / setup.tau;            // relaxation rate
   const real_t dx               = real_c(1);                        // lattice dx
   const real_t convergenceLimit = real_c(1e-5);                     // tolerance for relative change in torque
   const uint_t timesteps =
      funcTest ? 1 : (shortrun ? uint_c(150) : uint_c(5000)); // maximum number of time steps for the whole simulation

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

   ////////
   // PE //
   ////////

   auto mesapdDomain        = std::make_shared< mesa_pd::domain::BlockForestDomain >(blocks->getBlockForestPointer());
   auto ps                  = std::make_shared< mesa_pd::data::ParticleStorage >(1);
   auto ss                  = std::make_shared< mesa_pd::data::ShapeStorage >();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor            = walberla::make_shared< ParticleAccessor_T >(ps, ss);

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_c(1.5) * dx;

   if (setup.radius > real_c(setup.length) * real_c(0.5) - overlap)
   {
      std::cerr << "Periodic sphere is too large!" << std::endl;
      return EXIT_FAILURE;
   }

   // create the sphere in the middle of the domain
   Vector3< real_t > position(real_c(setup.length) * real_c(0.5));
   auto sphereShape = ss->create< mesa_pd::data::Sphere >(setup.radius);

   if (mesapdDomain->isContainedInProcessSubdomain(uint_c(walberla::mpi::MPIManager::instance()->rank()), position))
   {
      auto sphereParticle = ps->create();
      sphereParticle->setShapeID(sphereShape);
      sphereParticle->setType(0);
      sphereParticle->setPosition(position);
      sphereParticle->setAngularVelocity(Vector3(real_c(0), setup.angularVel, real_c(0)));
      sphereParticle->setOwner(walberla::MPIManager::instance()->rank());
      sphereParticle->setInteractionRadius(setup.radius);
   }

   // synchronize often enough for large particles
   std::function< void(void) > syncCall = [&]() {
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *mesapdDomain);
   };

   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add fields ( uInit = <0,0,0>, rhoInit = 1 )
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   BlockDataID pdfFieldID =
      field::addToStorage< PdfField_T >(blocks, "pdf field (fzyx)", real_c(std::nan("")), field::fzyx);
   BlockDataID pdfFieldCPUGPUID = gpu::addGPUFieldToStorage< PdfField_T >(blocks, pdfFieldID, "pdf field GPU", true);
#else
   BlockDataID pdfFieldCPUGPUID =
      field::addToStorage< PdfField_T >(blocks, "pdf field CPU", real_c(std::nan("")), field::fzyx);
#endif

   // add particle and volume fraction data structures
   ParticleAndVolumeFractionSoA_T< Weighting > particleAndVolumeFractionSoA(blocks, omega);
   // map particles and calculate solid volume fraction initially
   PSMSweepCollection psmSweepCollection(blocks, accessor, lbm_mesapd_coupling::RegularParticlesSelector(),
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

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   gpu::communication::UniformGPUScheme< Stencil_T > com(blocks, false, false);
#else
   walberla::blockforest::communication::UniformBufferedScheme< Stencil_T > com(blocks);
#endif
   com.addPackInfo(make_shared< PackInfo_T >(pdfFieldCPUGPUID));
   auto communication = std::function< void() >([&]() { com.communicate(); });

   using TorqueEval_T                    = TorqueEval< ParticleAccessor_T >;
   shared_ptr< TorqueEval_T > torqueEval = make_shared< TorqueEval_T >(&timeloop, &setup, blocks, accessor, fileIO);

   pystencils::PSMSweep PSMSweep(particleAndVolumeFractionSoA.BsFieldID, particleAndVolumeFractionSoA.BFieldID,
                                 particleAndVolumeFractionSoA.particleForcesFieldID,
                                 particleAndVolumeFractionSoA.particleVelocitiesFieldID, pdfFieldCPUGPUID, real_t(0.0),
                                 real_t(0.0), real_t(0.0), omega);

   // communication, streaming and force evaluation
   timeloop.add() << BeforeFunction(communication, "LBM Communication")
                  << Sweep(deviceSyncWrapper(psmSweepCollection.setParticleVelocitiesSweep),
                           "setParticleVelocitiesSweep");
   timeloop.add() << Sweep(deviceSyncWrapper(PSMSweep), "cell-wise LB sweep");
   timeloop.add() << Sweep(deviceSyncWrapper(psmSweepCollection.reduceParticleForcesSweep), "Reduce particle forces");
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   timeloop.add() << Sweep(gpu::fieldCpyFunctor< PdfField_T, gpu::GPUField< real_t > >(pdfFieldID, pdfFieldCPUGPUID),
                           "copy pdf from GPU to CPU")
#else
   struct emptySweep
   {
      void operator()(IBlock*) {}
   };
   timeloop.add() << Sweep(emptySweep(), "emptySweep")
#endif
                  << AfterFunction(SharedFunctor< TorqueEval_T >(torqueEval), "torque evaluation");

   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;

   timeloop.addFuncAfterTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   // time loop
   for (uint_t i = 0; i < timesteps; ++i)
   {
      // perform a single simulation step
      timeloop.singleStep(timeloopTiming);

      // resetting force
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor);

      // check if the relative change in the torque is below the specified convergence criterion
      if (i > setup.checkFrequency && torqueEval->getTorqueDiff() < convergenceLimit)
      {
         // if simulation has converged, terminate simulation
         break;
      }
   }

   timeloopTiming.logResultOnRoot();

   // check the result
   if (!funcTest && !shortrun)
   {
      real_t relErr = std::fabs((setup.analyticalTorque - torqueEval->getTorque()) / setup.analyticalTorque);
      if (fileIO)
      {
         WALBERLA_ROOT_SECTION()
         {
            std::cout << "Analytical torque: " << setup.analyticalTorque << "\n"
                      << "Simulated torque: " << torqueEval->getTorque() << "\n"
                      << "Relative error: " << relErr << "\n";
         }
      }
      // the relative error has to be below 10% (25% for SC2)
      WALBERLA_CHECK_LESS(relErr, (SC == 2) ? real_c(0.25) : real_c(0.1));
   }

   return 0;
}

} // namespace torque_sphere_psm

int main(int argc, char** argv) { torque_sphere_psm::main(argc, argv); }
