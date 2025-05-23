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
//! \file SettlingSpherePSMGPU.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Modification of momentum_exchange_method/SettlingSphere.cpp
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/all.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"

#include "field/AddToStorage.h"
#include "field/vtk/all.h"

#include "geometry/InitBoundaryHandling.h"

#include "gpu/AddGPUFieldToStorage.h"
#include "gpu/DeviceSelectMPI.h"
#include "gpu/communication/UniformGPUScheme.h"

#include "lbm/field/AddToStorage.h"
#include "lbm/vtk/all.h"

#include "lbm_mesapd_coupling/DataTypesCodegen.h"
#include "lbm_mesapd_coupling/partially_saturated_cells_method/codegen/PSMSweepCollection.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/SpringDashpot.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "vtk/all.h"

#include <functional>

#include "InitializeDomainForPSM.h"
#include "PSMPackInfo.h"
#include "PSMSweep.h"
#include "PSM_InfoHeader.h"
#include "PSM_MacroGetter.h"
#include "PSM_NoSlip.h"

namespace settling_sphere
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

const FlagUID Fluid_Flag("Fluid");
const FlagUID NoSlip_Flag("NoSlip");

//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the sphere
 *
 */
//*******************************************************************************************************************
template< typename ParticleAccessor_T >
class SpherePropertyLogger
{
 public:
   SpherePropertyLogger(const shared_ptr< ParticleAccessor_T >& ac, walberla::id_t sphereUid,
                        const std::string& fileName, bool fileIO, real_t dx_SI, real_t dt_SI, real_t diameter,
                        real_t gravitationalForceMag)
      : ac_(ac), sphereUid_(sphereUid), fileName_(fileName), fileIO_(fileIO), dx_SI_(dx_SI), dt_SI_(dt_SI),
        diameter_(diameter), gravitationalForceMag_(gravitationalForceMag), position_(real_t(0)),
        maxVelocity_(real_t(0))
   {
      if (fileIO_)
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open(fileName_.c_str());
            file << "#\t t\t posX\t posY\t gapZ\t velX\t velY\t velZ\n";
            file.close();
         }
      }
   }

   void operator()(const uint_t timestep)
   {
      Vector3< real_t > pos(real_t(0));
      Vector3< real_t > transVel(real_t(0));
      Vector3< real_t > hydForce(real_t(0));

      size_t idx = ac_->uidToIdx(sphereUid_);
      if (idx != ac_->getInvalidIdx())
      {
         if (!mesa_pd::data::particle_flags::isSet(ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            pos      = ac_->getPosition(idx);
            transVel = ac_->getLinearVelocity(idx);
            hydForce = ac_->getHydrodynamicForce(idx);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace(pos[0], mpi::SUM);
         mpi::allReduceInplace(pos[1], mpi::SUM);
         mpi::allReduceInplace(pos[2], mpi::SUM);

         mpi::allReduceInplace(transVel[0], mpi::SUM);
         mpi::allReduceInplace(transVel[1], mpi::SUM);
         mpi::allReduceInplace(transVel[2], mpi::SUM);

         mpi::allReduceInplace(hydForce[0], mpi::SUM);
         mpi::allReduceInplace(hydForce[1], mpi::SUM);
         mpi::allReduceInplace(hydForce[2], mpi::SUM);
      }

      position_    = pos[2];
      maxVelocity_ = std::max(maxVelocity_, -transVel[2]);

      if (fileIO_) writeToFile(timestep, pos, transVel, hydForce);
   }

   real_t getPosition() const { return position_; }

   real_t getMaxVelocity() const { return maxVelocity_; }

 private:
   void writeToFile(const uint_t timestep, const Vector3< real_t >& position, const Vector3< real_t >& velocity,
                    const Vector3< real_t >& hydForce)
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open(fileName_.c_str(), std::ofstream::app);

         auto scaledPosition     = position / diameter_;
         auto velocity_SI        = velocity * dx_SI_ / dt_SI_;
         auto normalizedHydForce = hydForce / gravitationalForceMag_;

         file << timestep << "\t" << real_c(timestep) * dt_SI_ << "\t" << "\t" << scaledPosition[0] << "\t"
              << scaledPosition[1] << "\t" << scaledPosition[2] - real_t(0.5) << "\t" << velocity_SI[0] << "\t"
              << velocity_SI[1] << "\t" << velocity_SI[2] << "\t" << normalizedHydForce[0] << "\t"
              << normalizedHydForce[1] << "\t" << normalizedHydForce[2] << "\n";
         file.close();
      }
   }

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;
   std::string fileName_;
   bool fileIO_;
   real_t dx_SI_, dt_SI_, diameter_, gravitationalForceMag_;

   real_t position_;
   real_t maxVelocity_;
};

void createPlaneSetup(const shared_ptr< mesa_pd::data::ParticleStorage >& ps,
                      const shared_ptr< mesa_pd::data::ShapeStorage >& ss, const math::AABB& simulationDomain)
{
   // create bounding planes
   mesa_pd::data::Particle p0 = *ps->create(true);
   p0.setPosition(simulationDomain.minCorner());
   p0.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p0.setShapeID(ss->create< mesa_pd::data::HalfSpace >(Vector3< real_t >(0, 0, 1)));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p1 = *ps->create(true);
   p1.setPosition(simulationDomain.maxCorner());
   p1.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p1.setShapeID(ss->create< mesa_pd::data::HalfSpace >(Vector3< real_t >(0, 0, -1)));
   p1.setOwner(mpi::MPIManager::instance()->rank());
   p1.setType(0);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p2 = *ps->create(true);
   p2.setPosition(simulationDomain.minCorner());
   p2.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p2.setShapeID(ss->create< mesa_pd::data::HalfSpace >(Vector3< real_t >(1, 0, 0)));
   p2.setOwner(mpi::MPIManager::instance()->rank());
   p2.setType(0);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p3 = *ps->create(true);
   p3.setPosition(simulationDomain.maxCorner());
   p3.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p3.setShapeID(ss->create< mesa_pd::data::HalfSpace >(Vector3< real_t >(-1, 0, 0)));
   p3.setOwner(mpi::MPIManager::instance()->rank());
   p3.setType(0);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p4 = *ps->create(true);
   p4.setPosition(simulationDomain.minCorner());
   p4.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p4.setShapeID(ss->create< mesa_pd::data::HalfSpace >(Vector3< real_t >(0, 1, 0)));
   p4.setOwner(mpi::MPIManager::instance()->rank());
   p4.setType(0);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p5 = *ps->create(true);
   p5.setPosition(simulationDomain.maxCorner());
   p5.setInteractionRadius(std::numeric_limits< real_t >::infinity());
   p5.setShapeID(ss->create< mesa_pd::data::HalfSpace >(Vector3< real_t >(0, -1, 0)));
   p5.setOwner(mpi::MPIManager::instance()->rank());
   p5.setType(0);
   mesa_pd::data::particle_flags::set(p5.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p5.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
}

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Testcase that simulates the settling of a sphere inside a rectangular column filled with viscous fluid
 *
 * see: ten Cate, Nieuwstad, Derksen, Van den Akker - "Particle imaging velocimetry experiments and lattice-Boltzmann
 * simulations on a single sphere settling under gravity" (2002), Physics of Fluids, doi: 10.1063/1.1512918
 */
//*******************************************************************************************************************

int main(int argc, char** argv)
{
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   logging::Logging::instance()->setLogLevel(logging::Logging::INFO);

   ///////////////////
   // Customization //
   ///////////////////

   // simulation control
   bool shortrun          = false;
   bool funcTest          = false;
   bool fileIO            = false;
   uint_t vtkIOFreq       = 0;
   std::string baseFolder = "vtk_out_SettlingSphere_CPU_GPU";

   // physical setup
   uint_t fluidType = 1;

   // numerical parameters
   uint_t numberOfCellsInHorizontalDirection = uint_t(135);
   bool averageForceTorqueOverTwoTimeSteps   = true;
   uint_t numRPDSubCycles                    = uint_t(1);
   bool useVelocityVerlet                    = false;

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
      if (std::strcmp(argv[i], "--vtkIOFreq") == 0)
      {
         vtkIOFreq = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--fluidType") == 0)
      {
         fluidType = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--numRPDSubCycles") == 0)
      {
         numRPDSubCycles = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--resolution") == 0)
      {
         numberOfCellsInHorizontalDirection = uint_c(std::atof(argv[++i]));
         continue;
      }
      if (std::strcmp(argv[i], "--noForceAveraging") == 0)
      {
         averageForceTorqueOverTwoTimeSteps = false;
         continue;
      }
      if (std::strcmp(argv[i], "--baseFolder") == 0)
      {
         baseFolder = argv[++i];
         continue;
      }
      if (std::strcmp(argv[i], "--useVV") == 0)
      {
         useVelocityVerlet = true;
         continue;
      }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if (funcTest) { walberla::logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::WARNING); }

   if (fileIO)
   {
      // create base directory if it does not yet exist
      filesystem::path tpath(baseFolder);
      if (!filesystem::exists(tpath)) filesystem::create_directory(tpath);
   }

   //////////////////////////////////////
   // SIMULATION PROPERTIES in SI units//
   //////////////////////////////////////

   // values are mainly taken from the reference paper
   const real_t diameter_SI      = real_t(15e-3);
   const real_t densitySphere_SI = real_t(1120);

   real_t densityFluid_SI, dynamicViscosityFluid_SI;
   real_t expectedSettlingVelocity_SI;
   switch (fluidType)
   {
   case 1:
      // Re_p around 1.5
      densityFluid_SI             = real_t(970);
      dynamicViscosityFluid_SI    = real_t(373e-3);
      expectedSettlingVelocity_SI = real_t(0.035986);
      break;
   case 2:
      // Re_p around 4.1
      densityFluid_SI             = real_t(965);
      dynamicViscosityFluid_SI    = real_t(212e-3);
      expectedSettlingVelocity_SI = real_t(0.05718);
      break;
   case 3:
      // Re_p around 11.6
      densityFluid_SI             = real_t(962);
      dynamicViscosityFluid_SI    = real_t(113e-3);
      expectedSettlingVelocity_SI = real_t(0.087269);
      break;
   case 4:
      // Re_p around 31.9
      densityFluid_SI             = real_t(960);
      dynamicViscosityFluid_SI    = real_t(58e-3);
      expectedSettlingVelocity_SI = real_t(0.12224);
      break;
   default:
      WALBERLA_ABORT("Only four different fluids are supported! Choose type between 1 and 4.");
   }
   const real_t kinematicViscosityFluid_SI = dynamicViscosityFluid_SI / densityFluid_SI;

   const real_t gravitationalAcceleration_SI = real_t(9.81);
   Vector3< real_t > domainSize_SI(real_t(100e-3), real_t(100e-3), real_t(160e-3));
   // shift starting gap a bit upwards to match the reported (plotted) values
   const real_t startingGapSize_SI = real_t(120e-3) + real_t(0.25) * diameter_SI;

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in SI units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid type = " << fluidType);
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter_SI << ", density = " << densitySphere_SI
                                                      << ", starting gap size = " << startingGapSize_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid_SI << ", dyn. visc = " << dynamicViscosityFluid_SI
                                                    << ", kin. visc = " << kinematicViscosityFluid_SI);
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = "
                             << expectedSettlingVelocity_SI << " --> Re_p = "
                             << expectedSettlingVelocity_SI * diameter_SI / kinematicViscosityFluid_SI);

   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   const real_t dx_SI = domainSize_SI[0] / real_c(numberOfCellsInHorizontalDirection);
   const Vector3< uint_t > domainSize(uint_c(floor(domainSize_SI[0] / dx_SI + real_t(0.5))),
                                      uint_c(floor(domainSize_SI[1] / dx_SI + real_t(0.5))),
                                      uint_c(floor(domainSize_SI[2] / dx_SI + real_t(0.5))));
   const real_t diameter     = diameter_SI / dx_SI;
   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;

   const real_t expectedSettlingVelocity = real_t(0.01);
   const real_t dt_SI                    = expectedSettlingVelocity / expectedSettlingVelocity_SI * dx_SI;

   const real_t viscosity      = kinematicViscosityFluid_SI * dt_SI / (dx_SI * dx_SI);
   const real_t relaxationTime = real_t(1) / lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t gravitationalAcceleration = gravitationalAcceleration_SI * dt_SI * dt_SI / dx_SI;

   const real_t densityFluid  = real_t(1);
   const real_t densitySphere = densityFluid * densitySphere_SI / densityFluid_SI;

   const real_t dx = real_t(1);

   const uint_t timesteps = funcTest ? 1 : (shortrun ? uint_t(200) : uint_t(250000));

   WALBERLA_LOG_INFO_ON_ROOT(" - dx_SI = " << dx_SI << ", dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter << ", density = " << densitySphere);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid << ", relaxation time (tau) = " << relaxationTime
                                                    << ", kin. visc = " << viscosity);
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration);
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedSettlingVelocity << " --> Re_p = "
                                                                << expectedSettlingVelocity * diameter / viscosity);
   WALBERLA_LOG_INFO_ON_ROOT(" - integrator = " << (useVelocityVerlet ? "Velocity Verlet" : "Explicit Euler"));

   if (vtkIOFreq > 0)
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder << "\" with frequency " << vtkIOFreq);
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3< uint_t > numberOfBlocksPerDirection(uint_t(1), uint_t(1), uint_t(MPIManager::instance()->numProcesses()));
   Vector3< uint_t > cellsPerBlockPerDirection(domainSize[0] / numberOfBlocksPerDirection[0],
                                               domainSize[1] / numberOfBlocksPerDirection[1],
                                               domainSize[2] / numberOfBlocksPerDirection[2]);
   WALBERLA_CHECK_EQUAL(
      numberOfBlocksPerDirection[0] * numberOfBlocksPerDirection[1] * numberOfBlocksPerDirection[2],
      uint_t(MPIManager::instance()->numProcesses()),
      "When using GPUs, the number of blocks ("
         << numberOfBlocksPerDirection[0] * numberOfBlocksPerDirection[1] * numberOfBlocksPerDirection[2]
         << ") has to match the number of MPI processes (" << uint_t(MPIManager::instance()->numProcesses()) << ")");

   for (uint_t i = 0; i < 3; ++i)
   {
      WALBERLA_CHECK_EQUAL(cellsPerBlockPerDirection[i] * numberOfBlocksPerDirection[i], domainSize[i],
                           "Unmatching domain decomposition in direction " << i << "!");
   }

   auto blocks = blockforest::createUniformBlockGrid(numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1],
                                                     numberOfBlocksPerDirection[2], cellsPerBlockPerDirection[0],
                                                     cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], dx, 0,
                                                     false, false, false, false, false, // periodicity
                                                     false);

   WALBERLA_LOG_INFO_ON_ROOT("Domain decomposition:");
   WALBERLA_LOG_INFO_ON_ROOT(" - blocks per direction = " << numberOfBlocksPerDirection);
   WALBERLA_LOG_INFO_ON_ROOT(" - cells per block = " << cellsPerBlockPerDirection);

   // write domain decomposition to file
   if (vtkIOFreq > 0) { vtk::writeDomainDecomposition(blocks, "initial_domain_decomposition", baseFolder); }

   //////////////////
   // RPD COUPLING //
   //////////////////

   auto rpdDomain = std::make_shared< mesa_pd::domain::BlockForestDomain >(blocks->getBlockForestPointer());

   // init data structures
   auto ps                  = walberla::make_shared< mesa_pd::data::ParticleStorage >(1);
   auto ss                  = walberla::make_shared< mesa_pd::data::ShapeStorage >();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor            = walberla::make_shared< ParticleAccessor_T >(ps, ss);

   // bounding planes
   createPlaneSetup(ps, ss, blocks->getDomain());

   // create sphere and store Uid
   Vector3< real_t > initialPosition(real_t(0.5) * real_c(domainSize[0]), real_t(0.5) * real_c(domainSize[1]),
                                     startingGapSize_SI / dx_SI + real_t(0.5) * diameter);
   auto sphereShape = ss->create< mesa_pd::data::Sphere >(diameter * real_t(0.5));
   ss->shapes[sphereShape]->updateMassAndInertia(densitySphere);

   walberla::id_t sphereUid = 0;
   if (rpdDomain->isContainedInProcessSubdomain(uint_c(mpi::MPIManager::instance()->rank()), initialPosition))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(initialPosition);
      p.setInteractionRadius(diameter * real_t(0.5));
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      sphereUid = p.getUid();
   }
   mpi::allReduceInplace(sphereUid, mpi::SUM);

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add PDF field
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   BlockDataID pdfFieldID =
      field::addToStorage< PdfField_T >(blocks, "pdf field (fzyx)", real_c(std::nan("")), field::fzyx);
   BlockDataID BFieldID         = field::addToStorage< BField_T >(blocks, "B field CPU", 0, field::fzyx, 1);
   BlockDataID BsFieldID        = field::addToStorage< BsField_T >(blocks, "Bs field CPU", 0, field::fzyx, 1);
   BlockDataID pdfFieldCPUGPUID = gpu::addGPUFieldToStorage< PdfField_T >(blocks, pdfFieldID, "pdf field GPU");
#else
   BlockDataID pdfFieldCPUGPUID =
      field::addToStorage< PdfField_T >(blocks, "pdf field CPU", real_c(std::nan("")), field::fzyx);
#endif

   BlockDataID densityFieldID = field::addToStorage< DensityField_T >(blocks, "Density", real_t(0), field::fzyx);
   BlockDataID velFieldID     = field::addToStorage< VelocityField_T >(blocks, "Velocity", real_t(0), field::fzyx);

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field");

   // set up RPD functionality
   std::function< void(void) > syncCall = [ps, rpdDomain]() {
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain);
   };

   syncCall();

   mesa_pd::kernel::ExplicitEuler explEulerIntegrator(real_t(1) / real_t(numRPDSubCycles));
   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(real_t(1) / real_t(numRPDSubCycles));
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(real_t(1) / real_t(numRPDSubCycles));

   mesa_pd::kernel::SpringDashpot collisionResponse(1);
   mesa_pd::mpi::ReduceProperty reduceProperty;

   // set up coupling functionality
   lbm_mesapd_coupling::RegularParticlesSelector sphereSelector;
   Vector3< real_t > gravitationalForce(real_t(0), real_t(0),
                                        -(densitySphere - densityFluid) * gravitationalAcceleration * sphereVolume);
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(
      viscosity, [](real_t r) { return real_t(0.0016) * r; });
   lbm::PSM_NoSlip noSlip(blocks, pdfFieldCPUGPUID);

   ///////////////
   // TIME LOOP //
   ///////////////

   // map no-slip boundaries into the LBM simulation
   std::string boundariesBlockString = " Boundaries"
                                       "{"
                                       "Border { direction T;    walldistance -1;  flag NoSlip; }"
                                       "Border { direction B;    walldistance -1;  flag NoSlip; }"
                                       "Border { direction N;    walldistance -1;  flag NoSlip; }"
                                       "Border { direction S;    walldistance -1;  flag NoSlip; }"
                                       "Border { direction W;    walldistance -1;  flag NoSlip; }"
                                       "Border { direction E;    walldistance -1;  flag NoSlip; }"
                                       "}";

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream boundariesFile("boundaries.prm");
      boundariesFile << boundariesBlockString;
      boundariesFile.close();
   }
   WALBERLA_MPI_BARRIER()

   auto boundariesCfgFile = Config();
   boundariesCfgFile.readParameterFile("boundaries.prm");
   auto boundariesConfig = boundariesCfgFile.getBlock("Boundaries");
   geometry::initBoundaryHandling< FlagField_T >(*blocks, flagFieldID, boundariesConfig);
   geometry::setNonBoundaryCellsToDomain< FlagField_T >(*blocks, flagFieldID, Fluid_Flag);
   noSlip.fillFromFlagField< FlagField_T >(blocks, flagFieldID, NoSlip_Flag, Fluid_Flag);

   // add particle and volume fraction data structures
   ParticleAndVolumeFractionSoA_T< 1 > particleAndVolumeFractionSoA(
      blocks, lbm::collision_model::omegaFromViscosity(viscosity));
   // map particles and calculate solid volume fraction initially
   PSMSweepCollection psmSweepCollection(blocks, accessor, sphereSelector, particleAndVolumeFractionSoA, Vector3(27));
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

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   gpu::communication::UniformGPUScheme< Stencil_T > com(blocks, false, false);
#else
   walberla::blockforest::communication::UniformBufferedScheme< Stencil_T > com(blocks);
#endif

   com.addPackInfo(make_shared< PackInfo_T >(pdfFieldCPUGPUID));
   auto communication = std::function< void() >([&]() { com.communicate(); });

   // create the timeloop
   SweepTimeloop timeloop(blocks->getBlockStorage(), timesteps);

   pystencils::PSMSweep PSMSweep(particleAndVolumeFractionSoA.BsFieldID, particleAndVolumeFractionSoA.BFieldID,
                                 particleAndVolumeFractionSoA.particleForcesFieldID,
                                 particleAndVolumeFractionSoA.particleVelocitiesFieldID, pdfFieldCPUGPUID, real_t(0.0),
                                 real_t(0.0), real_t(0.0), lbm::collision_model::omegaFromViscosity(viscosity));

   timeloop.addFuncBeforeTimeStep(RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger");

#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
   pystencils::PSM_MacroGetter getterSweep(BFieldID, densityFieldID, pdfFieldID, velFieldID, real_t(0.0), real_t(0.0),
                                           real_t(0.0));
#else
   pystencils::PSM_MacroGetter getterSweep(particleAndVolumeFractionSoA.BFieldID, densityFieldID, pdfFieldCPUGPUID,
                                           velFieldID, real_t(0.0), real_t(0.0), real_t(0.0));
#endif
   // vtk output
   if (vtkIOFreq != uint_t(0))
   {
      // spheres
      auto particleVtkOutput = make_shared< mesa_pd::vtk::ParticleVtkOutput >(ps);
      particleVtkOutput->setParticleSelector([sphereShape](const mesa_pd::data::ParticleStorage::iterator& pIt) {
         return pIt->getShapeID() == sphereShape;
      });
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleOwner >("owner");
      particleVtkOutput->addOutput< mesa_pd::data::SelectParticleLinearVelocity >("velocity");
      auto particleVtkWriter =
         vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", vtkIOFreq, baseFolder, "simulation_step");
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(particleVtkWriter), "VTK (sphere data)");

      // flag field (written only once in the first time step, ghost layers are also written)
      // auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, FieldGhostLayers, false,
      // baseFolder ); flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID,
      // "FlagField" ) ); timeloop.addFuncBeforeTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field", vtkIOFreq, 0, false, baseFolder);

      pdfFieldVTK->addBeforeFunction(communication);

      pdfFieldVTK->addBeforeFunction([&]() {
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
         gpu::fieldCpy< PdfField_T, gpu::GPUField< real_t > >(blocks, pdfFieldID, pdfFieldCPUGPUID);
         gpu::fieldCpy< BField_T, BFieldGPU_T >(blocks, BFieldID, particleAndVolumeFractionSoA.BFieldID);
         gpu::fieldCpy< BsField_T, BsFieldGPU_T >(blocks, BsFieldID, particleAndVolumeFractionSoA.BsFieldID);
#endif
         for (auto& block : *blocks)
            getterSweep(&block);
      });

      field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldID);
      fluidFilter.addFlag(Fluid_Flag);
      pdfFieldVTK->addCellInclusionFilter(fluidFilter);

      pdfFieldVTK->addCellDataWriter(make_shared< field::VTKWriter< VelocityField_T > >(velFieldID, "Velocity"));
      pdfFieldVTK->addCellDataWriter(make_shared< field::VTKWriter< DensityField_T > >(densityFieldID, "Density"));
#ifdef WALBERLA_BUILD_WITH_GPU_SUPPORT
      pdfFieldVTK->addCellDataWriter(make_shared< field::VTKWriter< BField_T > >(BFieldID, "Fraction mapping field B"));
      pdfFieldVTK->addCellDataWriter(
         make_shared< field::VTKWriter< BsField_T > >(BsFieldID, "Fraction mapping field Bs"));
#else
      pdfFieldVTK->addCellDataWriter(make_shared< field::VTKWriter< BField_T > >(particleAndVolumeFractionSoA.BFieldID,
                                                                                 "Fraction mapping field B"));
      pdfFieldVTK->addCellDataWriter(make_shared< field::VTKWriter< BsField_T > >(
         particleAndVolumeFractionSoA.BsFieldID, "Fraction mapping field Bs"));
#endif

      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip
   // treatment)
   timeloop.add() << BeforeFunction(communication, "LBM Communication")
                  << Sweep(noSlip.getSweep(), "Boundary Handling");

   // stream + collide LBM step
   addPSMSweepsToTimeloop(timeloop, psmSweepCollection, PSMSweep);

   // evaluation functionality
   std::string loggingFileName(baseFolder + "/LoggingSettlingSphereGPU_");
   loggingFileName += std::to_string(fluidType);
   loggingFileName += ".txt";
   if (fileIO) { WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\""); }
   SpherePropertyLogger< ParticleAccessor_T > logger(accessor, sphereUid, loggingFileName, fileIO, dx_SI, dt_SI,
                                                     diameter, -gravitationalForce[2]);

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   real_t terminationPosition = real_t(0.51) * diameter; // right before sphere touches the bottom wall

   const bool useOpenMP = false;

   // time loop
   for (uint_t i = 0; i < timesteps; ++i)
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      timeloop.singleStep(timeloopTiming);

      timeloopTiming["RPD"].start();

      // -> full hydrodynamic force/torque info is available on local particle
      reduceProperty.operator()< mesa_pd::HydrodynamicForceTorqueNotification >(*ps);

      if (averageForceTorqueOverTwoTimeSteps)
      {
         if (i == 0)
         {
            lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel
               initializeHydrodynamicForceTorqueForAveragingKernel;
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                                initializeHydrodynamicForceTorqueForAveragingKernel, *accessor);
         }
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque,
                             *accessor);
      }

      for (auto subCycle = uint_t(0); subCycle < numRPDSubCycles; ++subCycle)
      {
         if (useVelocityVerlet)
         {
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
            syncCall();
         }

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction,
                             *accessor);
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor);

         // lubrication correction
         ps->forEachParticlePairHalf(
            useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
            [&lubricationCorrectionKernel, rpdDomain](const size_t idx1, const size_t idx2, auto& ac) {
               // TODO change this to storing copy, not reference
               mesa_pd::collision_detection::AnalyticContactDetection acd;
               acd.getContactThreshold() = lubricationCorrectionKernel.getNormalCutOffDistance();
               mesa_pd::kernel::DoubleCast double_cast;
               mesa_pd::mpi::ContactFilter contact_filter;
               if (double_cast(idx1, idx2, ac, acd, ac))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                  {
                     double_cast(idx1, idx2, ac, lubricationCorrectionKernel, ac, acd.getContactNormal(),
                                 acd.getPenetrationDepth());
                  }
               }
            },
            *accessor);

         // one could add linked cells here

         // collision response
         ps->forEachParticlePairHalf(
            useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
            [collisionResponse, rpdDomain](const size_t idx1, const size_t idx2, auto& ac) {
               mesa_pd::collision_detection::AnalyticContactDetection acd;
               mesa_pd::kernel::DoubleCast double_cast;
               mesa_pd::mpi::ContactFilter contact_filter;
               if (double_cast(idx1, idx2, ac, acd, ac))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                  {
                     collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(),
                                       acd.getPenetrationDepth());
                  }
               }
            },
            *accessor);

         reduceProperty.operator()< mesa_pd::ForceTorqueNotification >(*ps);

         if (useVelocityVerlet)
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         else
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, explEulerIntegrator, *accessor);

         syncCall();
      }

      timeloopTiming["RPD"].end();

      // evaluation
      timeloopTiming["Logging"].start();
      logger(i);
      timeloopTiming["Logging"].end();

      // reset after logging here
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor);

      // check for termination
      if (logger.getPosition() < terminationPosition)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Sphere reached terminal position " << logger.getPosition() << " after " << i
                                                                       << " timesteps!");
         break;
      }
   }

   timeloopTiming.logResultOnRoot();

   // check the result
   if (!funcTest && !shortrun)
   {
      real_t relErr = std::fabs(expectedSettlingVelocity - logger.getMaxVelocity()) / expectedSettlingVelocity;
      WALBERLA_LOG_INFO_ON_ROOT("Expected maximum settling velocity: " << expectedSettlingVelocity);
      WALBERLA_LOG_INFO_ON_ROOT("Simulated maximum settling velocity: " << logger.getMaxVelocity());
      WALBERLA_LOG_INFO_ON_ROOT("Relative error: " << relErr);

      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS(relErr, real_t(0.1));
   }

   return EXIT_SUCCESS;
}

} // namespace settling_sphere

int main(int argc, char** argv) { settling_sphere::main(argc, argv); }
