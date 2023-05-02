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
//! \file LightRisingParticleInFluidAMR.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include "blockforest/communication/UniformToNonUniformPackInfoAdapter.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/all.h"
#include "core/math/Random.h"
#include "core/math/all.h"
#include "core/mpi/Broadcast.h"
#include "core/timing/RemainingTimeLogger.h"

#include "domain_decomposition/BlockSweepWrapper.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/adaptors/AdaptorCreators.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/all.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/field/QCriterionFieldWriter.h"
#include "lbm/field/VelocityFieldWriter.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/all.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/all.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"

#include <functional>

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/amr/InfoCollection.h"
#include "lbm_mesapd_coupling/amr/level_determination/ParticlePresenceLevelDetermination.h"
#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/ExtrapolationDirectionFinder.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/SubCyclingManager.h"
#include "lbm_mesapd_coupling/utility/virtualmass/AddVirtualForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/virtualmass/InitializeVirtualMassKernel.h"
#include "lbm_mesapd_coupling/utility/virtualmass/ParticleAccessorWithShapeVirtualMassWrapper.h"

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ClearGhostOwnerSync.h"
#include "mesa_pd/mpi/ClearNextNeighborSync.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/SyncGhostOwners.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

/**
 * \brief Showcase for a stabilized fluid - very light particle simulation using virtual mass
 * A light sphere, situated at the bottom of a cuboid domain, rises.
 * Depending on the density of the sphere (--sphereDensity) and galileo number (--galileoNumber), the sphere
 * then undergoes various movement patterns along its way to the top.
 * Employs adaptive mesh refinement and checkpointing.
 * Can output VTK files based on the q criterion to visualize vortices behind the rising sphere.
 * Most default parameters should be sufficient and do not need explicit values.
 * As used in the master thesis of Lukas Werner, 2020, and follow-up paper at http://doi.org/10.1002/fld.5034
 */

namespace light_rising_particle_amr {
using namespace walberla;
using walberla::uint_t;

#ifdef BUILD_WITH_MRT
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::D3Q19MRT>;
#else
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::TRT, false >;
#endif

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using ParticleField_T = lbm_mesapd_coupling::ParticleField_T;
using VelocityField_T = GhostLayerField<Vector3<real_t>, 1>;
using QCriterionField_T = GhostLayerField<real_t, 1>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

using ScalarField_T = GhostLayerField<real_t, 1>;

const uint_t FieldGhostLayers = 4;

using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;

using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t> ;
using MO_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
using BoundaryHandling_T = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T >;


//// flags

const FlagUID Fluid_Flag("fluid");
const FlagUID NoSlip_Flag("no slip");
const FlagUID MO_Flag("moving obstacle");
const FlagUID FormerMO_Flag("former moving obstacle");


//// block structure

static void refinementSelection(SetupBlockForest &forest, uint_t levels, AABB refinementBox) {
   real_t dx = real_t(1); // dx on finest level
   const uint_t finestLevel = levels - uint_t(1);
   for (auto block = forest.begin(); block != forest.end(); ++block) {
      uint_t blockLevel = block->getLevel();
      uint_t levelScalingFactor = (uint_t(1) << (finestLevel - blockLevel));
      real_t dxOnLevel = dx * real_c(levelScalingFactor);
      AABB blockAABB = block->getAABB();
      // extend block AABB by ghostlayers
      AABB extendedBlockAABB = blockAABB.getExtended(dxOnLevel * real_c(FieldGhostLayers));
      if (extendedBlockAABB.intersects(refinementBox))
         if (blockLevel < finestLevel) // only if the block is not on the finest level
            block->setMarker(true);
   }
}

static void workloadAndMemoryAssignment(SetupBlockForest &forest) {
   for (auto block = forest.begin(); block != forest.end(); ++block) {
      block->setWorkload(numeric_cast<workload_t>(uint_t(1) << block->getLevel()));
      block->setMemory(numeric_cast<memory_t>(1));
   }
}
static shared_ptr<StructuredBlockForest>
createBlockStructure(const AABB &domainAABB, Vector3<uint_t> blockSizeInCells, uint_t numberOfLevels, real_t diameter,
                     Vector3<real_t> spherePosition, Vector3<bool> periodicity,
                     bool readFromCheckPointFile, const std::string & forestCheckPointFileName,
                     bool keepGlobalBlockInformation = false) {
   if (readFromCheckPointFile) {
         WALBERLA_LOG_INFO_ON_ROOT("Reading block forest from file!");
         return blockforest::createUniformBlockGrid(forestCheckPointFileName,
                 blockSizeInCells[0], blockSizeInCells[1], blockSizeInCells[2],
                 false);
   } else {
      SetupBlockForest sforest;
      Vector3<uint_t> numberOfFineBlocksPerDirection(uint_c(domainAABB.size(0)) / blockSizeInCells[0],
                                                     uint_c(domainAABB.size(1)) / blockSizeInCells[1],
                                                     uint_c(domainAABB.size(2)) / blockSizeInCells[2]);
      for (uint_t i = 0; i < 3; ++i) {
         WALBERLA_CHECK_EQUAL(numberOfFineBlocksPerDirection[i] * blockSizeInCells[i], uint_c(domainAABB.size(i)),
                              "Domain can not be decomposed in direction " << i << " into fine blocks of size "
                                                                           << blockSizeInCells[i]);
      }

      uint_t levelScalingFactor = (uint_t(1) << (numberOfLevels - uint_t(1)));
      Vector3<uint_t> numberOfCoarseBlocksPerDirection(numberOfFineBlocksPerDirection / levelScalingFactor);
      for (uint_t i = 0; i < 3; ++i) {
         WALBERLA_CHECK_EQUAL(numberOfCoarseBlocksPerDirection[i] * levelScalingFactor,
                              numberOfFineBlocksPerDirection[i],
                              "Domain can not be refined in direction " << i
                                                                        << " according to the specified number of levels!");
      }
      AABB refinementBox( std::floor(spherePosition[0] - real_t(0.5) * diameter),
                          std::floor(spherePosition[1] - real_t(0.5) * diameter),
                          std::floor(spherePosition[2] - real_t(0.5) * diameter),
                          std::ceil( spherePosition[0] + real_t(0.5) * diameter),
                          std::ceil( spherePosition[1] + real_t(0.5) * diameter),
                          std::ceil( spherePosition[2] + real_t(0.5) * diameter) );
      WALBERLA_LOG_INFO_ON_ROOT(" - refinement box: " << refinementBox);

      sforest.addRefinementSelectionFunction(
              std::bind(refinementSelection, std::placeholders::_1, numberOfLevels, refinementBox));
      sforest.addWorkloadMemorySUIDAssignmentFunction(workloadAndMemoryAssignment);

      sforest.init(domainAABB, numberOfCoarseBlocksPerDirection[0], numberOfCoarseBlocksPerDirection[1],
                   numberOfCoarseBlocksPerDirection[2], periodicity[0], periodicity[1], periodicity[2]);

      // calculate process distribution
      const memory_t memoryLimit = math::Limits<memory_t>::inf();
      sforest.balanceLoad(blockforest::StaticLevelwiseCurveBalance(true),
                          uint_c(MPIManager::instance()->numProcesses()), real_t(0), memoryLimit, true);

      WALBERLA_LOG_INFO_ON_ROOT(sforest);

      MPIManager::instance()->useWorldComm();

      auto blockForest = make_shared<BlockForest>(uint_c(MPIManager::instance()->rank()), sforest,
                                                  keepGlobalBlockInformation);

      // create StructuredBlockForest (encapsulates a newly created BlockForest)
      shared_ptr<StructuredBlockForest> sbf =
              make_shared<StructuredBlockForest>(
                      blockForest,
                      blockSizeInCells[0], blockSizeInCells[1], blockSizeInCells[2]);
      sbf->createCellBoundingBoxes();

      return sbf;
   }
}

//// boundary handling

class MyBoundaryHandling : public blockforest::AlwaysInitializeBlockDataHandling< BoundaryHandling_T > {
public:
   MyBoundaryHandling(const weak_ptr< StructuredBlockStorage > & storage, const BlockDataID& flagFieldID,
         const BlockDataID& pdfFieldID, const BlockDataID& particleFieldID, const shared_ptr<ParticleAccessor_T>& ac) :
         storage_(storage), flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID), particleFieldID_(particleFieldID),
         ac_(ac) {

   }

   BoundaryHandling_T * initialize( IBlock * const block ) override;

private:
   weak_ptr< StructuredBlockStorage > storage_;

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;
};

BoundaryHandling_T * MyBoundaryHandling::initialize( IBlock * const block )
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   auto* flagField = block->getData<FlagField_T>(flagFieldID_);
   auto* pdfField = block->getData<PdfField_T>(pdfFieldID_);
   auto* particleField = block->getData<lbm_mesapd_coupling::ParticleField_T>(particleFieldID_);

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   auto storagePtr = storage_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( storagePtr );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storagePtr, *block ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}

void clearBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID ) {
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BoundaryHandling_T * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->clear( FieldGhostLayers );
   }
}

void recreateBoundaryHandling( BlockForest & forest, const BlockDataID & boundaryHandlingID ) {
   for( auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt )
   {
      BoundaryHandling_T * boundaryHandling = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID);
      boundaryHandling->fillWithDomain( FieldGhostLayers );
   }
}

void clearParticleField( BlockForest & forest, const BlockDataID & particleFieldID, const ParticleAccessor_T & accessor ) {
   for (auto blockIt = forest.begin(); blockIt != forest.end(); ++blockIt) {
      ParticleField_T *particleField = blockIt->getData<ParticleField_T>(particleFieldID);
      particleField->setWithGhostLayer(accessor.getInvalidUid());
   }
}

class SpherePropertyLogger {
public:
   SpherePropertyLogger(lbm_mesapd_coupling::SubCyclingManager* mesapdTimeStep, const shared_ptr<ParticleAccessor_T>& ac, walberla::id_t sphereUid,
                        const std::string& fileName, bool fileIO, bool writeHeader, real_t t_g, real_t u_g, real_t diameter) :
         mesapdTimeStep_(mesapdTimeStep), ac_(ac), sphereUid_(sphereUid), fileName_(fileName), fileIO_(fileIO), t_g_(t_g),
         u_g_(u_g), diameter_(diameter), position_(real_t(0)) {
      if (fileIO_ && writeHeader) {
         WALBERLA_ROOT_SECTION() {
            std::ofstream file;
            file.open(fileName_.c_str());
            file << "#\t posZ\t uZ\t u\t rotX\t rotY\t rotZ\t hydFZ\t hydT\t tN\t posXN\t posYN\t posZN\t uXN\t uYN\t uZN\t uN\t wX\t wY\t wZ\n";
            file.close();
         }
      }
   }

   void operator()() {
      const uint_t timestep (mesapdTimeStep_->getCurrentTimeStep());
      Vector3<real_t> transVel;
      Vector3<real_t> hydForce;
      Vector3<real_t> hydTorque;
      Vector3<real_t> angVel;
      Vector3<real_t> eulerRotation;

      size_t idx = ac_->uidToIdx(sphereUid_);
      if (idx != ac_->getInvalidIdx()) {
         if (!mesa_pd::data::particle_flags::isSet(ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST)) {
            transVel = ac_->getLinearVelocity(idx);
            hydForce = ac_->getHydrodynamicForce(idx);
            hydTorque = ac_->getHydrodynamicTorque(idx);
            angVel = ac_->getAngularVelocity(idx);
            eulerRotation = ac_->getRotation(idx).getMatrix().getEulerAnglesXYZ();
         }
      }

      WALBERLA_MPI_SECTION() {
         mpi::allReduceInplace(transVel, mpi::SUM);
         mpi::allReduceInplace(angVel, mpi::SUM);
         mpi::allReduceInplace(hydForce, mpi::SUM);
         mpi::allReduceInplace(hydTorque, mpi::SUM);
         mpi::allReduceInplace(eulerRotation, mpi::SUM);
      }

      syncPosition();

      if (fileIO_) {
         writeToFile(timestep, position_, eulerRotation, transVel, angVel, hydForce, hydTorque);
      }
   }

   Vector3<real_t> getPosition() const {
      return position_;
   }

   void syncPosition() {
      Vector3<real_t> pos(real_t(0));

      size_t idx = ac_->uidToIdx(sphereUid_);
      if (idx != ac_->getInvalidIdx()) {
         if (!mesa_pd::data::particle_flags::isSet(ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST)) {
            pos = ac_->getPosition(idx);
         }
      }

      WALBERLA_MPI_SECTION() {
         mpi::allReduceInplace(pos, mpi::SUM);
      }

      position_ = pos;
   }

private:
   void writeToFile(const uint_t timestep, const Vector3<real_t>& position, const Vector3<real_t>& eulerRotation,
         const Vector3<real_t>& velocity, const Vector3<real_t>& angularVelocity,
         const Vector3<real_t>& hydForce, const Vector3<real_t>& hydTorque) {
      WALBERLA_ROOT_SECTION() {
         std::ofstream file;
         file.open(fileName_.c_str(), std::ofstream::app);

         file << timestep << "\t"
              << position[2] << "\t"
              << velocity[2] << "\t"
              << velocity.length() << "\t"
              << eulerRotation[0] << "\t"
              << eulerRotation[1] << "\t"
              << eulerRotation[2] << "\t"
              << hydForce[2] << "\t"
              << hydTorque.length() << "\t"

              << real_t(timestep) / t_g_ << "\t"
              << position[0] / diameter_ << "\t"
              << position[1] / diameter_ << "\t"
              << position[2] / diameter_ << "\t"
              << velocity[0] / u_g_ << "\t"
              << velocity[1] / u_g_ << "\t"
              << velocity[2] / u_g_ << "\t"
              << velocity.length() / u_g_ << "\t"
              << angularVelocity[0] * diameter_ / u_g_ << "\t"
              << angularVelocity[1] * diameter_ / u_g_ << "\t"
              << angularVelocity[2] * diameter_ / u_g_ << "\n";

         file.close();
      }
   }

   lbm_mesapd_coupling::SubCyclingManager* mesapdTimeStep_;
   shared_ptr<ParticleAccessor_T> ac_;
   const walberla::id_t sphereUid_;
   std::string fileName_;
   bool fileIO_;

   real_t t_g_;
   real_t u_g_;
   real_t diameter_;

   Vector3<real_t> position_;
};

void createCheckpointFiles(const std::string& checkPointFileName, std::vector<std::string>& oldCheckpointFiles,
        uint_t lbmTimeStep, uint_t mesaTimeStep, const shared_ptr<BlockForest>& forest, const BlockDataID pdfFieldID,
        const BlockDataID particleStorageID) {
   WALBERLA_LOG_INFO_ON_ROOT("Writing checkpointing files in lbm time step " << lbmTimeStep << " (mesapd: " << mesaTimeStep << ").");

   const std::string timeStepTag = std::to_string(lbmTimeStep) + "_" + std::to_string(mesaTimeStep);

   const std::string lbmFileName = checkPointFileName + "_" + timeStepTag + "_lbm.txt";
   const std::string mesaFileName = checkPointFileName + "_" + timeStepTag + "_mesa.txt";
   const std::string forestFileName = checkPointFileName + "_" + timeStepTag + "_forest.txt";

   forest->saveBlockData(lbmFileName, pdfFieldID);
   forest->saveBlockData(mesaFileName, particleStorageID);
   forest->saveToFile(forestFileName);

   WALBERLA_ROOT_SECTION() {
      if (!oldCheckpointFiles.empty()) {
         for (const std::string& file: oldCheckpointFiles) {
            filesystem::path path(file);
            if (filesystem::exists(path)) {
               filesystem::remove(path);
            }
         }
         oldCheckpointFiles.clear();
         WALBERLA_LOG_INFO("Deleted old checkpoint files.");
      }

      oldCheckpointFiles.push_back(lbmFileName);
      oldCheckpointFiles.push_back(mesaFileName);
      oldCheckpointFiles.push_back(forestFileName);
   }
}

// main

int main(int argc, char** argv) {
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   math::seedRandomGenerator( static_cast<unsigned int>(1337 * mpi::MPIManager::instance()->worldRank()) );
   
   // general params
   bool fileIO = true;
   uint_t vtkIOFreq = 0;
   uint_t fullVtkIOFreq = 0;
   uint_t qCriterionVtkIOFreq = 0;
   std::string vtkIOSelection = "sliced"; // "none", "full", "thicksliced" (half the sphere still in), "trackingsliced" (sphere is followed in y direction)
   std::string baseFolder = "vtk_out_UnsettlingSphereDR";
   auto vtkLowerVelocityLimit = real_t(0.003);
   auto vtkLowerQCriterionLimit = real_t(1e-10);
   auto propertyLoggingFreq = uint_t(1);

   // numerical params
   auto blockSize = uint_t(32);
   auto XBlocks = uint_t(32);
   auto YBlocks = uint_t(32);
   auto ZBlocks = uint_t(64);
   auto refinementCheckFrequency = uint_t(0);
   auto initialZ = real_t(-1);
   auto diameter = real_t(16);
   bool offsetInitialPosition = true;
   bool averageForceTorqueOverTwoTimeSteps = true;
   auto numMESAPDSubCycles = uint_t(10);
   auto u_ref = real_t(0.01);
   auto densitySphere = real_t(0.1);
   auto timesteps = uint_t(2000);
   auto galileoNumber = real_t(300);
   auto Lambda_Bulk = real_t(100); // if = 1, normal TRT, else not. Up to 100.
   bool conserveMomentum = false;

   // virtual mass
   bool useVirtualMass = true;
   auto C_v = real_t(0.5);
   auto C_v_omega = real_t(0.5);

   // dynamic refinement
   auto numberOfLevels = uint_t(3);

   bool useVorticityCriterion = false;
   auto lowerFluidRefinementLimit = real_t(0);
   auto upperFluidRefinementLimit = std::numeric_limits<real_t>::infinity();

   bool useCurlCriterion = true;
   auto upperCurlLimit = real_t(0.001);
   auto curlCoarsenFactor = real_t(0.2);
   auto curlLengthScaleWeight = real_t(2);

   // checkpointing
   bool writeCheckPointFile = false;
   bool readFromCheckPointFile = false;
   uint_t checkPointingFreq = 0;
   uint_t checkPointLBMTimeStep = 0;
   uint_t checkPointMESAPDTimeStep = 0;

   for (int i = 1; i < argc; i++) {
      if (std::strcmp(argv[i], "--shortRun") == 0) { fileIO = false; densitySphere = real_t(0.02); galileoNumber = 100; timesteps = 2; continue; }

      if (std::strcmp(argv[i], "--noLogging")                == 0) { fileIO = false; continue; }
      if (std::strcmp(argv[i], "--vtkIOFreq")                == 0) { vtkIOFreq = uint_c( std::stod( argv[++i] ) ); continue; }
      if (std::strcmp(argv[i], "--fullVtkIOFreq")            == 0) { fullVtkIOFreq = uint_c( std::stod( argv[++i] ) ); continue; }
      if (std::strcmp(argv[i], "--qCriterionVtkIOFreq")      == 0) { qCriterionVtkIOFreq = uint_c( std::stod( argv[++i] ) ); continue; }
      if (std::strcmp(argv[i], "--vtkIOSelection")           == 0) { vtkIOSelection = argv[++i]; continue; }
      if (std::strcmp(argv[i], "--numMESAPDSubCycles")       == 0) { numMESAPDSubCycles = uint_c( std::stod( argv[++i] ) ); continue; }
      if (std::strcmp(argv[i], "--noForceAveraging")         == 0) { averageForceTorqueOverTwoTimeSteps = false; continue; }
      if (std::strcmp(argv[i], "--baseFolder")               == 0) { baseFolder = argv[++i]; continue; }
      if (std::strcmp(argv[i], "--vtkLowerVelocityLimit")    == 0) { vtkLowerVelocityLimit = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--vtkLowerQCriterionLimit")  == 0) { vtkLowerQCriterionLimit = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--propertyLoggingFreq")      == 0) { propertyLoggingFreq = uint_c( std::stod( argv[++i] ) ); continue; }

      if (std::strcmp(argv[i], "--dontOffsetInitialPos")     == 0) { offsetInitialPosition = false; continue; }
      if (std::strcmp(argv[i], "--sphereDensity")            == 0) { densitySphere = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--u_ref")                    == 0) { u_ref = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--timesteps")                == 0) { timesteps = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--blockSize")                == 0) { blockSize = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--XBlocks")                  == 0) { XBlocks = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--YBlocks")                  == 0) { YBlocks = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--ZBlocks")                  == 0) { ZBlocks = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--refinementCheckFrequency") == 0) { refinementCheckFrequency = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--initialZ")                 == 0) { initialZ = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--diameter")                 == 0) { diameter = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--galileoNumber")            == 0) { galileoNumber = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--lambdaBulk")               == 0) { Lambda_Bulk = real_c(std::stod(argv[++i])); continue; }

      if (std::strcmp(argv[i], "--noVirtualMass")            == 0) { useVirtualMass = false; continue; }
      if (std::strcmp(argv[i], "--c_v")                      == 0) { C_v = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--c_v_omega")                == 0) { C_v_omega = real_c(std::stod(argv[++i])); continue; }

      // dynamic refinement
      if (std::strcmp(argv[i], "--numberOfLevels")           == 0) { numberOfLevels = uint_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--useVorticityCriterion")    == 0) { useVorticityCriterion = true; useCurlCriterion = false; continue; }
      if (std::strcmp(argv[i], "--lowerLimit")               == 0) { lowerFluidRefinementLimit = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--upperLimit")               == 0) { upperFluidRefinementLimit = real_c(std::stod(argv[++i])); continue; }

      if (std::strcmp(argv[i], "--upperCurlLimit")           == 0) { upperCurlLimit = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--curlCoarsenFactor")        == 0) { curlCoarsenFactor = real_c(std::stod(argv[++i])); continue; }
      if (std::strcmp(argv[i], "--curlLengthScaleWeight")    == 0) { curlLengthScaleWeight = real_c(std::stod(argv[++i])); continue; }

      if (std::strcmp(argv[i], "--writeCheckPointFile")      == 0) { writeCheckPointFile = true; continue; }
      if (std::strcmp(argv[i], "--readFromCheckPointFile")   == 0) { readFromCheckPointFile = true; continue; }
      if (std::strcmp(argv[i], "--checkPointingFreq")        == 0) { checkPointingFreq = uint_c(std::stod( argv[++i]) ); continue; }
      if (std::strcmp(argv[i], "--checkPointLBMTimeStep")    == 0) { checkPointLBMTimeStep = uint_c(std::stod( argv[++i]) ); continue; }
      if (std::strcmp(argv[i], "--checkPointMESAPDTimeStep") == 0) { checkPointMESAPDTimeStep = uint_c(std::stod( argv[++i]) ); continue; }

      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   WALBERLA_CHECK(!realIsEqual(u_ref, real_t(0)), "u_ref has to be non-zero.")
   WALBERLA_CHECK(useCurlCriterion || useVorticityCriterion,
         "Use curl or vorticity (legacy) criterion for refinement.");
   WALBERLA_CHECK(!(useCurlCriterion && useVorticityCriterion),
         "Using curl and vorticity criterion together makes no sense.");

   // create base dir if it doesn't already exist
   filesystem::path bpath(baseFolder);
   if (!filesystem::exists(bpath)) {
      filesystem::create_directory(bpath);
   }
   filesystem::path qpath(baseFolder + "/q_criterion_mesh");
   if (filesystem::exists(qpath)) {
      filesystem::remove_all(qpath);
   }

   //// numerical parameters

   const auto densityFluid = real_t(1);

   // calculate relaxation time from u_ref
   real_t viscosity = u_ref * diameter / galileoNumber;
   real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   real_t relaxationTime = real_t(1) / omega;

   const auto d3 = diameter*diameter*diameter;
   const auto G2 = galileoNumber*galileoNumber;
   const auto v2 = viscosity*viscosity;
   const auto gravitationalAcceleration = G2*v2 / (std::fabs(densitySphere-real_t(1))*d3);

   const auto dx = real_t(1);
   const auto overlap = real_t(4.5) * dx;

   const Vector3<uint_t> domainSize( XBlocks * blockSize, YBlocks * blockSize, ZBlocks * blockSize );
   if (realIsIdentical(initialZ, real_t(-1))) {
      initialZ = real_c(blockSize) + real_t(0.5) * diameter;
   }
   Vector3<real_t> initialPosition(real_t(0.5) * real_c(domainSize[0]),
                                   real_t(0.5) * real_c(domainSize[1]),
                                   initialZ);
   if (offsetInitialPosition) {
      // offset sphere's initial position half a cell to introduce a slight disturbance.
      initialPosition += Vector3<real_t>(real_t(0.5), real_t(0.5), real_t(0));
   }

   if (useVorticityCriterion && floatIsEqual(lowerFluidRefinementLimit, real_t(0)) && std::isinf(upperFluidRefinementLimit)) {
      // use computed criterion instead of user input
      lowerFluidRefinementLimit = real_t(0.05) * u_ref;
      upperFluidRefinementLimit = real_t(0.1) * u_ref;
   }

   //const auto dx = real_t(1);

   const auto u_g = std::sqrt(std::fabs(densitySphere - real_t(1)) * gravitationalAcceleration * diameter);
   const auto t_g = std::sqrt(diameter / (std::fabs(densitySphere - real_t(1)) * gravitationalAcceleration));

   const bool rising = densitySphere < real_t(1);

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - number of processes = " << mpi::MPIManager::instance()->numProcesses());
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size / D = " << (Vector3<real_t>(domainSize) / diameter));
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter << ", density = " << densitySphere << ", initial position: " << initialPosition);
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid << ", relaxation time (tau) = " << relaxationTime << ", kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - lambda bulk = " << Lambda_Bulk );
   WALBERLA_LOG_INFO_ON_ROOT(" - galileo number = " << galileoNumber );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational velocity (u_g) = " << u_g);
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
#ifdef BUILD_WITH_MRT
   WALBERLA_LOG_INFO_ON_ROOT(" - lattice model: MRT");
#else
   WALBERLA_LOG_INFO_ON_ROOT(" - lattice model: TRT");
#endif
   WALBERLA_LOG_INFO_ON_ROOT(" - integrator = Velocity Verlet" );
   WALBERLA_LOG_INFO_ON_ROOT(" - use virtual mass = " << useVirtualMass );
   //WALBERLA_LOG_INFO_ON_ROOT(" - virtual mass acceleration estimation = " << vmAccelerationOption );
   WALBERLA_LOG_INFO_ON_ROOT(" - C_v = " << C_v );
   WALBERLA_LOG_INFO_ON_ROOT(" - C_v_omega = " << C_v_omega );
   WALBERLA_LOG_INFO_ON_ROOT(" - timesteps = " << timesteps);
   WALBERLA_LOG_INFO_ON_ROOT(" - averageForceTorqueOverTwoTimeSteps = " << averageForceTorqueOverTwoTimeSteps);
   WALBERLA_LOG_INFO_ON_ROOT(" - conserve momentum = " << conserveMomentum);
   if (useCurlCriterion) {
      WALBERLA_LOG_INFO_ON_ROOT(" - using curl criterion with upper limit = " << upperCurlLimit << ", coarsen factor = " << curlCoarsenFactor << " and length scale weight = " << curlLengthScaleWeight);
   }
   if (useVorticityCriterion) {
      WALBERLA_LOG_INFO_ON_ROOT(" - using vorticity criterion with lower limit = " << lowerFluidRefinementLimit << " and upper limit = " << upperFluidRefinementLimit );
   }
   if (vtkIOFreq > 0) {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder << "\" with frequency " << vtkIOFreq);
   }

   if(refinementCheckFrequency == 0 && numberOfLevels != 1){
      // determine check frequency automatically based on maximum admissible velocity and block sizes
      auto uMax = real_t(0.02); // this should never ever be reached
      const uint_t finestLevel = numberOfLevels - uint_t(1);
      const uint_t levelScalingFactor = ( uint_t(1) << finestLevel );
      const uint_t lbmTimeStepsPerTimeLoopIteration = levelScalingFactor;
      refinementCheckFrequency = uint_c((overlap + (real_c(blockSize) - real_t(2) * real_t(FieldGhostLayers)) * dx) / uMax) / lbmTimeStepsPerTimeLoopIteration;
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - refinement / load balancing check frequency (coarse time steps): " << refinementCheckFrequency);

   //// checkpointing setup

   WALBERLA_LOG_INFO_ON_ROOT(" - read from checkpoint file = " << readFromCheckPointFile);
   WALBERLA_LOG_INFO_ON_ROOT(" - write checkpoint file = " << writeCheckPointFile);

   std::string checkPointFileName = "checkpointing";
   if (readFromCheckPointFile || writeCheckPointFile) {
      WALBERLA_LOG_INFO_ON_ROOT(" - checkPointingFileName = " << checkPointFileName);
   }

   std::string readCheckPointFileName = checkPointFileName + "_" + std::to_string(checkPointLBMTimeStep) + "_" +
                                        std::to_string(checkPointMESAPDTimeStep);

   if (readFromCheckPointFile && writeCheckPointFile) {
      // decide which option to choose
      if (filesystem::exists(readCheckPointFileName+"_lbm.txt")) {
         WALBERLA_LOG_INFO_ON_ROOT("Checkpoint file already exists! Will skip writing check point file and start from this check point!");
         writeCheckPointFile = false;
      } else {
         WALBERLA_LOG_INFO_ON_ROOT("Checkpoint file does not exists yet! Will skip reading check point file and just regularly start the simulation from the beginning!");
         readFromCheckPointFile = false;
      }
   }

   // write parameters to file
   if (fileIO) {
      WALBERLA_ROOT_SECTION() {
         std::string paramsFileName(baseFolder + "/Params.txt");
         std::ofstream file;
         file.open(paramsFileName.c_str());
         file
               << "D\t "
               << "Ga\t "
               << "LBM\t "
               << "pi_p\t "
               << "tau\t "
               << "visc\t "
               << "u_g\t "
               << "t_g\t "
               << "acc_g\t "
               << "vm\t "
               << "C_v\t "
               << "C_vw\t "
               << "fa\t "
               << "curlCrit\t "
               << "curlMax\t "
               << "curlCoaF\t "
               << "curlLenS\t "
               << "refChk\t "
               << "numLvl\t "
               << "sX\t "
               << "sY\t "
               << "sZ\t "
               << "nProc\t "
               << "\n";
         file
               << diameter << "\t "
               << galileoNumber << "\t "
               #ifdef BUILD_WITH_MRT
               << "MRT\t "
               #else
               << "TRT\t "
               #endif
               << densitySphere << "\t "
               << relaxationTime << "\t "
               << viscosity << "\t "
               << u_g << "\t "
               << t_g << "\t "
               << gravitationalAcceleration << "\t "
               << (useVirtualMass ? "lt" : "no") << "\t "
               << C_v << "\t "
               << C_v_omega << "\t "
               << averageForceTorqueOverTwoTimeSteps << "\t "
               << useCurlCriterion << "\t "
               << upperCurlLimit << "\t "
               << curlCoarsenFactor << "\t "
               << curlLengthScaleWeight << "\t "
               << refinementCheckFrequency << "\t "
               << numberOfLevels << "\t "
               << domainSize[0] << "\t "
               << domainSize[1] << "\t "
               << domainSize[2] << "\t "
               << mpi::MPIManager::instance()->numProcesses() << "\t "
               << "\n";
         file.close();
      }
   }

   //

   uint_t currentLBMTimeStep = 0;
   uint_t currentMESAPDTimeStep = 0;

   if (readFromCheckPointFile) {
      currentLBMTimeStep = checkPointLBMTimeStep;
      currentMESAPDTimeStep = checkPointMESAPDTimeStep;
   }

   //// block structure setup

   const uint_t finestLevel = numberOfLevels - uint_t(1);

   Vector3<uint_t> blockSizeInCells(blockSize);

   AABB simulationDomain( real_t(0), real_t(0), real_t(0), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   Vector3<bool> periodicity(true, true, true);
   auto blocks = createBlockStructure(simulationDomain, blockSizeInCells, numberOfLevels, diameter, initialPosition,
           periodicity, readFromCheckPointFile, readCheckPointFileName+"_forest.txt");
   auto forest = blocks->getBlockForestPointer();

   if (vtkIOFreq > 0) {
      vtk::writeDomainDecomposition(blocks, "initial_domain_decomposition", baseFolder);
   }

   //// dynamic refinement I

   forest->recalculateBlockLevelsInRefresh(true);
   forest->alwaysRebalanceInRefresh(false);
   forest->reevaluateMinTargetLevelsAfterForcedRefinement(false);
   forest->allowRefreshChangingDepth(false);

   //// mesapd coupling

   auto mesapdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   // init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();

   auto accessor = walberla::make_shared<ParticleAccessor_T>(ps, ss);

   BlockDataID particleStorageID;
   if (readFromCheckPointFile) {
      WALBERLA_LOG_INFO_ON_ROOT("Initializing particles from checkpointing file.");
      particleStorageID = forest->loadBlockData( readCheckPointFileName + "_mesa.txt", mesa_pd::domain::createBlockForestDataHandling(ps), "Particle Storage" );
   } else {
      particleStorageID = forest->addBlockData(mesa_pd::domain::createBlockForestDataHandling(ps), "Particle Storage");
   }


   // bounding planes can be added here

   auto sphereShape = ss->create<mesa_pd::data::Sphere>(diameter * real_t(0.5));
   ss->shapes[sphereShape]->updateMassAndInertia(densitySphere);

   walberla::id_t sphereUid = 0;
   if (readFromCheckPointFile) {
      for (auto pIt = ps->begin(); pIt != ps->end(); ++pIt) {
         // find sphere in loaded data structure and store uid for later reference
         if (pIt->getShapeID() == sphereShape) {
            sphereUid = pIt->getUid();
            pIt->getNeighborStateRef().clear();
            pIt->getGhostOwnersRef().clear();
         }
      }
   } else {
      if (mesapdDomain->isContainedInProcessSubdomain(uint_c(mpi::MPIManager::instance()->rank()), initialPosition)) {
         mesa_pd::data::Particle &&p = *ps->create();
         p.setPosition(initialPosition);
         p.setInteractionRadius(diameter * real_t(0.5));
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         sphereUid = p.getUid();
      }
   }
   mpi::allReduceInplace(sphereUid, mpi::SUM);

   if(sphereUid == 0) {
      WALBERLA_ABORT("No sphere present - aborting!");
   }

   auto ghostOwnerSyncIterations = uint_t(real_t(0.5) * diameter/real_t(blockSize) + real_t(1));
   WALBERLA_LOG_INFO_ON_ROOT(" - Number of req'd ghost owner sync iterations: " << ghostOwnerSyncIterations << ", thus using " << (ghostOwnerSyncIterations == 1 ? "SNN" : "SGO") << " for synchronization.");

   bool useSyncNextNeighbors = (ghostOwnerSyncIterations == 1);

   std::function<void(void)> syncCall = [ps, mesapdDomain, useSyncNextNeighbors, overlap](){
      if (useSyncNextNeighbors) {
         mesa_pd::mpi::SyncNextNeighbors SNN; // alt. SyncGhostOwners
         SNN(*ps, *mesapdDomain, overlap);
      } else {
         // need to use syncghostowners since the particle is bigger than a block and thus, syncnextneighbors doesn't suffice.
         mesa_pd::mpi::SyncGhostOwners SGO;
         SGO(*ps, *mesapdDomain, overlap);
      }
   };

   std::function<void(void)> initializationSyncCall = [&syncCall, ghostOwnerSyncIterations](){
      for (uint_t j = 0; j < ghostOwnerSyncIterations; ++j) {
         syncCall();
      }
   };

   //initial sync
   initializationSyncCall();

   //// add data to blocks

   // create the lattice model
#ifdef BUILD_WITH_MRT
   const real_t Lambda_Magic = real_t(3) / real_t(16);
   const real_t lambda_e = lbm::collision_model::TRT::lambda_e( omega );
   const real_t lambda_d = lbm::collision_model::TRT::lambda_d( omega, Lambda_Magic );
   const real_t omegaBulk = real_t(1) / (Lambda_Bulk * ( real_t(1) / omega - real_t(1)/ real_t(2) ) + real_t(1)/ real_t(2) );

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT( omegaBulk, lambda_e, lambda_d, lambda_e, lambda_e, lambda_d, finestLevel ));
#else
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, lbm::collision_model::TRT::threeSixteenth, finestLevel ) );
#endif

   // add PDF field
   BlockDataID pdfFieldID;
   if (readFromCheckPointFile) {
      // add PDF field
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing PDF Field from checkpointing file!" );
      shared_ptr< lbm::internal::PdfFieldHandling< LatticeModel_T > > dataHandling =
              make_shared< lbm::internal::PdfFieldHandling< LatticeModel_T > >(blocks, latticeModel, false,
                      Vector3<real_t>(real_t(0)), real_t(1),
                      FieldGhostLayers, field::fzyx );

      pdfFieldID = blocks->loadBlockData( readCheckPointFileName+"_lbm.txt", dataHandling, "pdf field" );

   } else {
      // add PDF field
      pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >(blocks, "pdf field (fzyx)", latticeModel,
              Vector3<real_t>(real_t(0)), real_t(1),
              FieldGhostLayers, field::fzyx);
   }

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   using FluidFilter_T = typename field::FlagFieldEvaluationFilter<FlagField_T>;
   FluidFilter_T fluidFlagFieldEvaluationFilter(flagFieldID, Fluid_Flag);

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<ParticleField_T>(blocks, "particle field",
         accessor->getInvalidUid(),
         field::fzyx,
         FieldGhostLayers);

   // add velocity field and utility
   BlockDataID velocityFieldID = field::addToStorage<VelocityField_T>( blocks, "velocity field", Vector3<real_t>(real_t(0)), field::fzyx, uint_t(2) );

   typedef lbm::VelocityFieldWriter< PdfField_T, VelocityField_T > VelocityFieldWriter_T;
   BlockSweepWrapper< VelocityFieldWriter_T > velocityFieldWriter( blocks, VelocityFieldWriter_T( pdfFieldID, velocityFieldID ) );

   shared_ptr<blockforest::communication::NonUniformBufferedScheme<stencil::D3Q27> > velocityCommunicationScheme = make_shared<blockforest::communication::NonUniformBufferedScheme<stencil::D3Q27> >( blocks );
   velocityCommunicationScheme->addPackInfo( make_shared< field::refinement::PackInfo<VelocityField_T, stencil::D3Q27> >( velocityFieldID ) );

   // add q criterion field (only needed for mesh output)
   BlockDataID qCriterionFieldID = field::addToStorage<QCriterionField_T>(blocks, "q criterion field", real_t(0), field::fzyx, uint_t(1));

   typedef lbm::QCriterionFieldWriter<VelocityField_T, QCriterionField_T, FluidFilter_T> QCriterionFieldWriter_T;
   BlockSweepWrapper<QCriterionFieldWriter_T> qCriterionFieldWriter(blocks, QCriterionFieldWriter_T(blocks, velocityFieldID,
           qCriterionFieldID, fluidFlagFieldEvaluationFilter));

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addBlockData( make_shared< MyBoundaryHandling >( blocks, flagFieldID, pdfFieldID, particleFieldID, accessor ),
                                                          "boundary handling" );

   // map particles into simulation
   lbm_mesapd_coupling::RegularParticlesSelector sphereSelector;
   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks,
                                                                                                    boundaryHandlingID,
                                                                                                    particleFieldID);

   std::function<void(void)> mappingCall = [ps, sphereSelector, accessor, &movingParticleMappingKernel](){
      ps->forEachParticle(false, sphereSelector, *accessor, movingParticleMappingKernel, *accessor, MO_Flag);
   };

   mappingCall();

   syncCall();

   //// dynamic refinement II

   blockforest::CombinedMinTargetLevelDeterminationFunctions minTargetLevelDeterminationFunctions;

   auto couplingInfoCollection = walberla::make_shared<lbm_mesapd_coupling::amr::InfoCollection>();
   lbm_mesapd_coupling::amr::ParticlePresenceLevelDetermination particlePresenceRefinement(couplingInfoCollection, finestLevel);

   minTargetLevelDeterminationFunctions.add(particlePresenceRefinement);

   if (useVorticityCriterion) {
      // add refinement criterion based on vorticity magnitude
      // "legacy"
      lbm::refinement::VorticityBasedLevelDetermination<field::FlagFieldEvaluationFilter<FlagField_T>> vorticityRefinement(
            velocityFieldID, fluidFlagFieldEvaluationFilter, upperFluidRefinementLimit, lowerFluidRefinementLimit, finestLevel);

      minTargetLevelDeterminationFunctions.add(vorticityRefinement);
   }

   if (useCurlCriterion) {
      // add refinement criterion based on curl (=vorticity) magnitude
      // -> computes the same as vorticity, but has other params and uses correct length scale weighting
      const real_t lowerCurlLimit = curlCoarsenFactor*upperCurlLimit;

      lbm::refinement::CurlBasedLevelDetermination<field::FlagFieldEvaluationFilter<FlagField_T>> curlRefinement(
            velocityFieldID, *blocks, fluidFlagFieldEvaluationFilter, finestLevel, upperCurlLimit, lowerCurlLimit,
            curlLengthScaleWeight);

      minTargetLevelDeterminationFunctions.add(curlRefinement);
   }

   forest->setRefreshMinTargetLevelDeterminationFunction(minTargetLevelDeterminationFunctions);

   bool curveHilbert = false;
   bool curveAllGather = true;

   forest->setRefreshPhantomBlockMigrationPreparationFunction( blockforest::DynamicCurveBalance< blockforest::NoPhantomData >( curveHilbert, curveAllGather ) );


   //// time loop

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
   commFunction = scheme;

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );
   timeloop.setCurrentTimeStep(currentLBMTimeStep);

   shared_ptr<WcTimingPool> timeloopTiming = make_shared<WcTimingPool>();

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, sweep, pdfFieldID, boundaryHandlingID );

   //// MESAPD Time Step

   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;
   Vector3<real_t> gravitationalForce(real_t(0), real_t(0),
                                      -(densitySphere - densityFluid) * gravitationalAcceleration * sphereVolume);

   lbm_mesapd_coupling::SubCyclingManager mesapdTimeStep(numMESAPDSubCycles, timeloopTiming);
   mesapdTimeStep.setCurrentTimeStep(currentMESAPDTimeStep);
   bool useOpenMP = false;

   // before subcycling
   mesa_pd::mpi::ReduceProperty reduceProperty;
   lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueAveraging;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;

   mesapdTimeStep.addFuncBeforeSubCycles([ps, accessor, useOpenMP, &reduceProperty, averageForceTorqueOverTwoTimeSteps,
                                          &initializeHydrodynamicForceTorqueAveraging, &averageHydrodynamicForceTorque,
                                          &mesapdTimeStep](){
      reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);
      if (averageForceTorqueOverTwoTimeSteps) {
         if (mesapdTimeStep.getCurrentTimeStep() == 0) {
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueAveraging, *accessor);
         } else {
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor);
         }
      }
   }, "RPD Before Subcycling");

   // during subcycling
   using VirtualMass_ParticleAccessor_T = lbm_mesapd_coupling::ParticleAccessorWithShapeVirtualMassWrapper<ParticleAccessor_T>;
   auto virtualMassAccessor = walberla::make_shared<VirtualMass_ParticleAccessor_T>(ps, ss);

   mesa_pd::kernel::VelocityVerletPreForceUpdate vvIntegratorPreForce(real_t(1) / real_t(numMESAPDSubCycles));
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(real_t(1) / real_t(numMESAPDSubCycles));
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;

   mesapdTimeStep.addFuncDuringSubCycles([ps, accessor, virtualMassAccessor, useOpenMP, syncCall,
                                          useVirtualMass,
                                          vvIntegratorPreForce, addHydrodynamicInteraction, addGravitationalForce,
                                          vvIntegratorPostForce, mesapdDomain,
                                          C_v, C_v_omega, densityFluid](){
      if (useVirtualMass) {
         lbm_mesapd_coupling::InitializeVirtualMassKernel virtualMass;
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, virtualMass, *accessor, C_v, C_v_omega, densityFluid);

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *virtualMassAccessor, vvIntegratorPreForce, *virtualMassAccessor);
      } else {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
      }
      syncCall();

      // Hydrodynamic force and torque got reduced onto particle local process before subcycling,
      // thus only add them there.
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor);
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor);

      if (useVirtualMass) {
         lbm_mesapd_coupling::InitializeVirtualMassKernel virtualMass;
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, virtualMass, *accessor, C_v, C_v_omega, densityFluid);

         lbm_mesapd_coupling::AddVirtualForceTorqueKernel addVirtualForceAndTorque(ps);

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addVirtualForceAndTorque, *accessor);
      }

      if (useVirtualMass) {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *virtualMassAccessor, vvIntegratorPostForce, *virtualMassAccessor);
      } else {
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
      }

      syncCall();
   }, "RPD Subcycling");

   // after subcycling
   // evaluation
   std::string loggingFileName(baseFolder + "/LoggingUnsettlingSphereDynamicRefinement");
   if (useVirtualMass) {
      loggingFileName += "_vm";
   }
   loggingFileName += ".txt";
   if (fileIO) {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   SpherePropertyLogger logger(&mesapdTimeStep, accessor, sphereUid, loggingFileName, fileIO,
           !readFromCheckPointFile, t_g, u_g, diameter);

   auto periodicLogger = [&mesapdTimeStep, &logger, propertyLoggingFreq](){
      if (propertyLoggingFreq != 0 && mesapdTimeStep.getCurrentTimeStep() % propertyLoggingFreq == 0) {
         logger();
      }
   };

   mesapdTimeStep.addFuncAfterSubCycles(periodicLogger, "Sphere property logger" );

   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   mesapdTimeStep.addFuncAfterSubCycles([ps, accessor, useOpenMP, resetHydrodynamicForceTorque](){
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor);
   }, "RPD After Subcycling");

   // mesapdTimeStep seems to get copied or something by the functorwrapper
   auto wrapper = lbm::refinement::FunctorWrapper([&mesapdTimeStep](){
      mesapdTimeStep();
   });
   refinementTimestep->addPostStreamVoidFunction(wrapper, "rpd Time Step", finestLevel);


   //// Particle Mapping

   // sweep for updating the particle mapping into the LBM simulation
   auto movingParticleMapping = lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks,
                                                                                                               pdfFieldID,
                                                                                                               boundaryHandlingID,
                                                                                                               particleFieldID,
                                                                                                               accessor,
                                                                                                               MO_Flag,
                                                                                                               FormerMO_Flag,
                                                                                                               sphereSelector,
                                                                                                               conserveMomentum);
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::SweepAsFunctorWrapper(movingParticleMapping, blocks),
                                                 "Particle Mapping", finestLevel);

   //// PDF Restore

   blockforest::communication::NonUniformBufferedScheme<Stencil_T> fullPDFCommunicationScheme(blocks);
   auto uniformPackInfo = make_shared<field::communication::PackInfo<PdfField_T>>(pdfFieldID);
   auto packInfoAdapter = make_shared<blockforest::communication::UniformToNonUniformPackInfoAdapter>(uniformPackInfo);
   // using the pack info adapter is required, else the ghost layers are not correctly synced, resulting in invalid reconstructed cells (with NaN as value).
   fullPDFCommunicationScheme.addPackInfo(packInfoAdapter);
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::FunctorWrapper(fullPDFCommunicationScheme),
                                                 "PDF Communication", finestLevel);

   // add sweep for restoring PDFs in cells previously occupied by particles
   auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(
           blocks);
   auto extrapolationSphereNormalReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T, lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder, true>(
           blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder);
   //auto equilibriumRecon = lbm_mesapd_coupling::makeEquilibriumReconstructor<BoundaryHandling_T>(
   //        blocks, boundaryHandlingID);
   auto reconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T, BoundaryHandling_T>(
           blocks,
           pdfFieldID,
           boundaryHandlingID,
           particleFieldID,
           accessor,
           FormerMO_Flag,
           Fluid_Flag,
           extrapolationSphereNormalReconstructor,
           conserveMomentum);
   refinementTimestep->addPostStreamVoidFunction(lbm::refinement::SweepAsFunctorWrapper(*reconstructionManager, blocks),
                                                 "PDF Restore", finestLevel);

   // add LBM sweep with refinement
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( refinementTimestep ), "LBM refinement time step" );


   //// VTK output

   auto vtkInitialExecutionCount = uint_t(0);
   std::string vtkFileTag = "";
   if (readFromCheckPointFile) {
      vtkInitialExecutionCount = checkPointLBMTimeStep;
      vtkFileTag = "_" + std::to_string(vtkInitialExecutionCount);
   }

   // sphere(s)
   auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleOwner>("owner");
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleLinearVelocity>("velocity");
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleInteractionRadius>("interactionRadius");
   auto particleVtkWriter = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles"+vtkFileTag, vtkIOFreq,
         baseFolder, "simulation_step", false, true, true, true,
         vtkInitialExecutionCount);

   // flag field left out

   // pdf field
   auto pdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field"+vtkFileTag, vtkIOFreq, 0,
         false, baseFolder, "simulation_step", false, true, true, true,
         vtkInitialExecutionCount);

   vtk::ChainedFilter combinedFilter;
   if (vtkIOSelection == "sliced") {
      AABB sliceAABB( real_t(0), real_c(domainSize[1])*real_t(0.5), real_t(0),
                      real_c(domainSize[0]), real_c(domainSize[1])*real_t(0.5)+real_t(1), real_c(domainSize[2]) );
      vtk::AABBCellFilter aabbSliceFilter( sliceAABB );
      combinedFilter.addFilter( aabbSliceFilter );
   } else if (vtkIOSelection == "thicksliced") {
      AABB sliceAABB( real_t(0), real_c(domainSize[1])*real_t(0.5)-real_t(diameter)*real_t(0.5)-real_t(1), real_t(0),
                      real_c(domainSize[0]), real_c(domainSize[1])*real_t(0.5), real_c(domainSize[2]) );
      vtk::AABBCellFilter aabbSliceFilter( sliceAABB );
      combinedFilter.addFilter( aabbSliceFilter );
   } else if (vtkIOSelection == "trackingsliced") {
      auto filter = [logger, domainSize](CellSet& filteredCells, const IBlock& block, const StructuredBlockStorage& storage, const uint_t ghostLayers = uint_t(0)){
         auto yPosition = logger.getPosition()[1];
         real_t minY = std::floor(yPosition);
         real_t maxY = std::floor(yPosition)+real_t(1);
         AABB sliceAABB( real_t(0), minY, real_t(0),
                         real_c(domainSize[0]), maxY, real_c(domainSize[2]) );
         vtk::AABBCellFilter aabbSliceFilter( sliceAABB );
         aabbSliceFilter(filteredCells, block, storage, ghostLayers);
      };
      combinedFilter.addFilter( filter );
      logger.syncPosition(); // synchronize the position once to correctly supply the vtk output before the first timestep
   } else if (vtkIOSelection == "velocity") {
      field::VelocityCellFilter<PdfField_T, field::FlagFieldEvaluationFilter<FlagField_T>> filter(pdfFieldID,
            fluidFlagFieldEvaluationFilter, vtkLowerVelocityLimit);
      combinedFilter.addFilter(filter);
   } else if (vtkIOSelection != "none" && vtkIOSelection != "full") {
      WALBERLA_ABORT("Invalid vtk output selection.");
   }

   field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldID);
   fluidFilter.addFlag(Fluid_Flag);
   combinedFilter.addFilter(fluidFilter);

   pdfFieldVTK->addCellInclusionFilter(combinedFilter);

   pdfFieldVTK->addBeforeFunction(velocityFieldWriter);
   pdfFieldVTK->addBeforeFunction(*velocityCommunicationScheme);

   pdfFieldVTK->addCellDataWriter(make_shared<lbm::VelocityVTKWriter<LatticeModel_T, float>>(pdfFieldID, "VelocityFromPDF"));
   pdfFieldVTK->addCellDataWriter(make_shared<lbm::DensityVTKWriter<LatticeModel_T, float>>(pdfFieldID, "DensityFromPDF"));
   pdfFieldVTK->addCellDataWriter(make_shared<lbm::QCriterionVTKWriter<VelocityField_T, FluidFilter_T, float>>(blocks,
         fluidFlagFieldEvaluationFilter, velocityFieldID, "QCriterionFromPDF"));
   pdfFieldVTK->addCellDataWriter(make_shared<lbm::CurlMagnitudeVTKWriter<VelocityField_T, FluidFilter_T, float>>(blocks,
         fluidFlagFieldEvaluationFilter, velocityFieldID, "CurlSensor", real_t(2)));

   // full size vtk
   auto fullPdfFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field_full"+vtkFileTag, fullVtkIOFreq,
         0, false, baseFolder, "simulation_step", false, true, true, true,
         vtkInitialExecutionCount);
   fullPdfFieldVTK->addBeforeFunction(velocityFieldWriter);
   fullPdfFieldVTK->addBeforeFunction(*velocityCommunicationScheme);
   fullPdfFieldVTK->addCellInclusionFilter(fluidFilter);
   fullPdfFieldVTK->addCellDataWriter(make_shared<lbm::VelocityVTKWriter<LatticeModel_T, float>>(pdfFieldID, "VelocityFromPDF"));
   fullPdfFieldVTK->addCellDataWriter(make_shared<lbm::QCriterionVTKWriter<VelocityField_T , FluidFilter_T, float>>(blocks,
         fluidFlagFieldEvaluationFilter, velocityFieldID, "QCriterionFromPDF"));
   fullPdfFieldVTK->addCellDataWriter(make_shared<lbm::CurlMagnitudeVTKWriter<VelocityField_T, FluidFilter_T, float>>(blocks,
         fluidFlagFieldEvaluationFilter, velocityFieldID, "CurlSensor", curlLengthScaleWeight));

   // q value vtk
   auto qCriterionVTK = vtk::createVTKOutput_BlockData(blocks, "q_value"+vtkFileTag, qCriterionVtkIOFreq, 0,
         false, baseFolder, "simulation_step", false, true, true, true,
         vtkInitialExecutionCount);

   field::QCriterionCellFilter<VelocityField_T, field::FlagFieldEvaluationFilter<FlagField_T>> qFilter(velocityFieldID,
         fluidFlagFieldEvaluationFilter, vtkLowerQCriterionLimit);
   qCriterionVTK->addCellInclusionFilter(qFilter);
   qCriterionVTK->addBeforeFunction(velocityFieldWriter);
   qCriterionVTK->addBeforeFunction(*velocityCommunicationScheme);
   qCriterionVTK->addCellDataWriter(make_shared<lbm::QCriterionVTKWriter<VelocityField_T, FluidFilter_T, float>>(blocks,
         fluidFlagFieldEvaluationFilter, velocityFieldID, "QCriterionFromPDF"));
   qCriterionVTK->addCellDataWriter(make_shared<lbm::VorticityComponentVTKWriter<VelocityField_T, FluidFilter_T>>(blocks,
         fluidFlagFieldEvaluationFilter, velocityFieldID, uint_t(2), "VorticityZFromPDF", u_g/diameter));

   // domain decomposition
   auto domainDecompVTK = vtk::createVTKOutput_DomainDecomposition(blocks, "domain_decomposition"+vtkFileTag,
         vtkIOFreq, baseFolder, "simulation_step", false, true, true, true,
         vtkInitialExecutionCount);

   timeloop.addFuncBeforeTimeStep(vtk::writeFiles(particleVtkWriter), "VTK (sparse data)");
   if (vtkIOSelection != "none") {
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(pdfFieldVTK), "VTK (fluid field data)");
   }
   timeloop.addFuncAfterTimeStep(vtk::writeFiles(domainDecompVTK), "VTK (domain decomposition)");

   if (fullVtkIOFreq != uint_t(0)) {
      timeloop.addFuncBeforeTimeStep(vtk::writeFiles(fullPdfFieldVTK), "VTK (full fluid field data)");
   }

   timeloop.addFuncBeforeTimeStep(vtk::writeFiles(qCriterionVTK), "VTK (q value data)");

   timeloop.addFuncAfterTimeStep( RemainingTimeLogger(timeloop.getNrOfTimeSteps()), "Remaining Time Logger" );

   real_t terminationPosition;
   // when the sphere is only a diameter away from the initial position again
   if (rising) {
      terminationPosition = real_t(initialPosition[2]) - diameter;
   } else {
      terminationPosition = real_t(initialPosition[2]) + diameter;
   }
   if (terminationPosition < real_t(0)) {
      terminationPosition = real_t(domainSize[2]) + terminationPosition;
   } else if (terminationPosition > real_t(domainSize[2])) {
      terminationPosition = terminationPosition - real_t(domainSize[2]);
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - termination position = " << terminationPosition );

   std::vector<std::string> oldCheckpointFiles;

   //// execute simulation
   bool hasCrossedBoundary = false;

   for (uint_t i = 0; i < timesteps; ++i) {
      if (refinementCheckFrequency != 0 && i % refinementCheckFrequency == 0) {
         WALBERLA_LOG_INFO_ON_ROOT("Refinement: Refreshing...");

         lbm_mesapd_coupling::amr::updateAndSyncInfoCollection<BoundaryHandling_T, ParticleAccessor_T>(*forest, boundaryHandlingID,
               *accessor, numMESAPDSubCycles, *couplingInfoCollection);

         // for the fluid property based check
         if (useCurlCriterion || useVorticityCriterion) {
            velocityFieldWriter();
            (*velocityCommunicationScheme)();
         }

         // clear mesapd
         if (useSyncNextNeighbors) {
            mesa_pd::mpi::ClearNextNeighborSync CNNS;
            CNNS(*accessor);
         } else {
            mesa_pd::mpi::ClearGhostOwnerSync CGOS;
            CGOS(*accessor);
         }

         forest->refresh();
         mesapdDomain->refresh();
         initializationSyncCall();

         clearBoundaryHandling(*forest, boundaryHandlingID);
         clearParticleField(*forest, particleFieldID, *accessor);
         recreateBoundaryHandling(*forest, boundaryHandlingID);

         mappingCall();

         std::vector<uint_t> blocksPerLevel(numberOfLevels);
         for (uint_t l = 0; l < numberOfLevels; l++) {
            blocksPerLevel[l] = blocks->getNumberOfBlocks(l);
         }
         mpi::reduceInplace(blocksPerLevel, mpi::SUM);
         WALBERLA_ROOT_SECTION() {
            WALBERLA_LOG_INFO("Refinement: Number of blocks per level: ");
            for (uint_t l = 0; l < numberOfLevels; l++) {
               WALBERLA_LOG_INFO(" - Level " << l << ": " << blocksPerLevel[l] << " blocks");
            }
         }
      }

      // perform a single simulation step
      timeloop.singleStep( *timeloopTiming );

      // check for termination
      // check if the sphere crossed a boundary
      auto zPosition = logger.getPosition()[2];
      if (!hasCrossedBoundary && ((rising && zPosition < terminationPosition) ||
                                  (!rising && zPosition > terminationPosition))) {
         hasCrossedBoundary = true;
      }
      if (hasCrossedBoundary && ((rising && zPosition > terminationPosition)
                                 || (!rising && zPosition < terminationPosition))) {
         // end of simulation reached after crossing a boundary and then passing the termination position
         WALBERLA_LOG_INFO_ON_ROOT(
                 "Sphere reached terminal position " << logger.getPosition() << " after " << i << " timesteps!");
         break;
      }
      if (checkPointingFreq != 0 && i % checkPointingFreq == 0) {
         createCheckpointFiles(checkPointFileName, oldCheckpointFiles, timeloop.getCurrentTimeStep(),
                 mesapdTimeStep.getCurrentTimeStep(), forest, pdfFieldID, particleStorageID);
      }
   }

   size_t particleIdx = accessor->uidToIdx(sphereUid);
   if (particleIdx != accessor->getInvalidIdx()) {
      if (!mesa_pd::data::particle_flags::isSet(accessor->getFlags(particleIdx), mesa_pd::data::particle_flags::GHOST)) {
         Vector3<real_t> terminalVelocity = accessor->getLinearVelocity(particleIdx);
         real_t reynoldsNumber = terminalVelocity[2] * diameter / viscosity;
         WALBERLA_LOG_INFO("Terminal velocity: " << terminalVelocity);
         WALBERLA_LOG_INFO("Terminal Reynolds number: " << reynoldsNumber);
      }
   }

   timeloopTiming->logResultOnRoot();

   return EXIT_SUCCESS;
}

}

int main(int argc, char **argv) {
   light_rising_particle_amr::main(argc, argv);
}