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
//! \file FluidParticleWorkLoadEvaluation.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"
#include "core/mpi/Broadcast.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"
#include "lbm/BlockForestEvaluation.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <vector>
#include <iostream>

namespace fluid_particle_workload_evaluation
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

using LatticeModel_T = lbm::D3Q19< lbm::collision_model::TRT, false >;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID MO_Flag  ( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template <typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using MO_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, MO_T >;

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T>& ac,
                       bool useEntireFieldTraversal) :
         flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac ), useEntireFieldTraversal_(useEntireFieldTraversal) {}

   Type * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( storage );

      auto * flagField     = block->getData< FlagField_T >( flagFieldID_ );
      auto *  pdfField     = block->getData< PdfField_T > ( pdfFieldID_ );
      auto * particleField = block->getData< lbm_mesapd_coupling::ParticleField_T > ( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

      typename Type::Mode mode = (useEntireFieldTraversal_) ? Type::Mode::ENTIRE_FIELD_TRAVERSAL : Type::Mode::OPTIMIZED_SPARSE_TRAVERSAL ;

      Type * handling = new Type( "moving obstacle boundary handling", flagField, fluid,
                                  MO_T( "MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ),
                                  mode);

      handling->fillWithDomain( FieldGhostLayers );

      return handling;
   }

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;

   bool useEntireFieldTraversal_;
};

template< typename BoundaryHandling_T >
void evaluateFluidQuantities(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID boundaryHandlingID,
                             uint_t & numCells, uint_t & numFluidCells, uint_t & numNBCells )
{

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      auto * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID );
      auto xyzSize = boundaryHandling->getFlagField()->xyzSize();
      numCells += xyzSize.numCells();

      for(auto z = cell_idx_t(xyzSize.zMin()); z <= cell_idx_t(xyzSize.zMax()); ++z ){
         for(auto y = cell_idx_t(xyzSize.yMin()); y <= cell_idx_t(xyzSize.yMax()); ++y ){
            for(auto x = cell_idx_t(xyzSize.xMin()); x <= cell_idx_t(xyzSize.xMax()); ++x ) {
               if (boundaryHandling->isDomain(x, y, z)) {
                  ++numFluidCells;
               }
               if (boundaryHandling->isNearBoundary(x, y, z)) {
                  ++numNBCells;
               }
            }
         }
      }
   }
}

void evaluateRPDQuantities( const shared_ptr< mesa_pd::data::ParticleStorage > & ps,
                            uint_t & numLocalParticles, uint_t & numGhostParticles)
{

   for (auto pIt = ps->begin(); pIt != ps->end(); ++pIt)
   {
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(pIt->getFlags(), GHOST))
      {
         ++numGhostParticles;
      } else
      {
         //note: global particles are included here
         // use if (!isSet(pIt->getFlags(), GLOBAL)) if should be excluded
         ++numLocalParticles;
      }
   }
}

void evaluateTimers(WcTimingPool & timingPool,
                    const std::vector<std::vector<std::string> > & timerKeys,
                    std::vector<double> & timings )
{

   for (auto & timingsIt : timings)
   {
      timingsIt = 0.0;
   }

   timingPool.unifyRegisteredTimersAcrossProcesses();

   double scalingFactor = 1000.0; // milliseconds

   for (auto i = uint_t(0); i < timerKeys.size(); ++i )
   {
      auto keys = timerKeys[i];
      for (const auto &timerName : keys)
      {
         if(timingPool.timerExists(timerName))
         {
            timings[i] += timingPool[timerName].total() * scalingFactor;
         }
      }

   }
}


//*******************************************************************************************************************
/*! Application to evaluate the workload (time measurements) for a fluid-particle simulation
 *
 * This application is used in the paper
 *  Rettinger, Ruede - "Dynamic Load Balancing Techniques for Particulate Flow Simulations", 2019, Computation
 * in Section 3 to develop and calibrate the workload estimator.
 * The setup features settling particle inside a horizontally periodic box.
 * A comprehensive description is given in Sec. 3.3 of the paper.
 * It uses 4 x 4 x 5 blocks for domain partitioning.
 * For each block ( = each process), the block local quantities are evaluated as well as the timing infos of
 * the fluid-particle interaction algorithm. Those infos are then written to files that can be used later on
 * for function fitting to obtain a workload estimator.
 *
 * NOTE: Since this estimator relies on timing measurements, this evaluation procedure should be carried out everytime
 * a different implementation, hardware or algorithm is used.
 *
 */
//*******************************************************************************************************************
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto solidVolumeFraction = real_t(0.2);

   // LBM / numerical parameters
   auto blockSize  = uint_t(32);
   auto uSettling = real_t(0.1); // characteristic settling velocity
   auto diameter = real_t(10);

   auto Ga = real_t(30); //Galileo number
   auto numRPDSubCycles = uint_t(10);

   auto vtkIOFreq = uint_t(0);
   auto timestepsNonDim = real_t(2.5);
   auto numSamples = uint_t(2000);
   std::string baseFolder = "workload_files"; // folder for vtk and file output

   bool noFileOutput = false;
   bool useEntireFieldTraversal = true;
   bool useFusedStreamCollide = false;

   auto XBlocks = uint_t(4);
   auto YBlocks = uint_t(4);
   auto ZBlocks = uint_t(5);

   real_t topWallOffsetFactor = real_t(1.05);

   bool vtkDuringInit = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--vtkIOFreq"               ) == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noFileOutput"            ) == 0 ) { noFileOutput = true; continue; }
      if( std::strcmp( argv[i], "--vtkDuringInit"           ) == 0 ) { vtkDuringInit = true; continue; }
      if( std::strcmp( argv[i], "--basefolder"              ) == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--solidVolumeFraction"     ) == 0 ) { solidVolumeFraction = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--diameter"                ) == 0 ) { diameter = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--blockSize"               ) == 0 ) { blockSize = uint_c(std::atof( argv[++i]) ); continue; }
      if( std::strcmp( argv[i], "--XBlocks"                 ) == 0 ) { XBlocks = uint_c(std::atof( argv[++i]) ); continue; }
      if( std::strcmp( argv[i], "--YBlocks"                 ) == 0 ) { YBlocks = uint_c(std::atof( argv[++i]) ); continue; }
      if( std::strcmp( argv[i], "--ZBlocks"                 ) == 0 ) { ZBlocks = uint_c(std::atof( argv[++i]) ); continue; }
      if( std::strcmp( argv[i], "--uSettling"               ) == 0 ) { uSettling = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--topWallOffsetFactor"     ) == 0 ) { topWallOffsetFactor = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--Ga"                      ) == 0 ) { Ga = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--timestepsNonDim"         ) == 0 ) { timestepsNonDim = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--numRPDSubCycles"         ) == 0 ) { numRPDSubCycles = uint_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--useEntireFieldTraversal" ) == 0 ) { useEntireFieldTraversal = true; continue; }
      if( std::strcmp( argv[i], "--numSamples"              ) == 0 ) { numSamples = uint_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--useFusedStreamCollide"   ) == 0 ) { useFusedStreamCollide = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   WALBERLA_CHECK(diameter > real_t(1));
   WALBERLA_CHECK(uSettling > real_t(0));
   WALBERLA_CHECK(Ga > real_t(0));
   WALBERLA_CHECK(solidVolumeFraction > real_t(0));
   WALBERLA_CHECK(solidVolumeFraction < real_t(0.65));

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////


   if( MPIManager::instance()->numProcesses() != int(XBlocks * YBlocks * ZBlocks) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT("WARNING! You have specified less or more processes than number of blocks -> the time measurements are no longer blockwise!")
   }

   if( diameter > real_c(blockSize) )
   {
      WALBERLA_LOG_WARNING_ON_ROOT("The bodies might be too large to work with the currently used synchronization!");
   }

   WALBERLA_LOG_INFO_ON_ROOT("Using setup with sedimenting particles -> creating two planes and applying gravitational force")

   const uint_t XCells = blockSize * XBlocks;
   const uint_t YCells = blockSize * YBlocks;
   const uint_t ZCells = blockSize * ZBlocks;

   const real_t topWallOffset =  topWallOffsetFactor * real_t(blockSize); // move the top wall downwards to take away a certain portion of the overall domain

   // determine number of spheres to generate, if necessary scale diameter a bit to reach desired solid volume fraction
   real_t domainHeight = real_c(ZCells) - topWallOffset;
   real_t fluidVolume =  real_c( XCells * YCells ) * domainHeight;
   real_t solidVolume = solidVolumeFraction * fluidVolume;
   uint_t numberOfParticles = uint_c(std::ceil(solidVolume / ( math::pi / real_t(6) * diameter * diameter * diameter )));
   diameter = std::cbrt( solidVolume / ( real_c(numberOfParticles) * math::pi / real_t(6) ) );

   auto densityRatio = real_t(2.5);

   real_t viscosity = uSettling * diameter / Ga;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t gravitationalAcceleration = uSettling * uSettling / ( (densityRatio-real_t(1)) * diameter );

   real_t tref = diameter / uSettling;
   real_t Tref = domainHeight / uSettling;

   uint_t timesteps = uint_c(timestepsNonDim * Tref);

   const real_t dx = real_c(1);
   WALBERLA_LOG_INFO_ON_ROOT("viscosity = " << viscosity);
   WALBERLA_LOG_INFO_ON_ROOT("tau = " << real_t(1)/omega);
   WALBERLA_LOG_INFO_ON_ROOT("diameter = " << diameter);
   WALBERLA_LOG_INFO_ON_ROOT("solid volume fraction = " << solidVolumeFraction);
   WALBERLA_LOG_INFO_ON_ROOT("domain size (in cells) = " << XCells << " x " << YCells << " x " << ZCells);
   WALBERLA_LOG_INFO_ON_ROOT("number of bodies = " << numberOfParticles);
   WALBERLA_LOG_INFO_ON_ROOT("gravitational acceleration = " << gravitationalAcceleration);
   WALBERLA_LOG_INFO_ON_ROOT("Ga = " << Ga);
   WALBERLA_LOG_INFO_ON_ROOT("uSettling = " << uSettling);
   WALBERLA_LOG_INFO_ON_ROOT("tref = " << tref);
   WALBERLA_LOG_INFO_ON_ROOT("Tref = " << Tref);
   WALBERLA_LOG_INFO_ON_ROOT("timesteps = " << timesteps);
   WALBERLA_LOG_INFO_ON_ROOT("number of workload samples = " << numSamples);

   // create folder to store logging files
   WALBERLA_ROOT_SECTION()
   {
      walberla::filesystem::path path1( baseFolder );
      if( !walberla::filesystem::exists( path1 ) )
         walberla::filesystem::create_directory( path1 );
   }


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<bool> periodicity( true );
   periodicity[2] = false;

   // create domain
   shared_ptr< StructuredBlockForest > blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, blockSize, blockSize, blockSize, dx,
                                                    0, false, false, //one block per process!
                                                    periodicity[0], periodicity[1], periodicity[2], //periodicity
                                                    false );

   /////////
   // RPD //
   /////////

   const real_t restitutionCoeff = real_t(0.97);
   const real_t frictionCoeffStatic = real_t(0.8);
   const real_t frictionCoeffDynamic = real_t(0.15);
   const real_t collisionTime = real_t(4) * diameter; // from my paper
   const real_t poissonsRatio = real_t(0.22);
   const real_t kappa = real_t(2) * ( real_t(1) - poissonsRatio ) / ( real_t(2) - poissonsRatio ) ;

   auto rpdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   //init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   real_t timeStepSizeRPD = real_t(1)/real_t(numRPDSubCycles);
   mesa_pd::kernel::VelocityVerletPreForceUpdate  vvIntegratorPreForce(timeStepSizeRPD);
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(timeStepSizeRPD);

   // types: 0 = wall, 1: sphere
   mesa_pd::kernel::LinearSpringDashpot collisionResponse(2);
   collisionResponse.setFrictionCoefficientDynamic(0,1,frictionCoeffDynamic);
   collisionResponse.setFrictionCoefficientDynamic(1,1,frictionCoeffDynamic);
   collisionResponse.setFrictionCoefficientStatic(0,1,frictionCoeffStatic);
   collisionResponse.setFrictionCoefficientStatic(1,1,frictionCoeffStatic);

   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;
   const real_t particleMass = densityRatio * sphereVolume;
   const real_t effMass_sphereWall = particleMass;
   const real_t effMass_sphereSphere = particleMass * particleMass / ( real_t(2) * particleMass );
   collisionResponse.setStiffnessAndDamping(0,1,restitutionCoeff,collisionTime,kappa,effMass_sphereWall);
   collisionResponse.setStiffnessAndDamping(1,1,restitutionCoeff,collisionTime,kappa,effMass_sphereSphere);

   mesa_pd::mpi::ReduceProperty reduceProperty;
   mesa_pd::mpi::ReduceContactHistory reduceAndSwapContactHistory;

   //////////////
   // COUPLING //
   //////////////

   // connect to pe
   const real_t overlap = real_c( 1.5 ) * dx;

   std::function<void(void)> syncCall = [&ps,&rpdDomain,overlap](){
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   auto generationDomain = AABB( real_t(0), real_t(0), real_t(0), real_c(XCells), real_c(YCells), real_c(ZCells) - topWallOffset);

   // create plane at top and bottom

   mesa_pd::data::Particle&& p0 = *ps->create(true);
   p0.setPosition(generationDomain.minCorner());
   p0.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,1) ));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle&& p1 = *ps->create(true);
   p1.setPosition(generationDomain.maxCorner());
   p1.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,-1) ));
   p1.setOwner(mpi::MPIManager::instance()->rank());
   p1.setType(0);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   auto sphereShape = ss->create<mesa_pd::data::Sphere>( diameter * real_t(0.5) );
   ss->shapes[sphereShape]->updateMassAndInertia(densityRatio);

   auto xParticle = real_t(0);
   auto yParticle = real_t(0);
   auto zParticle = real_t(0);

   auto rank = mpi::MPIManager::instance()->rank();

   for( uint_t nPart = 0; nPart < numberOfParticles; ++nPart )
   {

      WALBERLA_ROOT_SECTION()
      {
         xParticle = math::realRandom<real_t>(generationDomain.xMin(), generationDomain.xMax());
         yParticle = math::realRandom<real_t>(generationDomain.yMin(), generationDomain.yMax());
         zParticle = math::realRandom<real_t>(generationDomain.zMin(), generationDomain.zMax());

      }

      WALBERLA_MPI_SECTION()
      {
         mpi::broadcastObject( xParticle );
         mpi::broadcastObject( yParticle );
         mpi::broadcastObject( zParticle );
      }

      auto position = Vector3<real_t>( xParticle, yParticle, zParticle );
      //WALBERLA_LOG_INFO_ON_ROOT(position);

      if (!rpdDomain->isContainedInProcessSubdomain(uint_c(rank), position)) continue;
      auto p                       = ps->create();
      p->setPosition(position);
      p->setInteractionRadius(diameter * real_t(0.5));
      p->setShapeID(sphereShape);
      p->setType(1);
      p->setOwner(rank);
   }

   syncCall();

   // resolve possible overlaps of the particles due to the random initialization via a particle-only simulation

   const bool useOpenMP = false;
   const real_t dt_RPD_Init = real_t(1);
   const auto initialParticleSimSteps = uint_t(20000);
   const real_t overlapLimit = real_t(0.001) * diameter;

   auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   particleVtkOutput->addOutput<mesa_pd::data::SelectParticleLinearVelocity>("velocity");
   particleVtkOutput->setParticleSelector( [sphereShape](const mesa_pd::data::ParticleStorage::iterator& pIt) {return !mesa_pd::data::particle_flags::isSet(pIt->getFlags(), mesa_pd::data::particle_flags::GHOST) && pIt->getShapeID() == sphereShape;} ); //limit output to local sphere
   auto particleVtkWriterInit = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles_init", 1, baseFolder, "simulation_step");

   for(auto pet = uint_t(0); pet <= initialParticleSimSteps; ++pet )
   {


      real_t maxPenetrationDepth = 0;
      ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                  [&collisionResponse, &rpdDomain, &maxPenetrationDepth, dt_RPD_Init]
                                        (const size_t idx1, const size_t idx2, auto& ac)
                                  {
                                     mesa_pd::collision_detection::AnalyticContactDetection acd;
                                     mesa_pd::kernel::DoubleCast double_cast;
                                     mesa_pd::mpi::ContactFilter contact_filter;
                                     if (double_cast(idx1, idx2, ac, acd, ac ))
                                     {
                                        if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                        {
                                           maxPenetrationDepth = std::max(maxPenetrationDepth, std::abs(acd.getPenetrationDepth()));
                                           collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(),
                                                 acd.getContactNormal(), acd.getPenetrationDepth(), dt_RPD_Init);
                                        }
                                     }
                                  },
                                  *accessor );

      reduceAndSwapContactHistory(*ps);

      mpi::allReduceInplace(maxPenetrationDepth, mpi::MAX);

      reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

      ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, mesa_pd::kernel::ExplicitEuler(dt_RPD_Init), *accessor);
      syncCall();

      if( pet % uint_t(20) == uint_t(0) )
      {
         if(vtkDuringInit)
         {
            particleVtkWriterInit->write();
         }
         WALBERLA_LOG_INFO_ON_ROOT(pet << " - current max overlap = " << maxPenetrationDepth / diameter * real_t(100) << "%");
      }

      if(maxPenetrationDepth < overlapLimit) break;

      // reset velocites to avoid too large ones

      ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
            [](const size_t idx, ParticleAccessor_T& ac){
               ac.setLinearVelocity(idx, ac.getLinearVelocity(idx) * real_t(0.5));
               ac.setAngularVelocity(idx, ac.getAngularVelocity(idx) * real_t(0.5));
               }, *accessor);


   }


   // reset all velocities to zero
   Vector3<real_t> initialSphereVelocity(real_t(0));

   ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, [](const size_t idx, ParticleAccessor_T& ac){
      ac.getNewContactHistoryRef(idx).clear();
      ac.getOldContactHistoryRef(idx).clear();
   }, *accessor);

   WALBERLA_LOG_INFO_ON_ROOT("Setting initial velocity " << initialSphereVelocity << " of all spheres");
   ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor,
                        [initialSphereVelocity](const size_t idx, ParticleAccessor_T& ac){
                           ac.setLinearVelocity(idx, initialSphereVelocity);
                           ac.setAngularVelocity(idx, Vector3<real_t>(real_t(0)));
                        }, *accessor);

   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         uint_t(1), field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::zyxf, FieldGhostLayers );

   // add boundary handling & initialize outer domain boundaries
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor, useEntireFieldTraversal), "boundary handling" );

   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densityRatio - real_t(1)) * gravitationalAcceleration * sphereVolume );
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(viscosity, [](real_t r){return (real_t(0.001 + real_t(0.00007)*r))*r;});
   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);

   WALBERLA_LOG_INFO_ON_ROOT(" - Lubrication correction:");
   WALBERLA_LOG_INFO_ON_ROOT("   - normal cut off distance = " << lubricationCorrectionKernel.getNormalCutOffDistance());
   WALBERLA_LOG_INFO_ON_ROOT("   - tangential translational cut off distance = " << lubricationCorrectionKernel.getTangentialTranslationalCutOffDistance());
   WALBERLA_LOG_INFO_ON_ROOT("   - tangential rotational cut off distance = " << lubricationCorrectionKernel.getTangentialRotationalCutOffDistance());
   const real_t maximumLubricationCutOffDistance = std::max(lubricationCorrectionKernel.getNormalCutOffDistance(), std::max(lubricationCorrectionKernel.getTangentialRotationalCutOffDistance(), lubricationCorrectionKernel.getTangentialTranslationalCutOffDistance()));


   // map planes into the LBM simulation -> act as no-slip boundaries
   ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

   // map particles into the LBM simulation
   ps->forEachParticle(false, lbm_mesapd_coupling::RegularParticlesSelector(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

   lbm::BlockForestEvaluation<FlagField_T> bfEval(blocks, flagFieldID, Fluid_Flag);

   WALBERLA_LOG_INFO_ON_ROOT(bfEval.loggingString());

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   timeloop.addFuncBeforeTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   if( vtkIOFreq != uint_t(0) )
   {
      auto particleVtkWriter = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", vtkIOFreq, baseFolder, "simulation_step");
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( particleVtkWriter ), "VTK (sphere data)" );

      // flag field
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", vtkIOFreq, 1, false, baseFolder );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", vtkIOFreq, 0, false, baseFolder );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );
      pdfFieldVTK->addCellInclusionFilter( fluidFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );

      auto domainDecompVTK = vtk::createVTKOutput_DomainDecomposition(blocks, "domain_decomposition", vtkIOFreq, baseFolder );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles(domainDecompVTK), "VTK (domain decomposition)");
   }


   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme( blocks );
   optimizedPDFCommunicationScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) ); // optimized sync

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

   if( !useFusedStreamCollide )
   {
      // Collide
      timeloop.add() << Sweep( makeCollideSweep(sweep), "Collide" );
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   if( useFusedStreamCollide )
   {
      // streaming & collide
      timeloop.add() << Sweep( makeSharedSweep(sweep), "Stream&Collide" );
   } else
   {
      // streaming
      timeloop.add() << Sweep( makeStreamSweep(sweep), "Stream" );

   }

   SweepTimeloop timeloopAfterParticles( blocks->getBlockStorage(), timesteps );

   bool conserveMomentum = false;
   // sweep for updating the particle mapping into the LBM simulation
   timeloopAfterParticles.add() << Sweep( lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID,boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag,
         lbm_mesapd_coupling::RegularParticlesSelector(), conserveMomentum), "Particle Mapping" );

   bool recomputeTargetDensity = false;
   auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, recomputeTargetDensity,true);


   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme( blocks );
   fullPDFCommunicationScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) ); // full sync

   timeloopAfterParticles.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                                << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag,
                                      gradReconstructor, conserveMomentum) ), "PDF Restore" );


   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   std::vector< std::vector<std::string> > timerKeys;
   std::vector<std::string> LBMTimer;
   LBMTimer.emplace_back("Stream&Collide");
   LBMTimer.emplace_back("Stream");
   LBMTimer.emplace_back("Collide");
   timerKeys.push_back(LBMTimer);

   std::vector<std::string> bhTimer;
   bhTimer.emplace_back("Boundary Handling");
   timerKeys.push_back(bhTimer);

   std::vector<std::string> couplingTimer1;
   couplingTimer1.emplace_back("Particle Mapping");
   std::vector<std::string> couplingTimer2;
   couplingTimer2.emplace_back("PDF Restore");
   timerKeys.push_back(couplingTimer1);
   timerKeys.push_back(couplingTimer2);

   std::vector<std::string> rpdTimer;
   rpdTimer.emplace_back("RPD Force");
   rpdTimer.emplace_back("RPD VV1");
   rpdTimer.emplace_back("RPD VV2");
   rpdTimer.emplace_back("RPD Lub");
   rpdTimer.emplace_back("RPD Collision");
   timerKeys.push_back(rpdTimer);

   uint_t numCells = uint_t(0);
   uint_t numFluidCells = uint_t(0);
   uint_t numNBCells = uint_t(0);
   uint_t numLocalParticles = uint_t(0);
   uint_t numGhostParticles = uint_t(0);
   uint_t numContacts = uint_t(0);

   std::vector<double> timings(timerKeys.size());

   // every rank writes its own file -> numProcesses number of samples!
   int myRank = MPIManager::instance()->rank();

   std::string logFileName = baseFolder + "/load";
   logFileName += "_settling";
   logFileName += "_spheres";
   logFileName += "_d" + std::to_string(int_c(std::ceil(diameter)));
   logFileName += "_bs" + std::to_string(blockSize);
   logFileName += "_" + std::to_string(myRank) + ".txt";


   std::ofstream file;

   if(!noFileOutput)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Writing load info to file " << logFileName);
      file.open( logFileName.c_str());
      file << "# svf = " << solidVolumeFraction << ", d = " << diameter << ", domain = " << XCells << "x" << YCells << "x" << ZCells << "\n";
   }


   auto timeStepOfFirstTiming = uint_t(50);

   if( timesteps - timeStepOfFirstTiming < numSamples )
   {
      WALBERLA_LOG_WARNING_ON_ROOT("Less actual time steps than number of required samples!");
   }

   uint_t nSample( 0 ); // number of current sample
   real_t samplingFrequency = real_c(timesteps - timeStepOfFirstTiming) / real_c(numSamples);

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {
      // perform a single simulation step

      timeloop.singleStep( timeloopTiming );

      reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);

      timeloopTiming["RPD Force"].start();
      if( i == 0 )
      {
         lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueForAveragingKernel;
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor );
      }
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor );
      timeloopTiming["RPD Force"].end();

      for(auto subCycle = uint_t(0); subCycle < numRPDSubCycles; ++subCycle )
      {

         timeloopTiming["RPD VV1"].start();
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
         timeloopTiming["RPD VV1"].end();

         syncCall();

         // lubrication correction
         timeloopTiming["RPD Lub"].start();
         ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                     [&lubricationCorrectionKernel,maximumLubricationCutOffDistance, &rpdDomain,&numContacts]
                                     (const size_t idx1, const size_t idx2, auto& ac)
                                     {
                                        mesa_pd::collision_detection::AnalyticContactDetection acd;
                                        acd.getContactThreshold() = maximumLubricationCutOffDistance;
                                        mesa_pd::kernel::DoubleCast double_cast;
                                        mesa_pd::mpi::ContactFilter contact_filter;
                                        if (double_cast(idx1, idx2, ac, acd, ac ))
                                        {
                                           if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                           {
                                              double_cast(acd.getIdx1(), acd.getIdx2(), ac, lubricationCorrectionKernel, ac, acd.getContactNormal(), acd.getPenetrationDepth());
                                              ++numContacts; // this then also includes all actual contacts since there the contact threshold is smaller
                                           }
                                        }
                                     },
                                     *accessor );
         timeloopTiming["RPD Lub"].end();

         // collision response
         timeloopTiming["RPD Collision"].start();
         ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                     [&collisionResponse, &rpdDomain, timeStepSizeRPD]
                                           (const size_t idx1, const size_t idx2, auto& ac)
                                     {
                                        mesa_pd::collision_detection::AnalyticContactDetection acd;
                                        mesa_pd::kernel::DoubleCast double_cast;
                                        mesa_pd::mpi::ContactFilter contact_filter;
                                        if (double_cast(idx1, idx2, ac, acd, ac ))
                                        {
                                           if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                           {
                                              collisionResponse(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), acd.getContactNormal(), acd.getPenetrationDepth(), timeStepSizeRPD);
                                           }
                                        }
                                     },
                                     *accessor );

         timeloopTiming["RPD Collision"].end();

         reduceAndSwapContactHistory(*ps);

         timeloopTiming["RPD Force"].start();
         lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor );
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );
         timeloopTiming["RPD Force"].end();

         reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

         // integration
         timeloopTiming["RPD VV2"].start();
         ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         timeloopTiming["RPD VV2"].end();

         syncCall();


      }

      timeloopTiming["RPD Force"].start();
      ps->forEachParticle( useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );
      timeloopTiming["RPD Force"].end();


      // update particle mapping
      timeloopAfterParticles.singleStep(timeloopTiming);

      //WALBERLA_LOG_INFO_ON_ROOT(timeloopTiming);

      // check if current time step should be included in sample
      if( i >= uint_c( samplingFrequency * real_c(nSample) ) + timeStepOfFirstTiming )
      {
         // include -> evaluate all timers and quantities

         evaluateFluidQuantities<BoundaryHandling_T>(blocks, boundaryHandlingID, numCells, numFluidCells, numNBCells);
         evaluateRPDQuantities(ps, numLocalParticles, numGhostParticles);

         evaluateTimers(timeloopTiming, timerKeys, timings);

         if(!noFileOutput)
         {
            auto totalTime = std::accumulate(timings.begin(), timings.end(), 0.0 );

            file << timeloop.getCurrentTimeStep() << " " << real_c(timeloop.getCurrentTimeStep()) / Tref << " "
                 << numCells << " " << numFluidCells << " " << numNBCells << " "
                 << numLocalParticles << " " << numGhostParticles << " "
                 << real_c(numContacts) / real_c(numRPDSubCycles) << " " << numRPDSubCycles;
            for (auto timing : timings) {
               file << " " << timing;
            }
            file << " " << totalTime << "\n";
         }

         ++nSample;
      }

      numCells = uint_t(0);
      numFluidCells = uint_t(0);
      numNBCells = uint_t(0);
      numLocalParticles = uint_t(0);
      numGhostParticles = uint_t(0);
      numContacts = uint_t(0);

      // reset timers to always include only a single time step in them
      timeloopTiming.clear();
   }

   if(!noFileOutput) {
      file.close();
   }

   WALBERLA_LOG_INFO_ON_ROOT("Simulation finished!");

   return 0;

}

}

int main( int argc, char **argv ){
   fluid_particle_workload_evaluation::main(argc, argv);
}
