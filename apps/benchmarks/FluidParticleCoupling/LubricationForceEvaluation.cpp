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
//! \file LubricationForceEvaluation.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/mpi/Broadcast.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/logging/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/LubricationCorrectionKernel.h"
#include "lbm_mesapd_coupling/utility/OmegaBulkAdaption.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/HalfSpace.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>

#ifdef WALBERLA_BUILD_WITH_CODEGEN
#include "GeneratedLBM.h"
#endif

namespace lubrication_force_evaluation
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

#ifdef WALBERLA_BUILD_WITH_CODEGEN
using LatticeModel_T = lbm::GeneratedLBM;
#else
using LatticeModel_T = lbm::D3Q19< lbm::collision_model::D3Q19MRT>;
#endif

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

using ScalarField_T = GhostLayerField< real_t, 1>;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID FreeSlip_Flag( "free slip" );
const FlagUID MO_SBB_Flag( "moving obstacle sbb" );
const FlagUID MO_CLI_Flag( "moving obstacle cli" );
const FlagUID FormerMO_Flag( "former moving obstacle" );

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template <typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using FreeSlip_T = lbm::FreeSlip< LatticeModel_T, FlagField_T>;
   using MO_SBB_T = lbm_mesapd_coupling::SimpleBB< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using MO_CLI_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, FreeSlip_T, MO_SBB_T, MO_CLI_T >;

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T>& ac) :
         flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac ) {}

   Type * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( storage );

      auto * flagField     = block->getData< FlagField_T >( flagFieldID_ );
      auto *  pdfField     = block->getData< PdfField_T > ( pdfFieldID_ );
      auto * particleField = block->getData< lbm_mesapd_coupling::ParticleField_T > ( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

      Type * handling = new Type( "moving obstacle boundary handling", flagField, fluid,
                                  FreeSlip_T( "FreeSlip", FreeSlip_Flag, pdfField, flagField, fluid ),
                                  MO_SBB_T( "SBB", MO_SBB_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ),
                                  MO_CLI_T( "CLI", MO_CLI_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block )  );

      const auto freeslip = flagField->getFlag( FreeSlip_Flag );

      CellInterval domainBB = storage->getDomainCellBB();
      domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
      domainBB.xMax() += cell_idx_c( FieldGhostLayers );

      domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
      domainBB.yMax() += cell_idx_c( FieldGhostLayers );

      // SOUTH
      CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
      storage->transformGlobalToBlockLocalCellInterval( south, *block );
      handling->forceBoundary( freeslip, south );

      // NORTH
      CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      storage->transformGlobalToBlockLocalCellInterval( north, *block );
      handling->forceBoundary( freeslip, north );

      domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
      domainBB.zMax() += cell_idx_c( FieldGhostLayers );

      // BOTTOM
      CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
      storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
      handling->forceBoundary( freeslip, bottom );

      // TOP
      CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
      storage->transformGlobalToBlockLocalCellInterval( top, *block );
      handling->forceBoundary( freeslip, top );

      handling->fillWithDomain( FieldGhostLayers );


      return handling;
   }

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;
};
//*******************************************************************************************************************

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Evaluates the hydrodynamic force for the lubrication case for sphere-sphere and sphere-wall case.
 *
 * 4 different setups are available that change the relative velocity to investigate the different components individually.
 * All particles are fixed but have a small velocity which is a valid assumption in Stokes flow.
 * The simulations are run until steady state is reached.
 *
 * The positions are always shifted by a random offset to account for the dependence of the resulting force/torque
 * on the exact location w.r.t. the underlying grid.
 * This can be eliminated by averaging over several realizations of the same physical setup.
 *
 * see also Rettinger, Ruede 2020 for the details
 */
//*******************************************************************************************************************

int main( int argc, char **argv )
{

   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   bool sphSphTest     = true;
   bool fileIO         = true;
   uint_t vtkIOFreq = 0;
   std::string fileNameEnding = "";
   std::string baseFolder = "vtk_out_Lubrication";

   real_t radius = real_t(5);
   real_t ReynoldsNumber = real_t(1e-2);
   real_t tau = real_t(1);
   real_t gapSize = real_t(0);
   real_t bulkViscRateFactor = real_t(1);
   real_t magicNumber = real_t(3)/real_t(16);
   bool useOmegaBulkAdaption = false;
   real_t adaptionLayerSize = real_t(2);
   bool useSBB = false;

   // 1: translation in normal direction -> normal Lubrication force
   // 2: translation in tangential direction -> tangential Lubrication force and torque
   // 3: rotation around tangential direction -> force & torque
   // 4: rotation around normal direction -> torque
   uint_t setup = 1;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--sphWallTest" )          == 0 ) { sphSphTest = false; continue; }
      if( std::strcmp( argv[i], "--noLogging" )            == 0 ) { fileIO = false; continue; }
      if( std::strcmp( argv[i], "--vtkIOFreq" )            == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--setup" )                == 0 ) { setup = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--baseFolder" )           == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--diameter" )             == 0 ) { radius = real_t(0.5)*real_c(std::atof(argv[++i])); continue; }
      if( std::strcmp( argv[i], "--gapSize" )              == 0 ) { gapSize = real_c(std::atof(argv[++i])); continue; }
      if( std::strcmp( argv[i], "--bulkViscRateFactor" )   == 0 ) { bulkViscRateFactor = real_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--tau" )                  == 0 ) { tau = real_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--fileName" )             == 0 ) { fileNameEnding = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--magicNumber" )          == 0 ) { magicNumber = real_c(std::atof(argv[++i])); continue; }
      if( std::strcmp( argv[i], "--useOmegaBulkAdaption" ) == 0 ) { useOmegaBulkAdaption = true; continue; }
      if( std::strcmp( argv[i], "--adaptionLayerSize" )    == 0 ) { adaptionLayerSize = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--useSBB" )               == 0 ) { useSBB = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   uint_t xSize = uint_c(real_t(24) * radius);
   uint_t ySize = uint_c(real_t(24) * radius);
   uint_t zSize = uint_c(real_t(24) * radius);

   uint_t xBlocks = uint_c(4); // number of blocks in x-direction
   uint_t yBlocks = uint_c(1); // number of blocks in y-direction
   uint_t zBlocks = uint_c(1); // number of blocks in z-direction

   uint_t xCells = xSize / xBlocks; // number of cells in x-direction on each block
   uint_t yCells = ySize / yBlocks; // number of cells in y-direction on each block
   uint_t zCells = zSize / zBlocks; // number of cells in z-direction on each block

   // Perform missing variable calculations
   real_t omega = real_t(1) / tau;
   real_t nu = walberla::lbm::collision_model::viscosityFromOmega(omega);
   real_t velocity = ReynoldsNumber * nu / (real_t(2) * radius);
   real_t omegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, bulkViscRateFactor);

   uint_t timesteps = uint_c( 10000 );

   real_t fStokes = real_t(6) * math::pi * nu * radius * velocity;
   real_t tStokes = real_t(8) * math::pi * nu * radius * radius * velocity;

   WALBERLA_LOG_INFO_ON_ROOT_SECTION()
   {
      std::stringstream ss;

      if ( sphSphTest )
      {
         ss << "-------------------------------------------------------\n"
            << "   Parameters for the sphere-sphere lubrication test \n"
            << "-------------------------------------------------------\n";
      }
      else
      {
         ss << "-------------------------------------------------------\n"
            << "   Parameters for the sphere-wall lubrication test \n"
            << "-------------------------------------------------------\n";
      }
      ss << " omega        = " << omega    << "\n"
         << " omegaBulk    = " << omegaBulk << " (bvrf " << bulkViscRateFactor << ")\n"
         << " use omega bulk adaption = " << useOmegaBulkAdaption << " (adaption layer size = " << adaptionLayerSize << ")\n"
         << " radius       = " << radius   << "\n"
         << " velocity     = " << velocity << "\n"
         << " Re           = " << ReynoldsNumber << "\n"
         << " gap size     = " << gapSize << "\n"
         << " time steps   = " << timesteps << "\n"
         << " fStokes      = " << fStokes << "\n"
         << " setup        = " << setup << "\n"
         << " useSBB       = " << useSBB << "\n"
         << "-------------------------------------------------------\n"
         << " domainSize = " << xSize << " x " << ySize << " x " << zSize  << "\n"
         << " blocks     = " << xBlocks << " x " << yBlocks << " x " << zBlocks  << "\n"
         << " blockSize  = " << xCells << " x " << yCells << " x " << zCells  << "\n"
         << "-------------------------------------------------------\n";
      WALBERLA_LOG_INFO( ss.str() );
   }


   auto blocks = blockforest::createUniformBlockGrid( xBlocks, yBlocks, zBlocks, xCells, yCells, zCells, real_t(1),
                                                      0, false, false,
                                                      sphSphTest, true, true,  //periodicity
                                                      false );

   //////////////////
   // RPD COUPLING //
   //////////////////

   auto rpdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   //init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );

   uint_t id1 (0);
   uint_t id2 (0);

   uint_t randomSeed = uint_c(std::chrono::system_clock::now().time_since_epoch().count());
   mpi::broadcastObject(randomSeed); // root process chooses seed and broadcasts it
   std::mt19937 randomNumberGenerator(static_cast<unsigned int>(randomSeed));

   Vector3<real_t> domainCenter( real_c(xSize) * real_t(0.5), real_c(ySize) * real_t(0.5), real_c(zSize) * real_t(0.5) );
   Vector3<real_t> offsetVector(math::realRandom<real_t>(real_t(0), real_t(1), randomNumberGenerator),
                                math::realRandom<real_t>(real_t(0), real_t(1), randomNumberGenerator),
                                math::realRandom<real_t>(real_t(0), real_t(1), randomNumberGenerator));

   if ( sphSphTest )
   {

      Vector3<real_t> pos1 = domainCenter + offsetVector;
      if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), pos1 ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos1);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         if(setup == 1)
         {
            // only normal translational vel
            p.setLinearVelocity(Vector3<real_t>( velocity, real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
         } else if (setup == 2)
         {
            // only tangential translational velocity
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), velocity));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
         } else if (setup == 3)
         {
            // only rotation around axis perpendicular to center line
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), velocity / radius, real_t(0)));
         } else if (setup == 4)
         {
            // only rotation around center line
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( velocity / radius, real_t(0) , real_t(0)));
         }
         id1 = p.getUid();
      }

      Vector3<real_t> pos2 = pos1 + Vector3<real_t>(real_t(2) * radius + gapSize, real_t(0), real_t(0));
      if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), pos2 ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos2);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         if(setup == 1)
         {
            // only normal translational vel
            p.setLinearVelocity(Vector3<real_t>( -velocity, real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
         } else if (setup == 2)
         {
            // only tangential translational velocity
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), -velocity));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
         } else if (setup == 3)
         {
            // only rotation around axis perpendicular to center line
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), velocity / radius, real_t(0)));
         } else if (setup == 4)
         {
            // only rotation around center line
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( -velocity / radius, real_t(0) , real_t(0)));
         }
         id2 = p.getUid();
      }

      mpi::allReduceInplace(id1, mpi::SUM);
      mpi::allReduceInplace(id2, mpi::SUM);

      WALBERLA_LOG_INFO_ON_ROOT("pos sphere 1 = " << pos1);
      WALBERLA_LOG_INFO_ON_ROOT("pos sphere 2 = " << pos2);
   } else
   {
      // sphere-wall test

      Vector3<real_t> referenceVector(offsetVector[0], domainCenter[1], domainCenter[2]);

      // create two planes
      mesa_pd::data::Particle&& p0 = *ps->create(true);
      p0.setPosition(referenceVector);
      p0.setInteractionRadius(std::numeric_limits<real_t>::infinity());
      p0.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(1,0,0) ));
      p0.setOwner(mpi::MPIManager::instance()->rank());
      mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
      mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
      id2 = p0.getUid();

      mesa_pd::data::Particle&& p1 = *ps->create(true);
      p1.setPosition(Vector3<real_t>(real_c(xSize),0,0));
      p1.setInteractionRadius(std::numeric_limits<real_t>::infinity());
      p1.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(-1,0,0) ));
      p1.setOwner(mpi::MPIManager::instance()->rank());
      mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
      mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

      Vector3<real_t> pos1 = referenceVector + Vector3<real_t>(radius + gapSize, real_t(0), real_t(0));
      if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), pos1 ))
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(pos1);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         if(setup == 1)
         {
            // only normal translational vel
            p.setLinearVelocity(Vector3<real_t>( -velocity, real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
         } else if (setup == 2)
         {
            // only tangential translational velocity
            p.setLinearVelocity(Vector3<real_t>( real_t(0), velocity, real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
         } else if (setup == 3)
         {
            // only rotation around axis perpendicular to center line
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( real_t(0), velocity / radius, real_t(0)));
         } else if (setup == 4)
         {
            // only rotation around center line
            p.setLinearVelocity(Vector3<real_t>( real_t(0), real_t(0), real_t(0)));
            p.setAngularVelocity(Vector3<real_t>( velocity / radius, real_t(0) , real_t(0)));
         }
         id1 = p.getUid();
      }

      mpi::allReduceInplace(id1, mpi::SUM);
      // id2 is globally known

      WALBERLA_LOG_INFO_ON_ROOT("pos plane = " << referenceVector);
      WALBERLA_LOG_INFO_ON_ROOT("pos sphere = " << pos1);
   }

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add omega bulk field
   BlockDataID omegaBulkFieldID = field::addToStorage<ScalarField_T>( blocks, "omega bulk field", omegaBulk, field::fzyx);

   // create the lattice model
   real_t lambda_e = lbm::collision_model::TRT::lambda_e( omega );
   real_t lambda_d = lbm::collision_model::TRT::lambda_d( omega, magicNumber );
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   WALBERLA_LOG_INFO_ON_ROOT("Using generated TRT-like lattice model!");
   LatticeModel_T latticeModel = LatticeModel_T(omegaBulkFieldID, lambda_d, lambda_e);
#else
   WALBERLA_LOG_INFO_ON_ROOT("Using waLBerla built-in MRT lattice model and ignoring omega bulk field since not supported!");
   LatticeModel_T latticeModel = LatticeModel_T(lbm::collision_model::D3Q19MRT( omegaBulk, omegaBulk, lambda_d, lambda_e, lambda_e, lambda_d ));
#endif

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         uint_t(1), field::fzyx );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );

   // set up RPD functionality
   std::function<void(void)> syncCall = [&ps,&rpdDomain](){
      const real_t overlap = real_t( 1.5 );
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();

   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;

   real_t lubricationCutOffDistanceNormal = real_t(2) / real_t(3);
   real_t lubricationCutOffDistanceTangentialTranslational = real_t(0.5);
   real_t lubricationCutOffDistanceTangentialRotational = real_t(0.5);
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(nu, [](real_t){return real_t(0);}, lubricationCutOffDistanceNormal, lubricationCutOffDistanceTangentialTranslational, lubricationCutOffDistanceTangentialRotational );
   real_t maximalCutOffDistance = std::max(lubricationCutOffDistanceNormal, std::max(lubricationCutOffDistanceTangentialTranslational,lubricationCutOffDistanceTangentialRotational ));

   lbm_mesapd_coupling::RegularParticlesSelector sphereSelector;

   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);

   ///////////////
   // TIME LOOP //
   ///////////////

   // map particles into the LBM simulation
   if(useSBB)
   {
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, movingParticleMappingKernel, *accessor, MO_SBB_Flag);
   } else
   {
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, movingParticleMappingKernel, *accessor, MO_CLI_Flag);
   }

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme( blocks );
   optimizedPDFCommunicationScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) ); // optimized sync

   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme( blocks );
   fullPDFCommunicationScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) ); // full sync

   timeloop.addFuncBeforeTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );


   if( vtkIOFreq != uint_t(0) )
   {
      // spheres
      auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleOwner>("owner");
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleLinearVelocity>("velocity");
      auto particleVtkWriter = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", vtkIOFreq, baseFolder, "simulation_step");
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( particleVtkWriter ), "VTK (sphere data)" );

      // flag field (written only once in the first time step, ghost layers are also written)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", vtkIOFreq, FieldGhostLayers, false, baseFolder );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      vtk::writeFiles( flagFieldVTK )();

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", vtkIOFreq, 0, false, baseFolder );

      pdfFieldVTK->addBeforeFunction( fullPDFCommunicationScheme );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );
      pdfFieldVTK->addCellInclusionFilter( fluidFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );

      // omega bulk field
      timeloop.addFuncBeforeTimeStep( field::createVTKOutput<ScalarField_T, float>( omegaBulkFieldID, *blocks, "omega_bulk_field", vtkIOFreq, uint_t(0), false, baseFolder ), "VTK (omega bulk field)" );

   }

   // update bulk omega in all cells to adapt to changed particle position
   if(useOmegaBulkAdaption)
   {
      using OmegaBulkAdapter_T = lbm_mesapd_coupling::OmegaBulkAdapter<ParticleAccessor_T, decltype(sphereSelector)>;
      real_t defaultOmegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, real_t(1));
      shared_ptr<OmegaBulkAdapter_T> omegaBulkAdapter = make_shared<OmegaBulkAdapter_T>(blocks, omegaBulkFieldID, accessor, defaultOmegaBulk, omegaBulk, adaptionLayerSize, sphereSelector);
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
      {
         (*omegaBulkAdapter)(&(*blockIt));
      }
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)

   timeloop.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // stream + collide LBM step
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   auto lbmSweep = LatticeModel_T::Sweep( pdfFieldID );
   timeloop.add() << Sweep( lbmSweep, "LB sweep" );
#else
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   timeloop.add() << Sweep( makeSharedSweep( lbmSweep ), "cell-wise LB sweep" );
#endif


   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   Vector3<real_t> hydForce(0.);
   Vector3<real_t> lubForce(0.);
   Vector3<real_t> hydTorque(0.);
   Vector3<real_t> lubTorque(0.);

   real_t curForceNorm = real_t(0);
   real_t oldForceNorm = real_t(0);
   real_t curTorqueNorm = real_t(0);
   real_t oldTorqueNorm = real_t(0);

   real_t convergenceLimit = real_t(1e-5);

   // time loop
   for (uint_t i = 1; i <= timesteps; ++i )
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      timeloop.singleStep( timeloopTiming );

      // lubrication correction
      mesa_pd::collision_detection::AnalyticContactDetection acd;
      acd.getContactThreshold() = maximalCutOffDistance;

      {
         auto idx1 = accessor->uidToIdx(id1);
         if( idx1 != accessor->getInvalidIdx() )
         {
            auto idx2 = accessor->uidToIdx(id2);
            if( idx2 != accessor->getInvalidIdx() )
            {
               mesa_pd::kernel::DoubleCast double_cast;
               mesa_pd::mpi::ContactFilter contact_filter;
               if (double_cast(idx1, idx2, *accessor, acd, *accessor ))
               {
                  if (contact_filter(acd.getIdx1(), acd.getIdx2(), *accessor, acd.getContactPoint(), *rpdDomain))
                  {
                     double_cast(acd.getIdx1(), acd.getIdx2(), *accessor, lubricationCorrectionKernel, *accessor, acd.getContactNormal(), acd.getPenetrationDepth());
                  }
               }
            }
         }
      }

      if( i% 100 == 0 && i > 1 )
      {
         oldForceNorm = curForceNorm;
         oldTorqueNorm = curTorqueNorm;

         hydForce.reset();
         lubForce.reset();
         hydTorque.reset();
         lubTorque.reset();

         auto idx1 = accessor->uidToIdx(id1);
         if( idx1 != accessor->getInvalidIdx() )
         {
            hydForce = accessor->getHydrodynamicForce(idx1);
            lubForce = accessor->getForce(idx1);
            hydTorque = accessor->getHydrodynamicTorque(idx1);
            lubTorque= accessor->getTorque(idx1);
         }

         WALBERLA_MPI_SECTION()
         {
            mpi::allReduceInplace( hydForce, mpi::SUM );
            mpi::reduceInplace( lubForce, mpi::SUM );
            mpi::allReduceInplace( hydTorque, mpi::SUM );
            mpi::reduceInplace( lubTorque, mpi::SUM );
         }

         curForceNorm = hydForce.length();
         curTorqueNorm = hydTorque.length();

         real_t forceDiff = std::fabs((curForceNorm - oldForceNorm) / oldForceNorm);
         real_t torqueDiff = std::fabs((curTorqueNorm - oldTorqueNorm) / oldTorqueNorm);

         WALBERLA_LOG_INFO_ON_ROOT( "F/Fs = " << hydForce/fStokes << " ( " << forceDiff << " ), T/Ts = " << hydTorque/tStokes << " ( " << torqueDiff << " )");

         if( i == 100 ) {
            WALBERLA_LOG_INFO_ON_ROOT("Flub = " << lubForce << ", Tlub = " << lubTorque);
            WALBERLA_LOG_INFO_ON_ROOT("Flub/Fs = " << lubForce/fStokes << ", Tlub/Ts = " << lubTorque/tStokes);
         }

         if( forceDiff < convergenceLimit && torqueDiff < convergenceLimit )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Force and torque norms converged - terminating simulation");
            break;
         }

      }

      // reset forces
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor,
                          [](const size_t idx, auto& ac){
                             ac.getForceRef(idx) = Vector3<real_t>(real_t(0));
                             ac.getTorqueRef(idx) = Vector3<real_t>(real_t(0));
                          }, *accessor );
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );


   }

   if( fileIO )
   {
      std::string loggingFileName( baseFolder + "/Logging" );
      if( sphSphTest ) loggingFileName += "_SphSph";
      else loggingFileName += "_SphPla";
      loggingFileName += "_Setup" + std::to_string(setup);
      loggingFileName += "_gapSize" + std::to_string(uint_c(gapSize*real_t(100)));
      loggingFileName += "_radius" + std::to_string(uint_c(radius));
      loggingFileName += "_bvrf" + std::to_string(uint_c(bulkViscRateFactor));
      loggingFileName += "_mn" + std::to_string(float(magicNumber));
      if( useOmegaBulkAdaption ) loggingFileName += "_uOBA" + std::to_string(uint_c(adaptionLayerSize));
      if( useSBB ) loggingFileName += "_SBB";
      if( !fileNameEnding.empty()) loggingFileName += "_" + fileNameEnding;
      loggingFileName += ".txt";

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file1;
         file1.open( loggingFileName.c_str(), std::ofstream::app );
         file1.setf( std::ios::unitbuf );
         file1.precision(15);
         file1 << radius << " " << gapSize << " " << fStokes << " "
               << hydForce[0] << " " << hydForce[1] << " " << hydForce[2] << " "
               << lubForce[0] << " " << lubForce[1] << " " << lubForce[2] << " "
               << hydTorque[0] << " " << hydTorque[1] << " " << hydTorque[2] << " "
               << lubTorque[0] << " " << lubTorque[1] << " " << lubTorque[2] << std::endl;
         file1.close();
      }
   }

   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

} // namespace lubrication_force_evaluation

int main( int argc, char **argv ){
   lubrication_force_evaluation::main(argc, argv);
}
