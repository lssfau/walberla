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
//! \file SphereMovingWithPrescribedVelocity.cpp
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
#include "core/timing/RemainingTimeLogger.h"
#include "core/logging/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/interpolators/FieldInterpolatorCreators.h"
#include "field/interpolators/TrilinearFieldInterpolator.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"
#include "field/adaptors/AdaptorCreators.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/Adaptors.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/Reconstructor.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/ExtrapolationDirectionFinder.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/reconstruction/PdfReconstructionManager.h"
#include "lbm_mesapd_coupling/utility/AddForceOnParticlesKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/AverageHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/AddHydrodynamicInteractionKernel.h"
#include "lbm_mesapd_coupling/utility/InitializeHydrodynamicForceTorqueForAveragingKernel.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
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
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/HydrodynamicForceTorqueNotification.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>

#ifdef WALBERLA_BUILD_WITH_CODEGEN
#include "GeneratedLBM.h"
#endif

namespace sphere_moving_with_prescribed_velocity
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
const FlagUID NoSlip_Flag( "no slip" );
const FlagUID MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template <typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using NoSlip_T = lbm::NoSlip< LatticeModel_T, flag_t >;
   using MO_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T >;

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
                                  NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                  MO_T( "MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ) );

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


//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the sphere
 *
 */
//*******************************************************************************************************************
template< typename ParticleAccessor_T>
class SpherePropertyLogger
{
public:
   SpherePropertyLogger( const shared_ptr< ParticleAccessor_T > & ac, walberla::id_t sphereUid,
         const std::string & fileName, bool fileIO,
         real_t diameter, real_t velocity, real_t forceMag) :
      ac_( ac ), sphereUid_( sphereUid ), fileName_( fileName ), fileIO_(fileIO),
      diameter_( diameter ), velocity_(velocity), forceMag_( forceMag ),
      position_( real_t(0) )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t t\t posZ\t posZ/D\t fZ\t fZ/fSt\n";
            file.close();
         }
      }
   }

   void operator()()
   {
      real_t pos(real_t(0));
      real_t hydForce(real_t(0));

      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != ac_->getInvalidIdx())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            pos = ac_->getPosition(idx)[2];
            hydForce = ac_->getHydrodynamicForce(idx)[2];
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( pos, mpi::SUM );
         mpi::allReduceInplace( hydForce, mpi::SUM );
      }

      position_ = pos;
      hydForce_ = hydForce;

   }

   real_t getPosition() const
   {
      return position_;
   }

   void writeToFile( const uint_t timestep, real_t density1, real_t density2, real_t totalMass )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         auto scaledPosition = position_ / diameter_;
         auto normalizedHydForce = hydForce_ / std::abs(forceMag_);

         file << std::setprecision(10)
              << timestep << "\t" << real_c(timestep) / (diameter_ / velocity_) << "\t"
              << "\t" << position_ << "\t" << scaledPosition
              << "\t" << hydForce_ << "\t" << normalizedHydForce
              << "\t" << density1 << "\t" << density2 << "\t" << totalMass
              << "\n";
         file.close();
      }
   }

private:
   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;
   std::string fileName_;
   bool fileIO_;
   real_t diameter_, velocity_, forceMag_;

   real_t position_, hydForce_;
};


void createPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss, const math::AABB & simulationDomain)
{
   mesa_pd::data::Particle p2 = *ps->create(true);
   p2.setPosition(simulationDomain.minCorner());
   p2.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(1,0,0) ));
   p2.setOwner(mpi::MPIManager::instance()->rank());
   p2.setType(0);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p3 = *ps->create(true);
   p3.setPosition(simulationDomain.maxCorner());
   p3.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(-1,0,0) ));
   p3.setOwner(mpi::MPIManager::instance()->rank());
   p3.setType(0);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p4 = *ps->create(true);
   p4.setPosition(simulationDomain.minCorner());
   p4.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,1,0) ));
   p4.setOwner(mpi::MPIManager::instance()->rank());
   p4.setType(0);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p5 = *ps->create(true);
   p5.setPosition(simulationDomain.maxCorner());
   p5.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,-1,0) ));
   p5.setOwner(mpi::MPIManager::instance()->rank());
   p5.setType(0);
   mesa_pd::data::particle_flags::set(p5.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p5.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

}

template< typename DensityInterpolator_T>
real_t getDensityAtPosition(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID densityInterpolatorID, Vector3<real_t> position)
{

   real_t density = real_t(0);
   for( auto & block: *blocks)
   {
      if(block.getAABB().contains(position))
      {
         auto densityInterpolator = block.getData<DensityInterpolator_T>(densityInterpolatorID);
         densityInterpolator->get(position, &density);
      }

   }
   mpi::reduceInplace(density, mpi::SUM);

   return density;
}

template< typename BoundaryHandling_T>
real_t getAverageDensityInSystem(const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfFieldID, BlockDataID boundaryHandlingID)
{
   real_t totalMass = real_t(0);
   uint_t count = uint_t(0);
   for( auto & block: *blocks)
   {
      auto pdfField = block.getData<PdfField_T>(pdfFieldID);
      auto boundaryHandling = block.getData<BoundaryHandling_T>(boundaryHandlingID);
      WALBERLA_FOR_ALL_CELLS_XYZ(pdfField,
            if(boundaryHandling->isDomain(x,y,z))
            {
               totalMass += pdfField->getDensity(x,y,z);
               ++count;
            }
      );
   }

   mpi::reduceInplace(totalMass, mpi::SUM);
   mpi::reduceInplace(count, mpi::SUM);

   return totalMass/real_c(count);
}

//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Test of a sphere that is moving with a prescribed velocity.
 *
 * It is used to investigate fluctuations in the density field around the particle and how the bulk viscosity
 * or others (PDF reconstructor, resolution, velocity) influence this behavior.
 *
 * With the option artificiallyAccelerateSphere, the sphere is gradually accelerated until the final velocity is reached.
 * This gets rid of the strong oscillations that would appear in an abrupt start.
 *
 * Results of force plot:
 * - mn=0.125 (2/16) reduces oscillations compared to 0.1875 (3/16) a bit, especially the peaks
 * - ug=0.005 leads to higher oscillations than ug=0.01 -> since relaxation time is smaller
 * - bvrf=100 reduces oscillations, especially below the mean. Peaks remain almost the same
 * - Grad vs Ext reconstructor -> almost same behavior
 * - D=25 decreases oscillations -> since relaxation time is bigger
 *
 * NOTE: Run at least with 3 processes due to periodicity!
 */
//*******************************************************************************************************************

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////
   // Customization //
   ///////////////////

   // simulation control
   bool funcTest = false;
   bool fileIO = true;
   uint_t vtkIOFreq = 0;
   std::string baseFolder = "vtk_out_MovingWithPrescribedVelocity";
   std::string fileNameEnding = "";
   bool logDensity = true;
   bool vtkOutputAtEnd = true;

   //numerical parameters
   bool averageForceTorqueOverTwoTimeSteps = true;
   bool conserveMomentum = false;
   uint_t numRPDSubCycles = uint_t(1);
   std::string reconstructorType = "Grad"; // Eq, EAN, Ext, Grad
   real_t bulkViscRateFactor = real_t(1);
   real_t magicNumber = real_t(0.1875);
   bool useOmegaBulkAdaption = false;
   real_t adaptionLayerSize = real_t(2);
   uint_t blocksInX = uint_t(1);
   uint_t blocksInY = uint_t(1);
   uint_t blocksInZ = uint_t(4);
   real_t initialSpherePosition = real_t(0.5); // in ratio to domain height

   // simulation setup
   real_t velocity = real_t(0.02);
   real_t Re = real_t(164);
   real_t diameter = real_t(20);
   real_t numberOfPasses = real_t(3);
   real_t domainWidthNonDim = real_t(4);
   real_t domainHeightNonDim = real_t(8);
   bool usePeriodicSetup = false;
   bool artificiallyAccelerateSphere = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--funcTest" )             == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--noLogging" )            == 0 ) { fileIO = false; continue; }
      if( std::strcmp( argv[i], "--noVtkOutputAtEnd" )     == 0 ) { vtkOutputAtEnd = false; continue; }
      if( std::strcmp( argv[i], "--logDensity" )           == 0 ) { logDensity = true; continue; }
      if( std::strcmp( argv[i], "--vtkIOFreq" )            == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--numRPDSubCycles" )      == 0 ) { numRPDSubCycles = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging" )     == 0 ) { averageForceTorqueOverTwoTimeSteps = false; continue; }
      if( std::strcmp( argv[i], "--conserveMomentum" )     == 0 ) { conserveMomentum = true; continue; }
      if( std::strcmp( argv[i], "--baseFolder" )           == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--reconstructorType" )    == 0 ) { reconstructorType = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--bulkViscRateFactor" )   == 0 ) { bulkViscRateFactor = real_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--magicNumber" )          == 0 ) { magicNumber = real_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--velocity" )             == 0 ) { velocity = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--Re" )                   == 0 ) { Re = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--diameter" )             == 0 ) { diameter = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--domainWidth" )          == 0 ) { domainWidthNonDim = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--domainHeight" )         == 0 ) { domainHeightNonDim = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--fileName" )             == 0 ) { fileNameEnding = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--useOmegaBulkAdaption" ) == 0 ) { useOmegaBulkAdaption = true; continue; }
      if( std::strcmp( argv[i], "--adaptionLayerSize" )    == 0 ) { adaptionLayerSize = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--numberOfPasses" )       == 0 ) { numberOfPasses = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--initialSpherePosition" )== 0 ) { initialSpherePosition = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--usePeriodicSetup" )     == 0 ) { usePeriodicSetup = true; continue; }
      if( std::strcmp( argv[i], "--artificiallyAccelerateSphere" )     == 0 ) { artificiallyAccelerateSphere = true; continue; }
      if( std::strcmp( argv[i], "--blocksInX" )            == 0 ) { blocksInX = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--blocksInY" )            == 0 ) { blocksInY = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--blocksInZ" )            == 0 ) { blocksInZ = uint_c( std::atof( argv[++i] ) ); continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if( funcTest )
   {
      walberla::logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::WARNING);
   }

   if( fileIO )
   {
      // create base directory if it does not yet exist
      filesystem::path tpath( baseFolder );
      if( !filesystem::exists( tpath ) )
         filesystem::create_directory( tpath );
   }

   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   real_t viscosity = velocity * diameter / Re;
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   const real_t relaxationTime = real_t(1) / omega;
   const real_t omegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, bulkViscRateFactor);

   const Vector3<uint_t> domainSize( uint_c(domainWidthNonDim * diameter), uint_c(domainWidthNonDim * diameter), uint_c(domainHeightNonDim * diameter) );

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - Re = " << Re);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter << ", constant velocity = " << velocity );
   WALBERLA_LOG_INFO_ON_ROOT(" - relaxation time (tau) = " << relaxationTime << ", omega = " << omega << " kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - magic number " << magicNumber);
   WALBERLA_LOG_INFO_ON_ROOT(" - omegaBulk = " << omegaBulk << ", bulk visc. = " << lbm_mesapd_coupling::bulkViscosityFromOmegaBulk(omegaBulk) << " (bvrf " << bulkViscRateFactor << ")");
   WALBERLA_LOG_INFO_ON_ROOT(" - use omega bulk adaption = " << useOmegaBulkAdaption << " (adaption layer size = " << adaptionLayerSize << ")");
   WALBERLA_LOG_INFO_ON_ROOT(" - reconstructor type = " << reconstructorType );

   if( vtkIOFreq > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder << "\" with frequency " << vtkIOFreq);
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<uint_t> numberOfBlocksPerDirection( blocksInX, blocksInY, blocksInZ );
   Vector3<uint_t> cellsPerBlockPerDirection( domainSize[0] / numberOfBlocksPerDirection[0],
                                              domainSize[1] / numberOfBlocksPerDirection[1],
                                              domainSize[2] / numberOfBlocksPerDirection[2] );
   for( uint_t i = 0; i < 3; ++i ) {
      WALBERLA_CHECK_EQUAL(cellsPerBlockPerDirection[i] * numberOfBlocksPerDirection[i], domainSize[i],
                           "Unmatching domain decomposition in direction " << i << "!");
   }

   real_t dx = real_t(1);
   auto blocks = blockforest::createUniformBlockGrid( numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1], numberOfBlocksPerDirection[2],
                                                      cellsPerBlockPerDirection[0], cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], dx,
                                                      0, false, false,
                                                      usePeriodicSetup, usePeriodicSetup, true, //periodicity
                                                      false );

   WALBERLA_LOG_INFO_ON_ROOT("Domain decomposition:");
   WALBERLA_LOG_INFO_ON_ROOT(" - blocks per direction = " << numberOfBlocksPerDirection );
   WALBERLA_LOG_INFO_ON_ROOT(" - cells per block = " << cellsPerBlockPerDirection );

   //write domain decomposition to file
   if( vtkIOFreq > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolder );
   }

   //////////////////
   // RPD COUPLING //
   //////////////////

   auto rpdDomain = std::make_shared<mesa_pd::domain::BlockForestDomain>(blocks->getBlockForestPointer());

   //init data structures
   auto ps = walberla::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = walberla::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = walberla::make_shared<ParticleAccessor_T >(ps, ss);

   // bounding planes
   if(!usePeriodicSetup) createPlaneSetup(ps,ss,blocks->getDomain());

   // create sphere and store Uid
   Vector3<real_t> initialPosition( real_t(0.5) * real_c(domainSize[0]), real_t(0.5) * real_c(domainSize[1]), initialSpherePosition * real_c(domainSize[2]));
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( diameter * real_t(0.5) );
   //ss->shapes[sphereShape]->updateMassAndInertia(densityRatio);

   walberla::id_t sphereUid = 0;
   if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), initialPosition ))
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(initialPosition);
      p.setInteractionRadius(diameter * real_t(0.5));
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);
      p.setLinearVelocity(Vector3<real_t>(0.0,0.0,velocity));
      sphereUid = p.getUid();
   }
   mpi::allReduceInplace(sphereUid, mpi::SUM);

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add omega bulk field
   BlockDataID omegaBulkFieldID = field::addToStorage<ScalarField_T>( blocks, "omega bulk field", omegaBulk, field::fzyx, uint_t(0) );

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
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         uint_t(1), field::fzyx );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::fzyx, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );

   // interpolation functionality
   using DensityAdaptor_T = typename lbm::Adaptor< LatticeModel_T >::Density;
   BlockDataID densityAdaptorID = field::addFieldAdaptor<  DensityAdaptor_T >( blocks, pdfFieldID, "density adaptor" );

   using DensityInterpolator_T = typename field::TrilinearFieldInterpolator<DensityAdaptor_T, FlagField_T>;
   BlockDataID densityInterpolatorID = field::addFieldInterpolator< DensityInterpolator_T, FlagField_T >( blocks, densityAdaptorID, flagFieldID, Fluid_Flag );

   // set up RPD functionality
   std::function<void(void)> syncCall = [ps,rpdDomain, adaptionLayerSize ](){
      const real_t overlap = real_t( 1.5 ) + std::max(real_t(0), adaptionLayerSize - real_t(2));
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();

   mesa_pd::kernel::ExplicitEuler explEulerIntegrator(real_t(1)/real_t(numRPDSubCycles));

   mesa_pd::mpi::ReduceProperty reduceProperty;

   // set up coupling functionality
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;

   lbm_mesapd_coupling::RegularParticlesSelector sphereSelector;

   lbm_mesapd_coupling::ParticleMappingKernel<BoundaryHandling_T> particleMappingKernel(blocks, boundaryHandlingID);
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);

   ///////////////
   // TIME LOOP //
   ///////////////

   // map planes into the LBM simulation -> act as no-slip boundaries
   ps->forEachParticle(false, lbm_mesapd_coupling::GlobalParticlesSelector(), *accessor, particleMappingKernel, *accessor, NoSlip_Flag);

   // map particles into the LBM simulation
   ps->forEachParticle(false, sphereSelector, *accessor, movingParticleMappingKernel, *accessor, MO_Flag);


   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks

   blockforest::communication::UniformBufferedScheme< Stencil_T > optimizedPDFCommunicationScheme( blocks );
   optimizedPDFCommunicationScheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) ); // optimized sync

   blockforest::communication::UniformBufferedScheme< Stencil_T > fullPDFCommunicationScheme( blocks );
   fullPDFCommunicationScheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) ); // full sync

   // create the timeloop

   const uint_t timesteps = uint_c(numberOfPasses * real_c(domainSize[2]) / velocity);
   WALBERLA_LOG_INFO_ON_ROOT("Running for " << timesteps << " time steps");
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );


   timeloop.addFuncBeforeTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   // vtk output
   if( vtkIOFreq != uint_t(0) )
   {
      // spheres
      auto particleVtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleOwner>("owner");
      particleVtkOutput->addOutput<mesa_pd::data::SelectParticleLinearVelocity>("velocity");
      auto particleVtkWriter = vtk::createVTKOutput_PointData(particleVtkOutput, "Particles", vtkIOFreq, baseFolder, "simulation_step");
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( particleVtkWriter ), "VTK (sphere data)" );

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

   // sweep for updating the particle mapping into the LBM simulation
   timeloop.add() << Sweep( lbm_mesapd_coupling::makeMovingParticleMapping<PdfField_T, BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, MO_Flag, FormerMO_Flag, sphereSelector, conserveMomentum), "Particle Mapping" );

   // sweep for restoring PDFs in cells previously occupied by particles
   if( reconstructorType == "EAN")
   {

      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto equilibriumAndNonEquilibriumSphereNormalReconstructor = lbm_mesapd_coupling::makeEquilibriumAndNonEquilibriumReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);
      auto reconstructionManager = lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, equilibriumAndNonEquilibriumSphereNormalReconstructor, conserveMomentum);

      timeloop.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                     << Sweep( makeSharedSweep(reconstructionManager), "PDF Restore" );

   } else if( reconstructorType == "Grad")
   {
      bool recomputeDensity = false;
      auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, recomputeDensity, true, true);

      timeloop.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                     << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, gradReconstructor, conserveMomentum) ), "PDF Restore" );
   } else if( reconstructorType == "Eq")
   {
      timeloop.add() << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, conserveMomentum) ), "PDF Restore" );
   } else if( reconstructorType == "Ext")
   {
      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto extrapolationReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T, lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder, true>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);
      timeloop.add() << BeforeFunction( fullPDFCommunicationScheme, "LBM Communication" )
                     << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, extrapolationReconstructor, conserveMomentum)), "PDF Restore" );
   } else {
      WALBERLA_ABORT("Unknown reconstructor type " << reconstructorType);
   }

   // update bulk omega in all cells to adapt to changed particle position
   if(useOmegaBulkAdaption)
   {
      using OmegaBulkAdapter_T = lbm_mesapd_coupling::OmegaBulkAdapter<ParticleAccessor_T, decltype(sphereSelector)>;
      real_t defaultOmegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, real_t(1));
      shared_ptr<OmegaBulkAdapter_T> omegaBulkAdapter = make_shared<OmegaBulkAdapter_T>(blocks, omegaBulkFieldID, accessor, defaultOmegaBulk, omegaBulk, adaptionLayerSize, sphereSelector);
      timeloop.add() << Sweep( makeSharedSweep(omegaBulkAdapter), "Omega Bulk Adapter");
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   auto bhSweep = BoundaryHandling_T::getBlockSweep( boundaryHandlingID );
   timeloop.add() << BeforeFunction( optimizedPDFCommunicationScheme, "LBM Communication" )
                  << Sweep(bhSweep, "Boundary Handling" );

   // stream + collide LBM step
#ifdef WALBERLA_BUILD_WITH_CODEGEN
   auto lbmSweep = LatticeModel_T::Sweep( pdfFieldID );
   timeloop.add() << Sweep( lbmSweep, "LB sweep" );
#else
   auto lbmSweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   timeloop.add() << Sweep( makeSharedSweep( lbmSweep ), "cell-wise LB sweep" );
#endif

   // evaluation functionality
   std::string loggingFileName( baseFolder + "/LoggingPrescribedMovingSphere");
   loggingFileName += "_Re" + std::to_string(uint_c(Re));
   loggingFileName += "_D" + std::to_string(uint_c(diameter));
   loggingFileName += "_vel" + std::to_string(float(velocity));
   loggingFileName += "_recon" + reconstructorType;
   loggingFileName += "_bvrf" + std::to_string(uint_c(bulkViscRateFactor));
   loggingFileName += "_mn" + std::to_string(float(magicNumber));
   if( useOmegaBulkAdaption ) loggingFileName += "_uOBA" + std::to_string(uint_c(adaptionLayerSize));
   if( conserveMomentum ) loggingFileName += "_conserveMomentum";
   if( artificiallyAccelerateSphere ) loggingFileName += "_acc";
   if( !fileNameEnding.empty()) loggingFileName += "_" + fileNameEnding;
   loggingFileName += ".txt";

   if( fileIO  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   SpherePropertyLogger<ParticleAccessor_T > logger( accessor, sphereUid, loggingFileName, fileIO, diameter, velocity,  real_t(3) * math::pi * diameter * viscosity * velocity);


   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   const bool useOpenMP = false;

   const real_t densityRatio = real_t(7800) / real_t(935);
   const real_t responseTime = densityRatio * diameter * diameter / ( real_t(18) * viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - response time " << responseTime);
   const real_t accelerationFactor = real_t(1) / (real_t(0.1) * responseTime);

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {
      // perform a single simulation step -> this contains LBM and setting of the hydrodynamic interactions
      timeloop.singleStep( timeloopTiming );

      timeloopTiming["RPD"].start();
      
      reduceProperty.operator()<mesa_pd::HydrodynamicForceTorqueNotification>(*ps);

      if( averageForceTorqueOverTwoTimeSteps )
      {
         if( i == 0 )
         {
            lbm_mesapd_coupling::InitializeHydrodynamicForceTorqueForAveragingKernel initializeHydrodynamicForceTorqueForAveragingKernel;
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, initializeHydrodynamicForceTorqueForAveragingKernel, *accessor );
         }
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, averageHydrodynamicForceTorque, *accessor );
      }

      reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);
         
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, explEulerIntegrator, *accessor);
      syncCall();

      timeloopTiming["RPD"].end();

      if( artificiallyAccelerateSphere )
      {
         // accelerating after reading from a checkpoint file does not make sense, as current actual runtime is not known
         real_t newSphereVel = velocity * (std::exp(-accelerationFactor * real_t(i) ) - real_t(1));
         ps->forEachParticle(useOpenMP, sphereSelector, *accessor, [newSphereVel](const size_t idx, ParticleAccessor_T& ac){ ac.setLinearVelocity(idx, Vector3<real_t>(real_t(0), real_t(0), newSphereVel));}, *accessor);
      }

      // evaluation
      timeloopTiming["Logging"].start();

      logger();

      real_t density1 = real_t(0);
      real_t density2 = real_t(0);
      real_t averageDensity = real_t(0);
      if(logDensity)
      {
         density1 = getDensityAtPosition<DensityInterpolator_T >(blocks, densityInterpolatorID, Vector3<real_t>(real_c(domainSize[0])*real_t(0.25), real_c(domainSize[1])*real_t(0.5), real_c(domainSize[2])*real_t(0.5)));
         density2 = getDensityAtPosition<DensityInterpolator_T >(blocks, densityInterpolatorID, Vector3<real_t>(real_c(domainSize[0])*real_t(0.75), real_c(domainSize[1])*real_t(0.5), logger.getPosition()));
         averageDensity = getAverageDensityInSystem<BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID);
      }
      logger.writeToFile(i, density1, density2, averageDensity);
      timeloopTiming["Logging"].end();

      // reset after logging here
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

   }

   timeloopTiming.logResultOnRoot();

   if(vtkOutputAtEnd)
   {
      std::string vtkFileName = "fluid_field";
      vtkFileName += "_bvrf" + std::to_string(uint_c(bulkViscRateFactor));
      if( useOmegaBulkAdaption ) vtkFileName += "_uOBA" + std::to_string(uint_c(adaptionLayerSize));

      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, vtkFileName, uint_t(1), 0, false, baseFolder );

      pdfFieldVTK->addBeforeFunction( fullPDFCommunicationScheme );

      AABB sliceAABB( real_t(0), real_c(domainSize[1])*real_t(0.5)-real_t(1), real_t(0),
                      real_c(domainSize[0]), real_c(domainSize[1])*real_t(0.5)+real_t(1), real_c(domainSize[2]) );
      vtk::AABBCellFilter aabbSliceFilter( sliceAABB );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );

      vtk::ChainedFilter combinedSliceFilter;
      combinedSliceFilter.addFilter( fluidFilter );
      combinedSliceFilter.addFilter( aabbSliceFilter );

      pdfFieldVTK->addCellInclusionFilter( combinedSliceFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      vtk::writeFiles(pdfFieldVTK)();
   }


   return EXIT_SUCCESS;
}

} // namespace sphere_moving_with_prescribed_velocity

int main( int argc, char **argv ){
   sphere_moving_with_prescribed_velocity::main(argc, argv);
}
