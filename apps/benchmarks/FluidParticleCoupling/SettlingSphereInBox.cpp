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
//! \file SettlingSphereInBox.cpp
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
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
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
#include "mesa_pd/kernel/ExplicitEuler.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
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

namespace settling_sphere_in_box
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
         real_t dx_SI, real_t dt_SI, real_t diameter, real_t gravitationalForceMag, real_t uref) :
      ac_( ac ), sphereUid_( sphereUid ), fileName_( fileName ), fileIO_(fileIO),
      dx_SI_( dx_SI ), dt_SI_( dt_SI ), diameter_( diameter ),
      gravitationalForceMag_( gravitationalForceMag ), uref_(uref), tref_(diameter / uref),
      position_( real_t(0) ), maxVelocity_( real_t(0) )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t t\t t/tref\t posZ\t gapZ/D\t velZ\t velZ_SI\t velZ/uref\t fZ\t fZ/fGravi\n";
            file.close();
         }
      }
   }

   void operator()(const uint_t timestep)
   {
      Vector3<real_t> pos(real_t(0));
      Vector3<real_t> transVel(real_t(0));
      Vector3<real_t> hydForce(real_t(0));

      size_t idx = ac_->uidToIdx(sphereUid_);
      if( idx != ac_->getInvalidIdx())
      {
         if(!mesa_pd::data::particle_flags::isSet( ac_->getFlags(idx), mesa_pd::data::particle_flags::GHOST))
         {
            pos = ac_->getPosition(idx);
            transVel = ac_->getLinearVelocity(idx);
            hydForce = ac_->getHydrodynamicForce(idx);
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( pos[0], mpi::SUM );
         mpi::allReduceInplace( pos[1], mpi::SUM );
         mpi::allReduceInplace( pos[2], mpi::SUM );

         mpi::allReduceInplace( transVel[0], mpi::SUM );
         mpi::allReduceInplace( transVel[1], mpi::SUM );
         mpi::allReduceInplace( transVel[2], mpi::SUM );

         mpi::allReduceInplace( hydForce[0], mpi::SUM );
         mpi::allReduceInplace( hydForce[1], mpi::SUM );
         mpi::allReduceInplace( hydForce[2], mpi::SUM );
      }

      position_ = pos[2];
      maxVelocity_ = std::max(maxVelocity_, -transVel[2]);

      if( fileIO_ )
         writeToFile( timestep, pos, transVel, hydForce);
   }

   real_t getPosition() const
   {
      return position_;
   }

   real_t getMaxVelocity() const
   {
      return maxVelocity_;
   }

private:
   void writeToFile( const uint_t timestep, const Vector3<real_t> & position, const Vector3<real_t> & velocity, const Vector3<real_t> & hydForce )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         auto scaledPosition = position / diameter_;
         auto velocity_SI = velocity * dx_SI_ / dt_SI_;
         auto normalizedHydForce = hydForce / gravitationalForceMag_;

         file << timestep << "\t" << real_c(timestep) * dt_SI_ << "\t" << real_c(timestep) / tref_
              << "\t" << position[2] << "\t" << scaledPosition[2] - real_t(0.5)
              << "\t" << velocity[2] << "\t" << velocity_SI[2] << "\t" << velocity[2] / uref_
              << "\t" << hydForce[2] << "\t" << normalizedHydForce[2]
              << "\n";
         file.close();
      }
   }

   shared_ptr< ParticleAccessor_T > ac_;
   const walberla::id_t sphereUid_;
   std::string fileName_;
   bool fileIO_;
   real_t dx_SI_, dt_SI_, diameter_, gravitationalForceMag_, uref_, tref_;

   real_t position_;
   real_t maxVelocity_;
};


void createPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss, const math::AABB & simulationDomain)
{
   // create bounding planes
   mesa_pd::data::Particle p0 = *ps->create(true);
   p0.setPosition(simulationDomain.minCorner());
   p0.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p0.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,1) ));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p1 = *ps->create(true);
   p1.setPosition(simulationDomain.maxCorner());
   p1.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p1.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,-1) ));
   p1.setOwner(mpi::MPIManager::instance()->rank());
   p1.setType(0);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p2 = *ps->create(true);
   p2.setPosition(simulationDomain.minCorner());
   p2.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p2.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(1,0,0) ));
   p2.setOwner(mpi::MPIManager::instance()->rank());
   p2.setType(0);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p3 = *ps->create(true);
   p3.setPosition(simulationDomain.maxCorner());
   p3.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p3.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(-1,0,0) ));
   p3.setOwner(mpi::MPIManager::instance()->rank());
   p3.setType(0);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p4 = *ps->create(true);
   p4.setPosition(simulationDomain.minCorner());
   p4.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p4.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,1,0) ));
   p4.setOwner(mpi::MPIManager::instance()->rank());
   p4.setType(0);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p5 = *ps->create(true);
   p5.setPosition(simulationDomain.maxCorner());
   p5.setInteractionRadius(std::numeric_limits<real_t>::infinity());
   p5.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,-1,0) ));
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

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////
   // Customization //
   ///////////////////

   // simulation control
   bool shortrun = false;
   bool funcTest = false;
   bool fileIO = true;
   uint_t vtkIOFreq = 0;
   std::string baseFolder = "vtk_out_SettlingSphere";
   std::string fileNameEnding = "";

   // physical setup
   uint_t fluidType = 1;

   //numerical parameters
   uint_t numberOfCellsInHorizontalDirection = uint_t(135);
   bool averageForceTorqueOverTwoTimeSteps = true;
   bool conserveMomentum = false;
   uint_t numRPDSubCycles = uint_t(1);
   bool useVelocityVerlet = true;
   std::string reconstructorType = "Grad"; // Eq, EAN, Ext, Grad
   real_t bulkViscRateFactor = real_t(1);
   real_t magicNumber = real_t(3)/real_t(16);
   real_t characteristicVelocity = real_t(0.02);
   bool useOmegaBulkAdaption = false;
   real_t adaptionLayerSize = real_t(2);
   bool useLubricationCorrection = true;

   bool useGalileoParameterization = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" )             == 0 ) { shortrun = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" )             == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--noLogging" )            == 0 ) { fileIO = false; continue; }
      if( std::strcmp( argv[i], "--vtkIOFreq" )            == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--fluidType" )            == 0 ) { fluidType = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--numRPDSubCycles" )      == 0 ) { numRPDSubCycles = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--resolution" )           == 0 ) { numberOfCellsInHorizontalDirection = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging" )     == 0 ) { averageForceTorqueOverTwoTimeSteps = false; continue; }
      if( std::strcmp( argv[i], "--conserveMomentum" )     == 0 ) { conserveMomentum = true; continue; }
      if( std::strcmp( argv[i], "--baseFolder" )           == 0 ) { baseFolder = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--useEulerIntegrator" )   == 0 ) { useVelocityVerlet = false; continue; }
      if( std::strcmp( argv[i], "--reconstructorType" )    == 0 ) { reconstructorType = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--bulkViscRateFactor" )   == 0 ) { bulkViscRateFactor = real_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--magicNumber" )          == 0 ) { magicNumber = real_c(std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--velocity" )             == 0 ) { characteristicVelocity = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--fileName" )             == 0 ) { fileNameEnding = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--useOmegaBulkAdaption" ) == 0 ) { useOmegaBulkAdaption = true; continue; }
      if( std::strcmp( argv[i], "--adaptionLayerSize" )    == 0 ) { adaptionLayerSize = real_c(std::atof( argv[++i] )); continue; }
      if( std::strcmp( argv[i], "--noLubricationCorrection" ) == 0 ) { useLubricationCorrection = false; continue; }
      if( std::strcmp( argv[i], "--useGalileoParameterization" ) == 0 ) { useGalileoParameterization = true; continue; }
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

   //////////////////////////////////////
   // SIMULATION PROPERTIES in SI units//
   //////////////////////////////////////

   // values are mainly taken from the reference paper
   const real_t diameter_SI = real_t(15e-3);
   const real_t densitySphere_SI = real_t(1120);

   real_t densityFluid_SI;
   real_t dynamicViscosityFluid_SI;
   real_t expectedSettlingVelocity_SI;

   // expected velocity given as uMax in experiments of ten Cate (Table 2, E1-E4), with uInfty from Table 1
   // > slightly different Re than in ten Cate's Table 1
   switch( fluidType )
   {
      case 1:
         // Re_p around 1.5
         densityFluid_SI = real_t(970);
         dynamicViscosityFluid_SI = real_t(373e-3);
         expectedSettlingVelocity_SI = real_t(0.035986);
         break;
      case 2:
         // Re_p around 4.1
         densityFluid_SI = real_t(965);
         dynamicViscosityFluid_SI = real_t(212e-3);
         expectedSettlingVelocity_SI = real_t(0.05718);
         break;
      case 3:
         // Re_p around 11.6
         densityFluid_SI = real_t(962);
         dynamicViscosityFluid_SI = real_t(113e-3);
         expectedSettlingVelocity_SI = real_t(0.087269);
         break;
      case 4:
         // Re_p around 31.9
         densityFluid_SI = real_t(960);
         dynamicViscosityFluid_SI = real_t(58e-3);
         expectedSettlingVelocity_SI = real_t(0.12224);
         break;
      default:
         WALBERLA_ABORT("Only four different fluids are supported! Choose type between 1 and 4.");
   }
   const real_t kinematicViscosityFluid_SI = dynamicViscosityFluid_SI / densityFluid_SI;

   const real_t gravitationalAcceleration_SI = real_t(9.81);

   const real_t ug_SI = std::sqrt((densitySphere_SI/densityFluid_SI-real_t(1))*diameter_SI*gravitationalAcceleration_SI);
   const real_t GalileoNumber = ug_SI * diameter_SI / (dynamicViscosityFluid_SI/densityFluid_SI);

   Vector3<real_t> domainSize_SI(real_t(100e-3), real_t(100e-3), real_t(160e-3));
   //shift starting gap a bit upwards to match the reported (plotted) values
   const real_t startingGapSize_SI = real_t(120e-3) + real_t(0.25) * diameter_SI;

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in SI units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid type = " << fluidType );
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter_SI << ", density = " << densitySphere_SI << ", starting gap size = " << startingGapSize_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid_SI << ", dyn. visc = " << dynamicViscosityFluid_SI << ", kin. visc = " << kinematicViscosityFluid_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedSettlingVelocity_SI << " --> Re_p = " << expectedSettlingVelocity_SI * diameter_SI / kinematicViscosityFluid_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational velocity = " << ug_SI << " --> Ga = " << GalileoNumber );


   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////


   const real_t dx_SI = domainSize_SI[0] / real_c(numberOfCellsInHorizontalDirection);
   const Vector3<uint_t> domainSize( uint_c(floor(domainSize_SI[0] / dx_SI + real_t(0.5)) ),
                                     uint_c(floor(domainSize_SI[1] / dx_SI + real_t(0.5)) ),
                                     uint_c(floor(domainSize_SI[2] / dx_SI + real_t(0.5)) ) );
   const real_t diameter = diameter_SI / dx_SI;
   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;


   const real_t dt_SI = (useGalileoParameterization) ? characteristicVelocity / ug_SI * dx_SI : // this uses Ga for parameterization, where ug is the characteristic velocity
                                                       characteristicVelocity / expectedSettlingVelocity_SI * dx_SI; // this uses Re for parameterization (only possible since settling velocity is known)

   const real_t viscosity =  kinematicViscosityFluid_SI * dt_SI / ( dx_SI * dx_SI );
   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   const real_t relaxationTime = real_t(1) / omega;
   const real_t omegaBulk = lbm_mesapd_coupling::omegaBulkFromOmega(omega, bulkViscRateFactor);

   const real_t gravitationalAcceleration = gravitationalAcceleration_SI * dt_SI * dt_SI / dx_SI;

   const real_t densityFluid = real_t(1);
   const real_t densitySphere = densityFluid * densitySphere_SI / densityFluid_SI;

   const real_t expectedSettlingVelocity = expectedSettlingVelocity_SI * dt_SI / dx_SI;
   const real_t ug = std::sqrt((densitySphere/densityFluid-real_t(1))*diameter*gravitationalAcceleration);
   const real_t GalileoNumber_LBM = ug * diameter / viscosity;

   const real_t dx = real_t(1);

   const uint_t timesteps = funcTest ? 1 : ( shortrun ? uint_t(200) : uint_t( 250000 ) );

   WALBERLA_LOG_INFO_ON_ROOT(" - dx_SI = " << dx_SI << ", dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter << ", density = " << densitySphere );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid << ", relaxation time (tau) = " << relaxationTime << ", omega = " << omega << " kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - magic number = " << magicNumber);
   WALBERLA_LOG_INFO_ON_ROOT(" - omegaBulk = " << omegaBulk << ", bulk visc. = " << lbm_mesapd_coupling::bulkViscosityFromOmegaBulk(omegaBulk) << " (bvrf " << bulkViscRateFactor << ")");
   WALBERLA_LOG_INFO_ON_ROOT(" - use omega bulk adaption = " << useOmegaBulkAdaption << " (adaption layer size = " << adaptionLayerSize << ")");
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedSettlingVelocity << " --> Re_p = " << expectedSettlingVelocity * diameter / viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational velocity = " << ug << " --> Ga = " << GalileoNumber_LBM );
   WALBERLA_LOG_INFO_ON_ROOT(" - integrator = " << (useVelocityVerlet ? "Velocity Verlet" : "Explicit Euler") );
   WALBERLA_LOG_INFO_ON_ROOT(" - conserve momentum = " << (conserveMomentum ? "yes" : "no") );
   WALBERLA_LOG_INFO_ON_ROOT(" - reconstructor type = " << reconstructorType );
   WALBERLA_LOG_INFO_ON_ROOT(" - lubrication correction = " << useLubricationCorrection );

   if( vtkIOFreq > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder << "\" with frequency " << vtkIOFreq);
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<uint_t> numberOfBlocksPerDirection( uint_t(1), uint_t(1), uint_t(4) );
   Vector3<uint_t> cellsPerBlockPerDirection( domainSize[0] / numberOfBlocksPerDirection[0],
                                              domainSize[1] / numberOfBlocksPerDirection[1],
                                              domainSize[2] / numberOfBlocksPerDirection[2] );
   for( uint_t i = 0; i < 3; ++i ) {
      WALBERLA_CHECK_EQUAL(cellsPerBlockPerDirection[i] * numberOfBlocksPerDirection[i], domainSize[i],
                           "Unmatching domain decomposition in direction " << i << "!");
   }

   auto blocks = blockforest::createUniformBlockGrid( numberOfBlocksPerDirection[0], numberOfBlocksPerDirection[1], numberOfBlocksPerDirection[2],
                                                      cellsPerBlockPerDirection[0], cellsPerBlockPerDirection[1], cellsPerBlockPerDirection[2], dx,
                                                      0, false, false,
                                                      false, false, false, //periodicity
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
   createPlaneSetup(ps,ss,blocks->getDomain());

   // create sphere and store Uid
   Vector3<real_t> initialPosition( real_t(0.5) * real_c(domainSize[0]), real_t(0.5) * real_c(domainSize[1]), startingGapSize_SI / dx_SI + real_t(0.5) * diameter);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( diameter * real_t(0.5) );
   ss->shapes[sphereShape]->updateMassAndInertia(densitySphere);

   walberla::id_t sphereUid = 0;
   if (rpdDomain->isContainedInProcessSubdomain( uint_c(mpi::MPIManager::instance()->rank()), initialPosition ))
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

   // set up RPD functionality
   std::function<void(void)> syncCall = [ps,rpdDomain](){
      const real_t overlap = real_t( 1.5 );
      mesa_pd::mpi::SyncNextNeighbors syncNextNeighborFunc;
      syncNextNeighborFunc(*ps, *rpdDomain, overlap);
   };

   syncCall();

   mesa_pd::kernel::ExplicitEuler explEulerIntegrator(real_t(1)/real_t(numRPDSubCycles));
   mesa_pd::kernel::VelocityVerletPreForceUpdate  vvIntegratorPreForce(real_t(1)/real_t(numRPDSubCycles));
   mesa_pd::kernel::VelocityVerletPostForceUpdate vvIntegratorPostForce(real_t(1)/real_t(numRPDSubCycles));

   mesa_pd::mpi::ReduceProperty reduceProperty;

   // set up coupling functionality
   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densitySphere - densityFluid) * gravitationalAcceleration * sphereVolume );
   lbm_mesapd_coupling::AddForceOnParticlesKernel addGravitationalForce(gravitationalForce);
   lbm_mesapd_coupling::AddHydrodynamicInteractionKernel addHydrodynamicInteraction;
   lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
   lbm_mesapd_coupling::AverageHydrodynamicForceTorqueKernel averageHydrodynamicForceTorque;
   lbm_mesapd_coupling::LubricationCorrectionKernel lubricationCorrectionKernel(viscosity, [](real_t r){ return real_t(0.0016) * r;});

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

      // flag field (written only once in the first time step, ghost layers are also written)
      //auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, FieldGhostLayers, false, baseFolder );
      //flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      //timeloop.addFuncBeforeTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

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

   } else if( reconstructorType == "Ext" )
   {
      auto sphereNormalExtrapolationDirectionFinder = make_shared<lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder>(blocks);
      auto extrapolationSphereNormalReconstructor = lbm_mesapd_coupling::makeExtrapolationReconstructor<BoundaryHandling_T,lbm_mesapd_coupling::SphereNormalExtrapolationDirectionFinder,true>(blocks, boundaryHandlingID, sphereNormalExtrapolationDirectionFinder, uint_t(3), true);

      timeloop.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                     << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, extrapolationSphereNormalReconstructor, conserveMomentum) ), "PDF Restore" );
   } else if( reconstructorType == "Grad")
   {
      auto gradReconstructor = lbm_mesapd_coupling::makeGradsMomentApproximationReconstructor<BoundaryHandling_T>(blocks, boundaryHandlingID, omega, false, true, true);

      timeloop.add() << BeforeFunction( fullPDFCommunicationScheme, "PDF Communication" )
                     << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, gradReconstructor, conserveMomentum) ), "PDF Restore" );
   } else if( reconstructorType == "Eq")
   {
      timeloop.add() << Sweep( makeSharedSweep(lbm_mesapd_coupling::makePdfReconstructionManager<PdfField_T,BoundaryHandling_T>(blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, FormerMO_Flag, Fluid_Flag, conserveMomentum) ), "PDF Restore" );
   } else {
      WALBERLA_ABORT("Unknown reconstructor type " << reconstructorType);
   }

   // update bulk omega in all cells to adapt to changed particle position
   if( useOmegaBulkAdaption )
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
   std::string loggingFileName( baseFolder + "/LoggingSettlingSphere_");
   loggingFileName += std::to_string(fluidType);
   loggingFileName += "_res" + std::to_string(numberOfCellsInHorizontalDirection);
   loggingFileName += "_recon" + reconstructorType;
   loggingFileName += "_bvrf" + std::to_string(uint_c(bulkViscRateFactor));
   loggingFileName += "_mn" + std::to_string(float(magicNumber));
   if( useOmegaBulkAdaption ) loggingFileName += "_uOBA" + std::to_string(uint_c(adaptionLayerSize));
   if( useGalileoParameterization ) loggingFileName += "_Ga";
   if( !fileNameEnding.empty()) loggingFileName += "_" + fileNameEnding;
   loggingFileName += ".txt";
   if( fileIO  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   SpherePropertyLogger<ParticleAccessor_T> logger( accessor, sphereUid, loggingFileName, fileIO, dx_SI, dt_SI, diameter, -gravitationalForce[2], characteristicVelocity );


   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   real_t terminationPosition = diameter * real_t(0.501); // right before sphere touches the bottom wall
   real_t oldPos = initialPosition[2];

   const bool useOpenMP = false;

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

      for(auto subCycle = uint_t(0); subCycle < numRPDSubCycles; ++subCycle )
      {

         if( useVelocityVerlet )
         {
            ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPreForce, *accessor);
            syncCall();
         }

         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addHydrodynamicInteraction, *accessor );
         ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, addGravitationalForce, *accessor );

         // lubrication correction
         if(useLubricationCorrection)
         {
            ps->forEachParticlePairHalf(useOpenMP, mesa_pd::kernel::ExcludeInfiniteInfinite(), *accessor,
                                        [&lubricationCorrectionKernel,rpdDomain](const size_t idx1, const size_t idx2, auto& ac)
                                        {
                                           //TODO change this to storing copy, not reference
                                           mesa_pd::collision_detection::AnalyticContactDetection acd;
                                           acd.getContactThreshold() = lubricationCorrectionKernel.getNormalCutOffDistance();
                                           mesa_pd::kernel::DoubleCast double_cast;
                                           mesa_pd::mpi::ContactFilter contact_filter;
                                           if (double_cast(idx1, idx2, ac, acd, ac ))
                                           {
                                              if (contact_filter(acd.getIdx1(), acd.getIdx2(), ac, acd.getContactPoint(), *rpdDomain))
                                              {
                                                 double_cast(idx1, idx2, ac, lubricationCorrectionKernel, ac, acd.getContactNormal(), acd.getPenetrationDepth());
                                              }
                                           }
                                        },
                                        *accessor );
         }

         reduceProperty.operator()<mesa_pd::ForceTorqueNotification>(*ps);

         if( useVelocityVerlet ) ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, vvIntegratorPostForce, *accessor);
         else ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectLocal(), *accessor, explEulerIntegrator, *accessor);

         syncCall();
      }

      timeloopTiming["RPD"].end();

      // evaluation
      timeloopTiming["Logging"].start();
      logger(i);
      timeloopTiming["Logging"].end();

      // reset after logging here
      ps->forEachParticle(useOpenMP, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

      // check for termination
      if ( logger.getPosition() < terminationPosition || oldPos < logger.getPosition() )
      {
         WALBERLA_LOG_INFO_ON_ROOT("Sphere reached terminal position " << logger.getPosition() << " after " << i << " timesteps!");
         break;
      }
      oldPos = logger.getPosition();
   }

   timeloopTiming.logResultOnRoot();

   // check the result
   if ( !funcTest && !shortrun )
   {
      real_t relErr = std::fabs( expectedSettlingVelocity - logger.getMaxVelocity()) / expectedSettlingVelocity;
      WALBERLA_LOG_INFO_ON_ROOT( "Expected maximum settling velocity: " << expectedSettlingVelocity );
      WALBERLA_LOG_INFO_ON_ROOT( "Simulated maximum settling velocity: " << logger.getMaxVelocity() );
      WALBERLA_LOG_INFO_ON_ROOT( "Relative error: " << relErr );

      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS( relErr, real_t(0.1) );
   }

   return EXIT_SUCCESS;
}

} // namespace settling_sphere_in_box

int main( int argc, char **argv ){
   settling_sphere_in_box::main(argc, argv);
}
