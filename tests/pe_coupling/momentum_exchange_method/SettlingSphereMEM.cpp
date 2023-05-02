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
//! \file SettlingSphereMEM.cpp
//! \ingroup pe_coupling
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

#include "pe/basic.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"
#include "pe/cr/ICR.h"
#include "pe/Types.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>

namespace settling_sphere_mem
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

//////////////
// TYPEDEFS //
//////////////

// PDF field, flag field & body field
typedef lbm::D3Q19< lbm::collision_model::TRT >  LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

const uint_t FieldGhostLayers = 1;

// boundary handling
typedef lbm::NoSlip< LatticeModel_T, flag_t > NoSlip_T;

typedef pe_coupling::CurvedLinear< LatticeModel_T, FlagField_T > MO_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, MO_T > BoundaryHandling_T;

typedef std::tuple< pe::Sphere, pe::Plane > BodyTypeTuple;

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
class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_ ( bodyFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;

}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );
   BodyField_T * bodyField       = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "moving obstacle boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   // boundary conditions (no-slip) are set by mapping the planes into the domain

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}
//*******************************************************************************************************************


//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the sphere
 *
 */
//*******************************************************************************************************************
class SpherePropertyLogger
{
public:
   SpherePropertyLogger( SweepTimeloop* timeloop, const shared_ptr< StructuredBlockStorage > & blocks,
                         const BlockDataID & pdfFieldID, const BlockDataID & boundaryHandlingID,
                         const BlockDataID & bodyStorageID, const std::string & fileName, bool fileIO,
                         real_t dx_SI, real_t dt_SI, real_t diameter, real_t mass) :
      timeloop_( timeloop ), blocks_( blocks ), pdfFieldID_( pdfFieldID ),
      boundaryHandlingID_( boundaryHandlingID ), bodyStorageID_( bodyStorageID ), fileName_( fileName ), fileIO_(fileIO),
      dx_SI_( dx_SI ), dt_SI_( dt_SI ), diameter_( diameter ), mass_( mass ),
      position_( real_t(0) ), maxVelocity_( real_t(0) )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t t\t posX\t posY\t gapZ\t velX\t velY\t velZ\n";
            file.close();
         }
      }
   }

   void operator()()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep());

      Vector3<real_t> pos(real_t(0));
      Vector3<real_t> transVel(real_t(0));
      Vector3<real_t> force(real_t(0));
      Vector3<real_t> fluidMomentum(real_t(0));

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         PdfField_T * pdfField = blockIt->getData<PdfField_T>(pdfFieldID_);
         BoundaryHandling_T * bh = blockIt->getData<BoundaryHandling_T>(boundaryHandlingID_);
         WALBERLA_FOR_ALL_CELLS_XYZ(pdfField,
               if( bh->isDomain(Cell(x,y,z)) ) fluidMomentum += pdfField->getMomentumDensity(x,y,z);)

         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            pos = bodyIt->getPosition();
            transVel = bodyIt->getLinearVel();
         }
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            force += bodyIt->getForce();
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( pos, mpi::SUM );
         mpi::allReduceInplace( transVel, mpi::SUM );
         mpi::allReduceInplace( force, mpi::SUM );
         mpi::allReduceInplace( fluidMomentum, mpi::SUM );
      }

      position_ = pos[2];
      maxVelocity_ = std::max(maxVelocity_, -transVel[2]);

      if( fileIO_ )
         writeToFile( timestep, pos, transVel, force, fluidMomentum);
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
   void writeToFile( uint_t timestep, const Vector3<real_t> & position, const Vector3<real_t> & velocity, const Vector3<real_t> & force, const Vector3<real_t> & fluidMomentum )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         auto scaledPosition = position / diameter_;
         auto velocity_SI = velocity * dx_SI_ / dt_SI_;


         file << timestep << "\t" << real_c(timestep) * dt_SI_ << "\t"
              << "\t" << scaledPosition[0] << "\t" << scaledPosition[1] << "\t" << scaledPosition[2] - real_t(0.5)
              << "\t" << velocity_SI[0] << "\t" << velocity_SI[1] << "\t" << velocity_SI[2]
              << "\t" << force[0] << "\t" << force[1] << "\t" << force[2]
              << "\t" << fluidMomentum[0] << "\t" << fluidMomentum[1] << "\t" << fluidMomentum[2]
              << "\t" << mass_ * velocity[0] << "\t" << mass_ * velocity[1] << "\t" << mass_ * velocity[2]
              << "\n";
         file.close();
      }
   }

   SweepTimeloop* timeloop_;
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID bodyStorageID_;
   std::string fileName_;
   bool fileIO_;
   real_t dx_SI_, dt_SI_, diameter_, mass_;

   real_t position_;
   real_t maxVelocity_;
};


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
   bool fileIO = false;
   uint_t vtkIOFreq = 0;
   std::string baseFolder = "vtk_out_SettlingSphere";

   // physical setup
   uint_t fluidType = 1;

   //numerical parameters
   uint_t numberOfCellsInHorizontalDirection = uint_t(100);
   bool averageForceTorqueOverTwoTimSteps = true;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" )         == 0 ) { shortrun = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" )         == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )           == 0 ) { fileIO = true; continue; }
      if( std::strcmp( argv[i], "--vtkIOFreq" )        == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--fluidType" )        == 0 ) { fluidType = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--resolution" )       == 0 ) { numberOfCellsInHorizontalDirection = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--noForceAveraging" ) == 0 ) { averageForceTorqueOverTwoTimSteps = false; continue; }
      if( std::strcmp( argv[i], "--baseFolder" )       == 0 ) { baseFolder = argv[++i]; continue; }
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

   real_t densityFluid_SI, dynamicViscosityFluid_SI;
   real_t expectedSettlingVelocity_SI;
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
   Vector3<real_t> domainSize_SI(real_t(100e-3), real_t(100e-3), real_t(160e-3));
   //shift starting gap a bit upwards to match the reported (plotted) values
   const real_t startingGapSize_SI = real_t(120e-3) + real_t(0.25) * diameter_SI;

   WALBERLA_LOG_INFO_ON_ROOT("Setup (in SI units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter_SI << ", density = " << densitySphere_SI << ", starting gap size = " << startingGapSize_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid_SI << ", dyn. visc = " << dynamicViscosityFluid_SI << ", kin. visc = " << kinematicViscosityFluid_SI );
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedSettlingVelocity_SI << " --> Re_p = " << expectedSettlingVelocity_SI * diameter_SI / kinematicViscosityFluid_SI );


   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////


   const real_t dx_SI = domainSize_SI[0] / real_c(numberOfCellsInHorizontalDirection);
   const Vector3<uint_t> domainSize( uint_c(floor(domainSize_SI[0] / dx_SI + real_t(0.5)) ),
                                     uint_c(floor(domainSize_SI[1] / dx_SI + real_t(0.5)) ),
                                     uint_c(floor(domainSize_SI[2] / dx_SI + real_t(0.5)) ) );
   const real_t diameter = diameter_SI / dx_SI;
   const real_t sphereVolume = math::pi / real_t(6) * diameter * diameter * diameter;

   const real_t expectedSettlingVelocity = real_t(0.01);
   const real_t dt_SI = expectedSettlingVelocity / expectedSettlingVelocity_SI * dx_SI;

   const real_t viscosity =  kinematicViscosityFluid_SI * dt_SI / ( dx_SI * dx_SI );
   const real_t relaxationTime = real_t(1) / lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t gravitationalAcceleration = gravitationalAcceleration_SI * dt_SI * dt_SI / dx_SI;

   const real_t densityFluid = real_t(1);
   const real_t densitySphere = densityFluid * densitySphere_SI / densityFluid_SI;

   const real_t dx = real_t(1);

   const uint_t timesteps = funcTest ? 1 : ( shortrun ? uint_t(200) : uint_t( 250000 ) );
   const uint_t numPeSubCycles = uint_t(1);

   WALBERLA_LOG_INFO_ON_ROOT(" - dx_SI = " << dx_SI << ", dt_SI = " << dt_SI);
   WALBERLA_LOG_INFO_ON_ROOT("Setup (in simulation, i.e. lattice, units):");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize);
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere: diameter = " << diameter << ", density = " << densitySphere );
   WALBERLA_LOG_INFO_ON_ROOT(" - fluid: density = " << densityFluid << ", relaxation time (tau) = " << relaxationTime << ", kin. visc = " << viscosity );
   WALBERLA_LOG_INFO_ON_ROOT(" - gravitational acceleration = " << gravitationalAcceleration );
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedSettlingVelocity << " --> Re_p = " << expectedSettlingVelocity * diameter / viscosity );

   if( vtkIOFreq > 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing vtk files to folder \"" << baseFolder << "\" with frequency " << vtkIOFreq);
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   Vector3<uint_t> numberOfBlocksPerDirection( uint_t(5), uint_t(5), uint_t(8) );
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
                                                      true, true, true, //periodicity
                                                      false );

   WALBERLA_LOG_INFO_ON_ROOT("Domain decomposition:");
   WALBERLA_LOG_INFO_ON_ROOT(" - blocks per direction = " << numberOfBlocksPerDirection );
   WALBERLA_LOG_INFO_ON_ROOT(" - cells per block = " << cellsPerBlockPerDirection );

   //write domain decomposition to file
   if( vtkIOFreq > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolder );
   }

   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage>  globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID          = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID          = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   // set up collision response, here DEM solver
   pe::cr::DEM cr(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, nullptr);

   // set up synchronization procedure
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );


   // create pe bodies

   // bounding planes (global)
   /*
   const auto planeMaterial = pe::createMaterial( "myPlaneMat", real_t(8920), real_t(0), real_t(1), real_t(1), real_t(0), real_t(1), real_t(1), real_t(0), real_t(0) );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(1,0,0), Vector3<real_t>(0,0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(-1,0,0), Vector3<real_t>(real_c(domainSize[0]),0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,1,0), Vector3<real_t>(0,0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,-1,0), Vector3<real_t>(0,real_c(domainSize[1]),0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(domainSize[2])), planeMaterial );
   */

   // add the sphere
   const auto sphereMaterial = pe::createMaterial( "mySphereMat", densitySphere , real_t(0.5), real_t(0.1), real_t(0.1), real_t(0.24), real_t(200), real_t(200), real_t(0), real_t(0) );
   Vector3<real_t> initialPosition( real_t(0.5) * real_c(domainSize[0]), real_t(0.5) * real_c(domainSize[1]), startingGapSize_SI / dx_SI + real_t(0.5) * diameter);
   pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, initialPosition, real_t(0.5) * diameter, sphereMaterial );

   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( real_t(1) / relaxationTime ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         uint_t(1), field::fzyx );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::fzyx );

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );

   // map planes into the LBM simulation -> act as no-slip boundaries
   //pe_coupling::mapBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, NoSlip_Flag, pe_coupling::selectGlobalBodies );

   // map pe bodies into the LBM simulation
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );

   ///////////////
   // TIME LOOP //
   ///////////////

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
   commFunction = scheme;

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // sweep for updating the pe body mapping into the LBM simulation
   timeloop.add()
         << Sweep( pe_coupling::BodyMapping< LatticeModel_T, BoundaryHandling_T >( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MO_Flag, FormerMO_Flag, pe_coupling::selectRegularBodies ), "Body Mapping" );

   // sweep for restoring PDFs in cells previously occupied by pe bodies
   pe_coupling::SphereNormalExtrapolationDirectionFinder extrapolationFinder( blocks, bodyFieldID );
   typedef pe_coupling::EquilibriumAndNonEquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T, pe_coupling::SphereNormalExtrapolationDirectionFinder > Reconstructor_T;
   Reconstructor_T reconstructor( blocks, boundaryHandlingID, bodyFieldID, extrapolationFinder );
   timeloop.add()
      << Sweep( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T >
                  ( blocks, pdfFieldID, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, reconstructor, FormerMO_Flag, Fluid_Flag ), "PDF Restore" );

   shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer1 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   std::function<void(void)> storeForceTorqueInCont1 = std::bind(&pe_coupling::BodiesForceTorqueContainer::store, bodiesFTContainer1);
   shared_ptr<pe_coupling::BodiesForceTorqueContainer> bodiesFTContainer2 = make_shared<pe_coupling::BodiesForceTorqueContainer>(blocks, bodyStorageID);
   std::function<void(void)> setForceTorqueOnBodiesFromCont2 = std::bind(&pe_coupling::BodiesForceTorqueContainer::setOnBodies, bodiesFTContainer2);
   shared_ptr<pe_coupling::ForceTorqueOnBodiesScaler> forceScaler = make_shared<pe_coupling::ForceTorqueOnBodiesScaler>(blocks, bodyStorageID, real_t(1));
   std::function<void(void)> setForceScalingFactorToHalf = std::bind(&pe_coupling::ForceTorqueOnBodiesScaler::resetScalingFactor,forceScaler,real_t(0.5));

   if( averageForceTorqueOverTwoTimSteps ) {
      bodiesFTContainer2->store();
   }

   // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // stream + collide LBM step
   timeloop.add()
      << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag ) ), "cell-wise LB sweep" );

   // Averaging the force/torque over two time steps is said to damp oscillations of the interaction force/torque.
   // See Ladd - " Numerical simulations of particulate suspensions via a discretized Boltzmann equation. Part 1. Theoretical foundation", 1994, p. 302
   if( averageForceTorqueOverTwoTimSteps ) {

      // store force/torque from hydrodynamic interactions in container1
      timeloop.addFuncAfterTimeStep(storeForceTorqueInCont1, "Force Storing");

      // set force/torque from previous time step (in container2)
      timeloop.addFuncAfterTimeStep(setForceTorqueOnBodiesFromCont2, "Force setting");

      // average the force/torque by scaling it with factor 1/2 (except in first timestep, there it is 1, which it is initially)
      timeloop.addFuncAfterTimeStep( SharedFunctor<pe_coupling::ForceTorqueOnBodiesScaler>(forceScaler),  "Force averaging");
      timeloop.addFuncAfterTimeStep( setForceScalingFactorToHalf, "Force scaling adjustment" );

      // swap containers
      timeloop.addFuncAfterTimeStep( pe_coupling::BodyContainerSwapper( bodiesFTContainer1, bodiesFTContainer2 ), "Swap FT container" );

   }

   // check for convergence of the particle position
   std::string loggingFileName( baseFolder + "/LoggingSettlingSphere_");
   loggingFileName += std::to_string(fluidType);
   loggingFileName += ".txt";
   if( fileIO  )
   {
      WALBERLA_LOG_INFO_ON_ROOT(" - writing logging output to file \"" << loggingFileName << "\"");
   }
   shared_ptr< SpherePropertyLogger > logger = walberla::make_shared< SpherePropertyLogger >( &timeloop, blocks, pdfFieldID, boundaryHandlingID, bodyStorageID,
                                                                                              loggingFileName, fileIO, dx_SI, dt_SI, diameter, sphereVolume * densitySphere );
   timeloop.addFuncAfterTimeStep( SharedFunctor< SpherePropertyLogger >( logger ), "Sphere property logger" );

   Vector3<real_t> gravitationalForce( real_t(0), real_t(0), -(densitySphere - densityFluid) * gravitationalAcceleration * sphereVolume );
   timeloop.addFuncAfterTimeStep( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, gravitationalForce ), "Gravitational force" );

   // add pe timesteps
   timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, real_t(1), numPeSubCycles ), "pe Time Step" );

   if( vtkIOFreq != uint_t(0) )
   {
      // spheres
      auto bodyVtkOutput   = make_shared<pe::SphereVtkOutput>( bodyStorageID, blocks->getBlockStorage() );
      auto bodyVTK   = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", vtkIOFreq, baseFolder );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( bodyVTK ), "VTK (sphere data)" );

      // flag field (written only once in the first time step, ghost layers are also written)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, FieldGhostLayers, false, baseFolder );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", vtkIOFreq, 0, false, baseFolder );

      blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync( blocks );
      pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
      pdfFieldVTK->addBeforeFunction( pdfGhostLayerSync );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );
      pdfFieldVTK->addCellInclusionFilter( fluidFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );
   }


   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   real_t terminationPosition = diameter; // right before sphere touches the bottom wall

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {
      // perform a single simulation step
      timeloop.singleStep( timeloopTiming );

      if ( logger->getPosition() < terminationPosition )
      {
         WALBERLA_LOG_INFO_ON_ROOT("Sphere reached terminal position " << logger->getPosition() << " after " << i << " timesteps!");
         break;
      }
   }

   timeloopTiming.logResultOnRoot();

   // check the result
   if ( !funcTest && !shortrun )
   {
      real_t relErr = std::fabs( expectedSettlingVelocity - logger->getMaxVelocity()) / expectedSettlingVelocity;
      WALBERLA_LOG_INFO_ON_ROOT( "Expected maximum settling velocity: " << expectedSettlingVelocity );
      WALBERLA_LOG_INFO_ON_ROOT( "Simulated maximum settling velocity: " << logger->getMaxVelocity() );
      WALBERLA_LOG_INFO_ON_ROOT( "Relative error: " << relErr );

      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS( relErr, real_t(0.1) );
   }

   return EXIT_SUCCESS;
}

} // namespace settling_sphere_mem

int main( int argc, char **argv ){
   settling_sphere_mem::main(argc, argv);
}
