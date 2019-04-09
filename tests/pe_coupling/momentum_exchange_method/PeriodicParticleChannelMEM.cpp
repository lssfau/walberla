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
//! \file PeriodicParticleChannelMEM.cpp
//! \ingroup pe_coupling
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
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

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/all.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/vtk/all.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"
#include "pe/cr/ICR.h"
#include "pe/Types.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "vtk/all.h"

#include <functional>
#include <vector>

namespace periodic_particle_channel_mem
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

//////////////
// TYPEDEFS //
//////////////

// pdf field & flag field
typedef lbm::D3Q19< lbm::collision_model::TRT >  LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

const uint_t FieldGhostLayers = 1;

typedef lbm::NoSlip< LatticeModel_T, flag_t >                 NoSlip_T;
typedef lbm::SimpleUBB< LatticeModel_T, flag_t >              UBB_T;
typedef pe_coupling::SimpleBB< LatticeModel_T, FlagField_T > MO_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, UBB_T, MO_T > BoundaryHandling_T;

typedef std::tuple< pe::Sphere, pe::Plane > BodyTypeTuple;

///////////
// FLAGS //
///////////

const FlagUID    Fluid_Flag( "fluid" );
const FlagUID      UBB_Flag( "velocity bounce back" );
const FlagUID   NoSlip_Flag( "no slip" );
const FlagUID       MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );



/////////////////////
// CHECK FUNCTIONS //
/////////////////////

template< typename LatticeModel_T, typename FlagField_T >
class VelocityCheck
{
public:

   using PdfField = lbm::PdfField<LatticeModel_T>;

   VelocityCheck( const ConstBlockDataID & pdfFieldID, const ConstBlockDataID & flagFieldID, const Set< FlagUID > & cellsToCheck,
                  const real_t uMax, const uint_t checkFrequency ) :
      executionCounter_( uint_c(0) ), checkFrequency_( ( checkFrequency > 0 ) ? checkFrequency : uint_c(1) ),
      cellsToCheck_( cellsToCheck ), flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), uMax_( uMax * uMax ) {}

   void operator()( const IBlock * const block )
   {
      ++executionCounter_;
      if( ( executionCounter_ - uint_c(1) ) % checkFrequency_ != 0 )
         return;

      const FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
      const PdfField    *  pdfField = block->getData< PdfField    >(  pdfFieldID_ );

      const CellInterval & cellBB = pdfField->xyzSize();
      WALBERLA_ASSERT_EQUAL( flagField->xyzSize(), cellBB );

      typename FlagField_T::flag_t mask = 0;
      for( auto flag = cellsToCheck_.begin(); flag != cellsToCheck_.end(); ++flag )
         mask = static_cast< typename FlagField_T::flag_t >( mask | flagField->getFlag( *flag ) );

      for( cell_idx_t z = cellBB.zMin(); z <= cellBB.zMax(); ++z ) {
         for( cell_idx_t y = cellBB.yMin(); y <= cellBB.yMax(); ++y ) {
            for( cell_idx_t x = cellBB.xMin(); x <= cellBB.xMax(); ++x )
            {
               if( flagField->isPartOfMaskSet(x,y,z,mask) )
               {
                   WALBERLA_CHECK_LESS_EQUAL( pdfField->getVelocity( x, y, z ).sqrLength(), uMax_ );
               }
            }
         }
      }
   }

private:

         uint_t executionCounter_;
   const uint_t checkFrequency_;

   const Set< FlagUID > cellsToCheck_;

   const ConstBlockDataID flagFieldID_;
   const ConstBlockDataID  pdfFieldID_;

   const real_t uMax_;

}; // VelocityCheck



class ObstacleLocationCheck
{
public:

   ObstacleLocationCheck( const shared_ptr< StructuredBlockStorage > & blocks,
                          const BlockDataID & bodyStorageID,
                          const AABB & aabb, const real_t additionalSpace )
   :  blocks_( blocks ), bodyStorageID_( bodyStorageID ),
      aabb_( aabb.xMin() - additionalSpace, aabb.yMin() - additionalSpace, aabb.zMin() - additionalSpace,
             aabb.xMax() + additionalSpace, aabb.yMax() + additionalSpace, aabb.zMax() + additionalSpace ) {}

   void operator()()
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            WALBERLA_CHECK( aabb_.contains( bodyIt->getPosition()[0], bodyIt->getPosition()[1], bodyIt->getPosition()[2] ) );
         }
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   const AABB aabb_;

}; // ObstacleLocationCheck



class ObstacleLinearVelocityCheck
{
public:

   ObstacleLinearVelocityCheck( const shared_ptr< StructuredBlockStorage > & blocks,
                                const BlockDataID & bodyStorageID, const real_t uMax )
   : blocks_( blocks ), bodyStorageID_( bodyStorageID ), uMax_( uMax * uMax )
   {}

   void operator()()
   {
      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            const auto & u = bodyIt->getLinearVel();
            WALBERLA_CHECK_LESS_EQUAL( (u[0]*u[0] + u[1]*u[1] + u[2]*u[2]), uMax_ );

         }
      }
   }

private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   const real_t uMax_;

}; // ObstacleLinearVelocityCheck



/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

// class MyBoundaryHandling is responsible for creating the boundary handling and setting up the geometry at the outer domain boundaries

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID, const real_t velocity ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_( bodyFieldID ), velocity_( velocity ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID  pdfFieldID_;
   const BlockDataID bodyFieldID_;

   const real_t velocity_;

}; // class MyBoundaryHandling



BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfFieldID_ );
   BodyField_T * bodyField = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "cf boundary handling", flagField, fluid,
                                    NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                       UBB_T( "UBB", UBB_Flag, pdfField, velocity_, real_c(0), real_c(0) ),
                                        MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   const auto ubb = flagField->getFlag( UBB_Flag );

   CellInterval domainBB = storage->getDomainCellBB();
   domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.xMax() += cell_idx_c( FieldGhostLayers );

   domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.yMax() += cell_idx_c( FieldGhostLayers );

   // SOUTH
   CellInterval south( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( south, *block );
   handling->forceBoundary( ubb, south );

   // NORTH
   CellInterval north( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( north, *block );
   handling->forceBoundary( ubb, north );

   domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.zMax() += cell_idx_c( FieldGhostLayers );

   // BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
   handling->forceBoundary( ubb, bottom );

   // TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( top, *block );
   handling->forceBoundary( ubb, top );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


//////////
// MAIN //
//////////

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto processes = MPIManager::instance()->numProcesses();
   if( processes != 1 && processes % 4 != 0  )
   {
      std::cerr << "Number of processes must be equal to either 1 or a multiple of 4!" << std::endl;
      return EXIT_FAILURE;
   }

   bool shortrun       = false;
   bool useFZYX        = false;
   bool useSplitKernel = false;
   bool fullPDFSync    = false;
   bool vtkIO          = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) { shortrun       = true; continue; }
      if( std::strcmp( argv[i], "--fzyx" )     == 0 ) { useFZYX        = true; continue; }
      if( std::strcmp( argv[i], "--split" )    == 0 ) { useSplitKernel = true; continue; }
      if( std::strcmp( argv[i], "--full" )     == 0 ) { fullPDFSync    = true; continue; }
      if( std::strcmp( argv[i], "--vtkIO" )    == 0 ) { vtkIO          = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const uint_t length = uint_c(200);
   const uint_t width  = uint_c(100);

   const real_t omega    = real_c(1.33);
   const real_t velocity = real_c(0.001);
   const real_t nu_SI    = real_c(0.00005); // kinematic viscosity - 50 times larger than water [ m^2 / s ]
   //const real_t rho_SI   = real_c(1000);    // density of water [ kg / m^3 ]
   const real_t L_SI     = real_c(0.04);    // length of channel [ m ]

   const real_t nu_L  = ( real_t(1) / omega - real_c(0.5) ) / real_t(3);
   const real_t dx_SI = L_SI / real_c( length );          // dx in [ m ]
   const real_t dt_SI = ( nu_L * dx_SI * dx_SI ) / nu_SI; // dt in [ s ]

   const uint_t timesteps = shortrun ? uint_c(10) : uint_c(150000); // -> ~10 seconds of real time

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   // (200x100x100) - periodic in x-direction!
   Vector3<uint_t> blockDist = getFactors3D( uint_c(processes), Vector3<real_t>( real_c(length), real_c(width), real_c(width) ) );
   if ( processes == 1 )
   {
      // due to periodicity, the pe needs at least three blocks in the periodic direction
      blockDist[0] = uint_t(4);
   }

   const uint_t xCells  = length / blockDist[0];
   const uint_t yCells  =  width / blockDist[1];
   const uint_t zCells  =  width / blockDist[2];

   const real_t dx = real_t(1);

   auto blocks = blockforest::createUniformBlockGrid( blockDist[0], blockDist[1], blockDist[2], xCells, yCells, zCells, dx,
                                                      ( processes != 1 ),
                                                      true, false, false );
   if( vtkIO )
   {
      vtk::writeDomainDecomposition( blocks );
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
   const real_t overlap = real_c( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // create pe bodies
   const auto material = pe::createMaterial( "granular", real_c(1.2), real_c(0.25), real_c(0.4), real_c(0.4), real_c(0.35), real_c(1.39e11), real_c(5.18e7), real_c(1.07e2), real_c(1.07e2) );

   // global bodies
   // bounding planes
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), material );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(width)), material );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,1,0), Vector3<real_t>(0,0,0), material );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,-1,0), Vector3<real_t>(0,real_c(width),0), material );

   // spheres as obstacles
   std::vector<pe::BodyID> globalBodiesToBeMapped;
   auto globalSphere1 = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>(real_c(length) / real_t(2), real_t(50), real_t(110)), real_t(60), material, true, false, true );
   globalBodiesToBeMapped.push_back(globalSphere1);
   auto globalSphere2 = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>(                 real_t(0), real_t(50), -real_t(60)), real_t(80), material, true, false, true );
   globalBodiesToBeMapped.push_back(globalSphere2);
   auto globalSphere3 = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>(            real_c(length), real_t(50), -real_t(60)), real_t(80), material, true, false, true );
   globalBodiesToBeMapped.push_back(globalSphere3);

   // local bodies: moving spheres
   const real_t radius = real_t(10);

   auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(15), real_t(50), real_t(35) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(15), real_t(35), real_t(50) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(15), real_t(65), real_t(50) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(15), real_t(50), real_t(65) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(35), real_t(35), real_t(35) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(35), real_t(65), real_t(35) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(35), real_t(35), real_t(65) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(35), real_t(65), real_t(65) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( real_t(50), real_t(50), real_t(50) ), radius, material );
   if( sphere != nullptr ) sphere->setLinearVel( velocity, real_t(0), real_t(0) );

   //synchronize the pe set up on all processes
   syncCall();

   ////////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add pdf field
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   BlockDataID pdfFieldID = useFZYX ? lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel,
                                                                 Vector3< real_t >( velocity, real_t(0), real_t(0) ), real_t(1),
                                                                 uint_t(1), field::fzyx ) :
                                      lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel,
                                                                 Vector3< real_t >( velocity, real_t(0), real_t(0) ), real_t(1),
                                                                 uint_t(1), field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf );

   // add boundary handling & initialize outer domain boundaries (moving walls on the front, back, top, and bottom plane)
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
              MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID, real_c(6) * velocity ), "boundary handling" );

   // map pe bodies into the LBM simulation
   // global bodies act as no-slip obstacles and are not treated by the fluid-solid coupling
   // special care has to be taken here that only the global spheres (not the planes) are mapped since the planes would overwrite the already set boundary flags
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt)
   {
      for( auto globalBodyIt = globalBodiesToBeMapped.begin(); globalBodyIt != globalBodiesToBeMapped.end(); ++globalBodyIt )
      {
         pe_coupling::mapBody< BoundaryHandling_T >( *globalBodyIt, *blockIt, *blocks, boundaryHandlingID, NoSlip_Flag );
      }
   }

   // moving bodies are handled by the momentum exchange method
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag, pe_coupling::selectRegularBodies );

   ///////////////
   // TIME LOOP //
   ///////////////

   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // sweep for updating the pe body mapping into the LBM simulation
   timeloop.add()
      << Sweep( pe_coupling::BodyMapping< BoundaryHandling_T >( blocks, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID, MO_Flag, FormerMO_Flag, pe_coupling::selectRegularBodies ),
                "Body Mapping" );

   // sweep for restoring PDFs in cells previously occupied by pe bodies
   typedef pe_coupling::EquilibriumReconstructor< LatticeModel_T, BoundaryHandling_T > Reconstructor_T;
   Reconstructor_T reconstructor( blocks, boundaryHandlingID, pdfFieldID, bodyFieldID );
   timeloop.add()
      << Sweep( pe_coupling::PDFReconstruction< LatticeModel_T, BoundaryHandling_T, Reconstructor_T >( blocks, boundaryHandlingID, bodyStorageID, globalBodyStorage, bodyFieldID,
                                                                                                       reconstructor, FormerMO_Flag, Fluid_Flag  ), "PDF Restore" );

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   if( fullPDFSync )
   {
      blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > scheme( blocks );
      scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
      commFunction = scheme;
   }
   else
   {
      blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
      scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
      commFunction = scheme;
   }

   // add LBM communication function and boundary handling sweep
   timeloop.add()
      << BeforeFunction( commFunction, "LBM Communication" )
      << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // LBM stream collide sweep
   if( useSplitKernel )
      timeloop.add() << Sweep( lbm::SplitPureSweep< LatticeModel_T >( pdfFieldID ), "LBM SPLIT PURE" );
   else
      timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag ) ), "LBM DEFAULT" );

   // some checks if everything is still fine/stable (checks run every 100th time step)
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( field::makeStabilityChecker< lbm::PdfField< LatticeModel_T >, FlagField_T >( blocks, pdfFieldID, flagFieldID, Fluid_Flag,
                                                                                                                                  uint_t(100), false, true ) ),
                                  "LBM stability check" );
   timeloop.add()
      << Sweep( VelocityCheck< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag, real_c(0.5), uint_c(100) ), "LBM Velocity Check" );

   // add lubrication force correction sweep
   timeloop.addFuncAfterTimeStep( pe_coupling::LubricationCorrection( blocks, globalBodyStorage, bodyStorageID, nu_L ), "Lubrication Force" );

   // add pe timesteps
   timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, real_c(1.0), uint_c(10) ), "pe Time Step" );

   // some checks if everything is still fine in the pe simulation world
   timeloop.addFuncAfterTimeStep( ObstacleLocationCheck( blocks, bodyStorageID, blocks->getDomain(), radius + real_c(2) * overlap ), "Obstacle Location Check" );
   timeloop.addFuncAfterTimeStep( ObstacleLinearVelocityCheck( blocks, bodyStorageID, real_c(0.5) ), "Obstacle Velocity Check" );

   ////////////////
   // VTK OUTPUT //
   ////////////////

   if( vtkIO )
   {
      const uint_t writeFrequency = vtk::determineWriteFrequency( dt_SI, uint_c(30) );

      // spheres
      auto bodyVtkOutput   = make_shared<pe::DefaultBodyVTKOutput>( bodyStorageID, blocks->getBlockStorage() );
      auto bodyVTK   = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", writeFrequency );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( bodyVTK ), "VTK (sphere data)" );

      // flag field (written only once in the first time step, ghost layers are also written)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, FieldGhostLayers );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field (ghost layers cannot be written because re-sampling/coarsening is applied)
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", writeFrequency );
      pdfFieldVTK->setSamplingResolution( real_c(2.5) );

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
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();

   return 0;
}

} //namespace periodic_particle_channel_mem

int main( int argc, char **argv ){
   periodic_particle_channel_mem::main(argc, argv);
}
