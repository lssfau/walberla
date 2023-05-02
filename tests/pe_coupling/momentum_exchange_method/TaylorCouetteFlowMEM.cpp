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
//! \file TaylorCouetteFlowMEM.cpp
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
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
#include "field/communication/PackInfo.h"
#include "field/vtk/all.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/vtk/all.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"
#include "pe/Types.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"

#include "vtk/all.h"

#include <functional>

namespace taylor_coette_flow_mem
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

typedef pe_coupling::CurvedLinear< LatticeModel_T, FlagField_T > MO_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, MO_T > BoundaryHandling_T;

typedef std::tuple< pe::Capsule, pe::CylindricalBoundary > BodyTypeTuple;

///////////
// FLAGS //
///////////

const FlagUID    Fluid_Flag( "fluid" );
const FlagUID       MO_Flag( "moving obstacle" );
const FlagUID FormerMO_Flag( "former moving obstacle" );



/////////////////////
// CHECK FUNCTIONS //
/////////////////////

class VelocityChecker
{
public:
   VelocityChecker( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & pdfFieldID,
                    real_t radius1, real_t radius2, real_t angularVel1, real_t angularVel2, real_t domainLength, real_t domainWidth ) :
   blocks_( blocks ), pdfFieldID_( pdfFieldID ), radius1_( radius1 ), radius2_( radius2 ), angularVel1_( angularVel1 ), angularVel2_( angularVel2 ),
   domainLength_( domainLength ), domainWidth_( domainWidth )
   { }
   real_t getMaximumRelativeError()
   {
      Vector3<real_t> midPoint( domainLength_ * real_t(0.5), domainWidth_ * real_t(0.5), domainWidth_ * real_t(0.5) );
      Cell midCell = blocks_->getCell( midPoint );
      CellInterval evaluationCells( midCell.x(), midCell.y(), midCell.z(), midCell.x(), blocks_->getDomainCellBB().yMax(), midCell.z() );

      real_t maxError = real_t(0);

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         auto pdfField = blockIt->getData< PdfField_T > ( pdfFieldID_ );
         CellInterval intersection = blocks_->getBlockCellBB( *blockIt );
         intersection.intersect(  evaluationCells );
         blocks_->transformGlobalToBlockLocalCellInterval( intersection, *blockIt );
         for( auto fieldIt = pdfField->beginSliceXYZ(intersection); fieldIt != pdfField->end(); ++fieldIt )
         {
            const Vector3< real_t > cellCenter = blocks_->getBlockLocalCellCenter( *blockIt, fieldIt.cell() );
            auto diff = cellCenter - midPoint;
            real_t radius = diff.length();
            if( radius > radius1_ && radius < radius2_ )
            {
               real_t velSimulation = pdfField->getVelocity( fieldIt.cell() )[2];
               real_t velAnalytical = getAnalyticalVelocity( radius );
               real_t relVelDiff = std::fabs(velSimulation - velAnalytical) / std::fabs(velAnalytical);
               maxError = std::max( maxError, relVelDiff );
            }
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( maxError, mpi::MAX );
      }
      return maxError;
   }
private:
   real_t getAnalyticalVelocity( real_t r )
   {
      real_t radius1sqr = radius1_ * radius1_;
      real_t radius2sqr = radius2_ * radius2_;
      real_t a = ( angularVel2_ * radius2sqr - angularVel1_ * radius1sqr ) / ( radius2sqr - radius1sqr );
      real_t b = ( angularVel1_ - angularVel2_ ) * radius1sqr * radius2sqr / ( radius2sqr - radius1sqr );
      real_t velInAgularDir = a * r + b / r;
      return velInAgularDir;
   }

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID pdfFieldID_;
   real_t radius1_, radius2_, angularVel1_, angularVel2_, domainLength_, domainWidth_;
};

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

// class MyBoundaryHandling is responsible for creating the boundary handling and setting up the geometry at the outer domain boundaries

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyFieldID ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), bodyFieldID_( bodyFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID  pdfFieldID_;
   const BlockDataID bodyFieldID_;

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
                                                           MO_T( "MO", MO_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   CellInterval domainBB = storage->getDomainCellBB();
   domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.xMax() += cell_idx_c( FieldGhostLayers );

   domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.yMax() += cell_idx_c( FieldGhostLayers );

   domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.zMax() += cell_idx_c( FieldGhostLayers );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


//////////
// MAIN //
//////////
/*
 * This test features a so-called Taylor-Couette flow, i.e. a Couette flow between two rotating, coaxial cylinders
 * For such a setup, the velocity in angular direction can be calculated analytically.
 */
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

   bool shortrun = false;
   uint_t vtkIOFreq = 0;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) { shortrun = true; continue; }
      if( std::strcmp( argv[i], "--vtkIOFreq" ) == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   bool vtkIO =  vtkIOFreq != 0;

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const uint_t length = uint_c(20);
   const uint_t width = uint_c(101);
   if( width % 2 == 0 ) WALBERLA_ABORT("Width has to be odd, since otherwise the current implementation of the velocity evaluation does not work correctly!");

   const real_t omega = real_c(1.33);

   const uint_t timesteps = shortrun ? uint_c(10) : uint_c(5000);

   const real_t radius1 = real_c(width) * real_t(0.25); // radius of internal cylinder
   const real_t radius2 = real_c(width) * real_t(0.45); // radius of external cylinder
   const real_t angularVel1 = real_t(0.001); // angular velocity of internal cylinder
   const real_t angularVel2 = real_t(-0.001); // angular velocity of external cylinder

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   // (length x width x width) - periodic in x-direction!
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
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");

   // set up synchronization procedure
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // create pe bodies
   const auto material = pe::createMaterial( "granular", real_t(1.2), real_t(0.25), real_t(0.4), real_t(0.4), real_t(0.35), real_t(1.39e11), real_t(5.18e7), real_t(1.07e2), real_t(1.07e2) );

   // instead of a cylinder, we use the capsule and make sure it extends the computational domain in x-direction to effectively get a cylinder
   auto cap = pe::createCapsule( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, pe::Vec3(real_c(length),real_c(width),real_c(width)) * real_t(0.5), radius1, real_c(length)*real_t(2), material, true, false, true );
   cap->setAngularVel( pe::Vec3(1,0,0) * angularVel1 );

   auto cb = pe::createCylindricalBoundary( *globalBodyStorage, 0, pe::Vec3(real_c(length),real_c(width),real_c(width))*real_t(0.5), radius2 );
   cb->setAngularVel( pe::Vec3(1,0,0) * angularVel2 );

   //synchronize the pe set up on all processes
   syncCall();

   ////////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // add pdf field
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage( blocks, "pdf field (fzyx)", latticeModel, Vector3< real_t >( real_t(0) ), real_t(1), uint_t(1), field::fzyx );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::fzyx );

   // add boundary handling & initialize outer domain boundaries (moving walls on the front, back, top, and bottom plane)
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
              MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );

   // map pe bodies into the LBM simulation
   // moving bodies are handled by the momentum exchange method, thus act here as velocity boundary condition
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_Flag );

   ///////////////
   // TIME LOOP //
   ///////////////

   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
   commFunction = scheme;

   // add LBM communication function and boundary handling sweep
   timeloop.add()
      << BeforeFunction( commFunction, "LBM Communication" )
      << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // LBM stream collide sweep
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag ) ), "LBM DEFAULT" );

   ////////////////
   // VTK OUTPUT //
   ////////////////

   if( vtkIO )
   {
      const uint_t writeFrequency = 10; //vtk::determineWriteFrequency( dt_SI, uint_c(30) );

      // flag field (written only once in the first time step, ghost layers are also written)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, FieldGhostLayers );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field (ghost layers cannot be written because re-sampling/coarsening is applied)
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", writeFrequency );
      pdfFieldVTK->setSamplingResolution( real_c(1) );

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

   VelocityChecker velChecker(blocks, pdfFieldID, radius1, radius2, angularVel1, angularVel2, real_c(length), real_c(width) );
   real_t maxError = velChecker.getMaximumRelativeError();
   if( !shortrun )
   {
      timeloopTiming.logResultOnRoot();
      WALBERLA_LOG_RESULT_ON_ROOT("Maximum relative error in velocity in angular direction: " << maxError );
      // error has to be below 10%
      WALBERLA_CHECK_LESS( maxError, real_t(0.1) );
   }

   return 0;
}

} //namespace taylor_coette_flow_mem

int main( int argc, char **argv ){
   taylor_coette_flow_mem::main(argc, argv);
}
