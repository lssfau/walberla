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
//! \file ForcesOnSphereNearPlaneInShearFlow.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"
#include <blockforest/loadbalancing/StaticCurve.h>
#include <blockforest/SetupBlockForest.h>

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

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/refinement/all.h"
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

namespace forces_on_sphere_near_plane_in_shear_flow
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
typedef lbm::D3Q19< lbm::collision_model::TRT, false >  LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

const uint_t FieldGhostLayers = 4;

// boundary handling
typedef pe_coupling::SimpleBB< LatticeModel_T, FlagField_T > MO_SBB_T;
typedef pe_coupling::CurvedLinear< LatticeModel_T, FlagField_T > MO_CLI_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, MO_SBB_T, MO_CLI_T > BoundaryHandling_T;

typedef std::tuple< pe::Sphere, pe::Plane > BodyTypeTuple;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID MO_SBB_Flag( "moving obstacle SBB" );
const FlagUID MO_CLI_Flag( "moving obstacle CLI" );


/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, uint_t levels, const AABB& refinementBox )
{
   real_t dx = real_t(1); // dx on finest level
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      uint_t blockLevel = block->getLevel();
      uint_t levelScalingFactor = ( uint_t(1) << (levels - uint_t(1) - blockLevel) );
      real_t dxOnLevel = dx * real_c(levelScalingFactor);
      AABB blockAABB = block->getAABB();

      // extend block AABB by ghostlayers
      AABB extendedBlockAABB = blockAABB.getExtended( dxOnLevel * real_c(FieldGhostLayers) );

      if( extendedBlockAABB.intersects( refinementBox ) )
         if( blockLevel < ( levels - uint_t(1) ) )
            block->setMarker( true );
   }
}

static void workloadAndMemoryAssignment( SetupBlockForest& forest )
{
   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      block->setWorkload( numeric_cast< workload_t >( uint_t(1) << block->getLevel() ) );
      block->setMemory( numeric_cast< memory_t >(1) );
   }
}

static shared_ptr< StructuredBlockForest > createBlockStructure( const AABB & domainAABB, Vector3<uint_t> blockSizeInCells,
                                                                 uint_t numberOfLevels, real_t diameter, Vector3<real_t> initialPosition,
                                                                 bool useLargeRefinementRegion,
                                                                 bool keepGlobalBlockInformation = false )
{
   SetupBlockForest sforest;

   Vector3<uint_t> numberOfFineBlocksPerDirection( uint_c(domainAABB.size(0)) / blockSizeInCells[0],
                                                   uint_c(domainAABB.size(1)) / blockSizeInCells[1],
                                                   uint_c(domainAABB.size(2)) / blockSizeInCells[2] );

   for(uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_CHECK_EQUAL( numberOfFineBlocksPerDirection[i] * blockSizeInCells[i], uint_c(domainAABB.size(i)),
                            "Domain can not be decomposed in direction " << i << " into fine blocks of size " << blockSizeInCells[i] );
   }

   uint_t levelScalingFactor = ( uint_t(1) << ( numberOfLevels - uint_t(1) ) );
   Vector3<uint_t> numberOfCoarseBlocksPerDirection( numberOfFineBlocksPerDirection / levelScalingFactor );

   for(uint_t i = 0; i < 3; ++i )
   {
      WALBERLA_CHECK_EQUAL(numberOfCoarseBlocksPerDirection[i] * levelScalingFactor, numberOfFineBlocksPerDirection[i],
                            "Domain can not be refined in direction " << i << " according to the specified number of levels!" );
   }

   AABB refinementBox;
   if( useLargeRefinementRegion ) {
      // refinement area is the complete area along the bottom plane, containing the sphere
      // this avoids refinement borders in flow direction
      refinementBox = AABB(domainAABB.xMin(), domainAABB.yMin(), domainAABB.zMin(),
                           domainAABB.xMax(), domainAABB.yMax(), initialPosition[2] + real_t(0.5) * diameter);
   } else{
      // refinement area is just around the sphere
      refinementBox = AABB(initialPosition[0] - real_t(0.5) * diameter, initialPosition[1] - real_t(0.5) * diameter, domainAABB.zMin(),
                           initialPosition[0] + real_t(0.5) * diameter, initialPosition[1] + real_t(0.5) * diameter, initialPosition[2] + real_t(0.5) * diameter);
   }

   WALBERLA_LOG_INFO_ON_ROOT(" - refinement box: " << refinementBox);

   sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, numberOfLevels, refinementBox ) );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   sforest.init( domainAABB, numberOfCoarseBlocksPerDirection[0], numberOfCoarseBlocksPerDirection[1], numberOfCoarseBlocksPerDirection[2], true, true, false );

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

   WALBERLA_LOG_INFO_ON_ROOT( sforest );

   MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   shared_ptr< StructuredBlockForest > sbf =
         make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation ),
                                               blockSizeInCells[0], blockSizeInCells[1], blockSizeInCells[2]);
   sbf->createCellBoundingBoxes();

   return sbf;
}

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
                                                           MO_SBB_T( "MO_SBB", MO_SBB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                           MO_CLI_T( "MO_CLI", MO_CLI_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   // boundary conditions are set by mapping the (moving) planes into the domain

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}
//*******************************************************************************************************************


class SpherePropertyLogger
{
public:
   SpherePropertyLogger( const shared_ptr< StructuredBlockStorage > & blocks,
                         const BlockDataID & bodyStorageID,
                         const std::string & fileName, bool fileIO,
                         real_t dragNormalizationFactor, real_t liftNormalizationFactor, real_t physicalTimeScale ) :
      blocks_( blocks ), bodyStorageID_( bodyStorageID ),
      fileName_( fileName ), fileIO_(fileIO),
      dragNormalizationFactor_( dragNormalizationFactor ), liftNormalizationFactor_( liftNormalizationFactor ),
      physicalTimeScale_( physicalTimeScale ),
      counter_( 0 )
   {
      if ( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( fileName_.c_str() );
            file << "#\t Cd\t Cl\t fX\t fY\t fZ\t tX\t tY\t tZ\n";
            file.close();
         }
      }
   }

   void operator()()
   {

      Vector3<real_t> force(real_t(0));
      Vector3<real_t> torque(real_t(0));

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end(); ++bodyIt )
         {
            force += bodyIt->getForce();
            torque += bodyIt->getTorque();
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( force[0], mpi::SUM );
         mpi::allReduceInplace( force[1], mpi::SUM );
         mpi::allReduceInplace( force[2], mpi::SUM );

         mpi::allReduceInplace( torque[0], mpi::SUM );
         mpi::allReduceInplace( torque[1], mpi::SUM );
         mpi::allReduceInplace( torque[2], mpi::SUM );
      }

      if( fileIO_ )
         writeToFile( counter_, force, torque);

      dragForce_ = force[0];
      liftForce_ = force[2];

      ++counter_;

   }

   real_t getDragForce()
   {
      return dragForce_;
   }

   real_t getLiftForce()
   {
      return liftForce_;
   }

   real_t getDragCoefficient()
   {
      return dragForce_ / dragNormalizationFactor_;
   }

   real_t getLiftCoefficient()
   {
      return liftForce_ / liftNormalizationFactor_;
   }


private:
   void writeToFile( uint_t timestep, const Vector3<real_t> & force, const Vector3<real_t> & torque )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName_.c_str(), std::ofstream::app );

         file << real_c(timestep) / physicalTimeScale_
              << "\t" << force[0] / dragNormalizationFactor_ << "\t" << force[2] / liftNormalizationFactor_
              << "\t" << force[0] << "\t" << force[1] << "\t" << force[2]
              << "\t" << torque[0] << "\t" << torque[1] << "\t" << torque[2]
              << "\n";
         file.close();
      }
   }

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID bodyStorageID_;
   std::string fileName_;
   bool fileIO_;
   real_t dragForce_, liftForce_;
   real_t dragNormalizationFactor_, liftNormalizationFactor_;
   real_t physicalTimeScale_;
   uint_t counter_;

};

void initializeCouetteProfile( const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & pdfFieldID, const BlockDataID & boundaryHandlingID,
                               const real_t & domainHeight, const real_t wallVelocity )
{
   const real_t rho = real_c(1);

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto pdfField = blockIt->getData< PdfField_T > ( pdfFieldID );
      BoundaryHandling_T * boundaryHandling = blockIt->getData< BoundaryHandling_T >( boundaryHandlingID );

      WALBERLA_FOR_ALL_CELLS_XYZ(pdfField,
                                 if( !boundaryHandling->isDomain(x,y,z) ) continue;

                                 const Vector3< real_t > coord = blocks->getBlockLocalCellCenter( *blockIt, Cell(x,y,z) );

                                 Vector3< real_t > velocity( real_c(0) );

                                 velocity[0] =  wallVelocity * coord[2] / domainHeight;

                                 pdfField->setToEquilibrium( x, y, z, velocity, rho );
      )
   }
}

void logFinalResult( const std::string & fileName, real_t Re, real_t wallDistance, real_t diameter, uint_t numLevels,
                     uint_t domainLength, uint_t domainWidth, uint_t domainHeight,
                     real_t dragCoeff, real_t dragCoeffRef,
                     real_t liftCoeff, real_t liftCoeffRef,
                     uint_t timestep, real_t nonDimTimestep )
{
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open( fileName.c_str(), std::ofstream::app );

      file << Re << "\t " << wallDistance << "\t " << diameter << "\t " << numLevels << "\t "
           << domainLength << "\t " << domainWidth << "\t " << domainHeight << "\t "
           << dragCoeff << "\t " << dragCoeffRef << "\t "<< std::abs(dragCoeff-dragCoeffRef) / dragCoeffRef << "\t "
           << liftCoeff << "\t " << liftCoeffRef << "\t "<< std::abs(liftCoeff-liftCoeffRef) / liftCoeffRef << "\t "
           << timestep << "\t " << nonDimTimestep
           << "\n";
      file.close();
   }
}


//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Testcase that evaluates the drag and lift force on a sphere that is close to the bottom plane in shear flow
 *
 * see overview paper:
 * Zeng, Najjar, Balachandar, Fischer - "Forces on a finite-sized particle located close to a wall in a linear shear flow", 2009
 *
 * contains references to:
 * Leighton, Acrivos - "The lift on a small sphere touching a plane in the presence of a simple shear flow", 1985
 * Zeng, Balachandar, Fischer - "Wall-induced forces on a rigid sphere at finite Reynolds number", 2005
 *
 * CFD-IBM simulations in:
 * Lee, Balachandar - "Drag and lift forces on a spherical particle moving on a wall in a shear flow at finite Re", 2010
 *
 * Description:
 *  - Domain size [x, y, z] = [48 x 16 x 8 ] * diameter = [L(ength), W(idth), H(eight)]
 *  - horizontally periodic, bounded by two planes in z-direction
 *  - top plane is moving with constant wall velocity -> shear flow
 *  - sphere is placed in the vicinity of the bottom plane at [ L/2 + xOffset, W/2 + yOffset, dist * diameter]
 *  - distance of sphere center to the bottom plane is crucial parameter
 *  - viscosity is adjusted to match specified Reynolds number = shearRate * diameter * wallDistance / viscosity
 *  - dimensionless drag and lift forces are evaluated and written to logging file
 *  - different variants of grid refinement can be chosen (number of levels, refinement region)
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
   bool fileIO = false;
   bool zeroShearTest = false; // shear rate is zero such that a quiescent fluid is obtained -> drag and lift have to be and remain zero -> ignores Reynolds number
   uint_t vtkIOFreq = 0;
   std::string baseFolderVTK = "vtk_out_ForcesNearPlane";
   std::string baseFolderLogging = ".";

   // physical setup
   real_t diameter = real_t(20); // cells per diameter -> determines overall resolution
   real_t normalizedWallDistance = real_t(1); // distance of the sphere center to the bottom wall, normalized by the diameter
   real_t ReynoldsNumberShear = real_t(1); // = shearRate * wallDistance * diameter / viscosity

   //numerical parameters
   real_t minimumNonDimTimesteps = real_t(100); // minimum number of non-dimensional time steps before simulation can be terminated by convergence
   uint_t numberOfLevels = uint_t(4); // number of grid levels for static refinement ( 1 = no refinement)
   real_t xOffsetOfSpherePosition = real_t(0); // offset in x-direction of sphere position
   real_t yOffsetOfSpherePosition = real_t(0); // offset in y-direction of sphere position
   bool useSBB = false; // use simple bounce-back boundary condition for sphere
   bool useLargeRefinementRegion = false; // uses the whole area near the bottom plane as the finest grid, else refinement is only applied around the sphere


   // command line arguments
   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--zeroShearTest" )    == 0 ) { zeroShearTest = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )           == 0 ) { fileIO = true; continue; }
      if( std::strcmp( argv[i], "--vtkIOFreq" )        == 0 ) { vtkIOFreq = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--diameter" )         == 0 ) { diameter = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--numLevels" )        == 0 ) { numberOfLevels = uint_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--timesteps" )        == 0 ) { minimumNonDimTimesteps = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--wallDistance" )     == 0 ) { normalizedWallDistance = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--Re" )               == 0 ) { ReynoldsNumberShear = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--xOffset" )          == 0 ) { xOffsetOfSpherePosition = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--yOffset" )          == 0 ) { yOffsetOfSpherePosition = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--useSBB" )           == 0 ) { useSBB = true; continue; }
      if( std::strcmp( argv[i], "--useLargeRefinementRegion" ) == 0 ) { useLargeRefinementRegion = true; continue; }
      if( std::strcmp( argv[i], "--baseFolderVTK" ) == 0 ) { baseFolderVTK = argv[++i]; continue; }
      if( std::strcmp( argv[i], "--baseFolderLogging" ) == 0 ) { baseFolderVTK = argv[++i]; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   WALBERLA_CHECK_GREATER_EQUAL(normalizedWallDistance, real_t(0.5));
   WALBERLA_CHECK_GREATER_EQUAL(ReynoldsNumberShear, real_t(0));
   WALBERLA_CHECK_GREATER_EQUAL(diameter, real_t(0));

   //////////////////////////
   // NUMERICAL PARAMETERS //
   //////////////////////////

   const real_t domainLength = real_t(48) * diameter; //x
   const real_t domainWidth  = real_t(16) * diameter; //y
   const real_t domainHeight = real_t( 8) * diameter; //z

   Vector3<uint_t> domainSize( uint_c( std::ceil(domainLength)), uint_c( std::ceil(domainWidth)), uint_c( std::ceil(domainHeight)) );

   real_t wallVelocity = real_t(0.01);
   if( zeroShearTest ) wallVelocity = real_t(0);

   const real_t wallDistance = diameter * normalizedWallDistance;
   const real_t shearRate = wallVelocity / domainHeight;
   const real_t velAtSpherePosition = shearRate * wallDistance;
   const real_t viscosity = ( zeroShearTest ) ? real_t(0.015) : ( velAtSpherePosition * diameter / ReynoldsNumberShear );

   const real_t relaxationTime = real_t(1) / lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t densityFluid = real_t(1);

   const real_t dx = real_t(1);

   const real_t physicalTimeScale = ( zeroShearTest ) ? real_t(10) : (diameter / velAtSpherePosition);
   const uint_t minimumLBMtimesteps = uint_c(minimumNonDimTimesteps * physicalTimeScale);

   const real_t omega = real_t(1) / relaxationTime;
   const uint_t finestLevel = numberOfLevels - uint_t(1);
   Vector3<real_t> initialPosition( domainLength * real_t(0.5) + xOffsetOfSpherePosition,
                                    domainWidth * real_t(0.5) + yOffsetOfSpherePosition,
                                    wallDistance );

   WALBERLA_LOG_INFO_ON_ROOT("Setup:");
   WALBERLA_LOG_INFO_ON_ROOT(" - domain size = " << domainSize );
   WALBERLA_LOG_INFO_ON_ROOT(" - normalized wall distance = " << normalizedWallDistance );
   WALBERLA_LOG_INFO_ON_ROOT(" - shear rate = " << shearRate );
   WALBERLA_LOG_INFO_ON_ROOT(" - wall velocity = " << wallVelocity );
   WALBERLA_LOG_INFO_ON_ROOT(" - Reynolds number (shear rate based) = " << ReynoldsNumberShear << ", vel at sphere pos = " << velAtSpherePosition);
   WALBERLA_LOG_INFO_ON_ROOT(" - density = " << densityFluid  );
   std::stringstream omega_msg;
   for( uint_t i = 0; i < numberOfLevels; ++i )
   {
      omega_msg << lbm::collision_model::levelDependentRelaxationParameter( i, omega, finestLevel ) << " ( on level " << i << " ), ";
   }
   WALBERLA_LOG_INFO_ON_ROOT(" - viscosity = " << viscosity << " -> omega = " << omega_msg.str());
   WALBERLA_LOG_INFO_ON_ROOT(" - sphere diameter = " << diameter << ", position = " << initialPosition << " ( xOffset = " << xOffsetOfSpherePosition << ", yOffset = " << yOffsetOfSpherePosition << " )");
   WALBERLA_LOG_INFO_ON_ROOT(" - using " << ((useSBB) ? "SBB" : "CLI") << " boundary condition along sphere");
   WALBERLA_LOG_INFO_ON_ROOT(" - base folder VTK = " << baseFolderVTK << ", base folder logging = " << baseFolderLogging );
   if( zeroShearTest ) WALBERLA_LOG_INFO_ON_ROOT(" === ATTENTION: USING zeroShearTest SETUP -> MANY PHYSICAL PARAMETERS ARE IGNORED AND THUS MIGHT BE REPORTED INCORRECTLY! === ") ;

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t levelScalingFactor = ( uint_t(1) << finestLevel );
   const uint_t lbmTimeStepsPerTimeLoopIteration = levelScalingFactor;

   const uint_t timesteps = uint_c( minimumLBMtimesteps / lbmTimeStepsPerTimeLoopIteration );

   Vector3<uint_t> coarseBlocksPerDirection( 6, 2, 1 );
   if( numberOfLevels == 1 || numberOfLevels == 2 )
   {
      // have enough coarse blocks to run on cluster
      coarseBlocksPerDirection = Vector3<uint_t>(12,4,2);
   }
   WALBERLA_CHECK( domainSize[0] % coarseBlocksPerDirection[0] == 0 &&
                   domainSize[1] % coarseBlocksPerDirection[1] == 0 &&
                   domainSize[2] % coarseBlocksPerDirection[2] == 0 );
   Vector3<uint_t> blockSizeInCells(domainSize[0] / ( coarseBlocksPerDirection[0] * levelScalingFactor ),
                                    domainSize[1] / ( coarseBlocksPerDirection[1] * levelScalingFactor ),
                                    domainSize[2] / ( coarseBlocksPerDirection[2] * levelScalingFactor ) );

   AABB simulationDomain( real_t(0), real_t(0), real_t(0), real_c(domainSize[0]), real_c(domainSize[1]), real_c(domainSize[2]) );
   auto blocks = createBlockStructure( simulationDomain, blockSizeInCells, numberOfLevels, diameter, initialPosition, useLargeRefinementRegion );

   //write domain decomposition to file
   if( vtkIOFreq > 0 )
   {
      vtk::writeDomainDecomposition( blocks, "initial_domain_decomposition", baseFolderVTK );
   }


   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage>  globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");

   // set up synchronization procedure
   const real_t overlap = real_t( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // create pe bodies

   // bounding planes (global)
   const auto planeMaterial = pe::createMaterial( "myPlaneMat", real_t(8920), real_t(0), real_t(1), real_t(1), real_t(0), real_t(1), real_t(1), real_t(0), real_t(0) );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), planeMaterial );
   auto topPlane = pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,domainHeight), planeMaterial );
   topPlane->setLinearVel(wallVelocity, real_t(0), real_t(0));

   // add the sphere
   pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, initialPosition, real_t(0.5) * diameter, planeMaterial );

   uint_t minBlockSizeInCells = blockSizeInCells.min();
   for( uint_t i = 0; i < uint_c(diameter / real_c(minBlockSizeInCells)) + 1; ++i)
      syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega, lbm::collision_model::TRT::threeSixteenth, finestLevel ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         FieldGhostLayers, field::zyxf );
   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf, FieldGhostLayers );

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );


   // map planes into the LBM simulation -> act as no-slip boundaries
   // since we have refinement boundaries along the bottom plane, we use SBB here
   // (see test/pe_coupling/momentum_exchange_method/GlobalBodyAsBoundaryMEMStaticRefinement test for more infos)
   pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_SBB_Flag, pe_coupling::selectGlobalBodies );

   // map movable sphere into the LBM simulation
   // the whole sphere resides on the same refinement level, so we can also use CLI instead of SBB
   if( useSBB )
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_SBB_Flag, pe_coupling::selectRegularBodies );
   else
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_CLI_Flag, pe_coupling::selectRegularBodies );


   // initialize Couette velocity profile in whole domain
   if( !zeroShearTest ) initializeCouetteProfile(blocks, pdfFieldID, boundaryHandlingID, domainHeight, wallVelocity);

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

   if( vtkIOFreq != uint_t(0) )
   {
      // flag field (written only once in the first time step, ghost layers are also written)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", vtkIOFreq, uint_t(0), false, baseFolderVTK );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", vtkIOFreq, uint_t(0), false, baseFolderVTK );

      field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
      fluidFilter.addFlag( Fluid_Flag );
      pdfFieldVTK->addCellInclusionFilter( fluidFilter );

      pdfFieldVTK->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "VelocityFromPDF" ) );
      pdfFieldVTK->addCellDataWriter( make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "DensityFromPDF" ) );

      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTK ), "VTK (fluid field data)" );
   }

   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, sweep, pdfFieldID, boundaryHandlingID );

   // add force evaluation and logging
   real_t normalizationFactor = ( zeroShearTest ) ? real_t(1) : ( math::M_PI / real_t(8) * densityFluid * shearRate * shearRate * wallDistance * wallDistance * diameter * diameter );
   std::string loggingFileName( baseFolderLogging + "/LoggingForcesNearPlane");
   loggingFileName += "_lvl" + std::to_string(numberOfLevels);
   loggingFileName += "_D" + std::to_string(uint_c(diameter));
   loggingFileName += "_Re" + std::to_string(uint_c(ReynoldsNumberShear));
   loggingFileName += "_WD" + std::to_string(uint_c(normalizedWallDistance * real_t(1000)));
   loggingFileName += ".txt";
   shared_ptr< SpherePropertyLogger > logger = walberla::make_shared< SpherePropertyLogger >( blocks, bodyStorageID,
                                                                                              loggingFileName, fileIO,
                                                                                              normalizationFactor, normalizationFactor, physicalTimeScale);
   refinementTimestep->addPostStreamVoidFunction( lbm::refinement::FunctorWrapper( SharedFunctor< SpherePropertyLogger >( logger ) ), "Sphere property logger", finestLevel );

   // add pe timesteps
   refinementTimestep->addPostStreamVoidFunction( lbm::refinement::FunctorWrapper(pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID )),
                                                  "force reset", finestLevel );

   // add LBM sweep with refinement
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( refinementTimestep ), "LBM refinement time step" );


   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   // compute reference values from literature

   const real_t normalizedGapSize = normalizedWallDistance - real_t(0.5);

   // drag correlation for the drag coefficient
   const real_t standardDragCorrelation = real_t(24) / ReynoldsNumberShear * (real_t(1) + real_t(0.15) * std::pow(ReynoldsNumberShear, real_t(0.687))); // Schiller-Naumann correlation
   const real_t dragCorrelationWithGapSizeStokes = real_t(24) / ReynoldsNumberShear * (real_t(1) + real_t(0.138) * std::exp(real_t(-2) * normalizedGapSize) + real_t(9)/( real_t(16) * (real_t(1) + real_t(2) * normalizedGapSize) ) ); // Goldman et al. (1967)
   const real_t alphaDragS = real_t(0.15) - real_t(0.046) * ( real_t(1) - real_t(0.16) * normalizedGapSize * normalizedGapSize ) * std::exp( -real_t(0.7) *  normalizedGapSize);
   const real_t betaDragS = real_t(0.687) + real_t(0.066)*(real_t(1) - real_t(0.76) * normalizedGapSize * normalizedGapSize) * std::exp( - std::pow( normalizedGapSize, real_t(0.9) ) );
   const real_t dragCorrelationZeng = dragCorrelationWithGapSizeStokes * ( real_t(1) + alphaDragS * std::pow( ReynoldsNumberShear, betaDragS ) ); // Zeng et al. (2009) - Eqs. (13) and (14)

   // lift correlations for the lift coefficient
   const real_t liftCorrelationZeroGapStokes = real_t(5.87); // Leighton, Acrivos (1985)
   const real_t liftCorrelationZeroGap = real_t(3.663) / std::pow( ReynoldsNumberShear * ReynoldsNumberShear + real_t(0.1173), real_t(0.22) ); //  Zeng et al. (2009) - Eq. (19)
   const real_t alphaLiftS = - std::exp( -real_t(0.3) + real_t(0.025) * ReynoldsNumberShear);
   const real_t betaLiftS = real_t(0.8) + real_t(0.01) * ReynoldsNumberShear;
   const real_t lambdaLiftS = ( real_t(1) - std::exp(-normalizedGapSize)) * std::pow( ReynoldsNumberShear / real_t(250), real_t(5) / real_t(2) );
   const real_t liftCorrelationZeng = liftCorrelationZeroGap * std::exp( - real_t(0.5) * normalizedGapSize * std::pow( ReynoldsNumberShear / real_t(250), real_t(4)/real_t(3))) *
                                      ( std::exp( alphaLiftS * std::pow( normalizedGapSize, betaLiftS ) ) - lambdaLiftS ); // Zeng et al. (2009) - Eqs. (28) and (29)

   ////////////////////////
   // EXECUTE SIMULATION //
   ////////////////////////

   WcTimingPool timeloopTiming;

   const real_t relativeChangeConvergenceEps = real_t(1e-3);
   const real_t physicalCheckingFrequency = real_t(0.00625);
   const uint_t checkingFrequency = (zeroShearTest) ? uint_t(1) : uint_c( physicalCheckingFrequency * physicalTimeScale );

   WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with at least " << timesteps << " (coarse) time steps");
   WALBERLA_LOG_INFO_ON_ROOT("Convergence checking frequency = " << checkingFrequency << " (physically = " << physicalCheckingFrequency << ")");

   real_t maxDragCurrentCheckingPeriod = -math::Limits<real_t>::inf();
   real_t minDragCurrentCheckingPeriod = math::Limits<real_t>::inf();
   real_t maxLiftCurrentCheckingPeriod = -math::Limits<real_t>::inf();
   real_t minLiftCurrentCheckingPeriod = math::Limits<real_t>::inf();
   uint_t timestep = 0;

   // time loop
   while( true )
   {
      // perform a single simulation step
      timeloop.singleStep( timeloopTiming );

      real_t curDrag = logger->getDragCoefficient();
      real_t curLift = logger->getLiftCoefficient();

      maxDragCurrentCheckingPeriod = std::max( maxDragCurrentCheckingPeriod, curDrag);
      minDragCurrentCheckingPeriod = std::min( minDragCurrentCheckingPeriod, curDrag);
      maxLiftCurrentCheckingPeriod = std::max( maxLiftCurrentCheckingPeriod, curLift);
      minLiftCurrentCheckingPeriod = std::min( minLiftCurrentCheckingPeriod, curLift);

      real_t dragDiffCurrentCheckingPeriod = std::fabs(maxDragCurrentCheckingPeriod - minDragCurrentCheckingPeriod)/std::fabs(maxDragCurrentCheckingPeriod);
      real_t liftDiffCurrentCheckingPeriod = std::fabs(maxLiftCurrentCheckingPeriod - minLiftCurrentCheckingPeriod)/std::fabs(maxLiftCurrentCheckingPeriod);

      // continuous output during simulation
      if( timestep % (checkingFrequency * uint_t(10)) == 0)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Drag: current C_D = " << curDrag );
         if( !zeroShearTest )
         {
            WALBERLA_LOG_INFO_ON_ROOT(" - standard C_D = " <<  standardDragCorrelation );
            WALBERLA_LOG_INFO_ON_ROOT(" - C_D ( Stokes fit ) = " << dragCorrelationWithGapSizeStokes );
            WALBERLA_LOG_INFO_ON_ROOT(" - C_D ( Zeng ) = " << dragCorrelationZeng );
         }

         WALBERLA_LOG_INFO_ON_ROOT("Lift: current C_L = " << curLift );
         if( !zeroShearTest )
         {
            WALBERLA_LOG_INFO_ON_ROOT(" - C_L ( Stokes, zero gap ) = " << liftCorrelationZeroGapStokes);
            WALBERLA_LOG_INFO_ON_ROOT(" - C_L ( zero gap ) = " << liftCorrelationZeroGap);
            WALBERLA_LOG_INFO_ON_ROOT(" - C_L ( Zeng ) = " << liftCorrelationZeng);
            WALBERLA_LOG_INFO_ON_ROOT( "Drag difference [(max-min)/max] = " << dragDiffCurrentCheckingPeriod << ", lift difference = " << liftDiffCurrentCheckingPeriod);
         }
      }

      // check for convergence ( = difference between min and max values in current checking period is below limit)
      if( timestep % checkingFrequency == 0 )
      {
         if( timestep > timesteps &&
             dragDiffCurrentCheckingPeriod < relativeChangeConvergenceEps &&
             liftDiffCurrentCheckingPeriod < relativeChangeConvergenceEps )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Forces converged with an eps of " << relativeChangeConvergenceEps );
            WALBERLA_LOG_INFO_ON_ROOT(" - drag min = " << minDragCurrentCheckingPeriod << " , max = " << maxDragCurrentCheckingPeriod);
            WALBERLA_LOG_INFO_ON_ROOT(" - lift min = " << minLiftCurrentCheckingPeriod << " , max = " << maxLiftCurrentCheckingPeriod);
            break;
         }

         //reset min and max values for new checking period
         maxDragCurrentCheckingPeriod = -math::Limits<real_t>::inf();
         minDragCurrentCheckingPeriod = math::Limits<real_t>::inf();
         maxLiftCurrentCheckingPeriod = -math::Limits<real_t>::inf();
         minLiftCurrentCheckingPeriod = math::Limits<real_t>::inf();
      }

      if( zeroShearTest && timestep > timesteps ) break;

      ++timestep;

   }

   timeloopTiming.logResultOnRoot();

   if ( !zeroShearTest )
   {
      std::string resultFileName( baseFolderLogging + "/ResultForcesNearPlane.txt");
      logFinalResult(resultFileName, ReynoldsNumberShear, normalizedWallDistance, diameter, numberOfLevels, domainSize[0], domainSize[1], domainSize[2],
                     logger->getDragCoefficient(), dragCorrelationZeng, logger->getLiftCoefficient(), liftCorrelationZeng, timestep, real_c(timestep) / physicalTimeScale );
   }

   return EXIT_SUCCESS;
}

} // namespace forces_on_sphere_near_plane_in_shear_flow

int main( int argc, char **argv ){
   forces_on_sphere_near_plane_in_shear_flow::main(argc, argv);
}
