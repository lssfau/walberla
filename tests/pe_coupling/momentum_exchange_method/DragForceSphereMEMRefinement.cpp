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
//! \file DragForceSphereMEMRefinement.cpp
//! \ingroup pe_coupling
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
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/all.h"
#include "core/SharedFunctor.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/MacroscopicValueCalculation.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/refinement/all.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>
#include <vector>
#include <iomanip>
#include <iostream>

namespace drag_force_sphere_mem_refinement
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
typedef lbm::D3Q19< lbm::collision_model::TRT, false, lbm::force_model::SimpleConstant >  LatticeModel_T;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

const uint_t FieldGhostLayers = 4;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

// boundary handling
typedef pe_coupling::SimpleBB       < LatticeModel_T, FlagField_T >  MO_BB_T;
typedef pe_coupling::CurvedLinear   < LatticeModel_T, FlagField_T > MO_CLI_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, MO_BB_T, MO_CLI_T > BoundaryHandling_T;

using BodyTypeTuple = std::tuple<pe::Sphere> ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID MO_BB_Flag   ( "moving obstacle BB" );
const FlagUID MO_CLI_Flag  ( "moving obstacle CLI" );


////////////////
// PARAMETERS //
////////////////

struct Setup{
   uint_t checkFrequency;
   real_t visc;
   real_t tau;
   real_t radius;
   real_t length;
   real_t chi;
   real_t extForce;
   real_t analyticalDrag;
   uint_t levels;
};


/////////////////////
// BLOCK STRUCTURE //
/////////////////////

static void refinementSelection( SetupBlockForest& forest, const uint_t levels, const real_t diameter, const real_t domainLength )
{
   const real_t lowerCorner = real_c(0.5) * ( domainLength - diameter );

   AABB refinementBox( lowerCorner, lowerCorner, lowerCorner,
                       lowerCorner + diameter, lowerCorner + diameter, lowerCorner + diameter );

   for( auto block = forest.begin(); block != forest.end(); ++block )
   {
      if( block->getAABB().intersects( refinementBox ) )
         if( block->getLevel() < ( levels - uint_t(1) ) )
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

static shared_ptr< StructuredBlockForest > createBlockStructure( const Setup & setup,
                                                                 const uint_t numberOfXBlocks, const uint_t numberOfYBlocks, const uint_t numberOfZBlocks,
                                                                 const bool xPeriodic = false, const bool yPeriodic = false, const bool zPeriodic = false,
                                                                 const bool keepGlobalBlockInformation = false )
{
   // initialize SetupBlockForest = determine domain decomposition
   SetupBlockForest sforest;

   sforest.addRefinementSelectionFunction( std::bind( refinementSelection, std::placeholders::_1, setup.levels, setup.radius * real_c(2), setup.length ) );
   sforest.addWorkloadMemorySUIDAssignmentFunction( workloadAndMemoryAssignment );

   sforest.init( AABB( real_c(0), real_c(0), real_c(0), setup.length, setup.length, setup.length ),
                 numberOfXBlocks, numberOfYBlocks, numberOfZBlocks, xPeriodic, yPeriodic, zPeriodic );

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   sforest.balanceLoad( blockforest::StaticLevelwiseCurveBalance(true), uint_c( MPIManager::instance()->numProcesses() ), real_t(0), memoryLimit, true );

   WALBERLA_LOG_INFO_ON_ROOT( sforest );

   MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   uint_t numberOfXCellsPerBlock = uint_c( setup.length / real_c(numberOfXBlocks) );
   uint_t numberOfYCellsPerBlock = uint_c( setup.length / real_c(numberOfYBlocks) );
   uint_t numberOfZCellsPerBlock = uint_c( setup.length / real_c(numberOfZBlocks) );
   shared_ptr< StructuredBlockForest > sbf =
         make_shared< StructuredBlockForest >( make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ), sforest, keepGlobalBlockInformation ),
                                               numberOfXCellsPerBlock, numberOfYCellsPerBlock, numberOfZCellsPerBlock );
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

   BoundaryHandling_T * handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                           MO_BB_T (  "MO_BB",  MO_BB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                           MO_CLI_T ( "MO_CLI", MO_CLI_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


//*******************************************************************************************************************
/*!\brief Evaluating the drag force on a fixed sphere
 */
//*******************************************************************************************************************
class ForceEval
{
   public:
      ForceEval( SweepTimeloop* timeloop, Setup* setup, const shared_ptr< StructuredBlockStorage > & blocks,
                 const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & bodyStorageID,
                 bool fileIO ) :
         timeloop_( timeloop ), setup_( setup ), blocks_( blocks ),
         flagFieldID_ (flagFieldID), pdfFieldID_ (pdfFieldID), bodyStorageID_( bodyStorageID ),
         fileIO_(fileIO), normalizedDragOld_(0.0), normalizedDragNew_(0.0)
      {
         // calculate the analytical drag force value based on the series expansion of chi
         // see also Sangani - Slow flow through a periodic array of spheres, IJMF 1982. Eq. 60 and Table 1
         real_t analyticalDrag = real_c(0);
         real_t tempChiPowS = real_c(1);

         // coefficients to calculate the drag in a series expansion
         real_t dragCoefficients[31] = {  real_c(1.000000),  real_c(1.418649),  real_c(2.012564),  real_c(2.331523),  real_c(2.564809),   real_c(2.584787),
                                          real_c(2.873609),  real_c(3.340163),  real_c(3.536763),  real_c(3.504092),  real_c(3.253622),   real_c(2.689757),
                                          real_c(2.037769),  real_c(1.809341),  real_c(1.877347),  real_c(1.534685), real_c(0.9034708),  real_c(0.2857896),
                                        real_c(-0.5512626), real_c(-1.278724),  real_c(1.013350),  real_c(5.492491),  real_c(4.615388), real_c(-0.5736023),
                                         real_c(-2.865924), real_c(-4.709215), real_c(-6.870076), real_c(0.1455304),  real_c(12.51891),   real_c(9.742811),
                                         real_c(-5.566269)};

         for(uint_t s = 0; s <= uint_t(30); ++s){
            analyticalDrag += dragCoefficients[s] * tempChiPowS;
            tempChiPowS *= setup->chi;
         }
         setup_->analyticalDrag = analyticalDrag;

         if (fileIO_)
         {
            std::ofstream file;
            filename_ = "DragForceSphereRefinement.txt";
            WALBERLA_ROOT_SECTION()
            {
               file.open( filename_.c_str() );
               file << "#\t forceX\t uBar\t dragSim\t dragAnaly\n";
               file.close();
            }
         }

         // when evaluating the force after the coarse level step, it has already been applied several times on the finest levels
         forceScalingFactor_ = real_c(1) / real_c( uint_t(1) << ( setup->levels - uint_t(1) ) );
      }

      // evaluate the acting drag force
      void operator()()
      {
         const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

         if( timestep % setup_->checkFrequency != 0) return;

         // get force in x-direction acting on the sphere
         real_t forceX = forceScalingFactor_ * getForce();
         // get average volumetric flowrate in the domain
         real_t uBar = getAverageVel();

         // f_total = f_drag + f_buoyancy
         real_t totalForce = forceX  + real_c(4.0/3.0) * math::pi * setup_->radius * setup_->radius * setup_->radius * setup_->extForce ;

         real_t normalizedDragForce = totalForce / real_c( 6.0 * math::pi * setup_->visc * setup_->radius * uBar );

         // update drag force values
         normalizedDragOld_ = normalizedDragNew_;
         normalizedDragNew_ = normalizedDragForce;

         // write to file if desired
         WALBERLA_ROOT_SECTION()
         {
            if (fileIO_)
            {
               std::ofstream file;
               file.open( filename_.c_str(), std::ofstream::app );
               file.setf( std::ios::unitbuf );
               file.precision(15);
               file << timestep << " " << forceX << " " << uBar << " " << normalizedDragNew_ << " " << setup_->analyticalDrag << "\n";
               file.close();
            }
         }
      }

      // obtain the drag force acting on the sphere by summing up all the process local parts of fX
      real_t getForce()
      {
         real_t forceX = real_c( 0.0 );

         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_ ); bodyIt != pe::BodyIterator::end(); ++bodyIt)
            {
               forceX += bodyIt->getForce()[0];
            }
         }
         WALBERLA_MPI_SECTION()
         {
            mpi::allReduceInplace( forceX, mpi::SUM );
         }

         return forceX;
      }

      // calculate the average velocity in forcing direction (here: x) inside the domain
      real_t getAverageVel()
      {
         real_t velocity_sum = real_t(0);
         // iterate all blocks stored locally on this process
         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            real_t dx = blocks_->dx(blocks_->getLevel(*blockIt));
            real_t volume = dx * dx * dx;
            // retrieve the pdf field and the flag field from the block
            PdfField_T  *  pdfField = blockIt->getData< PdfField_T >( pdfFieldID_ );
            FlagField_T * flagField = blockIt->getData< FlagField_T >( flagFieldID_ );

            // get the flag that marks a cell as being fluid
            auto fluid = flagField->getFlag( Fluid_Flag );

            auto xyzFieldSize = pdfField->xyzSize();
            for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
            {
               Cell cell = *fieldIt;
               if( flagField->isFlagSet( cell.x(), cell.y(), cell.z(), fluid ) )
               {
                velocity_sum += volume * pdfField->getVelocity(cell)[0];
               }
            }
         }

         WALBERLA_MPI_SECTION()
         {
            mpi::allReduceInplace( velocity_sum, mpi::SUM );
         }

         return velocity_sum / ( setup_->length * setup_->length * setup_->length );
      }

      // return the relative temporal change in the normalized drag
      real_t getDragForceDiff() const
      {
         return std::fabs( ( normalizedDragNew_ - normalizedDragOld_ ) / normalizedDragNew_ );
      }

      // return the drag force
      real_t getDragForce() const
      {
         return normalizedDragNew_;
      }

   private:
      SweepTimeloop* timeloop_;

      Setup* setup_;

      shared_ptr< StructuredBlockStorage > blocks_;
      const BlockDataID flagFieldID_;
      const BlockDataID pdfFieldID_;
      const BlockDataID bodyStorageID_;

      bool fileIO_;
      std::string filename_;
      real_t forceScalingFactor_;

      // drag coefficient
      real_t normalizedDragOld_;
      real_t normalizedDragNew_;

};


//////////
// MAIN //
//////////


//*******************************************************************************************************************
/*!\brief Testcase that checks the drag force acting on a fixed sphere in the center of a cubic domain in Stokes flow
 *
 * The drag force for this problem (often denoted as Simple Cubic setup) is given by a semi-analytical series expansion.
 * The cubic domain is periodic in all directions, making it a physically infinite periodic array of spheres.
   \verbatim
           _______________
        ->|               |->
        ->|      ___      |->
      W ->|     /   \     |-> E
      E ->|    |  x  |    |-> A
      S ->|     \___/     |-> S
      T ->|               |-> T
        ->|_______________|->

   \endverbatim
 *
 * The collision model used for the LBM is TRT with a relaxation parameter tau=1.5 and the magic parameter 3/16.
 * The Stokes approximation of the equilibrium PDFs is used.
 * The flow is driven by a constant body force of 1e-5.
 * The domain is 40x40x40, and the sphere diameter is given by chi * domainlength
 * Refinement is applied with the specified amount of levels.
 * The simulation is run until the relative change in the dragforce between 100 time steps is less than 1e-5.
 * The pe is not used since the sphere is kept fixed and the force is explicitly reset after each time step.
 *
 */
//*******************************************************************************************************************

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   ///////////////////
   // Customization //
   ///////////////////

   bool shortrun       = false;
   bool funcTest       = false;
   bool fileIO         = false;
   bool vtkIO          = false;
   bool MO_CLI         = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) { shortrun  = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" ) == 0 ) { funcTest  = true; continue; }
      if( std::strcmp( argv[i], "--fileIO"   ) == 0 ) { fileIO    = true; continue; }
      if( std::strcmp( argv[i], "--vtkIO"    ) == 0 ) { vtkIO     = true; continue; }
      if( std::strcmp( argv[i], "--MO_CLI"   ) == 0 ) { MO_CLI    = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }


   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   Setup setup;


   setup.length         = real_c( 40 );                      // length of the cubic domain in lattice cells
   setup.chi            = real_c( 0.1 );                     // porosity parameter: diameter / length
   setup.levels         = uint_t( 3 );                       // number of refinement levels (1: no refinement)
   setup.tau            = real_c( 1.5 );                     // relaxation time
   setup.extForce       = real_c( 1e-5 );                    // constant body force in lattice units
   setup.checkFrequency = uint_t( 100 );                     // evaluate the dragforce only every checkFrequency time steps
   setup.radius = real_c(0.5) * setup.chi * setup.length;    // sphere radius
   setup.visc   = ( setup.tau - real_c(0.5) ) / real_c(3);   // viscosity in lattice units

   const real_t omega      = real_c(1) / setup.tau;          // relaxation rate
   const real_t dx         = real_c(1);                      // lattice dx
   const real_t convergenceLimit = real_c( 1e-5 );           // tolerance for relative change in drag force
   const uint_t timesteps  =  funcTest ? 1 : ( shortrun ? uint_c(150) : uint_c( 10000 ) );  // maximum number of time steps for the whole simulation

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = uint_t( 4 );
   const uint_t YBlocks = XBlocks;
   const uint_t ZBlocks = XBlocks;

   // create fully periodic domain with refined blocks
   auto blocks = createBlockStructure( setup, XBlocks, YBlocks, ZBlocks, true, true, true );

   if( vtkIO )
   {
      vtk::writeDomainDecomposition( blocks );
   }

   ////////
   // PE //
   ////////

   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe, overlap region is defined by the finest level as the bodies act there
   const real_t overlap = real_c( 1.5 ) * dx / real_c( uint_t(1) << ( setup.levels - uint_t(1) ) );

   // create the sphere in the middle of the domain
   Vector3<real_t> position ( setup.length * real_c(0.5));
   pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, position, setup.radius );

   // synchronize often enough for large bodies
   for( uint_t i = 0; i < std::max( uint_c(1), XBlocks / uint_c(2) ) * ( uint_t(1) << ( setup.levels - uint_t(1) ) ); ++i)
      pe::syncShadowOwners<BodyTypeTuple>( blocks->getBlockForest(), bodyStorageID, nullptr, overlap);

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ),
                                                 lbm::force_model::SimpleConstant( Vector3<real_t> ( setup.extForce, 0, 0 ) ) );

   // add PDF field ( uInit = <0,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_c(1),
                                                                         FieldGhostLayers, field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field", FieldGhostLayers );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf, FieldGhostLayers );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
                                    MyBoundaryHandling( flagFieldID, pdfFieldID, bodyFieldID ), "boundary handling" );

   // initially map pe bodies into the LBM simulation
   if( MO_CLI )
   {
      // uses a higher order boundary condition (CLI)
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_CLI_Flag );
   }else{
      // uses standard bounce back boundary conditions
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID,  MO_BB_Flag );
   }


   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // add LBM sweep with refinement
   auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T >( blocks, sweep, pdfFieldID, boundaryHandlingID );
   timeloop.addFuncBeforeTimeStep( makeSharedFunctor( refinementTimestep ), "LBM refinement time step" );

   // evaluate force
   shared_ptr< ForceEval > forceEval = make_shared< ForceEval >( &timeloop, &setup, blocks, flagFieldID, pdfFieldID, bodyStorageID, fileIO );
   timeloop.addFuncAfterTimeStep( makeSharedFunctor( forceEval ), "drag force evaluation" );

   // resetting force
   timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID ), "reset force on sphere");

   if( vtkIO )
   {
      const uint_t writeFrequency = setup.checkFrequency;

      // spheres
      auto bodyVtkOutput   = make_shared< pe::SphereVtkOutput >( bodyStorageID, blocks->getBlockStorage() );
      auto bodyVTK   = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", writeFrequency );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( bodyVTK ), "VTK (sphere data)" );

      // flag field (written only once in the first time step)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", timesteps, 0 );
      flagFieldVTK->addCellDataWriter( make_shared< field::VTKWriter< FlagField_T > >( flagFieldID, "FlagField" ) );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( flagFieldVTK ), "VTK (flag field data)" );

      // pdf field
      auto pdfFieldVTK = vtk::createVTKOutput_BlockData( blocks, "fluid_field", writeFrequency );

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

   // time loop
   for (uint_t i = 0; i < timesteps; ++i )
   {
      // perform a single simulation step
      timeloop.singleStep( timeloopTiming );

      // check if the relative change in the normalized drag force is below the specified convergence criterion
      if ( i > setup.checkFrequency && forceEval->getDragForceDiff() < convergenceLimit )
      {
         // if simulation has converged, terminate simulation
         break;
      }

   }

   timeloopTiming.logResultOnRoot();

   if ( !funcTest && !shortrun ){
      // check the result
      real_t relErr = std::fabs( ( setup.analyticalDrag - forceEval->getDragForce() ) / setup.analyticalDrag );
      if ( fileIO )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::cout << "Analytical drag: " << setup.analyticalDrag << "\n"
                      << "Simulated drag: " << forceEval->getDragForce() << "\n"
                      << "Relative error: " << relErr << "\n";
         }
      }
      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS( relErr, real_c(0.1) );
   }

   return 0;

}

} //namespace drag_force_sphere_mem_refinement

int main( int argc, char **argv ){
   drag_force_sphere_mem_refinement::main(argc, argv);
}
