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
//! \file TorqueSphereMEM.cpp
//! \ingroup pe_coupling
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
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "pe/basic.h"

#include "pe_coupling/mapping/all.h"
#include "pe_coupling/momentum_exchange_method/all.h"
#include "pe_coupling/utility/all.h"

#include <vector>
#include <iomanip>
#include <iostream>

namespace torque_sphere_mem
{


///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
typedef lbm::D3Q19< lbm::collision_model::TRT >  LatticeModel_T;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

const uint_t FieldGhostLayers = 1;

// boundary handling
typedef pe_coupling::SimpleBB       < LatticeModel_T, FlagField_T >  MO_BB_T;
typedef pe_coupling::CurvedLinear   < LatticeModel_T, FlagField_T > MO_CLI_T;
typedef pe_coupling::CurvedQuadratic< LatticeModel_T, FlagField_T >  MO_MR_T;

typedef BoundaryHandling< FlagField_T, Stencil_T, MO_BB_T, MO_CLI_T, MO_MR_T > BoundaryHandling_T;

using BodyTypeTuple = std::tuple<pe::Sphere> ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID MO_BB_Flag   ( "moving obstacle BB" );
const FlagUID MO_CLI_Flag  ( "moving obstacle CLI" );
const FlagUID MO_MR_Flag   ( "moving obstacle MR" );


////////////////
// PARAMETERS //
////////////////

struct Setup{
   uint_t checkFrequency;
   real_t visc;
   real_t tau;
   real_t radius;
   uint_t length;
   real_t phi;
   real_t angularVel;
   real_t analyticalTorque;
};


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, const BlockDataID & pdfFieldPreColID, const BlockDataID & bodyFieldID ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), pdfFieldPreColID_ ( pdfFieldPreColID ), bodyFieldID_ ( bodyFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID pdfFieldPreColID_;
   const BlockDataID bodyFieldID_;

}; // class MyBoundaryHandling

BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );
   PdfField_T *  pdfFieldPreCol  = block->getData< PdfField_T > ( pdfFieldPreColID_ );
   BodyField_T * bodyField       = block->getData< BodyField_T >( bodyFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "fixed obstacle boundary handling", flagField, fluid,
                                                            MO_BB_T (  "MO_BB",  MO_BB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                           MO_CLI_T ( "MO_CLI", MO_CLI_Flag, pdfField, flagField, bodyField, fluid, *storage, *block ),
                                                            MO_MR_T (  "MO_MR",  MO_MR_Flag, pdfField, flagField, bodyField, fluid, *storage, *block, pdfFieldPreCol ) );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}


//*******************************************************************************************************************
/*!\brief Evaluating the torque on a sphere, rotating with a constant angular velocity
 */
//*******************************************************************************************************************
class TorqueEval
{
   public:
      TorqueEval( SweepTimeloop* timeloop, Setup* setup,
                  const shared_ptr< StructuredBlockStorage > & blocks,
                  const BlockDataID & bodyStorageID,
                  bool fileIO ) :
         timeloop_( timeloop ), setup_( setup ),
         blocks_( blocks ),
         bodyStorageID_ ( bodyStorageID ),
         fileIO_(fileIO), torqueOld_(0.0), torqueNew_(0.0)
      {
         // calculate the (semi)analytical torque value
         // see also Hofmann et al. - Hydrodynamic interactions in colloidal crystals:(II). Application to dense cubic and tetragonal arrays (1999), Eqs. 5.1 and 5.5
         const real_t S = real_c( 1.95708 );
         setup_->analyticalTorque = - setup_->visc * ( real_c(6) * setup_->phi /( real_c(1) - setup_->phi - S * std::pow(setup_->phi, real_c(10./3.) ) ) ) *
                                    setup_->angularVel * real_c( setup_->length * setup_->length * setup_->length );

         if (fileIO_)
         {
            std::ofstream file;
            filename_ = "TorqueSphere.txt";
            WALBERLA_ROOT_SECTION()
            {
               file.open( filename_.c_str() );
               file << "#\t torqueSim\t torqueAnaly\n";
               file.close();
            }
         }
      }

      // evaluate the acting torque
      void operator()()
      {
         const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

         if( timestep % setup_->checkFrequency != 0) return;

         // update torque values
         torqueOld_ = torqueNew_;
         torqueNew_ = calculateTorque();

         // write to file if desired
         WALBERLA_ROOT_SECTION()
         {
            if (fileIO_)
            {
               std::ofstream file;
               file.open( filename_.c_str(), std::ofstream::app );
               file.setf( std::ios::unitbuf );
               file.precision(15);
               file << timestep << " " << torqueNew_ << " " << setup_->analyticalTorque << "\n";
               file.close();
            }
         }
      }

      // obtain the torque acting on the sphere by summing up all the process local parts
      real_t calculateTorque()
      {
         real_t torque = real_c( 0 );
         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_ ); bodyIt != pe::BodyIterator::end(); ++bodyIt)
            {
               torque += bodyIt->getTorque()[1];
            }
         }

         WALBERLA_MPI_SECTION()
         {
            mpi::allReduceInplace( torque, mpi::SUM );
         }
         return torque;
      }

      // return the relative temporal change in the torque
      real_t getTorqueDiff() const
      {
         return std::fabs( ( torqueNew_ - torqueOld_ ) / torqueNew_ );
      }

      // return the torque
      real_t getTorque() const
      {
         return torqueNew_;
      }

   private:
      SweepTimeloop* timeloop_;

      Setup* setup_;


      shared_ptr< StructuredBlockStorage > blocks_;
      const BlockDataID bodyStorageID_;

      bool fileIO_;
      std::string filename_;

      real_t torqueOld_;
      real_t torqueNew_;

};

class PDFFieldCopy
{
public:
   PDFFieldCopy( const BlockDataID & pdfFieldID, const BlockDataID & pdfFieldPreColID )
      : pdfFieldID_( pdfFieldID ), pdfFieldPreColID_( pdfFieldPreColID )
   {}
   void operator()( IBlock * block )
   {
      auto pdfField = block->getData< PdfField_T>( pdfFieldID_ );
      auto pdfFieldPreCol = block->getData< PdfField_T>( pdfFieldPreColID_ );
      std::copy( pdfField->data(), pdfField->data() + pdfField->allocSize(), pdfFieldPreCol->data() );
      pdfFieldPreCol->resetLatticeModel( pdfField->latticeModel() );

   }
private:
   BlockDataID pdfFieldID_;
   BlockDataID pdfFieldPreColID_;
};


//////////
// MAIN //
//////////


//*******************************************************************************************************************
/*!\brief Testcase that checks the torque acting on a constantly rotating sphere in the center of a cubic domain
 *
 * The torque for this problem (often denoted as Simple Cubic setup) is given by a semi-analytical formula.
 * The cubic domain is periodic in all directions, making it a physically infinite periodic array of spheres.
   \verbatim
         _______________
        |       <-      |
        |      ___      |
        |     /   \     |
        |    |  x  |    |
        |     \___/     |
        |      ->       |
        |_______________|

   \endverbatim
 *
 * The collision model used for the LBM is TRT with a relaxation parameter tau=1.5 and the magic parameter 3/16.
 * The Stokes approximation of the equilibrium PDFs is used.
 * The sphere rotates with a angular velocity of 1e-5.
 * The domain is 32x32x32, and the sphere has a diameter of around 27 cells ( chi * domainlength )
 * The simulation is run until the relative change in the torque between 100 time steps is less than 1e-5.
 * The pe is not used since the angular velocity is kept constant and the force is explicitly reset after each time step.
 *
 */
//*******************************************************************************************************************


int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto processes = MPIManager::instance()->numProcesses();

   if( processes != 1 && processes != 2 && processes != 4 && processes != 8)
   {
      std::cerr << "Number of processes must be equal to either 1, 2, 4, or 8!" << std::endl;
      return EXIT_FAILURE;
   }

   ///////////////////
   // Customization //
   ///////////////////

   bool shortrun       = false;
   bool funcTest       = false;
   bool fileIO         = false;
   bool MO_CLI         = false;
   bool MO_MR          = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) { shortrun  = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" ) == 0 ) { funcTest  = true; continue; }
      if( std::strcmp( argv[i], "--fileIO" )   == 0 ) { fileIO    = true; continue; }
      if( std::strcmp( argv[i], "--MO_CLI" )   == 0 ) { MO_CLI    = true; continue; }
      if( std::strcmp( argv[i], "--MO_MR"  )   == 0 ) { MO_MR     = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }


   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   Setup setup;

   setup.length         = uint_t( 32 );                       // length of the cubic domain in lattice cells
   const real_t chi     = real_c( 0.85 );                     // porosity parameter: diameter / length
   setup.tau            = real_c( 1.5 );                      // relaxation time
   setup.angularVel     = real_c( 1e-5 );                     // angular velocity of the sphere
   setup.checkFrequency = uint_t( 100 );                      // evaluate the torque only every checkFrequency time steps
   setup.radius = real_c(0.5) * chi * real_c( setup.length ); // sphere radius
   setup.visc   = ( setup.tau - real_c(0.5) ) / real_c(3);    // viscosity in lattice units
   setup.phi    = real_c(4.0/3.0) * math::pi * setup.radius * setup.radius * setup.radius /
                  ( real_c( setup.length * setup.length * setup.length ) ); // solid volume fraction
   const real_t omega      = real_c(1) / setup.tau;           // relaxation rate
   const real_t dx         = real_c(1);                       // lattice dx
   const real_t convergenceLimit = real_c( 1e-5 );            // tolerance for relative change in torque
   const uint_t timesteps  =  funcTest ? 1 : ( shortrun ? uint_c(150) : uint_c( 5000 ) );  // maximum number of time steps for the whole simulation

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = (processes >= 2) ? uint_t( 2 ) : uint_t( 1 );
   const uint_t YBlocks = (processes >= 4) ? uint_t( 2 ) : uint_t( 1 );
   const uint_t ZBlocks = (processes == 8) ? uint_t( 2 ) : uint_t( 1 );
   const uint_t XCells = setup.length / XBlocks;
   const uint_t YCells = setup.length / YBlocks;
   const uint_t ZCells = setup.length / ZBlocks;

   // create fully periodic domain
   auto blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells, dx, true,
                                                      true, true, true );

   ////////
   // PE //
   ////////

   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_c( 1.5 ) * dx;

   if( setup.radius > real_c( setup.length ) * real_c(0.5) - overlap )
   {
      std::cerr << "Periodic sphere is too large and would lead to invalid PE state!" << std::endl;
      return EXIT_FAILURE;
   }

   // create the sphere in the middle of the domain
   Vector3<real_t> position (real_c(setup.length) * real_c(0.5));
   auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, position, setup.radius );
   if ( sphere != nullptr )
   {
      sphere->setAngularVel( real_c(0), setup.angularVel, real_c(0) );
   }

   // synchronize often enough for large bodies
   for( uint_t i = 0; i < XBlocks / 2; ++i)
      pe::syncShadowOwners<BodyTypeTuple>( blocks->getBlockForest(), bodyStorageID, nullptr, overlap);


   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   // add PDF field ( uInit = <0,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_c(1),
                                                                         uint_t(1), field::zyxf );

   // add PDF field, for MR boundary condition only
   BlockDataID pdfFieldPreColID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, " pre collision pdf field", latticeModel,
                                                                               Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_c(1),
                                                                               uint_t(1), field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add body field
   BlockDataID bodyFieldID = field::addToStorage<BodyField_T>( blocks, "body field", nullptr, field::zyxf );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
                                    MyBoundaryHandling( flagFieldID, pdfFieldID, pdfFieldPreColID, bodyFieldID ), "boundary handling" );

   ///////////////
   // TIME LOOP //
   ///////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;

   blockforest::communication::UniformBufferedScheme< stencil::D3Q27 > scheme( blocks );
   scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
   commFunction = scheme;

   // initially map pe bodies into the LBM simulation
   if( MO_CLI )
   {
      // uses a higher order boundary condition (CLI)
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID, MO_CLI_Flag );
   }else if ( MO_MR ){
      // uses a higher order boundary condition (MR)
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID,  MO_MR_Flag );
   }else{
      // uses standard bounce back boundary conditions
      pe_coupling::mapMovingBodies< BoundaryHandling_T >( *blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage, bodyFieldID,  MO_BB_Flag );
   }

   shared_ptr< TorqueEval > torqueEval = make_shared< TorqueEval >( &timeloop, &setup, blocks, bodyStorageID, fileIO );

   if( MO_MR )
   {
      // copy pdfs to have the pre collision PDFs available, needed for MR boundary treatment
      timeloop.add() << Sweep( PDFFieldCopy( pdfFieldID, pdfFieldPreColID ), "pdf field copy" );

      auto sweep = lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );

      // collision sweep
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      // stream and torque evaluation
      timeloop.add() << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   }
   else
   {
      // add LBM communication function and boundary handling sweep (does the hydro force calculations and the no-slip treatment)
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      // LBM stream + collide and torque evaluation
      timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LatticeModel_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag ) ), "cell-wise LB sweep" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   }
   // resetting force
   timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID ), "reset force on sphere");

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

      // check if the relative change in the torque is below the specified convergence criterion
      if ( i > setup.checkFrequency && torqueEval->getTorqueDiff() < convergenceLimit )
      {
         // if simulation has converged, terminate simulation
         break;
      }

   }

   timeloopTiming.logResultOnRoot();

   // check the result
   if ( !funcTest && !shortrun ){
      real_t relErr = std::fabs( ( setup.analyticalTorque - torqueEval->getTorque() ) / setup.analyticalTorque );
      if ( fileIO )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::cout << "Analytical torque: " << setup.analyticalTorque << "\n"
                      << "Simulated torque: " << torqueEval->getTorque() << "\n"
                      << "Relative error: " << relErr << "\n";
         }
      }
      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS( relErr, real_c(0.1) );
   }

   return 0;

}

} //namespace torque_sphere_mem

int main( int argc, char **argv ){
   torque_sphere_mem::main(argc, argv);
}
