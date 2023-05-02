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
//! \file TorqueSpherePSM.cpp
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

#include "pe_coupling/partially_saturated_cells_method/all.h"
#include "pe_coupling/utility/all.h"

#include <vector>
#include <iomanip>
#include <iostream>

namespace torque_sphere_psm
{


///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
typedef lbm::D3Q19< lbm::collision_model::SRT, false >  LatticeModel_T;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

typedef std::pair< pe::BodyID, real_t >                              BodyAndVolumeFraction_T;
typedef GhostLayerField< std::vector< BodyAndVolumeFraction_T >, 1 > BodyAndVolumeFractionField_T;

using BodyTypeTuple = std::tuple<pe::Sphere> ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );

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
            filename_ = "TorqueSpherePSM.txt";
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

   bool shortrun  = false;
   bool funcTest  = false;
   bool fileIO    = false;
   bool SC1W1     = false;
   bool SC2W1     = false;
   bool SC3W1     = false;
   bool SC1W2     = false;
   bool SC2W2     = false;
   bool SC3W2     = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) { shortrun  = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" ) == 0 ) { funcTest  = true; continue; }
      if( std::strcmp( argv[i], "--fileIO"   ) == 0 ) { fileIO    = true; continue; }
      if( std::strcmp( argv[i], "--SC1W1"    ) == 0 ) { SC1W1     = true; continue; }
      if( std::strcmp( argv[i], "--SC2W1"    ) == 0 ) { SC2W1     = true; continue; }
      if( std::strcmp( argv[i], "--SC3W1"    ) == 0 ) { SC3W1     = true; continue; }
      if( std::strcmp( argv[i], "--SC1W2"    ) == 0 ) { SC1W2     = true; continue; }
      if( std::strcmp( argv[i], "--SC2W2"    ) == 0 ) { SC2W2     = true; continue; }
      if( std::strcmp( argv[i], "--SC3W2"    ) == 0 ) { SC3W2     = true; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   if( !SC1W1 && !SC2W1 && !SC3W1 && !SC1W2 && !SC2W2 && !SC3W2 )
   {
      std::cerr << "Specify the model (--SC_W_) you want to use for the partially saturated cells method!" << std::endl;
      return EXIT_FAILURE;
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
   LatticeModel_T latticeModel = LatticeModel_T( omega );

   // add PDF field ( uInit = <0,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_c(1),
                                                                         uint_t(1), field::fzyx );

   // add body and volume fraction field
   BlockDataID bodyAndVolumeFractionFieldID = field::addToStorage< BodyAndVolumeFractionField_T >( blocks, "body and volume fraction field",
                                                                                                   std::vector<BodyAndVolumeFraction_T>(), field::fzyx, 0 );
   // map bodies and calculate solid volume fraction initially
   pe_coupling::BodyAndVolumeFractionMapping bodyMapping( blocks, globalBodyStorage, bodyStorageID, bodyAndVolumeFractionFieldID );
   bodyMapping();

   // initialize the PDF field for PSM
   if( SC1W1 || SC2W1 || SC3W1 )
   {
      pe_coupling::initializeDomainForPSM< LatticeModel_T, 1 >( *blocks, pdfFieldID, bodyAndVolumeFractionFieldID );
   }
   else
   {
      pe_coupling::initializeDomainForPSM< LatticeModel_T, 2 >( *blocks, pdfFieldID, bodyAndVolumeFractionFieldID );
   }

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

   shared_ptr< TorqueEval > torqueEval = make_shared< TorqueEval >( &timeloop, &setup, blocks, bodyStorageID, fileIO );

   if( SC1W1 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,1,1>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   } else if( SC2W1 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,2,1>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   } else if( SC3W1 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,3,1>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   } else if( SC1W2 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,1,2>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   } else if( SC2W2 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,2,2>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" )
                     << AfterFunction( SharedFunctor< TorqueEval >( torqueEval ), "torque evaluation" );
   } else if( SC3W2 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,3,2>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" )
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
      // the relative error has to be below 10% (25% for SC2)
      WALBERLA_CHECK_LESS( relErr, (SC2W1 || SC2W2) ? real_c(0.25) : real_c(0.1) );
   }

   return 0;

}

} //namespace torque_sphere_psm

int main( int argc, char **argv ){
   torque_sphere_psm::main(argc, argv);
}
