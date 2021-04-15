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
//! \file DragForceSpherePSM.cpp
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
#include "lbm/lattice_model/CollisionModel.h"
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

namespace drag_force_sphere_psm
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

// PDF field, flag field & body field
using ForceModel_T = lbm::force_model::LuoConstant;
typedef lbm::D3Q19< lbm::collision_model::SRT, false, ForceModel_T>  LatticeModel_T;

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
   real_t chi;
   real_t extForce;
   real_t analyticalDrag;
};

enum PSMVariant { SC1W1, SC2W1, SC3W1, SC1W2, SC2W2, SC3W2 };

PSMVariant to_PSMVariant( const std::string& s )
{
   if( s == "SC1W1"  ) return PSMVariant::SC1W1;
   if( s == "SC2W1"  ) return PSMVariant::SC2W1;
   if( s == "SC3W1"  ) return PSMVariant::SC3W1;
   if( s == "SC1W2"  ) return PSMVariant::SC1W2;
   if( s == "SC2W2"  ) return PSMVariant::SC2W2;
   if( s == "SC3W2"  ) return PSMVariant::SC3W2;
   throw std::runtime_error("invalid conversion from text to PSMVariant");
}

std::string PSMVariant_to_string ( const PSMVariant& m )
{
   if( m == PSMVariant::SC1W1 ) return "SC1W1";
   if( m == PSMVariant::SC2W1 ) return "SC2W1";
   if( m == PSMVariant::SC3W1 ) return "SC3W1";
   if( m == PSMVariant::SC1W2 ) return "SC1W2";
   if( m == PSMVariant::SC2W2 ) return "SC2W2";
   if( m == PSMVariant::SC3W2 ) return "SC3W2";
   throw std::runtime_error("invalid conversion from PSMVariant to string");
}


class DragForceEvaluator
{
public:
   DragForceEvaluator( SweepTimeloop* timeloop, Setup* setup, const shared_ptr< StructuredBlockStorage > & blocks,
                       const BlockDataID & pdfFieldID, const BlockDataID & bodyStorageID, const BlockDataID & bodyAndVolumeFractionFieldID,
                       PSMVariant method ) :
      timeloop_( timeloop ), setup_( setup ), blocks_( blocks ),
      pdfFieldID_( pdfFieldID ), bodyStorageID_( bodyStorageID ), bodyAndVolumeFractionFieldID_( bodyAndVolumeFractionFieldID ),
      method_( method ), normalizedDragOld_( 0.0 ), normalizedDragNew_( 0.0 )
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
   }

   // evaluate the acting drag force
   void operator()()
   {
      const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

      if( timestep % setup_->checkFrequency != 0) return;

      // get force in x-direction acting on the sphere
      real_t forceX = computeDragForce();
      // get average volumetric flowrate in the domain
      real_t uBar = computeAverageVel();

      // f_total = f_drag + f_buoyancy
      real_t totalForce = forceX  + real_c(4.0/3.0) * math::pi * setup_->radius * setup_->radius * setup_->radius * setup_->extForce ;

      real_t normalizedDragForce = totalForce / real_c( 6.0 * math::pi * setup_->visc * setup_->radius * uBar );

      // update drag force values
      normalizedDragOld_ = normalizedDragNew_;
      normalizedDragNew_ = normalizedDragForce;
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

   void logResultToFile( const std::string & filename ) const
   {
      // write to file if desired
      // format: length tau viscosity simulatedDrag analyticalDrag\n
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( filename.c_str(), std::ofstream::app );
         file.precision(8);
         file << setup_->length << " " << setup_->tau << " " << setup_->visc << " " << normalizedDragNew_ << " " << setup_->analyticalDrag << "\n";
         file.close();
      }
   }

private:

   // obtain the drag force acting on the sphere by summing up all the process local parts of fX
   real_t computeDragForce()
   {
      real_t force = real_c( 0.0 );

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin( *blockIt, bodyStorageID_ ); bodyIt != pe::BodyIterator::end(); ++bodyIt)
         {
            force += bodyIt->getForce()[0];
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( force, mpi::SUM );
      }

      return force;
   }

   // calculate the average velocity in forcing direction (here: x) inside the domain (assuming dx=1)
   real_t computeAverageVel()
   {
      real_t velocity_sum = real_t(0);

      for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
      {
         PdfField_T* pdfField = blockIt->getData< PdfField_T >( pdfFieldID_ );
         BodyAndVolumeFractionField_T * bodyAndVolumeFractionField = blockIt->getData< BodyAndVolumeFractionField_T >( bodyAndVolumeFractionFieldID_ );

         auto xyzFieldSize = pdfField->xyzSize();
         if( method_ == PSMVariant::SC1W1 || method_ == PSMVariant::SC2W1 || method_ == PSMVariant::SC3W1 )
         {
            for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
            {
               velocity_sum += pe_coupling::getPSMMacroscopicVelocity<LatticeModel_T, 1>( *blockIt, pdfField, bodyAndVolumeFractionField, *blocks_, *fieldIt )[0];
            }
         }else{
            for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
            {
               velocity_sum += pe_coupling::getPSMMacroscopicVelocity<LatticeModel_T, 2>( *blockIt, pdfField, bodyAndVolumeFractionField, *blocks_, *fieldIt )[0];
            }
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( velocity_sum, mpi::SUM );
      }

      return velocity_sum / real_c( setup_->length * setup_->length * setup_->length );
   }

   SweepTimeloop* timeloop_;

   Setup* setup_;

   shared_ptr< StructuredBlockStorage > blocks_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyStorageID_;
   const BlockDataID bodyAndVolumeFractionFieldID_;

   PSMVariant method_;

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
 * The domain is 32x32x32, and the sphere has a diameter of 16 cells ( chi * domainlength )
 * The simulation is run until the relative change in the dragforce between 100 time steps is less than 1e-5.
 * The pe is not used since the sphere is kept fixed and the force is explicitly reset after each time step.
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
   bool logging   = false;
   PSMVariant method = PSMVariant::SC3W2;
   real_t tau     = real_c( 1.5 );
   uint_t length  = uint_c( 32 );

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun"  ) == 0 ) { shortrun = true; continue; }
      if( std::strcmp( argv[i], "--funcTest"  ) == 0 ) { funcTest = true; continue; }
      if( std::strcmp( argv[i], "--logging"   ) == 0 ) { logging  = true; continue; }
      if( std::strcmp( argv[i], "--PSMVariant") == 0 ) { method   = to_PSMVariant( argv[++i] ); continue; }
      if( std::strcmp( argv[i], "--tau"       ) == 0 ) { tau      = real_c( std::atof( argv[++i] ) ); continue; }
      if( std::strcmp( argv[i], "--length"    ) == 0 ) { length   = uint_c( std::atof( argv[++i] ) ); continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   Setup setup;

   setup.length         = length;                            // length of the cubic domain in lattice cells
   setup.chi            = real_c( 0.5 );                     // porosity parameter: diameter / length
   setup.tau            = tau;                               // relaxation time
   setup.extForce       = real_c( 1e-5 );                    // constant body force in lattice units
   setup.checkFrequency = uint_t( 100 );                     // evaluate the dragforce only every checkFrequency time steps
   const real_t convergenceLimit = real_c( 1e-5 );           // tolerance for relative change in drag force

   setup.radius            = real_c(0.5) * setup.chi * real_c( setup.length );
   const real_t omega      = real_c(1) / setup.tau;
   setup.visc              = lbm::collision_model::viscosityFromOmega( omega );
   const real_t dx         = real_c(1);
   const uint_t timesteps  = funcTest ? 2 : ( shortrun ? uint_c(150) : uint_c( 1000000 ) );  // maximum number of time steps for the whole simulation

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
   pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, position, setup.radius );

   // synchronize often enough for large bodies
   for( uint_t i = 0; i < XBlocks / 2; ++i)
      pe::syncShadowOwners<BodyTypeTuple>( blocks->getBlockForest(), bodyStorageID, nullptr, overlap);

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omega, ForceModel_T( Vector3<real_t> ( setup.extForce, 0, 0 ) ) );

   // add PDF field ( uInit = <0,0,0>, rhoInit = 1 )
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_c(1),
                                                                         uint_t(1), field::zyxf );

   // add body and volume fraction field
   BlockDataID bodyAndVolumeFractionFieldID = field::addToStorage< BodyAndVolumeFractionField_T >( blocks, "body and volume fraction field",
                                                                                                   std::vector<BodyAndVolumeFraction_T>(), field::zyxf, 0 );

   // map bodies and calculate solid volume fraction initially
   pe_coupling::BodyAndVolumeFractionMapping bodyMapping( blocks, globalBodyStorage, bodyStorageID, bodyAndVolumeFractionFieldID );
   bodyMapping();

   // initialize the PDF field for PSM
   if( method == PSMVariant::SC1W1 || method == PSMVariant::SC2W1 || method == PSMVariant::SC3W1 )
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

   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< field::communication::PackInfo< PdfField_T > >( pdfFieldID ) );
   commFunction = scheme;

   // since external forcing is applied, the evaluation of the velocity has to be carried out directly after the streaming step
   // however, the default sweep is a  stream - collide step, i.e. after the sweep, the velocity evaluation is not correct
   // solution: split the sweep explicitly into collide and stream

   // define partially saturated cells sweep with the templates SC and W for different variants
   // here special case, as we do not have a flag field
   if( method == PSMVariant::SC1W1 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,1,1>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // collision sweep & calculation of drag force
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
   } else if( method == PSMVariant::SC2W1 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,2,1>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // collision sweep & calculation of drag force
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
   } else if( method == PSMVariant::SC3W1 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,3,1>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // collision sweep & calculation of drag force
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
   } else if( method == PSMVariant::SC1W2 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,1,2>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // collision sweep & calculation of drag force
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
   } else if( method == PSMVariant::SC2W2 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,2,2>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // collision sweep & calculation of drag force
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
   } else if( method == PSMVariant::SC3W2 )
   {
      auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,3,2>(pdfFieldID, bodyAndVolumeFractionFieldID, blocks);
      // collision sweep & calculation of drag force
      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "cell-wise LB sweep (collide)" );

      // communication, streaming and force evaluation
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( lbm::makeStreamSweep( sweep ), "cell-wise LB sweep (stream)" );
   } else
   {
      throw std::runtime_error("Invalid PSM method!");
   }

   shared_ptr< DragForceEvaluator > forceEval = make_shared< DragForceEvaluator >( &timeloop, &setup, blocks, pdfFieldID, bodyStorageID, bodyAndVolumeFractionFieldID, method );
   timeloop.addFuncAfterTimeStep( SharedFunctor< DragForceEvaluator >(forceEval), "drag force evaluation" );

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
      if ( logging )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::cout << "Analytical drag: " << setup.analyticalDrag << "\n"
                      << "Simulated drag: " << forceEval->getDragForce() << "\n"
                      << "Relative error: " << relErr << "\n";
         }
         forceEval->logResultToFile( "log_DragForceSpherePSM_"+PSMVariant_to_string(method)+".txt" );
      }
      // the relative error has to be below 10% (25% for SC2)
      WALBERLA_CHECK_LESS( relErr, (method == PSMVariant::SC2W1 || method == PSMVariant::SC2W2) ? real_c(0.25) : real_c(0.1) );
   }

   return 0;
}

} //namespace drag_force_sphere_psm

int main( int argc, char **argv ){
   drag_force_sphere_psm::main(argc, argv);
}
