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
//! \file SegreSilberbergPSM.cpp
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

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/StabilityChecker.h"
#include "field/communication/PackInfo.h"

#include "lbm/boundary/all.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/ForceModel.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "pe/basic.h"
#include "pe/vtk/BodyVtkOutput.h"
#include "pe/vtk/SphereVtkOutput.h"
#include "pe/cr/ICR.h"
#include "pe/Types.h"

#include "pe_coupling/partially_saturated_cells_method/all.h"
#include "pe_coupling/utility/all.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/all.h"
#include "field/vtk/all.h"
#include "lbm/vtk/all.h"

#include <functional>
#include <memory>

namespace segre_silberberg_psm
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
typedef lbm::D3Q19< lbm::collision_model::SRT, false, lbm::force_model::SimpleConstant >  LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;
typedef GhostLayerField< pe::BodyID, 1 >  BodyField_T;

const uint_t FieldGhostLayers = 1;

typedef std::pair< pe::BodyID, real_t >                              BodyAndVolumeFraction_T;
typedef GhostLayerField< std::vector< BodyAndVolumeFraction_T >, 1 > BodyAndVolumeFractionField_T;

// boundary handling
typedef lbm::NoSlip< LatticeModel_T, flag_t > NoSlip_T;
typedef BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T > BoundaryHandling_T;

typedef std::tuple< pe::Sphere, pe::Plane > BodyTypeTuple;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID NoSlip_Flag  ( "no slip" );

////////////////
// Parameters //
////////////////

struct Setup
{
   // domain size in x,y,z direction
   uint_t xlength;
   uint_t ylength;
   uint_t zlength;

   // number of block (= processes) in x,y,z, direction
   Vector3<uint_t> nBlocks;

   // cells inside each block in x,y,z direction (uniform distribution)
   Vector3<uint_t> cellsPerBlock;

   real_t particleDiam; // diameter of particle
   real_t nu;           // kin. viscosity
   real_t rho_p;        // particle density
   real_t tau;          // relaxation time
   real_t forcing;      // forcing to drive flow
   real_t Re_p;         // particle Reynoldsnumber
   real_t Re;           // bulk Reynoldsnumber

   uint_t checkFrequency;
   uint_t timesteps;    // maximum number of timesteps to simulate

};

//*******************************************************************************************************************

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

//*******************************************************************************************************************
/*!\brief Responsible for creating the boundary handling and setting up the geometry at the outer domain boundaries
 */
//*******************************************************************************************************************
class MyBoundaryHandling
{
public:

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;

}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField       = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField        = block->getData< PdfField_T > ( pdfFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ) );

   const auto noslip = flagField->getFlag( NoSlip_Flag );

   CellInterval domainBB = storage->getDomainCellBB();

   domainBB.xMin() -= cell_idx_c( FieldGhostLayers ); // -1
   domainBB.xMax() += cell_idx_c( FieldGhostLayers ); // cellsX

   domainBB.yMin() -= cell_idx_c( FieldGhostLayers ); // 0
   domainBB.yMax() += cell_idx_c( FieldGhostLayers ); // cellsY+1

   domainBB.zMin() -= cell_idx_c( FieldGhostLayers ); // 0
   domainBB.zMax() += cell_idx_c( FieldGhostLayers ); // cellsZ+1

   // BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
   handling->forceBoundary( noslip, bottom );

   // TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( top, *block );
   handling->forceBoundary( noslip, top );

   handling->fillWithDomain( FieldGhostLayers );

   return handling;
}
//*******************************************************************************************************************


// initialises the PDF field with a parabolic profile
void initPDF(const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & pdfFieldID, Setup & setup)
{
   // iterate all blocks with an iterator 'block'
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      // get the field data out of the block
      auto pdf = block->getData< PdfField_T > ( pdfFieldID );

      // obtain a CellInterval object that holds information about the number of cells in x,y,z direction of the field excluding ghost layers
      CellInterval xyz = pdf->xyzSize();

      // iterate all (inner) cells in the field
      for( auto cell = xyz.begin(); cell != xyz.end(); ++cell ){

         // obtain the physical coordinate of the center of the current cell
         const Vector3< real_t > p = blocks->getBlockLocalCellCenter( *block, *cell );

         const real_t rho = real_c(1);
         const real_t height = real_c(setup.zlength) * real_c(0.5);

         // parabolic Poiseuille profile for a given external force
         Vector3< real_t > velocity ( setup.forcing / ( real_c(2) * setup.nu) * (height * height - ( p[2] - height ) * ( p[2] - height ) ), real_c(0), real_c(0) );
         pdf->setToEquilibrium( *cell, velocity, rho );
      }
   }
}

//*******************************************************************************************************************
/*!\brief Evaluating the position and velocity of the body
 *
 */
//*******************************************************************************************************************
class SteadyStateCheck
{
   public:
      SteadyStateCheck( SweepTimeloop* timeloop, Setup * setup,
                        const shared_ptr< StructuredBlockStorage > & blocks, const BlockDataID & bodyStorageID,
                        bool fileIO, bool SC1W1, bool SC2W1, bool SC3W1, bool SC1W2, bool SC2W2, bool SC3W2 )
      :  timeloop_( timeloop ), setup_(setup),
         blocks_( blocks ), bodyStorageID_( bodyStorageID ),
         fileIO_(fileIO), SC1W1_(SC1W1), SC2W1_(SC2W1), SC3W1_(SC3W1), SC1W2_(SC1W2), SC2W2_(SC2W2), SC3W2_(SC3W2),
         posOld_(0.0), posNew_(0.0)
      {
         if ( fileIO_ )
         {
            std::ofstream file;
            filename_ = "SegreSilberbergTest_";
            if( SC1W1_ )      filename_ += "SC1W1";
            else if (SC2W1_ ) filename_ += "SC2W1";
            else if (SC3W1_ ) filename_ += "SC3W1";
            else if (SC1W2_ ) filename_ += "SC1W2";
            else if (SC2W2_ ) filename_ += "SC2W2";
            else if (SC3W2_ ) filename_ += "SC3W2";
            filename_ += ".txt";
            WALBERLA_ROOT_SECTION()
            {
               file.open( filename_.c_str() );
               file << "#\t posX\t posY\t posZ\t velX\t velY\t velZ\n";
               file.close();
            }
         }
      }

      void operator()()
      {
         const uint_t timestep (timeloop_->getCurrentTimeStep()+1);

         if( timestep % setup_->checkFrequency != 0) return;

         pe::Vec3 pos      = pe::Vec3(0.0);
         pe::Vec3 transVel = pe::Vec3(0.0);

         for( auto blockIt = blocks_->begin(); blockIt != blocks_->end(); ++blockIt )
         {
            for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
            {
               pos = bodyIt->getPosition();
               transVel = bodyIt->getLinearVel();
            }
         }

         WALBERLA_MPI_SECTION()
         {
            mpi::allReduceInplace( pos, mpi::SUM );
            mpi::allReduceInplace( transVel, mpi::SUM );
         }

         // update position values
         posOld_ = posNew_;
         posNew_ = pos[2];

         if( fileIO_ )
            writeToFile( timestep, pos, transVel);
      }

      // return the relative change in the sphere position ( z coordinate )
      real_t getPositionDiff() const
      {
         return std::fabs( ( posNew_ - posOld_ ) / posNew_ );
      }

      // return the position (z coordinate)
      real_t getPosition() const
      {
         return posNew_;
      }

   private:
      void writeToFile( uint_t timestep, const pe::Vec3 & position, const pe::Vec3 & velocity )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            file.open( filename_.c_str(), std::ofstream::app );
            file.setf( std::ios::unitbuf );
            file.precision(15);

            file << timestep << "\t" << position[0] << "\t" << position[1] << "\t" << position[2]
                             << "\t" << velocity[0] << "\t" << velocity[1] << "\t" << velocity[2]
                             << "\n";
            file.close();
         }
      }

      SweepTimeloop* timeloop_;

      Setup * setup_;

      shared_ptr< StructuredBlockStorage > blocks_;
      const BlockDataID bodyStorageID_;

      bool fileIO_;
      bool SC1W1_, SC2W1_, SC3W1_, SC1W2_, SC2W2_, SC3W2_;

      std::string filename_;

      real_t posOld_;
      real_t posNew_;
};


//////////
// MAIN //
//////////

//*******************************************************************************************************************
/*!\brief Testcase that simulates a neutrally buoyant sphere in a Poiseuille flow
 *
 *   The Segre Silberberg effect denotes that a particle in a Poiseuille flow will eventually arrive
 *   at a distinct equilibrium position, depending on the particle size and its velocity.
   \verbatim
        __________________
   ->
   -->        ___
   --->      /   \
   --->     |  x  |-\
   --->      \___/   \-\
   -->                  \------ equilibrium position
   ->  ___________________

   \endverbatim
 *
 * The collision model used for the LBM is TRT with a relaxation parameter tau=1.4 and the magic parameter 3/16.
 * The particle is a sphere with diameter d=12 cells and density 1 (LU).
 * The domain is a channel with 4d x 4d x 4d, with walls in z-direction and periodicity else.
 * The simulation is run until the relative change in the z-coordinate between 100 time steps is less than 1e-5.
 *
 */
//*******************************************************************************************************************

int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   auto processes = MPIManager::instance()->numProcesses();

   if( processes != 1 && processes != 9 && processes != 18 )
   {
      std::cerr << "Number of processes must be equal to either 1, 9, or 18!" << std::endl;
      return EXIT_FAILURE;
   }

   ///////////////////
   // Customization //
   ///////////////////

   bool shortrun  = false;
   bool functest  = false;
   bool fileIO    = false;
   bool vtkIO     = false;
   bool SC1W1     = false;
   bool SC2W1     = false;
   bool SC3W1     = false;
   bool SC1W2     = false;
   bool SC2W2     = false;
   bool SC3W2     = false;

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--shortrun" ) == 0 ) { shortrun  = true; continue; }
      if( std::strcmp( argv[i], "--funcTest" ) == 0 ) { functest  = true; continue; }
      if( std::strcmp( argv[i], "--fileIO"   ) == 0 ) { fileIO    = true; continue; }
      if( std::strcmp( argv[i], "--vtkIO"    ) == 0 ) { vtkIO     = true; continue; }
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

   setup.particleDiam =  real_c( 12 );                          // particle diameter
   setup.rho_p = real_c( 1 );                                   // particle density
   setup.tau = real_c( 1.4 );                                   // relaxation time
   setup.nu = (real_c(2) * setup.tau - real_c(1) ) / real_c(6); // viscosity in lattice units
   setup.timesteps = functest ? 1 : ( shortrun ? uint_c(100) : uint_c( 250000 ) ); // maximum number of time steps
   setup.zlength = uint_c( 4.0 * setup.particleDiam );          // zlength in lattice cells of the channel
   setup.xlength = setup.zlength;                               // xlength in lattice cells of the channel
   setup.ylength = setup.zlength;                               // ylength in lattice cells of the channel
   setup.Re = real_c( 13 );                                     // bulk Reynoldsnumber = uMean * width / nu

   const real_t particleRadius = real_c(0.5) * setup.particleDiam; // particle radius
   const real_t dx = real_c(1);                                    // lattice grid spacing
   const uint_t numLbmSubCycles = uint_c(2);                       // number of LBM subcycles, i.e. number of LBM steps until one pe step is carried out
   const real_t dt_pe       = real_c(numLbmSubCycles);             // time step size for PE (LU): here 2 times as large as for LBM
   const uint_t pe_interval = uint_c(1);                           // subcycling interval for PE: only 1 PE step
   const real_t omega = real_c(1) / setup.tau;                     // relaxation rate
   const real_t convergenceLimit = real_c( 1e-5 );                 // tolerance for relative change in position
   setup.checkFrequency = uint_c( real_c(100) / dt_pe );           // evaluate the position only every checkFrequency LBM time steps

   const real_t uMean = setup.Re * setup.nu / real_c( setup.zlength );
   const real_t uMax = real_c(2) * uMean;
   setup.forcing = uMax * ( real_c(2) * setup.nu ) / real_c( real_c(3./4.) * real_c(setup.zlength) * real_c(setup.zlength) ); // forcing to drive the flow
   setup.Re_p = uMean * particleRadius * particleRadius / ( setup.nu * real_c( setup.zlength ) );  // particle Reynolds number = uMean * r * r / ( nu * width )

   if( fileIO || vtkIO ){
      WALBERLA_ROOT_SECTION()
      {
         std::cout << "size: <" << setup.xlength << ", " << setup.ylength << ", " << setup.zlength << ">\n"
                   << "Re_p: " << setup.Re_p << ", Re_mean = " << setup.Re << "\n"
                   << "u_max: " << uMax << ", visc = " << setup.nu << "\n";
      }
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   setup.nBlocks[0] = uint_c(3);
   setup.nBlocks[1] = uint_c(3);
   setup.nBlocks[2] = (processes == 18) ? uint_c(2) : uint_c(1);
   setup.cellsPerBlock[0] = setup.xlength / setup.nBlocks[0];
   setup.cellsPerBlock[1] = setup.ylength / setup.nBlocks[1];
   setup.cellsPerBlock[2] = setup.zlength / setup.nBlocks[2];

   auto blocks = blockforest::createUniformBlockGrid( setup.nBlocks[0], setup.nBlocks[1], setup.nBlocks[2],
                                                      setup.cellsPerBlock[0], setup.cellsPerBlock[1], setup.cellsPerBlock[2], dx,
                                                      (processes != 1),
                                                      true, true, false );

   /////////////////
   // PE COUPLING //
   /////////////////

   // set up pe functionality
   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();

   auto bodyStorageID  = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID          = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID          = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   // set up collision response, here DEM solver
   pe::cr::DEM cr( globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID, nullptr );

   // set up synchronization procedure
   const real_t overlap = real_c( 1.5 ) * dx;
   std::function<void(void)> syncCall = std::bind( pe::syncShadowOwners<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // create pe bodies

   // bounding planes (global)
   const auto planeMaterial = pe::createMaterial( "myPlaneMat", real_c(8920), real_c(0), real_c(1), real_c(1), real_c(0), real_c(1), real_c(1), real_c(0), real_c(0) );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,1), Vector3<real_t>(0,0,0), planeMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>(0,0,-1), Vector3<real_t>(0,0,real_c(setup.zlength)), planeMaterial );

   // add the sphere
   const auto sphereMaterial = pe::createMaterial( "mySphereMat", setup.rho_p , real_c(0.5), real_c(0.1), real_c(0.1), real_c(0.24), real_c(200), real_c(200), real_c(0), real_c(0) );
   Vector3<real_t> position( real_c(setup.xlength) * real_c(0.5), real_c(setup.ylength) * real_c(0.5), real_c(setup.zlength) * real_c(0.5) - real_c(1) );
   auto sphere = pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, position, particleRadius, sphereMaterial );
   if( sphere != nullptr )
   {
      real_t height = real_c( 0.5 ) * real_c( setup.zlength );
      // set sphere velocity to undisturbed fluid velocity at sphere center
      sphere->setLinearVel( setup.forcing/( real_c(2) * setup.nu) * ( height * height - ( position[2] - height ) * ( position[2] - height ) ), real_c(0), real_c(0) );
   }

   syncCall();

   ///////////////////////
   // ADD DATA TO BLOCKS //
   ////////////////////////

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omega, lbm::force_model::SimpleConstant( Vector3<real_t> ( setup.forcing, real_c(0), real_c(0) ) ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (fzyx)", latticeModel,
                                                                         Vector3< real_t >( real_c(0), real_c(0), real_c(0) ), real_c(1),
                                                                         uint_t(1), field::fzyx );

   // initialize already with the Poiseuille flow profile
   initPDF( blocks, pdfFieldID, setup);

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
                                    MyBoundaryHandling( flagFieldID, pdfFieldID ), "boundary handling" );

   // add body and volume fraction field
   BlockDataID bodyAndVolumeFractionFieldID = field::addToStorage< BodyAndVolumeFractionField_T >( blocks, "body and volume fraction field",
                                                                                                   std::vector<BodyAndVolumeFraction_T>(), field::fzyx, 0 );
   // map bodies and calculate solid volume fraction initially
   pe_coupling::BodyAndVolumeFractionMapping bodyMapping( blocks, globalBodyStorage, bodyStorageID, bodyAndVolumeFractionFieldID, pe_coupling::selectRegularBodies );
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

   // setup of the LBM communication for synchronizing the pdf field between neighboring blocks
   std::function< void () > commFunction;
   blockforest::communication::UniformBufferedScheme< Stencil_T > scheme( blocks );
   scheme.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldID ) );
   commFunction = scheme;

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), setup.timesteps );

   // compute the volume fractions of the cells
   timeloop.addFuncBeforeTimeStep( bodyMapping, "volume fraction mapping" );

   for( uint_t lbmSubCycle = uint_c(0); lbmSubCycle < numLbmSubCycles; ++lbmSubCycle )
   {
      // add LBM communication function and boundary handling sweep
      timeloop.add() << BeforeFunction( commFunction, "LBM Communication" )
                     << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      // stream + collide LBM step
      if( SC1W1 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,1,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( SC2W1 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,2,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( SC3W1 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,3,1>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( SC1W2 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,1,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( SC2W2 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,2,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      } else if( SC3W2 )
      {
         auto sweep = pe_coupling::makePSMSweep<LatticeModel_T,FlagField_T,3,2>( pdfFieldID, bodyAndVolumeFractionFieldID, blocks, flagFieldID, Fluid_Flag );
         timeloop.add() << Sweep( makeSharedSweep( sweep ), "cell-wise LB sweep" );
      }
   }

   // average the forces acting on the bodies over the number of LBM steps
   timeloop.addFuncAfterTimeStep( pe_coupling::ForceTorqueOnBodiesScaler( blocks, bodyStorageID, real_t(1)/real_c(numLbmSubCycles) ), "Force averaging for several LBM steps" );

   // add pressure force contribution
   // The flow in the channel is here driven by a body force and is periodic
   // Thus, the pressure gradient, usually encountered in a channel flow, is not present here
   // The buoyancy force on the body due to this pressure gradient has to be added 'artificially'
   // F_{buoyancy} = - V_{body} * grad ( p ) = V_{body} * \rho_{fluid} *  a
   // ( V_{body} := volume of body,  a := acceleration driving the flow )
   Vector3<real_t> buoyancyForce(math::pi / real_t(6) * setup.forcing * setup.particleDiam * setup.particleDiam * setup.particleDiam ,
                                 real_t(0), real_t(0));
   timeloop.addFuncAfterTimeStep( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, buoyancyForce ), "Buoyancy force" );

   // add pe timesteps
   timeloop.addFuncAfterTimeStep( pe_coupling::TimeStep( blocks, bodyStorageID, cr, syncCall, dt_pe, pe_interval ), "pe Time Step" );

   // check for convergence of the particle position
   shared_ptr< SteadyStateCheck > check = std::make_shared< SteadyStateCheck >( &timeloop, &setup, blocks, bodyStorageID,
                                                                                                fileIO, SC1W1, SC2W1, SC3W1, SC1W2, SC2W2, SC3W2 );
   timeloop.addFuncAfterTimeStep( SharedFunctor< SteadyStateCheck >( check ), "steady state check" );

   if( vtkIO )
   {
      const uint_t writeFrequency = setup.checkFrequency;

      // spheres
      auto bodyVtkOutput   = make_shared<pe::DefaultBodyVTKOutput>( bodyStorageID, blocks->getBlockStorage() );
      auto bodyVTK   = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", writeFrequency );
      timeloop.addFuncAfterTimeStep( vtk::writeFiles( bodyVTK ), "VTK (sphere data)" );

      // flag field (written only once in the first time step, ghost layers are also written)
      auto flagFieldVTK = vtk::createVTKOutput_BlockData( blocks, "flag_field", setup.timesteps, FieldGhostLayers );
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

   // time loop
   for (uint_t i = 0; i < setup.timesteps; ++i )
   {
      // perform a single simulation step
      timeloop.singleStep( timeloopTiming );

      // check if the relative change in the position is below the specified convergence criterion
      if ( i > setup.checkFrequency && check->getPositionDiff() < convergenceLimit )
      {
         // if simulation has converged, terminate simulation
         break;
      }
   }

   timeloopTiming.logResultOnRoot();

   // check the result
   if ( !functest && !shortrun )
   {
      // reference value from Inamuro et al. - "Flow between parallel walls containing the lines of neutrally buoyant circular cylinders" (2000), Table 2
      // note that this reference value provides only a rough estimate
      real_t referencePos = real_c( 0.2745 ) * real_c( setup.zlength );
      real_t relErr = std::fabs( referencePos - check->getPosition() ) / referencePos;
      if ( fileIO )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::cout << "Expected position: " << referencePos << "\n"
                      << "Simulated position: " << check->getPosition() << "\n"
                      << "Relative error: " << relErr << "\n";
         }
      }
      // the relative error has to be below 10%
      WALBERLA_CHECK_LESS( relErr, real_c(0.1) );
   }

   return EXIT_SUCCESS;
}

} // namespace segre_silberberg_psm

int main( int argc, char **argv ){
   segre_silberberg_psm::main(argc, argv);
}
