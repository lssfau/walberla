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
//! \file SphereWallCollisionBehaviorDPM.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/all.h"
#include "boundary/all.h"
#include "lbm/all.h"
#include "pe_coupling/all.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/all.h"

#include "core/grid_generator/SCIterator.h"

#include <domain_decomposition/SharedSweep.h>

#include "field/AddToStorage.h"
#include "field/communication/PackInfo.h"
#include "field/distributors/all.h"
#include "field/interpolators/all.h"

#include "pe/basic.h"
#include "pe/vtk/SphereVtkOutput.h"
#include "pe/Types.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <vector>
#include <functional>
#include <iomanip>
#include <iostream>
#include <random>

namespace sphere_wall_collision_behavior_dpm
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

///////////////
// CONSTANTS //
///////////////

const uint_t FieldGhostLayers( 1 );

//////////////
// TYPEDEFS //
//////////////

// PDF field, flag field & body field
typedef GhostLayerField< Matrix3<real_t>, 1 >                          TensorField_T;
typedef GhostLayerField< Vector3<real_t>, 1 >                          Vec3Field_T;
typedef GhostLayerField< real_t, 1 >                                   ScalarField_T;
using ForceModel_T = lbm::force_model::GuoField<Vec3Field_T>;

typedef lbm::D3Q19< lbm::collision_model::SRTField<ScalarField_T>, false, ForceModel_T >   LatticeModel_T;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

// boundary handling
typedef lbm::NoSlip< LatticeModel_T, flag_t >                          NoSlip_T;

typedef BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T> BoundaryHandling_T;

typedef std::tuple<pe::Plane, pe::Sphere> BodyTypeTuple ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag  ( "fluid" );
const FlagUID NoSlip_Flag ( "no slip flag" );


// coupling methods
enum DPMethod { GNS };

// interpolation (for PDFs, velocity and/or solid volume fraction)
enum Interpolation { INearestNeighbor, IKernel, ITrilinear };

// force distribution
enum Distribution { DNearestNeighbor, DKernel };

//drag correlation
enum DragCorrelation { ErgunWenYu, Tang, Stokes, Felice, Tenneti, NoDrag };

//lift correlation
enum LiftCorrelation { NoLift, Saffman };

//added mass correlation
enum AddedMassCorrelation { NoAM, Finn };

//effective viscosity
enum EffectiveViscosity { None, Rescaled, Eilers };


/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
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

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField  = block->getData< PdfField_T > ( pdfFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "Boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ) );

   const auto noslip = flagField->getFlag( NoSlip_Flag );

   CellInterval domainBB = storage->getDomainCellBB();

   domainBB.xMin() -= cell_idx_c( FieldGhostLayers ); // -1
   domainBB.xMax() += cell_idx_c( FieldGhostLayers ); // cellsX

   // LEFT
   CellInterval left( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMin(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( left, *block );
   handling->forceBoundary( noslip, left );

   // RIGHT
   CellInterval right( domainBB.xMax(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( right, *block );
   handling->forceBoundary( noslip, right );

   domainBB.yMin() -= cell_idx_c( FieldGhostLayers ); // 0
   domainBB.yMax() += cell_idx_c( FieldGhostLayers ); // cellsY+1

   // FRONT
   CellInterval front( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMin(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( front, *block );
   handling->forceBoundary( noslip, front );

   // BACK
   CellInterval back( domainBB.xMin(), domainBB.yMax(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( back, *block );
   handling->forceBoundary( noslip, back );

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


DPMethod to_DPMethod( const std::string& s )
{
   if( s == "GNS"  ) return DPMethod::GNS;
   throw std::runtime_error("invalid conversion from " + s + " to DPMethod");
}

Interpolation to_interpolation( const std::string& s )
{
   if( s == "nearest" ) return Interpolation::INearestNeighbor;
   if( s == "kernel"  ) return Interpolation::IKernel;
   if( s == "trilinear"  ) return Interpolation::ITrilinear;
   throw std::runtime_error("invalid conversion from " + s + " to Interpolation");
}

Distribution to_distribution( const std::string& s )
{
   if( s == "nearest" ) return Distribution::DNearestNeighbor;
   if( s == "kernel"  ) return Distribution::DKernel;
   throw std::runtime_error("invalid conversion from " + s + " to Distribution");
}

DragCorrelation to_dragCorrelation( const std::string& s )
{
   if( s == "ergun" ) return DragCorrelation::ErgunWenYu;
   if( s == "tang" ) return DragCorrelation::Tang;
   if( s == "stokes" ) return DragCorrelation::Stokes;
   if( s == "felice" ) return DragCorrelation::Felice;
   if( s == "tenneti" ) return DragCorrelation::Tenneti;
   if( s == "none" ) return DragCorrelation::NoDrag;
   throw std::runtime_error("invalid conversion from " + s + " to DragCorrelation");
}

LiftCorrelation to_liftCorrelation( const std::string& s )
{
   if( s == "none" ) return LiftCorrelation::NoLift;
   if( s == "saffman" ) return LiftCorrelation::Saffman;
   throw std::runtime_error("invalid conversion from " + s + " to LiftCorrelation");
}

AddedMassCorrelation to_addedMassCorrelation( const std::string& s )
{
   if( s == "none" ) return AddedMassCorrelation::NoAM;
   if( s == "finn" ) return AddedMassCorrelation::Finn;
   throw std::runtime_error("invalid conversion from " + s + " to AddedMassCorrelation");
}

EffectiveViscosity to_effvisc( const std::string& s )
{
   if( s == "none" ) return EffectiveViscosity::None;
   if( s == "rescaled"  ) return EffectiveViscosity::Rescaled;
   if( s == "eilers"  ) return EffectiveViscosity::Eilers;
   throw std::runtime_error("invalid conversion from " + s + " to effective viscosity");
}

std::string dpm_to_string( const DPMethod& dpm )
{
   if( dpm == DPMethod::GNS ) return "GNS";
   throw std::runtime_error("invalid conversion from DPMethod to string");
}

class QuantityEvaluator
{
public:

   QuantityEvaluator( SweepTimeloop* timeloop,  StructuredBlockStorage & blocks,
                      const BlockDataID & bodyStorageID, real_t diameter, real_t densityRatio, real_t Galileo,
                      const std::string & baseFolder, bool writeLogging ) :
      timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ), diameter_( diameter ),
      terminalVelocity_(0), evaluateTerminalVelocity_( true ),
      writeLogging_( writeLogging )
   {
      if( writeLogging_ )
      {
         fileName_ = baseFolder+"/evalCollisionBehaviorDPM_"+std::to_string(uint_c(real_t(10) * densityRatio))+"_"+std::to_string(uint_c(real_t(10) * Galileo))+"_"+std::to_string(uint_c(real_t(10) * diameter))+".txt";
         std::ofstream file;
         file.open( fileName_.c_str() );
         file << "#t z velz x y velx vely fx fy fz\n";
         file.close();
      }
   }

   void operator()()
   {
      // get particle position
      Vector3<real_t> particlePosition = getParticlePosition();
      // get particle velocity
      Vector3<real_t> particleVelocity = getMeanParticleVelocity();
      // get force on particles
      Vector3<real_t> particleForce = getParticleForce();

      if( evaluateTerminalVelocity_ ) terminalVelocity_ = std::min(particleVelocity[2], terminalVelocity_);

      WALBERLA_ROOT_SECTION()
      {
         positionsOverTime_.push_back(particlePosition[2]);
         velocitiesOverTime_.push_back(particleVelocity[2]);
         if( writeLogging_ ) writeToFile( particlePosition, particleVelocity, particleForce );
      }
   }

   real_t getVelocity()
   {
      return getMeanParticleVelocity()[2];
   }

   real_t getTerminalVelocity()
   {
      return terminalVelocity_;
   }

   real_t getPosition()
   {
      return getParticlePosition()[2];
   }

   void stopTerminalVelocityEvaluation()
   {
      evaluateTerminalVelocity_ = false;
   }

   void writeScaledOutputToFile(uint_t tImpact, real_t terminalVelocity)
   {
      std::ofstream file;
      file.open( fileName_.c_str() );
      file.precision(8);
      file << "# t position velocity\n";
      for(uint_t t = 0; t < positionsOverTime_.size(); ++t)
      {
         file << (real_c(t) - real_c(tImpact)) * terminalVelocity / diameter_ << " " << ( positionsOverTime_[t] - real_t(0.5) * diameter_ ) / diameter_ << " " << velocitiesOverTime_[t] / terminalVelocity << "\n";

      }
      file.close();

   }

private:

   Vector3<real_t> getParticleForce()
   {
      Vector3<real_t> force(0.);
      for( auto blockIt = blocks_.begin(); blockIt != blocks_.end(); ++blockIt )
      {
         for( auto bodyIt = pe::BodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID_); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            force += bodyIt->getForce();
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( force[0], mpi::SUM );
         mpi::allReduceInplace( force[1], mpi::SUM );
         mpi::allReduceInplace( force[2], mpi::SUM );

      }
      return force;
   }

   Vector3<real_t> getMeanParticleVelocity()
   {
      Vector3<real_t> velocity(0.);
      for( auto blockIt = blocks_.begin(); blockIt != blocks_.end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            velocity += bodyIt->getLinearVel();
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( velocity[0], mpi::SUM );
         mpi::allReduceInplace( velocity[1], mpi::SUM );
         mpi::allReduceInplace( velocity[2], mpi::SUM );

      }

      return velocity;
   }

   Vector3<real_t> getParticlePosition()
   {
      Vector3<real_t> position(0.);
      uint_t counter( 0 );
      for( auto blockIt = blocks_.begin(); blockIt != blocks_.end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            position = bodyIt->getPosition();
            ++counter;
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( position[0], mpi::SUM );
         mpi::allReduceInplace( position[1], mpi::SUM );
         mpi::allReduceInplace( position[2], mpi::SUM );
         mpi::allReduceInplace( counter, mpi::SUM );
      }

      WALBERLA_CHECK_EQUAL(counter, uint_t(1), "No more sphere in domain!");

      return position;
   }

   void writeToFile( const Vector3<real_t> & position, const Vector3<real_t> & vel, const Vector3<real_t> & force )
   {
      std::ofstream file;
      file.open( fileName_.c_str(), std::ofstream::app );
      file.precision(8);

      file << timeloop_->getCurrentTimeStep() << " \t"  << position[2] << " \t" << vel[2] << " \t "
           << position[0] << " \t " << position[1] << " \t "
           << vel[0] << " \t " << vel[1] << " \t "
           << force[0] << " \t " << force[1] << " \t " << force[2] << "\n";
      file.close();
   }

   SweepTimeloop* timeloop_;
   StructuredBlockStorage & blocks_;
   const BlockDataID bodyStorageID_;
   real_t diameter_;
   real_t terminalVelocity_;
   std::string fileName_;
   std::vector<real_t> positionsOverTime_;
   std::vector<real_t> velocitiesOverTime_;
   bool evaluateTerminalVelocity_;
   bool writeLogging_;

};


class CollisionPropertiesEvaluator
{
public:
   CollisionPropertiesEvaluator( pe::cr::ICR & collisionResponse ) : collisionResponse_( collisionResponse ), maximumPenetration_(real_t(0))
   {}

   void operator()()
   {
      real_t maxPen = collisionResponse_.getMaximumPenetration();
      maximumPenetration_ = std::max( maximumPenetration_, maxPen );
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( maximumPenetration_, mpi::MAX );
      }
   }
   real_t getMaximumPenetrationInSimulation()
   {
      return maximumPenetration_;
   }
private:
   pe::cr::ICR & collisionResponse_;
   real_t maximumPenetration_;
};

void emptyFunction(){}

//*******************************************************************************************************************
/*!\brief Testcase that evaluates the collision behavior of a settling sphere and a resting wall in a viscous fluid
 *
 * A glass sphere of radius D is settling under gravity towards the bottom plane.
 * The domain size is [32 x 32 x 512] * D, to ensure that the terminal settling velocity is reached before it hits the wall.
 *     _______________
 *    |               |
 *    |       o       | |
 *    |       |       | | gravity (z-direction)
 *    |       |       | v
 *    |       |       |
 *    |       |       |
 *    |_______|_______|
 *
 * The sphere properties are: dry coeff of restitution = 0.97, coeff of friction = 0.1
 * Depending on the Stokes number ( St = ( rho_s / rho_f ) * Re / 9 ), the actual coefficient of resitution changes
 * due to the presence of the viscous fluid.
 * The Reynolds number is evaluated as Re = u_T * D / visc, with the terminal settling velocity u_T.
 * Thus, this quantity is not a-priori available.
 * The main input parameter is the Galileo number, which, together with the density ratio, determines the fluid viscosity.
 *
 * Many parts of the discrete particle method can be modified via command line arguments:
 *  - discrete particle method
 *  - interpolation method
 *  - distribution method
 *  - number of interaction subcycles
 *  - number of pe steps per interaction subcycle
 *  - correlations for drag, lift, added mass
 *  - effective viscosity
 *  - turbulence model
 *  - lubrication force and its cut off distance
 *
 * An extensive description of the different components and this test case is given in Rettinger, Ruede (2017).
 *
 * The actual coeff of restitution can be compared to experimental findings from Joseph et al and Gondret et al.
 *
 *
 * References:
 * Experimental:
 * - G. G. Joseph, R. Zenit, M. L. Hunt, A. M. Rosenwinkel - "Particle–wall collisions in a viscous fluid", Journal of Fluid
 *   Mechanics 433 (2001) 329–346. doi:10.1017/S0022112001003470.
 * - P. Gondret, M. Lance, L. Petit - "Bouncing motion of spherical particles in fluids", Physics of Fluids
 *   14 (2) (2002) 643–652. doi:10.1063/1.1427920.
 *
 * Fully resolved simulations:
 * - A. Kidanemariam, M. Uhlmann - "Direct numerical simulation of pattern formation in subaqueous sediment", Journal of
 *   Fluid Mechanics 750. doi:10.1017/jfm.2014.284.
 *
 * Discrete particle simulations:
 * - C. Rettinger, U. Ruede - "A Coupled Lattice Boltzmann Method and Discrete Element Method for Discrete Particle
 *   Simulations of Particulate Flows". arXiv preprint arXiv:1711.00336 (2017)
 * - R. Sun, H. Xiao - "SediFoam: A general–purpose, open–source CFD-–DEM solver for particle–laden flow with emphasis
 *   on sediment transport", Computers & Geosciences 89 (2016) 207–219. doi:10.1016/j.cageo.2016.01.011.
 * - J. R. Finn, M. Li, S. V. Apte - "Particle based modelling and simulation of natural sand dynamics in the wave bottom
 *   boundary layer", Journal of Fluid Mechanics 796 (2016) 340–385. doi:10.1017/jfm.2016.246.
 *
 */
//*******************************************************************************************************************
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   /////////////////////////////////////////////
   //                                         //
   //   Command Line Argument Customization   //
   //                                         //
   /////////////////////////////////////////////

   bool funcTest = false;
   bool fileIO   = false;
   std::string baseFolder = "vtk_out_SphereWallDPM";

   real_t gravity = real_t(1e-4);
   real_t densityRatio = real_t(2.0);
   real_t diameter = real_t(0.5);
   real_t GalileoNumber = real_t(30.9);
   uint_t interactionSubCycles = uint_t(1); // number of subcycles that involve evaluation of the interaction force
   uint_t peSubSteps = uint_t(1); // number of pe only calls in each subcycle
   real_t collisionTime = real_t(1);

   DPMethod dpm = DPMethod::GNS;
   Interpolation interpol = Interpolation::IKernel;
   Distribution dist = Distribution::DKernel;
   DragCorrelation dragCorr = DragCorrelation::Tenneti;
   LiftCorrelation liftCorr = LiftCorrelation ::NoLift;
   AddedMassCorrelation addedMassCorr = AddedMassCorrelation::NoAM;
   EffectiveViscosity effVisc = EffectiveViscosity::None;
   bool useTurbulenceModel = false;
   real_t lubricationCutOffDistance = real_t(0); //0 switches it off, should be <= diameter for sphere-wall collision, and <= diameter/2 for sphere-sphere collision
   bool useLubricationCorrection = false; // false: use full lubrication force, true: use only correction part
   const real_t smagorinskyConstant = real_t(0.1); //for turbulence model

   uint_t dimlessTimesteps = uint_t(500);

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--funcTest" )                      == 0 ) funcTest = true;
      else if( std::strcmp( argv[i], "--fileIO" )                   == 0 ) fileIO = true;
      else if( std::strcmp( argv[i], "--baseFolder" )               == 0 ) baseFolder = argv[++i];
      else if( std::strcmp( argv[i], "--gravity" )                  == 0 ) gravity = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--Galileo" )                  == 0 ) GalileoNumber = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--densityRatio" )             == 0 ) densityRatio = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--diameter" )                 == 0 ) diameter = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--collisionTime" )            == 0 ) collisionTime = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--interactionSubCycles" )     == 0 ) interactionSubCycles = uint_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--peSubSteps" )               == 0 ) peSubSteps = uint_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--DPMethod" )                 == 0 ) dpm = to_DPMethod( argv[++i] );
      else if( std::strcmp( argv[i], "--interpolation" )            == 0 ) interpol = to_interpolation( argv[++i] );
      else if( std::strcmp( argv[i], "--distribution" )             == 0 ) dist = to_distribution( argv[++i] );
      else if( std::strcmp( argv[i], "--dragCorrelation" )          == 0 ) dragCorr = to_dragCorrelation( argv[++i] );
      else if( std::strcmp( argv[i], "--liftCorrelation" )          == 0 ) liftCorr = to_liftCorrelation( argv[++i] );
      else if( std::strcmp( argv[i], "--addedMassCorrelation" )     == 0 ) addedMassCorr = to_addedMassCorrelation( argv[++i] );
      else if( std::strcmp( argv[i], "--effectiveViscosity" )       == 0 ) effVisc = to_effvisc( argv[++i] );
      else if( std::strcmp( argv[i], "--useTurbulenceModel" )       == 0 ) useTurbulenceModel = true;
      else if( std::strcmp( argv[i], "--useLubricationCorrection" ) == 0 ) useLubricationCorrection = true;
      else if( std::strcmp( argv[i], "--lubricationCutOff" )        == 0 ) lubricationCutOffDistance = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--timesteps" )                == 0 ) dimlessTimesteps = uint_c( std::atof( argv[++i] ) );
      else WALBERLA_ABORT("Found invalid command line argument: \"" << argv[i] << "\" - aborting...");
   }

   WALBERLA_CHECK( diameter <= real_t(1), "Diameter is not allowed to be > 1!" );
   WALBERLA_CHECK( interactionSubCycles > uint_t(0), "Number of interaction sub cycles has to be at least 1!");
   WALBERLA_CHECK( peSubSteps > uint_t(0), "Number of pe sub steps has to be at least 1!");
   WALBERLA_CHECK( lubricationCutOffDistance >= real_t(0), "Lubrication cut off distance has to be non-negative!");

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

   ///////////////////////////////
   //                           //
   //   SIMULATION PROPERTIES   //
   //                           //
   ///////////////////////////////

   // roughly resembles the experimental setup from Gondret et al (2002)
   const uint_t xlength = uint_t( real_t(32) * diameter);
   const uint_t ylength = uint_t( real_t(32) * diameter);
   const uint_t zlength = uint_t(real_t(512) * diameter);
   const real_t dx = real_t(1);

   const real_t ug = std::sqrt(( densityRatio - real_t(1)) * gravity * diameter );
   const real_t viscosity = ug * diameter / GalileoNumber;
   const real_t tau = real_t(1) / lbm::collision_model::omegaFromViscosity(viscosity);

   const real_t sphereVolume = math::pi * diameter * diameter * diameter / real_t(6);
   Vector3<real_t> gravitationalForce ( real_t(0), real_t(0), ( densityRatio - real_t(1) ) * sphereVolume * gravity );

   if( !funcTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT("Lx x Ly x Lz = " << xlength << " x " << ylength << " x " << zlength);
      WALBERLA_LOG_INFO_ON_ROOT("density ratio = " << densityRatio );
      WALBERLA_LOG_INFO_ON_ROOT("Ga = " << GalileoNumber );
      WALBERLA_LOG_INFO_ON_ROOT("diameter = " << diameter );
      WALBERLA_LOG_INFO_ON_ROOT("tau = " << tau );
      WALBERLA_LOG_INFO_ON_ROOT("viscosity = " << viscosity );
      WALBERLA_LOG_INFO_ON_ROOT("u_g = " << ug );
      WALBERLA_LOG_INFO_ON_ROOT("gravity = " << gravity );
      WALBERLA_LOG_INFO_ON_ROOT("total external (gravity & buoyancy) force on sphere = " << gravitationalForce );
      WALBERLA_LOG_INFO_ON_ROOT("-------------------------------------------------------------------------------" );
      WALBERLA_LOG_INFO_ON_ROOT("discrete particle method = " << dpm_to_string(dpm) );
      WALBERLA_LOG_INFO_ON_ROOT("interpolator = " << interpol );
      WALBERLA_LOG_INFO_ON_ROOT("distribution = " << dist );
      WALBERLA_LOG_INFO_ON_ROOT("dragCorrelation = " << dragCorr );
      WALBERLA_LOG_INFO_ON_ROOT("addedMassCorrelation = " << addedMassCorr );
      WALBERLA_LOG_INFO_ON_ROOT("liftCorrelation = " << liftCorr );
      WALBERLA_LOG_INFO_ON_ROOT("effective viscosity = " << effVisc );
      WALBERLA_LOG_INFO_ON_ROOT("turbulence model = " << ( useTurbulenceModel ? "yes" : "no" ) );
      WALBERLA_LOG_INFO_ON_ROOT("interaction sub cycles = " << interactionSubCycles );
      WALBERLA_LOG_INFO_ON_ROOT("pe sub steps = " << peSubSteps );
      WALBERLA_LOG_INFO_ON_ROOT("lubrication cut off distance = " << lubricationCutOffDistance );
      WALBERLA_LOG_INFO_ON_ROOT("use lubrication correction term instead of full formula = " << ( useLubricationCorrection ? "yes" : "no" ) );
      WALBERLA_LOG_INFO_ON_ROOT("dt_DEM = " << real_t(1) / real_c(interactionSubCycles * peSubSteps) );
   }


   const real_t dt = real_t(1);
   const real_t dtInteractionSubCycle = dt / real_c(interactionSubCycles);
   const real_t dtBodyVelocityTimeDerivativeEvaluation = dtInteractionSubCycle;


   const real_t tStokes = densityRatio * diameter * diameter / ( real_t(18) * viscosity );
   const uint_t timesteps = (funcTest) ? uint_t(3) : dimlessTimesteps * uint_c(tStokes); // total number of time steps for the whole simulation


   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = uint_t(1);
   const uint_t YBlocks = uint_t(1);
   const uint_t ZBlocks = uint_t(32);

   const uint_t XCells  = xlength / XBlocks;
   const uint_t YCells  = ylength / YBlocks;
   const uint_t ZCells  = zlength / ZBlocks;

   if( XBlocks * XCells != xlength ||
       YBlocks * YCells != ylength ||
       ZBlocks * ZCells != zlength ) WALBERLA_ABORT("Domain decomposition failed!");

   auto blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells,
                                                      dx, 0, false, false,
                                                      false, false, false,
                                                      false );


   ////////
   // PE //
   ////////

   shared_ptr<pe::BodyStorage>  globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID         = blocks->addBlockData(pe::ccd::createHashGridsDataHandling( globalBodyStorage, bodyStorageID ), "CCD");
   auto fcdID         = blocks->addBlockData(pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   pe::cr::ICR* cr;
   pe::cr::DEM cr_dem(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID );
   cr = &cr_dem;

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_t( 1.5 ) * dx;
   auto syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );

   // create the sphere
   const real_t restitutionCoeff = real_t(0.97);
   const real_t frictionCoeff = real_t(0.1);

   const real_t scaledCollisionTime = collisionTime * diameter / dx; // smaller diameter -> smaller collision time to keep maximum penetration small -> more pe substeps to resolve it

   const real_t particleMass = densityRatio * sphereVolume;
   const real_t Mij = particleMass; // * particleMass / ( real_t(2) * particleMass ); // Mij = M for sphere-wall collision
   const real_t lnDryResCoeff = std::log(restitutionCoeff);
   const real_t stiffnessCoeff = math::pi * math::pi * Mij / ( scaledCollisionTime * scaledCollisionTime * ( real_t(1) - lnDryResCoeff * lnDryResCoeff / ( math::pi * math::pi + lnDryResCoeff* lnDryResCoeff ) ) );
   const real_t normalizedStiffnessCoeff = stiffnessCoeff / ( ( densityRatio - real_t(1) ) * gravity * sphereVolume / diameter );

   const real_t dampingCoeff = - real_t(2) * std::sqrt( Mij * stiffnessCoeff ) *
                               ( std::log(restitutionCoeff) / std::sqrt( math::pi * math::pi + (std::log(restitutionCoeff) * std::log(restitutionCoeff) ) ) );
   const real_t contactDuration = real_t(2) * math::pi * Mij / ( std::sqrt( real_t(4) * Mij * stiffnessCoeff - dampingCoeff * dampingCoeff )); //formula from Uhlman
   const real_t contactDuration2 = std::sqrt(( math::pi * math::pi + std::log(restitutionCoeff) * std::log(restitutionCoeff)) / ( stiffnessCoeff / Mij)); //formula from Finn

   if( !funcTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Created material with:\n"
                                       << " - coefficient of restitution = " << restitutionCoeff << "\n"
                                       << " - coefficient of friction = " << frictionCoeff << "\n"
                                       << " - normalized stiffness coefficient kn* = " << normalizedStiffnessCoeff << "\n"
                                       << " - stiffness coefficient kn = " << stiffnessCoeff << "\n"
                                       << " - damping coefficient cdn = " << dampingCoeff << "\n"
                                       << " - contact time Tc = " << contactDuration << "\n"
                                       << " - contact time Tc2 = " << contactDuration2);
   }

   auto peMaterial = pe::createMaterial( "sedimentMat", densityRatio, restitutionCoeff, frictionCoeff, frictionCoeff, real_t(0), real_t(200), stiffnessCoeff, dampingCoeff, dampingCoeff );
   real_t zPosition = real_c(zlength) - real_t(3) * diameter;
   pe::createSphere( *globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, Vector3<real_t>(real_c(xlength) * real_t(0.5), real_c(ylength) * real_t(0.5), zPosition), diameter * real_t(0.5), peMaterial );
   pe::createPlane( *globalBodyStorage, 0, Vector3<real_t>( 0, 0, 1 ), Vector3<real_t>( 0, 0, 0 ), peMaterial );


   //////////////////////
   //                  //
   //    BLOCK DATA    //
   //                  //
   //////////////////////

   // create force field
   BlockDataID forceFieldID = field::addToStorage< Vec3Field_T >( blocks, "force field", Vector3<real_t>(real_t(0)), field::zyxf, FieldGhostLayers );

   BlockDataID dragForceFieldID = field::addToStorage< Vec3Field_T >( blocks, "drag force field", Vector3<real_t>(real_t(0)), field::zyxf, FieldGhostLayers );
   BlockDataID amForceFieldID = field::addToStorage< Vec3Field_T >( blocks, "am force field", Vector3<real_t>(real_t(0)), field::zyxf, FieldGhostLayers );
   BlockDataID liftForceFieldID = field::addToStorage< Vec3Field_T >( blocks, "lift force field", Vector3<real_t>(real_t(0)), field::zyxf, FieldGhostLayers );

   // create omega field
   BlockDataID omegaFieldID = field::addToStorage< ScalarField_T >( blocks, "omega field", real_t(0), field::zyxf, FieldGhostLayers );

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omegaFieldID, ForceModel_T( forceFieldID ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel, Vector3<real_t>(0), real_t(1), FieldGhostLayers, field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldID, pdfFieldID ), "boundary handling" );

   // field to store fluid velolcity
   BlockDataID velocityFieldID = field::addToStorage< Vec3Field_T >( blocks, "velocity field", Vector3<real_t>(0), field::zyxf, FieldGhostLayers );
   BlockDataID swappedOldVelocityFieldID = field::addToStorage< Vec3Field_T >( blocks, "swapped old velocity field", Vector3<real_t>(0), field::zyxf, FieldGhostLayers );

   // create pressure field
   BlockDataID pressureFieldID = field::addToStorage< ScalarField_T >( blocks, "pressure field", real_t(0), field::zyxf, FieldGhostLayers );

   // create solid volume fraction field
   BlockDataID svfFieldID = field::addToStorage< ScalarField_T >( blocks, "svf field", real_t(0), field::zyxf, FieldGhostLayers );

   // field to store pressure gradient
   BlockDataID pressureGradientFieldID = field::addToStorage< Vec3Field_T >( blocks, "pressure gradient field", Vector3<real_t>(real_c(0)), field::zyxf, FieldGhostLayers );

   // field to store curl of fluid velocity
   BlockDataID velocityCurlFieldID = field::addToStorage< Vec3Field_T >( blocks, "velocity curl field", Vector3<real_t>(real_c(0)), field::zyxf, FieldGhostLayers );

   // field to store velocity gradient
   BlockDataID velocityGradientFieldID = field::addToStorage< TensorField_T >( blocks, "velocity gradient field", Matrix3<real_t>(real_c(0)), field::zyxf, FieldGhostLayers );

   // field to store gradient of stress tensor
   BlockDataID stressTensorGradientFieldID = field::addToStorage< Vec3Field_T >( blocks, "stress tensor gradient field", Vector3<real_t>(real_c(0)), field::zyxf, FieldGhostLayers );

   // field to store time derivative of fluid velocity
   BlockDataID timeDerivativeVelocityFieldID = field::addToStorage< Vec3Field_T >( blocks, "time derivative velocity field", Vector3<real_t>(real_c(0)), field::zyxf, FieldGhostLayers );

   // communication schemes
   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<Vec3Field_T> forceComm( blocks, forceFieldID );

   blockforest::communication::UniformBufferedScheme<stencil::D3Q27> pdfScheme( blocks );
   pdfScheme.addPackInfo( make_shared< field::communication::PackInfo<PdfField_T> >( pdfFieldID ) );

   // setup of the communication for synchronizing the velocity field between neighboring blocks
   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > velocityCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   velocityCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<Vec3Field_T> >( velocityFieldID ) );

   // setup of the communication for synchronizing the solid volume fraction field between neighboring blocks which takes into account the svf values inside the ghost layers
   shared_ptr<pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<ScalarField_T> >svfCommunicationScheme = make_shared<pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<ScalarField_T> >( blocks, svfFieldID );

   // setup of the communication for synchronizing the pressure field between neighboring blocks
   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > pressureCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   pressureCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<ScalarField_T> >( pressureFieldID ) );

   // setup of the communication for synchronizing the pressure gradient field between neighboring blocks
   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > pressureGradientCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   pressureGradientCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<Vec3Field_T> >( pressureGradientFieldID ) );

   // setup of the communication for synchronizing the velocity curl field between neighboring blocks
   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > velocityCurlCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   velocityCurlCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<Vec3Field_T> >( velocityCurlFieldID ) );

   // communication for synchronizing the velocity gradient values
   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > velocityGradientCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   velocityGradientCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<TensorField_T> >( velocityGradientFieldID ) );

   // communication for synchronizing the stress tensor gradient values
   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > stressTensorGradientCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   stressTensorGradientCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<Vec3Field_T> >( stressTensorGradientFieldID ) );

   // communication for synchronizing the interaction force field
   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<Vec3Field_T> dragForceComm( blocks, dragForceFieldID );
   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<Vec3Field_T> amForceComm( blocks, amForceFieldID );
   pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<Vec3Field_T> liftForceComm( blocks, liftForceFieldID );

   shared_ptr<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> dragForceFieldToForceFieldAdder =
            make_shared<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder>(blocks, forceFieldID, dragForceFieldID, uint_t(1) );
   shared_ptr<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> amForceFieldToForceFieldAdder =
         make_shared<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder>(blocks, forceFieldID, amForceFieldID, uint_t(1) );
   shared_ptr<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> liftForceFieldToForceFieldAdder =
         make_shared<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder>(blocks, forceFieldID, liftForceFieldID, uint_t(1) );

   /////////////////////////////////
   //                             //
   //    CORRELATION FUNCTIONS    //
   //                             //
   /////////////////////////////////

   // drag correlation function
   std::function<Vector3<real_t> ( const Vector3<real_t>&, const Vector3<real_t> &, real_t, real_t, real_t, real_t)> dragCorrelationFunction;
   if( dragCorr == DragCorrelation::ErgunWenYu )
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceErgunWenYu;
   }
   else if( dragCorr == DragCorrelation::Tang )
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceTang;
   }
   else if( dragCorr == DragCorrelation::Stokes )
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceStokes;
   }
   else if( dragCorr == DragCorrelation::Felice )
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceFelice;
   }
   else if( dragCorr == DragCorrelation::Tenneti )
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::dragForceTenneti;
   }
   else if( dragCorr == DragCorrelation::NoDrag )
   {
      dragCorrelationFunction = pe_coupling::discrete_particle_methods::noDragForce;
   }
   else
   {
      WALBERLA_ABORT("Drag correlation not yet implemented!");
   }

   // lift correlation function
   std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t, real_t )> liftCorrelationFunction;
   if( liftCorr == LiftCorrelation::NoLift )
   {
      liftCorrelationFunction = pe_coupling::discrete_particle_methods::noLiftForce;
   }
   else if( liftCorr == LiftCorrelation::Saffman )
   {
      liftCorrelationFunction = pe_coupling::discrete_particle_methods::liftForceSaffman;
   }
   else
   {
      WALBERLA_ABORT("Lift correlation not yet implemented!");
   }

   // added mass correlation function
   std::function<Vector3<real_t> ( const Vector3<real_t> &, const Vector3<real_t> &, real_t, real_t )> addedMassCorrelationFunction;
   if( addedMassCorr == AddedMassCorrelation::NoAM )
   {
      addedMassCorrelationFunction = pe_coupling::discrete_particle_methods::noAddedMassForce;
   }
   else if( addedMassCorr == AddedMassCorrelation::Finn )
   {
      addedMassCorrelationFunction = pe_coupling::discrete_particle_methods::addedMassForceFinn;
   }
   else
   {
      WALBERLA_ABORT("Added mass correlation not yet implemented!");
   }

   // set up effective viscosity calculation
   std::function<real_t ( real_t, real_t)> effectiveViscosityFunction;
   if( effVisc == EffectiveViscosity::None )
   {
      effectiveViscosityFunction = pe_coupling::discrete_particle_methods::calculateUnchangedEffectiveViscosity;
   }
   else if( effVisc == EffectiveViscosity::Rescaled )
   {
      effectiveViscosityFunction = pe_coupling::discrete_particle_methods::calculateRescaledEffectiveViscosity;
   }
   else if( effVisc == EffectiveViscosity::Eilers )
   {
      effectiveViscosityFunction = pe_coupling::discrete_particle_methods::calculateEilersEffectiveViscosity;
   }
   else
   {
      WALBERLA_ABORT("Effective viscosity not yet implemented!");
   }
   shared_ptr<pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator> effectiveViscosityEvaluator =
         walberla::make_shared<pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator>(omegaFieldID, svfFieldID, viscosity, effectiveViscosityFunction );


   ////////////////////////////////
   //                            //
   //    EVALUATION FUNCTIONS    //
   //                            //
   ////////////////////////////////

   // evaluator for bodies' velocity time derivative
   shared_ptr<pe_coupling::discrete_particle_methods::BodyVelocityTimeDerivativeEvaluator> bodyVelocityTimeDerivativeEvaluator = make_shared<pe_coupling::discrete_particle_methods::BodyVelocityTimeDerivativeEvaluator>( blocks, bodyStorageID, dtBodyVelocityTimeDerivativeEvaluation );
   (*bodyVelocityTimeDerivativeEvaluator)();

   // function used to evaluate the interaction force between fluid and particles
   std::function<void(void)> dragAndPressureForceEvaluationFunction;
   if( dpm == DPMethod::GNS ) {
      if (interpol == Interpolation::INearestNeighbor) {
         if (dist == Distribution::DNearestNeighbor) {
            typedef pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor> IFE_T;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         } else if (dist == Distribution::DKernel) {
            typedef pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor> IFE_T;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         }
      } else if (interpol == Interpolation::IKernel) {
         if (dist == Distribution::DNearestNeighbor) {
            typedef pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor> IFE_T;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         } else if (dist == Distribution::DKernel) {
            typedef pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor> IFE_T;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         }
      } else if (interpol == Interpolation::ITrilinear) {
         if (dist == Distribution::DNearestNeighbor) {
            typedef pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::NearestNeighborDistributor> IFE_T;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         } else if (dist == Distribution::DKernel) {
            typedef pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::KernelDistributor> IFE_T;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         }
      }
   }
   else
   {
      WALBERLA_ABORT("Discrete particle method not yet implemented!");
   }


   // function to evaluate the lift force contribution
   std::function<void(void)> liftForceEvaluationFunction;
   if( interpol == Interpolation::INearestNeighbor )
   {
      if( dist == Distribution::DNearestNeighbor )
      {
         typedef pe_coupling::discrete_particle_methods::LiftForceEvaluator< FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         typedef pe_coupling::discrete_particle_methods::LiftForceEvaluator< FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if( interpol == Interpolation::IKernel )
   {
      if( dist == Distribution::DNearestNeighbor )
      {
         typedef pe_coupling::discrete_particle_methods::LiftForceEvaluator< FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         typedef pe_coupling::discrete_particle_methods::LiftForceEvaluator< FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if( interpol == Interpolation::ITrilinear )
   {
      if( dist == Distribution::DNearestNeighbor )
      {
         typedef pe_coupling::discrete_particle_methods::LiftForceEvaluator< FlagField_T, field::TrilinearFieldInterpolator, field::NearestNeighborDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         typedef pe_coupling::discrete_particle_methods::LiftForceEvaluator< FlagField_T, field::TrilinearFieldInterpolator, field::KernelDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }


   // function to evaluate the added mass contribution
   std::function<void(void)> addedMassEvaluationFunction;
   if( interpol == Interpolation::INearestNeighbor )
   {
      if( dist == Distribution::DNearestNeighbor )
      {
         typedef pe_coupling::discrete_particle_methods::AddedMassForceEvaluator< FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         typedef pe_coupling::discrete_particle_methods::AddedMassForceEvaluator< FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if( interpol == Interpolation::IKernel )
   {
      if( dist == Distribution::DNearestNeighbor )
      {
         typedef pe_coupling::discrete_particle_methods::AddedMassForceEvaluator< FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         typedef pe_coupling::discrete_particle_methods::AddedMassForceEvaluator< FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }
   else if( interpol == Interpolation::ITrilinear )
   {
      if( dist == Distribution::DNearestNeighbor )
      {
         typedef pe_coupling::discrete_particle_methods::AddedMassForceEvaluator< FlagField_T, field::TrilinearFieldInterpolator, field::NearestNeighborDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         typedef pe_coupling::discrete_particle_methods::AddedMassForceEvaluator< FlagField_T, field::TrilinearFieldInterpolator, field::KernelDistributor > IFE_T;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
   }

   // function to evaluate lubrication forces
   std::function<void(void)> lubricationEvaluationFunction;
   if( lubricationCutOffDistance > real_t(0) )
   {
      if( useLubricationCorrection )
      {
         using LE_T = pe_coupling::LubricationCorrection;
         shared_ptr<LE_T> lubEval = make_shared<LE_T>( blocks, globalBodyStorage, bodyStorageID, viscosity, lubricationCutOffDistance );
         lubricationEvaluationFunction = std::bind(&LE_T::operator(), lubEval);
      }
      else
      {
         using LE_T = pe_coupling::discrete_particle_methods::LubricationForceEvaluator;
         shared_ptr<LE_T> lubEval = make_shared<LE_T>( blocks, globalBodyStorage, bodyStorageID, viscosity, lubricationCutOffDistance );
         lubricationEvaluationFunction = std::bind(&LE_T::operator(), lubEval);
      }
   } else {
      lubricationEvaluationFunction = emptyFunction;
   }


   ///////////////////////////
   //                       //
   //    CREATE TIMELOOP    //
   //                       //
   ///////////////////////////

   // initial SVF calculation
   if( dist == Distribution::DNearestNeighbor )
   {
      pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::NearestNeighborDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );
   } else {
      pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::KernelDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );
   }
   (*svfCommunicationScheme)();

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // evaluation functions for simulation
   shared_ptr<QuantityEvaluator> quantityEvaluator = walberla::make_shared< QuantityEvaluator >( &timeloop, *blocks, bodyStorageID, diameter, densityRatio, GalileoNumber, baseFolder, fileIO );
   shared_ptr<CollisionPropertiesEvaluator> collisionPropertiesEvaluator = walberla::make_shared<CollisionPropertiesEvaluator>( *cr );

   if( addedMassCorr != AddedMassCorrelation::NoAM )
   {
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::FieldDataSwapper<Vec3Field_T>( velocityFieldID, swappedOldVelocityFieldID ), "Velocity Field Swap" );
   }

   // subcycling loop begin
   for( uint_t subcycle = 1; subcycle <= interactionSubCycles; ++subcycle )
   {

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (liftForceFieldToForceFieldAdder ), "Lift Force Field To Force Field Adder" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (amForceFieldToForceFieldAdder ), "AM Force Field To Force Field Adder" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Old Velocity Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCommunicationScheme), "Old Velocity Field Communication" );

      if( liftCorr != LiftCorrelation::NoLift )
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::VelocityCurlFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( velocityCurlFieldID, velocityFieldID, boundaryHandlingID ), "Velocity Curl Field Evaluation" )
                        << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCurlCommunicationScheme), "Velocity Curl Field Communication" );
      }

      if( addedMassCorr != AddedMassCorrelation::NoAM )
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::VelocityGradientFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( velocityGradientFieldID, velocityFieldID, boundaryHandlingID ), "Velocity Gradient Field Evaluation" )
                        << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityGradientCommunicationScheme), "Velocity Gradient Field Communication" );
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::VelocityTotalTimeDerivativeFieldEvaluator( timeDerivativeVelocityFieldID, velocityFieldID, swappedOldVelocityFieldID, velocityGradientFieldID, dt ), "Velocity Time Derivative Field Evaluation" );
      }

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSPressureFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( pressureFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Pressure Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(pressureCommunicationScheme), "Pressure Field Communication" );

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::PressureGradientFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( pressureGradientFieldID, pressureFieldID, boundaryHandlingID ), "Pressure Gradient Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(pressureGradientCommunicationScheme), "Pressure Gradient Field Communication" );


      if( liftCorr != LiftCorrelation::NoLift )
      {
         // evaluate Flift
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( liftForceFieldID ), "Lift Force Field Reset" )
                        << AfterFunction( liftForceEvaluationFunction, "Lift Force Evaluation")
                        << AfterFunction( liftForceComm, "Lift Force Field Communication" );
      }

      if( addedMassCorr != AddedMassCorrelation::NoAM )
      {

         // evaluate Fam
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( amForceFieldID ), "AM Force Field Reset" )
                        << AfterFunction( addedMassEvaluationFunction, "Added Mass Force Evaluation")
                        << AfterFunction( amForceComm, "Force Field Communication" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::BodyVelocityTimeDerivativeEvaluator>( bodyVelocityTimeDerivativeEvaluator ), "Body Velocity Time Derivative Evaluation" );
      }

      if( liftCorr != LiftCorrelation::NoLift || addedMassCorr != AddedMassCorrelation::NoAM )
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (liftForceFieldToForceFieldAdder ), "Lift Force Field To Force Field Adder" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (amForceFieldToForceFieldAdder ), "AM Force Field To Force Field Adder" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Velocity Field Evaluation" )
                        << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCommunicationScheme), "Velocity Field Communication" );
      }

      // evaluate Fdrag
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( dragForceFieldID ), "Drag Force Field Reset" )
                     << AfterFunction( dragAndPressureForceEvaluationFunction, "Fluid-Particle Interaction Force Evaluation" )
                     << AfterFunction( dragForceComm, "Drag Force Field Communication" )
                     << AfterFunction( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, -gravitationalForce  ), "Gravitational and buoyancy force" )
                     << AfterFunction( pe_coupling::TimeStep( blocks, bodyStorageID, *cr, syncCall, dtInteractionSubCycle, peSubSteps, lubricationEvaluationFunction ), "Pe Time Step" )
                     << AfterFunction( SharedFunctor<CollisionPropertiesEvaluator>( collisionPropertiesEvaluator ), "Collision properties evaluator");

      if( dist == Distribution::DNearestNeighbor )
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::NearestNeighborDistributor> ( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag ), "Solid Volume Fraction Field Evaluation" )
                        << AfterFunction( SharedFunctor< pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<ScalarField_T> >(svfCommunicationScheme), "Solid Volume Fraction Field Communication" );
      }
      else
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::KernelDistributor> ( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag ), "Solid Volume Fraction Field Evaluation" )
                        << AfterFunction( SharedFunctor< pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<ScalarField_T> >(svfCommunicationScheme), "Solid Volume Fraction Field Communication" );
      }

   }

   timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" )
                  << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (liftForceFieldToForceFieldAdder ), "Lift Force Field To Force Field Adder" )
                  << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (amForceFieldToForceFieldAdder ), "AM Force Field To Force Field Adder" )
                  << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

   timeloop.add() << Sweep( makeSharedSweep<pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator>( effectiveViscosityEvaluator ), "Effective Viscosity Evaluation");
   if( useTurbulenceModel ) timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSSmagorinskyLESField<LatticeModel_T>(blocks, omegaFieldID, pdfFieldID, svfFieldID, smagorinskyConstant), "Turbulence Model" );

   // execute GNS-LBM sweep, boundary handling, PDF communication
   auto sweep = pe_coupling::discrete_particle_methods::makeGNSSweep< LatticeModel_T, FlagField_T >( pdfFieldID, svfFieldID, flagFieldID, Fluid_Flag );

   timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "GNS-LBM sweep (collide)" );

   timeloop.add() << BeforeFunction( pdfScheme, "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   timeloop.add() << Sweep( lbm::makeStreamSweep( sweep ), "GNS-LBM sweep (stream)" );

   timeloop.addFuncAfterTimeStep( SharedFunctor<QuantityEvaluator>( quantityEvaluator ), "Quantity Evaluator" );

   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );


   //////////////////////
   //                  //
   //    SIMULATION    //
   //                  //
   //////////////////////

   WcTimingPool timeloopTiming;

   real_t currentVelocity(0);
   real_t oldVelocity(0);
   uint_t numberOfBouncesUntilTermination = uint_t(4);
   uint_t countedBounces = uint_t(0);
   uint_t tImpact = uint_t(0);
   real_t terminalVelocity(0);
   real_t maxPositionAfterBounce = real_t(0);

   // three different evaluation times for the rebound velocity (see Kidanemariam, Uhlmann (2014) )
   real_t evalOffsetFactorS = real_t(0.05);
   real_t evalOffsetFactorM = real_t(0.1);
   real_t evalOffsetFactorL = real_t(0.15);
   real_t reboundVelocityS( 0. );
   real_t reboundVelocityM( 0. );
   real_t reboundVelocityL( 0. );
   uint_t tEvalOffsetS = 0;
   uint_t tEvalOffsetM = 0;
   uint_t tEvalOffsetL = 0;

   for( size_t t = 0; t < timesteps; ++t)
   {
      timeloop.singleStep(timeloopTiming);
      currentVelocity = quantityEvaluator->getVelocity();
      if( currentVelocity > real_t(0) && oldVelocity < real_t(0) )
      {
         ++countedBounces;
         if( countedBounces == 1 )
         {
            tImpact = t;
            terminalVelocity = std::fabs(quantityEvaluator->getTerminalVelocity());
            quantityEvaluator->stopTerminalVelocityEvaluation();
            tEvalOffsetS = uint_c(std::ceil(evalOffsetFactorS * diameter / terminalVelocity));
            tEvalOffsetM = uint_c(std::ceil(evalOffsetFactorM * diameter / terminalVelocity));
            tEvalOffsetL = uint_c(std::ceil(evalOffsetFactorL * diameter / terminalVelocity));
         }
         WALBERLA_LOG_INFO_ON_ROOT("----------- " << countedBounces << ". bounce detected ----------------")
      }
      if( currentVelocity > oldVelocity && oldVelocity > real_t(0) && countedBounces == 1)
      {
         // impact was wrongly detected in an intermediate step, but should be the time when the collision is fully resolved and maximal velocity is reached
         ++tImpact;
      }
      if( t == tImpact + tEvalOffsetS ) reboundVelocityS = currentVelocity;
      if( t == tImpact + tEvalOffsetM ) reboundVelocityM = currentVelocity;
      if( t == tImpact + tEvalOffsetL ) reboundVelocityL = currentVelocity;
      if( countedBounces >= numberOfBouncesUntilTermination) break;
      if( countedBounces >= 1 )
      {
         maxPositionAfterBounce = std::max(maxPositionAfterBounce, quantityEvaluator->getPosition());
      }
      if( countedBounces >= 1 && ( real_c(t-tImpact) * terminalVelocity / diameter ) > real_t(100) ) break;

      oldVelocity = currentVelocity;
   }

   maxPositionAfterBounce -= diameter / real_t(2);

   timeloopTiming.logResultOnRoot();

   if( !funcTest )
   {
      WALBERLA_LOG_INFO_ON_ROOT("v_T = " << terminalVelocity);
      real_t Re = terminalVelocity * diameter / viscosity;
      WALBERLA_LOG_INFO_ON_ROOT("Re_T = " << Re);
      real_t St = densityRatio * Re / real_t(9);
      WALBERLA_LOG_INFO_ON_ROOT("density ratio = " << densityRatio);
      WALBERLA_LOG_INFO_ON_ROOT("St = " << St);


      WALBERLA_LOG_INFO_ON_ROOT("tI = " << tImpact )

      WALBERLA_LOG_INFO_ON_ROOT(" - tRS = " <<  evalOffsetFactorS * diameter / terminalVelocity << " = "  << tEvalOffsetS  );
      real_t effectiveCoefficientOfRestitutionS = reboundVelocityS / terminalVelocity;
      WALBERLA_LOG_INFO_ON_ROOT("   e_eff = " << effectiveCoefficientOfRestitutionS << ", e/e_dry = " << effectiveCoefficientOfRestitutionS / restitutionCoeff);

      WALBERLA_LOG_INFO_ON_ROOT(" - tRM = " <<  evalOffsetFactorM * diameter / terminalVelocity << " = "  << tEvalOffsetM  );
      real_t effectiveCoefficientOfRestitutionM = reboundVelocityM / terminalVelocity;
      WALBERLA_LOG_INFO_ON_ROOT("   e_eff = " << effectiveCoefficientOfRestitutionM << ", e/e_dry = " << effectiveCoefficientOfRestitutionM / restitutionCoeff);

      WALBERLA_LOG_INFO_ON_ROOT(" - tRL = " <<  evalOffsetFactorL * diameter / terminalVelocity << " = "  << tEvalOffsetL  );
      real_t effectiveCoefficientOfRestitutionL = reboundVelocityL / terminalVelocity;
      WALBERLA_LOG_INFO_ON_ROOT("   e_eff = " << effectiveCoefficientOfRestitutionL << ", e/e_dry = " << effectiveCoefficientOfRestitutionL / restitutionCoeff);

      WALBERLA_LOG_INFO_ON_ROOT("maximum position after first bounce (zMax) = " << maxPositionAfterBounce << ", zMax / D = " << maxPositionAfterBounce / diameter );
      real_t maxPen = collisionPropertiesEvaluator->getMaximumPenetrationInSimulation();
      WALBERLA_LOG_INFO_ON_ROOT("maximum penetration (maxPen) = " <<  maxPen << ", maxPen / D = " << real_t(100) * maxPen / diameter << "%");

      if( fileIO ) quantityEvaluator->writeScaledOutputToFile(tImpact, terminalVelocity);

      //log to file
      if( fileIO ) {
         WALBERLA_ROOT_SECTION() {
            std::string fileName = baseFolder + "/logCollisionBehavior.txt";
            std::ofstream file;
            file.open(fileName.c_str(), std::ofstream::app);
            file.precision(8);

            file << densityRatio << " \t" << GalileoNumber << " \t" << diameter << " \t"
                 << tau << " \t" << gravity << " \t" << normalizedStiffnessCoeff << " \t" << interpol << " \t" << dist << " \t"
                 << addedMassCorr << " \t" << liftCorr << " \t" << effVisc << " \t" << ((useTurbulenceModel) ? 1 : 0) << " \t"
                 << interactionSubCycles << " \t" << peSubSteps << " \t" << lubricationCutOffDistance << " \t"
                 << terminalVelocity << " \t" << Re << " \t" << St << " \t"
                 << effectiveCoefficientOfRestitutionM << " \t" << restitutionCoeff << " \t"
                 << maxPositionAfterBounce << " \t " << maxPen << "\n";
         }
      }

   }

   return 0;
}

} //namespace sphere_wall_collision_behavior_dpm

int main( int argc, char **argv ){
   sphere_wall_collision_behavior_dpm::main(argc, argv);
}

