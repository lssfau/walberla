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
//! \file HinderedSettlingDynamicsDPM.cpp
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
#include "core/logging/all.h"

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

#include <functional>
#include <vector>
#include <iomanip>
#include <iostream>
#include <random>

namespace hindered_settling_dynamics_dpm
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

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const/* storage */ ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;

}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const /*storage*/ ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField  = block->getData< PdfField_T > ( pdfFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( "Boundary handling", flagField, fluid,
                                                           NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ) );

   handling->fillWithDomain( FieldGhostLayers );
   return handling;
}
//*******************************************************************************************************************


// VTK Output
shared_ptr<vtk::VTKOutput> createFluidFieldVTKWriter( shared_ptr< StructuredBlockForest > & blocks,
                                                      const BlockDataID & pdfFieldID, const BlockDataID & flagFieldID, uint_t writeFrequency,
                                                      const std::string & vtkBaseFolder )
{
    auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, "fluid_field", writeFrequency, uint_t(1), false, vtkBaseFolder );

    blockforest::communication::UniformBufferedScheme<stencil::D3Q27> pdfGhostLayerSync( blocks );
    pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo<PdfField_T> >( pdfFieldID ) );
    pdfFieldVTKWriter->addBeforeFunction( pdfGhostLayerSync );

    field::FlagFieldCellFilter< FlagField_T > fluidFilter( flagFieldID );
    fluidFilter.addFlag( Fluid_Flag );
    pdfFieldVTKWriter->addCellInclusionFilter( fluidFilter );

    auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "Velocity (Lattice)" );
    auto  densityWriter = make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "Density (Lattice)" );
    pdfFieldVTKWriter->addCellDataWriter( velocityWriter );
    pdfFieldVTKWriter->addCellDataWriter( densityWriter );

    return pdfFieldVTKWriter;
}

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

uint_t createSpheresRandomly( StructuredBlockForest & forest, pe::BodyStorage & globalBodyStorage,
                              const BlockDataID & bodyStorageID, const AABB & generationDomain,
                              real_t diameter, real_t solidVolumeFraction, pe::MaterialID & material, real_t initialZVelocity )
{
   real_t domainVolume = generationDomain.volume();
   real_t totalSphereVolume = domainVolume * solidVolumeFraction;
   real_t sphereVolume = diameter * diameter * diameter * math::pi / real_t(6);
   uint_t numberOfSpheres = uint_c( totalSphereVolume / sphereVolume );

   real_t xParticle = real_t(0);
   real_t yParticle = real_t(0);
   real_t zParticle = real_t(0);

   for( uint_t nSphere = 0; nSphere < numberOfSpheres; ++nSphere )
   {

      WALBERLA_ROOT_SECTION()
      {
         xParticle = math::realRandom<real_t>(generationDomain.xMin(), generationDomain.xMax());
         yParticle = math::realRandom<real_t>(generationDomain.yMin(), generationDomain.yMax());
         zParticle = math::realRandom<real_t>(generationDomain.zMin(), generationDomain.zMax());
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::broadcastObject( xParticle );
         mpi::broadcastObject( yParticle );
         mpi::broadcastObject( zParticle );
      }


      pe::SphereID sp = pe::createSphere( globalBodyStorage, forest.getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( xParticle, yParticle, zParticle ), diameter * real_t(0.5), material );
      if( sp != nullptr )
      {
         sp->setLinearVel(Vector3<real_t>(real_t(0),real_t(0),initialZVelocity));
      }
   }

   return numberOfSpheres;
}

uint_t createSphereLattice( StructuredBlockForest & forest, pe::BodyStorage & globalBodyStorage,
                            const BlockDataID & bodyStorageID, const AABB & generationDomain,
                            real_t diameter, real_t solidVolumeFraction, pe::MaterialID & material, real_t initialZVelocity )
{
   real_t sphereVolume = math::pi * diameter * diameter * diameter / real_t(6);
   real_t numSpheresDesired = solidVolumeFraction * generationDomain.volume() / sphereVolume;
   uint_t spheresPerDirection = uint_c(std::cbrt(numSpheresDesired) );

   real_t spacing = generationDomain.xSize() / real_c(spheresPerDirection);

   WALBERLA_ASSERT( spacing >= diameter );

   Vector3<real_t> generationOrigin( generationDomain.xMin() + spacing * real_t(0.5), generationDomain.yMin() + spacing * real_t(0.5), generationDomain.zMin() + spacing * real_t(0.5));

   uint_t numSpheres( 0 );

   for( auto it = grid_generator::SCIterator(generationDomain, generationOrigin, spacing); it != grid_generator::SCIterator(); ++it )
   {
      pe::SphereID sp = pe::createSphere( globalBodyStorage, forest.getBlockStorage(), bodyStorageID, 0, *it, diameter * real_t(0.5), material );

      if( sp != nullptr )
      {
         sp->setLinearVel(Vector3<real_t>(real_t(0),real_t(0),initialZVelocity));
         ++numSpheres;
      }
   }

   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( numSpheres, mpi::SUM );
   }

   WALBERLA_LOG_INFO_ON_ROOT("Created spheres in lattice arrangement with a center-to-center spacing of " << spacing << " ( " << real_t(100)*spacing/diameter << "% of diameter )" );

   return numSpheres;
}

void resetSphereVelocities( const shared_ptr<StructuredBlockForest> & blocks, const BlockDataID & bodyStorageID, Vector3<real_t> vel )
{
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      for( auto bodyIt = pe::BodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID); bodyIt != pe::BodyIterator::end<pe::Sphere>(); ++bodyIt )
      {
         bodyIt->setAngularVel(real_t(0),real_t(0),real_t(0));
         bodyIt->setLinearVel(vel);
      }
   }
}

class QuantityEvaluator
{
public:

   QuantityEvaluator( SweepTimeloop* timeloop,  StructuredBlockStorage & blocks,
                      const BlockDataID & bodyStorageID, const BlockDataID & forceFieldID, const BlockDataID & pdfFieldID,
                      const std::string & fileName, bool fileIO, const uint_t & numSpheres, const real_t & velExpected,
                      const std::function<Vector3<real_t> ()>  & evaluateMeanFluidVelocity ) :
      timeloop_( timeloop ), blocks_( blocks ), bodyStorageID_( bodyStorageID ), forceFieldID_( forceFieldID ),
      pdfFieldID_( pdfFieldID ), fileName_( fileName ), fileIO_( fileIO ),
      numSpheres_( numSpheres ), velExpected_( velExpected ), evaluateMeanFluidVelocity_( evaluateMeanFluidVelocity )
   {
      if( fileIO ) {
         std::ofstream file;
         file.open(fileName_.c_str());
         file
               << "#t velExpected fX_particle fY_particle fZ_particle fX_fluid fY_fluid fZ_fluid velX_particle velY_particle velZ_particle velX_fluid velY_fluid velZ_fluid velX_rel velY_rel velZ_rel\n";
         file.close();
      }
   }

   void operator()() {
      // get mean particle velocity
      Vector3<real_t> particleVelocity = getMeanParticleVelocity();
      // get force on particles
      Vector3<real_t> particleForce = getParticleForce();
      // get force on fluid
      Vector3<real_t> fluidForce = getFluidForce();
      //get mean fluid velocity
      Vector3<real_t> fluidVelocity = evaluateMeanFluidVelocity_();

      WALBERLA_LOG_INFO_ON_ROOT(
            "mean particle force = " << particleForce[2] << ", mean fluid velocity = " << fluidVelocity[2]
                                     << ", mean fluid force = " << fluidForce[2]);

      if( fileIO_ )
      {
         WALBERLA_ROOT_SECTION()
         {
            writeToFile(particleVelocity, particleForce, fluidForce, fluidVelocity);
         }
      }
   }

   real_t getSettlingVelocity()
   {
      Vector3<real_t> particleVelocity = getMeanParticleVelocity();
      Vector3<real_t> fluidVelocity = evaluateMeanFluidVelocity_();
      return particleVelocity[2] - fluidVelocity[2];
   }

   Vector3<real_t> getMeanFluidVelocity()
   {
      return evaluateMeanFluidVelocity_();
   }

   Vector3<real_t> getMeanParticleVelocity()
   {
      Vector3<real_t> velocity(0.);
      uint_t counter = 0;
      for( auto blockIt = blocks_.begin(); blockIt != blocks_.end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            velocity += bodyIt->getLinearVel();
            ++counter;
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( velocity[0], mpi::SUM );
         mpi::allReduceInplace( velocity[1], mpi::SUM );
         mpi::allReduceInplace( velocity[2], mpi::SUM );
         mpi::allReduceInplace( counter, mpi::SUM );
      }

      WALBERLA_CHECK_EQUAL( counter, numSpheres_, "Total number of spheres in the domain has changed!" );

      return velocity / real_c(counter);
   }

   Vector3<real_t> getMeanForceOnParticles()
   {
      Vector3<real_t> force = getParticleForce();
      return force;
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
   Vector3<real_t> getFluidForce()
   {

      Vector3<real_t> force(0.);
      for( auto blockIt = blocks_.begin(); blockIt != blocks_.end(); ++blockIt )
      {
         Vec3Field_T* forceField = blockIt->getData< Vec3Field_T >( forceFieldID_ );
         WALBERLA_FOR_ALL_CELLS_XYZ( forceField,
            force += forceField->get(x,y,z);
         );
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( force[0], mpi::SUM );
         mpi::allReduceInplace( force[1], mpi::SUM );
         mpi::allReduceInplace( force[2], mpi::SUM );
      }
      return force;
   }

   void writeToFile( const Vector3<real_t> & meanParticleVel, const Vector3<real_t> & particleForce, const Vector3<real_t> & fluidForce,
                     const Vector3<real_t> & meanFluidVel )
   {
      std::ofstream file;
      file.open( fileName_.c_str(), std::ofstream::app );
      file.precision(8);

      file << timeloop_->getCurrentTimeStep() << "\t " << velExpected_ << "\t "
           << particleForce[0] << "\t " << particleForce[1] << "\t " << particleForce[2] << "\t "
           << fluidForce[0] << "\t " << fluidForce[1] << "\t " << fluidForce[2] << "\t "
           << meanParticleVel[0] << "\t " << meanParticleVel[1] << "\t " << meanParticleVel[2] << "\t "
           << meanFluidVel[0] << "\t " << meanFluidVel[1] << "\t " << meanFluidVel[2] << "\t"
           << meanParticleVel[0] - meanFluidVel[0] << "\t " << meanParticleVel[1] - meanFluidVel[1] << "\t " << meanParticleVel[2] - meanFluidVel[2] << "\n";
      file.close();
   }

   SweepTimeloop* timeloop_;
   StructuredBlockStorage & blocks_;
   const BlockDataID bodyStorageID_;
   const BlockDataID forceFieldID_;
   const BlockDataID pdfFieldID_;
   std::string fileName_;
   bool fileIO_;
   const uint_t numSpheres_;
   const real_t velExpected_;

   std::function<Vector3<real_t> ()> evaluateMeanFluidVelocity_;


};

Vector3<real_t> getGNSMeanFluidVelocity( const shared_ptr<StructuredBlockStorage> & blocks, BlockDataID pdfFieldID, BlockDataID svfFieldID, real_t domainVolume )
{
   Vector3<real_t> velocity( real_t(0) );
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      PdfField_T* pdfField = blockIt->getData< PdfField_T >( pdfFieldID );
      ScalarField_T* svfField = blockIt->getData< ScalarField_T >( svfFieldID );
      WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,
         real_t svf = svfField->get(x,y,z);
         velocity += pdfField->getVelocity(x,y,z) / (real_t(1) - svf );
      );
   }
   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( velocity[0], mpi::SUM );
      mpi::allReduceInplace( velocity[1], mpi::SUM );
      mpi::allReduceInplace( velocity[2], mpi::SUM );
   }
   return velocity / domainVolume;
}

void logSingleResultToFile( const std::string & fileName, const real_t & solidVolumeFraction, const real_t & settlingVel,
                            const real_t & expectedVel, const real_t & velUnhindered, const real_t & ReynoldsNumber, const real_t & tau )
{
   WALBERLA_ROOT_SECTION()
   {
      std::ofstream file;
      file.open( fileName.c_str(), std::ofstream::app );
      file.precision(8);

      file << solidVolumeFraction << " \t" << settlingVel << " \t" << expectedVel << " \t"
           << settlingVel/velUnhindered << " \t" << expectedVel/velUnhindered << " \t" << ReynoldsNumber << " \t" << tau << "\n";
   }
}


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
   void resetMaximumPenetration()
   {
      maximumPenetration_ = real_t(0);
   }
private:
   pe::cr::ICR & collisionResponse_;
   real_t maximumPenetration_;
};

class DummySweep
{
public:
   DummySweep( )
   = default;

   void operator()(IBlock * const /*block*/)
   {}
};

void emptyFunction(){}

//*******************************************************************************************************************
/*!\brief Testcase that evaluates the hindered settling velocity of several spheres in a periodic box filled with fluid
 *
 * Sphere of radius D are randomly generated inside the domain to obtain a specific global solid volume fraction
 * and then start to settle under gravity.
 * The domain size is [32 x 32 x 32] * D and fully periodic.
 *     _______________
 *    | o  o   o   o  |
 *    |   o   o       | |
 *    | o        o    | | gravity (z-direction)
 *    |    o   o     o| v
 *    |o   o    o   o |
 *    |  o   o    o   |
 *    |_______o_____o_|
 *
 * An external force is applied onto the fluid in opposite gravitational direction to counteract the settling motion and
 * avoid infinitely large settling velocity (since no walls are present).
 *
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
 *
 * References:
 * Experimental:
 * - J. Richardson, W. Zaki - "The sedimentation of a suspension of uniform spheres under conditions of viscous flow",
 *   Chemical Engineering Science 3 (2) (1954) 65–73. doi:10.1016/0009-2509(54)85015-9.
 * - T. Baldock, M. Tomkins, P. Nielsen, M. Hughes - "Settling velocity of sediments at high concentrations",
 *   Coastal Engineering 51 (1) (2004) 91–100. doi:10.1016/j.coastaleng.2003.12.004.
 *
 * Discrete particle simulations:
 * - C. Rettinger, U. Ruede - "A Coupled Lattice Boltzmann Method and Discrete Element Method for Discrete Particle
 *   Simulations of Particulate Flows". arXiv preprint arXiv:1711.00336 (2017)
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

   bool funcTest  = false;
   bool vtkIO     = false;
   bool fileIO    = false;
   uint_t vtkWriteFrequency = 0;
   std::string baseFolder = "vtk_out_HinderedSettlingDPM";

   real_t densityRatio = real_t(2500) / real_t(1000);
   real_t gravity = real_t(0.002); // has to be small enough to keep settling velocities small
   real_t diameter = real_t(0.5);
   uint_t interactionSubCycles = uint_t(1); // number of subcycles that involve evaluation of the interaction force -> around 3 for stability with added mass
   uint_t peSubSteps = uint_t(1); // number of pe only calls in each subcycle
   real_t solidVolumeFraction = real_t(0.05);

   DPMethod dpm = DPMethod::GNS;
   Interpolation interpol = Interpolation::IKernel;
   Distribution dist = Distribution::DKernel;
   DragCorrelation dragCorr = DragCorrelation::Tenneti;
   LiftCorrelation liftCorr = LiftCorrelation ::NoLift;
   AddedMassCorrelation addedMassCorr = AddedMassCorrelation::NoAM;
   EffectiveViscosity effVisc = EffectiveViscosity::None;

   bool useTurbulenceModel = false;
   const real_t smagorinskyConstant = real_t(0.1); //for turbulence model

   real_t lubricationCutOffDistance = real_t(0); //0 switches it off
   bool useLubricationCorrection = false; // false: use full lubrication force, true: use only correction part

   bool createSpheresInLattice = false;
   bool initialSimulationToAdjustFluidForcing = false;

   real_t initialSphereVelocityZ( real_t(0) );

   real_t relativeVelDiffLimit( real_t(1e-5) ); //set negative to avoid convergence before 30 dimensionless timesteps

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--funcTest" )                      == 0 ) funcTest = true;
      else if( std::strcmp( argv[i], "--vtkIOFreq" )                == 0 ) vtkWriteFrequency = uint_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--fileIO" )                   == 0 ) fileIO = true;
      else if( std::strcmp( argv[i], "--baseFolder" )               == 0 ) baseFolder = argv[++i];
      else if( std::strcmp( argv[i], "--gravity" )                  == 0 ) gravity = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--diameter" )                 == 0 ) diameter = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--solidVolumeFraction" )      == 0 ) solidVolumeFraction = real_c( std::atof( argv[++i] ) );
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
      else if( std::strcmp( argv[i], "--createSpheresInLattice" )   == 0 ) createSpheresInLattice = true;
      else if( std::strcmp( argv[i], "--initialSimulation" )        == 0 ) initialSimulationToAdjustFluidForcing = true;
      else if( std::strcmp( argv[i], "--convergenceLimit" )         == 0 ) relativeVelDiffLimit = real_c( std::atof( argv[++i] ) );
      else if( std::strcmp( argv[i], "--initialSphereVelocityZ" )   == 0 ) initialSphereVelocityZ = real_c( std::atof( argv[++i] ) );
      else WALBERLA_ABORT("Found invalid command line argument \"" << argv[i] << "\" - aborting...");
   }

   if( vtkWriteFrequency > 0 ) vtkIO = true;

   WALBERLA_CHECK( diameter <= real_t(1), "Diameter is not allowed to be > 1!" );
   WALBERLA_CHECK( solidVolumeFraction <= real_t(0.65), "Solid volume fraction is not allowed to be > 0.65!" );
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

   const uint_t xlength = uint_c( real_t(32) * diameter )  ;
   const uint_t ylength = xlength;
   const uint_t zlength = xlength;

   if( solidVolumeFraction < real_t(1e-10) )
   {
      // create only a single sphere
      solidVolumeFraction = math::pi * diameter * diameter * diameter / ( real_t(6) * real_c( xlength * ylength * zlength));
   }

   const real_t diameter_SI = real_t(0.00035); // m, Finn et al, Tab 5
   const real_t gravity_SI = real_t(9.81); // m/s^2

   const real_t dx_SI = diameter_SI / diameter;
   const real_t dx = real_t(1);
   const real_t viscosity_SI = real_t(1e-3); // kg/(ms)
   const real_t densityFluid_SI = real_t(1e3); // kg/m^3

   const real_t dt_SI = std::sqrt(gravity * dx_SI / gravity_SI);
   const real_t viscosity = ( viscosity_SI/densityFluid_SI ) * dt_SI / ( dx_SI * dx_SI );
   const real_t tau = real_t(1) / lbm::collision_model::omegaFromViscosity( viscosity );

   real_t gravitationalForce = - gravity * ( densityRatio - real_t(1) ) * diameter * diameter * diameter * math::pi / real_t(6);

   // unhindered settling velocity of a single sphere in infinite fluid, would come from experiments or DNS, NOT Stokes settling velocity, from Finn et al, Tab 5
   const real_t velUnhindered_SI = real_t(-0.048); // m/s

   const real_t velStokes = -( densityRatio - real_t(1) ) * diameter * diameter * gravity / ( real_t(18) * viscosity );
   const real_t velUnhindered = velUnhindered_SI * dt_SI / dx_SI;

   const real_t dt_DEM = real_t(1) / real_c(interactionSubCycles * peSubSteps);

   const real_t dt = real_t(1);
   const real_t dtInteractionSubCycle = dt / real_c(interactionSubCycles);
   const real_t dtBodyVelocityTimeDerivativeEvaluation = dtInteractionSubCycle;

   const uint_t tStokes = uint_c(densityRatio * diameter * diameter / ( real_t(18) * viscosity ));

   const uint_t timesteps = (funcTest) ? uint_t(10) : uint_t(30) * tStokes; // total number of time steps for the whole simulation

   const Vector3<real_t> initialFluidVelocity( real_t(0) );

   const std::string fileNameLoggingInit = baseFolder+"/evalHinderedSettlingDynamicsSubgrid_eps"+std::to_string(uint_c(real_t(100) * solidVolumeFraction))+"_d"+std::to_string(uint_c(real_t(100) * diameter))+"_init.txt";
   const std::string fileNameLogging = baseFolder+"/evalHinderedSettlingDynamicsSubgrid_eps"+std::to_string(uint_c(real_t(100) * solidVolumeFraction))+"_d"+std::to_string(uint_c(real_t(100) * diameter))+".txt";

   if( !funcTest ) {
      WALBERLA_LOG_INFO_ON_ROOT("Lx x Ly x Lz = " << xlength << " x " << ylength << " x " << zlength);
      WALBERLA_LOG_INFO_ON_ROOT("tau = " << tau);
      WALBERLA_LOG_INFO_ON_ROOT("viscosity = " << viscosity);
      WALBERLA_LOG_INFO_ON_ROOT("gravity = " << gravity);
      WALBERLA_LOG_INFO_ON_ROOT("total external (gravity & buoyancy) force on single sphere = " << gravitationalForce);
      WALBERLA_LOG_INFO_ON_ROOT("diameter = " << diameter);
      WALBERLA_LOG_INFO_ON_ROOT("input solid volume fraction = " << solidVolumeFraction);
      WALBERLA_LOG_INFO_ON_ROOT("discrete particle method = " << dpm_to_string(dpm));
      WALBERLA_LOG_INFO_ON_ROOT("interpolator = " << interpol);
      WALBERLA_LOG_INFO_ON_ROOT("distribution = " << dist);
      WALBERLA_LOG_INFO_ON_ROOT("dragCorrelation = " << dragCorr);
      WALBERLA_LOG_INFO_ON_ROOT("addedMassCorrelation = " << addedMassCorr);
      WALBERLA_LOG_INFO_ON_ROOT("liftCorrelation = " << liftCorr);
      WALBERLA_LOG_INFO_ON_ROOT("effective viscosity = " << effVisc);
      WALBERLA_LOG_INFO_ON_ROOT("turbulence model = " << (useTurbulenceModel ? "yes" : "no"));
      WALBERLA_LOG_INFO_ON_ROOT("interaction sub cycles = " << interactionSubCycles);
      WALBERLA_LOG_INFO_ON_ROOT("pe sub steps = " << peSubSteps);
      WALBERLA_LOG_INFO_ON_ROOT("lubrication cut off distance = " << lubricationCutOffDistance);
      WALBERLA_LOG_INFO_ON_ROOT("use lubrication correction term instead of full formula = " << (useLubricationCorrection ? "yes" : "no"));
      WALBERLA_LOG_INFO_ON_ROOT("Ga = " << std::sqrt((densityRatio - real_t(1)) * gravity * diameter * diameter * diameter) / viscosity);
      WALBERLA_LOG_INFO_ON_ROOT("dx_SI = " << dx_SI << ", dt_SI = " << dt_SI);
      WALBERLA_LOG_INFO_ON_ROOT("dt_DEM = " << dt_DEM);
      WALBERLA_LOG_INFO_ON_ROOT("t_ref = " << tStokes << " simulation steps");
      WALBERLA_LOG_INFO_ON_ROOT("logging is written to " + fileNameLogging );
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = uint_t(4);
   const uint_t YBlocks = uint_t(4);
   const uint_t ZBlocks = uint_t(4);

   const uint_t XCells  = xlength / XBlocks;
   const uint_t YCells  = ylength / YBlocks;
   const uint_t ZCells  = zlength / ZBlocks;

   if( XBlocks * XCells != xlength ||
       YBlocks * YCells != ylength ||
       ZBlocks * ZCells != zlength ) WALBERLA_ABORT("Domain decomposition failed!");

   auto blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells,
                                                      dx, 0, false, false,
                                                      true, true, true,
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

   const real_t restitutionCoeff = real_t(0.88);
   const real_t frictionCoeff = real_t(0.25);

   real_t sphereVolume = diameter * diameter * diameter * math::pi / real_t(6);
   const real_t particleMass = densityRatio * sphereVolume;
   const real_t Mij = particleMass * particleMass / ( real_t(2) * particleMass );
   const real_t lnDryResCoeff = std::log(restitutionCoeff);
   const real_t collisionTime = real_t(0.5);
   const real_t stiffnessCoeff = math::pi * math::pi * Mij / ( collisionTime * collisionTime * ( real_t(1) - lnDryResCoeff * lnDryResCoeff / ( math::pi * math::pi + lnDryResCoeff* lnDryResCoeff ) ) );
   const real_t dampingCoeff = - real_t(2) * std::sqrt( Mij * stiffnessCoeff ) *
                               ( std::log(restitutionCoeff) / std::sqrt( math::pi * math::pi + (std::log(restitutionCoeff) * std::log(restitutionCoeff) ) ) );
   const real_t contactDuration = real_t(2) * math::pi * Mij / ( std::sqrt( real_t(4) * Mij * stiffnessCoeff - dampingCoeff * dampingCoeff )); //formula from Uhlmann

   if( !funcTest ) {
      WALBERLA_LOG_INFO_ON_ROOT("Created sediment material with:\n"
                                      << " - coefficient of restitution = " << restitutionCoeff << "\n"
                                      << " - coefficient of friction = " << frictionCoeff << "\n"
                                      << " - stiffness coefficient kn = " << stiffnessCoeff << "\n"
                                      << " - damping coefficient cdn = " << dampingCoeff << "\n"
                                      << " - contact time Tc = " << contactDuration);
   }
   auto peMaterial = pe::createMaterial( "sedimentMat", densityRatio, restitutionCoeff, frictionCoeff, frictionCoeff, real_t(0), real_t(200), stiffnessCoeff, dampingCoeff, dampingCoeff );

   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_t( 1.5 ) * dx;
   auto syncCall = std::bind( pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()), bodyStorageID, static_cast<WcTimingTree*>(nullptr), overlap, false );
   shared_ptr<CollisionPropertiesEvaluator> collisionPropertiesEvaluator = walberla::make_shared<CollisionPropertiesEvaluator>( *cr );

   // create the spheres
   uint_t numSpheres = 0;

   if ( createSpheresInLattice )
   {
      numSpheres = createSphereLattice( *blocks, *globalBodyStorage, bodyStorageID, AABB( real_t(0), real_t(0), real_t(0), real_c(xlength), real_c(ylength), real_c(zlength) ),
                                        diameter, solidVolumeFraction, peMaterial, real_t(0) );
      syncCall();

   } else {
      numSpheres = createSpheresRandomly( *blocks, *globalBodyStorage, bodyStorageID, AABB( real_t(0), real_t(0), real_t(0), real_c(xlength), real_c(ylength), real_c(zlength) ),
                                          diameter, solidVolumeFraction, peMaterial, real_t(0) );
      syncCall();

      const uint_t initialPeSteps = uint_t(50000);
      const real_t dt_DEM_init = collisionTime / real_t(10);
      const real_t overlapLimit = real_t(0.05) * diameter;
      
      WALBERLA_LOG_INFO_ON_ROOT("Sphere creation done --- resolving overlaps with goal all < " << overlapLimit / diameter * real_t(100) << "%");
      
      for( uint_t pet = uint_t(1); pet <= initialPeSteps; ++pet )
      {
         cr->timestep( dt_DEM_init );
         syncCall();
         (*collisionPropertiesEvaluator)();
         real_t maxPen = collisionPropertiesEvaluator->getMaximumPenetrationInSimulation();
         if( maxPen < overlapLimit )
         {
            WALBERLA_LOG_INFO_ON_ROOT("Carried out " << pet << " DEM-only time steps to resolve initial overlaps");
            break;
         }else{
            if( pet % uint_t(200) == uint_t(0) )
            {
               WALBERLA_LOG_INFO_ON_ROOT(pet << " - current max overlap = " << maxPen / diameter * real_t(100) << "%");
            }
         }
         collisionPropertiesEvaluator->resetMaximumPenetration();
      }
   }

   Vector3<real_t> initialSphereVelocity(real_t(0), real_t(0), initialSphereVelocityZ);
   WALBERLA_LOG_INFO_ON_ROOT("Resetting sphere velocity to " << initialSphereVelocity );
   resetSphereVelocities( blocks, bodyStorageID, initialSphereVelocity );


   const real_t domainVolume = real_c( xlength * ylength * zlength );
   real_t actualSolidVolumeFraction = real_c( numSpheres ) * diameter * diameter * diameter * math::pi / ( real_t(6) * domainVolume );
   real_t ReynoldsNumber = std::fabs(velUnhindered) * diameter / viscosity;

   // apply external forcing on fluid to approximately balance the force from the settling particles to avoid too large fluid or particle velocities
   const real_t extForceZ = actualSolidVolumeFraction * gravity * (densityRatio - real_t(1));
   Vector3<real_t> extForce = Vector3<real_t>(real_t(0), real_t(0), extForceZ );

   // apply estimate by Richardson & Zaki (1954) for unbounded flow (DomainLength->infty)
   real_t n = 0;
   if( ReynoldsNumber < real_t(0.2) )
   {
      n = real_t(4.65);
   } else if ( ReynoldsNumber < real_t(1) )
   {
      n = real_t(4.35) * std::pow( ReynoldsNumber, real_t(-0.03) );
   } else if ( ReynoldsNumber < real_t(500) )
   {
      n = real_t(4.45) * std::pow( ReynoldsNumber, real_t(-0.1) );
   } else
   {
      n = real_t(2.39);
   }

   real_t expectedVelocity = velUnhindered * std::pow( ( real_t(1) - actualSolidVolumeFraction ), n );

   WALBERLA_LOG_INFO_ON_ROOT("solid volume fraction = " << actualSolidVolumeFraction );
   WALBERLA_LOG_INFO_ON_ROOT("number of spheres = " << numSpheres );
   WALBERLA_LOG_INFO_ON_ROOT("velocity_Stokes = " << velStokes );
   WALBERLA_LOG_INFO_ON_ROOT("velocity_Unhindered = " << velUnhindered );
   WALBERLA_LOG_INFO_ON_ROOT("Re = " << ReynoldsNumber );
   WALBERLA_LOG_INFO_ON_ROOT("expected settling velocity = " << expectedVelocity <<" = " << expectedVelocity * dx_SI / dt_SI << " m/s" );
   WALBERLA_LOG_INFO_ON_ROOT("external forcing on fluid = " << extForce );
   WALBERLA_LOG_INFO_ON_ROOT("total external forcing applied on all fluid cells = " << extForce[2] * (real_t(1) - actualSolidVolumeFraction) * real_c( xlength * ylength * zlength ) );
   WALBERLA_LOG_INFO_ON_ROOT("total external (gravity & buoyancy) force on all spheres = " << gravitationalForce * real_c(numSpheres) );

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
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel, initialFluidVelocity, real_t(1), FieldGhostLayers, field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldID, pdfFieldID ), "boundary handling" );

   // field to store fluid velocity
   BlockDataID velocityFieldID = field::addToStorage< Vec3Field_T >( blocks, "velocity field", initialFluidVelocity, field::zyxf, FieldGhostLayers );
   BlockDataID oldVelocityFieldID = field::addToStorage< Vec3Field_T >( blocks, "old velocity field", initialFluidVelocity, field::zyxf, FieldGhostLayers );
   BlockDataID swappedOldVelocityFieldID = field::addToStorage< Vec3Field_T >( blocks, "swapped old velocity field", initialFluidVelocity, field::zyxf, FieldGhostLayers );

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

   shared_ptr<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> > oldVelocityCommunicationScheme = make_shared<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >( blocks );
   oldVelocityCommunicationScheme->addPackInfo( make_shared< field::communication::PackInfo<Vec3Field_T> >( oldVelocityFieldID ) );

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


   // function to evaluate the mean fluid velocity
   std::function<Vector3<real_t> ()> evaluateMeanFluidVelocity = std::bind(getGNSMeanFluidVelocity, blocks, pdfFieldID, svfFieldID, domainVolume);


   //////////////////////////////
   //                          //
   //    INITIAL SIMULATION    //
   //                          //
   //////////////////////////////

   // initial simulation: keep particles fix, evaluate acting force on particles
   // if this force is approximately converged, see if it matches the gravitational force
   // if not, change the external force accordingly

   if( solidVolumeFraction > real_t(0.01) && initialSimulationToAdjustFluidForcing )
   {
      WALBERLA_LOG_INFO_ON_ROOT("===================================================================================" );
      WALBERLA_LOG_INFO_ON_ROOT("Starting initial simulation to equilibrate fluid forcing and interaction force");

      // create the timeloop
      SweepTimeloop timeloopInit( blocks->getBlockStorage(), timesteps );


      // evaluation function
      shared_ptr<QuantityEvaluator> quantityEvaluator = walberla::make_shared< QuantityEvaluator >( &timeloopInit, *blocks, bodyStorageID, forceFieldID, pdfFieldID, fileNameLoggingInit, fileIO, numSpheres, expectedVelocity, evaluateMeanFluidVelocity );

      shared_ptr<pe_coupling::discrete_particle_methods::GNSExternalForceToForceFieldAdder> gnsExternalForceOnForceFieldAdder =
               walberla::make_shared<pe_coupling::discrete_particle_methods::GNSExternalForceToForceFieldAdder>( forceFieldID, svfFieldID, extForce );


      //initial setup
      {
         // init solid volume fraction field
         if( dist == Distribution::DNearestNeighbor )
         {
            pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::NearestNeighborDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
            for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );
         } else {
            pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::KernelDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
            for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );
         }

         (*svfCommunicationScheme)();

         pe_coupling::discrete_particle_methods::ForceFieldResetter forceFieldResetter( forceFieldID );
         for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  forceFieldResetter( &(*blockIt) );

         for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  (*gnsExternalForceOnForceFieldAdder)( &(*blockIt) );

      }

      (*quantityEvaluator)();

      ///////// begin of GNS timeloop ///////////////////////

      timeloopInit.addFuncBeforeTimeStep( pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID), "Resetting force on bodies" );


      timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::GNSPressureFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( pressureFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Pressure Field Evaluation" )
                         << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(pressureCommunicationScheme), "Pressure Field Communication" );

      timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::PressureGradientFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( pressureGradientFieldID, pressureFieldID, boundaryHandlingID ), "Pressure Gradient Field Evaluation" )
                         << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(pressureGradientCommunicationScheme), "Pressure Gradient Field Communication" );

      // subcycling loop begin
      for( uint_t subcycle = 1; subcycle <= interactionSubCycles; ++subcycle )
      {
         timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" );
         timeloopInit.add() << Sweep( makeSharedSweep(gnsExternalForceOnForceFieldAdder), "Force Field Add" )
                            << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

         timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Velocity Field Evaluation" )
                            << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCommunicationScheme), "Velocity Field Communication" );

         // evaluate Fdrag
         timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( dragForceFieldID ), "Drag Force Field Reset" )
                            << AfterFunction( dragAndPressureForceEvaluationFunction, "Fluid-Particle Interaction Force Evaluation" )
                            << AfterFunction( dragForceComm, "Drag Force Field Communication" );

         if( subcycle != interactionSubCycles )
         {
            timeloopInit.add() << Sweep( DummySweep(), "DummySweep")
                               << AfterFunction( pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID), "Resetting force on bodies" );
         }
      }

      timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" );
      timeloopInit.add() << Sweep( makeSharedSweep(gnsExternalForceOnForceFieldAdder), "Force Field Add" )
                         << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

      timeloopInit.add() << Sweep( makeSharedSweep<pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator>( effectiveViscosityEvaluator ), "Effective Viscosity Evaluation");
      if( useTurbulenceModel ) timeloopInit.add() << Sweep( pe_coupling::discrete_particle_methods::GNSSmagorinskyLESField<LatticeModel_T>(blocks, omegaFieldID, pdfFieldID, svfFieldID, smagorinskyConstant), "Turbulence Model" );

      // execute GNS-LBM sweep, boundary handling, PDF communication
      auto sweep = pe_coupling::discrete_particle_methods::makeGNSSweep< LatticeModel_T, FlagField_T >( pdfFieldID, svfFieldID, flagFieldID, Fluid_Flag );

      timeloopInit.add() << Sweep( lbm::makeCollideSweep( sweep ), "GNS-LBM sweep (collide)" );

      timeloopInit.add() << BeforeFunction( pdfScheme, "LBM Communication" )
                     << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      timeloopInit.add() << Sweep( lbm::makeStreamSweep( sweep ), "GNS-LBM sweep (stream)" );


      // evaluate the forces / velocities
      timeloopInit.addFuncAfterTimeStep( SharedFunctor<QuantityEvaluator>( quantityEvaluator ), "Quantity Evaluator" );

      // configure vtk output
      if( vtkIO )
      {
         auto fluidVTKWriter = vtk::createVTKOutput_BlockData( blocks, "fluid_field_init", vtkWriteFrequency, uint_t(1), false, baseFolder );
         blockforest::communication::UniformBufferedScheme<stencil::D3Q27> pdfGhostLayerSync( blocks );
         pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo<PdfField_T> >( pdfFieldID ) );
         fluidVTKWriter->addBeforeFunction( pdfGhostLayerSync );
         auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldID, "Velocity (Lattice)" );
         auto  densityWriter = make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldID, "Density (Lattice)" );
         fluidVTKWriter->addCellDataWriter( velocityWriter );
         fluidVTKWriter->addCellDataWriter( densityWriter );
         timeloopInit.addFuncAfterTimeStep( vtk::writeFiles( fluidVTKWriter ), "VTK (fluid field data)" );

         timeloopInit.addFuncAfterTimeStep( field::createVTKOutput<Vec3Field_T, float>( forceFieldID, *blocks, "force_field_init", vtkWriteFrequency, 1, false, baseFolder ), "VTK (force field init)" );
      }

      // execute simulation
      WcTimingPool timeloopInitTiming;

      real_t oldInteractionForce( real_t(0) );
      real_t curInteractionForce( real_t(0) );
      real_t relativeForceDiffLimit( real_t(1e-4) );
      real_t relativeForceConvergenceLimit( real_t( 1e-3 ) );
      for( uint_t t = 0; t <= timesteps; ++t )
      {
         timeloopInit.singleStep( timeloopInitTiming );
         curInteractionForce = quantityEvaluator->getMeanForceOnParticles()[2];

         if( std::isnan(curInteractionForce) ) WALBERLA_ABORT( "Nan found in interaction force!" );

         if( t > tStokes ) // leave enough timesteps to let the velocity field adapt to the new forcing
         {
            if( std::fabs((curInteractionForce - oldInteractionForce ) / curInteractionForce) < relativeForceDiffLimit || t == timesteps)
            {
               WALBERLA_LOG_INFO_ON_ROOT("initial simulation ended with relative difference of interaction forces of " << relativeForceDiffLimit << " after " << t << " time steps.");


               real_t actingExternalForceOnSpheres = real_c(numSpheres) * ( ( - gravity * densityRatio * diameter * diameter * diameter * math::pi / real_t(6)  ) +
                                                                            ( gravity * real_t(1) * diameter * diameter * diameter * math::pi / real_t(6) ) +
                                                                            ( extForce[2] * real_t(1) * diameter * diameter * diameter * math::pi / real_t(6) ) );
               WALBERLA_LOG_INFO_ON_ROOT("f_interaction_z = " << curInteractionForce << ", f_ext_z = " << actingExternalForceOnSpheres );
               if( std::fabs( ( std::fabs( curInteractionForce ) - std::fabs( actingExternalForceOnSpheres ) )/ std::fabs( curInteractionForce ) ) < relativeForceConvergenceLimit )
               {
                  WALBERLA_LOG_INFO_ON_ROOT("initial simulation converged");
                  break;
               }
               else
               {
                  t = 0u;
                  //timeloopInit.setCurrentTimeStepToZero();
                  extForce[2] = extForce[2] * ( std::fabs( actingExternalForceOnSpheres ) / std::fabs( curInteractionForce ) );
                  gnsExternalForceOnForceFieldAdder->reset(extForce);
                  WALBERLA_LOG_INFO_ON_ROOT("restarting initial simulation with new external force = " << extForce[2]);
                  curInteractionForce = real_t(0);
               }
            }
            oldInteractionForce = curInteractionForce;
         }

      }

   }

   // reset remaining force on all bodies
   pe_coupling::ForceTorqueOnBodiesResetter forceResetter( blocks, bodyStorageID);
   forceResetter();

   WALBERLA_LOG_INFO_ON_ROOT("===================================================================================" );
   WALBERLA_LOG_INFO_ON_ROOT("Starting simulation with:" );
   WALBERLA_LOG_INFO_ON_ROOT("external forcing on fluid = " << extForce );
   WALBERLA_LOG_INFO_ON_ROOT("total external forces on all particles = " << real_c(numSpheres) * ( - gravity * ( densityRatio - real_t(1) ) + extForce[2] ) * diameter * diameter * diameter * math::pi / real_t(6) );
   WALBERLA_LOG_INFO_ON_ROOT("simulating " << timesteps << " time steps" );


   /////////////////////////////
   //                         //
   //    ACTUAL SIMULATION    //
   //                         //
   /////////////////////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   // evaluation function
   shared_ptr<QuantityEvaluator> quantityEvaluator = walberla::make_shared< QuantityEvaluator >( &timeloop, *blocks, bodyStorageID, forceFieldID, pdfFieldID, fileNameLogging, fileIO, numSpheres, expectedVelocity, evaluateMeanFluidVelocity );
   {
      (*collisionPropertiesEvaluator)();
      real_t maxPen = collisionPropertiesEvaluator->getMaximumPenetrationInSimulation();
      WALBERLA_LOG_INFO_ON_ROOT("maximum penetration before the simulation (maxPen) = " <<  maxPen << ", maxPen / D = " << real_t(100) * maxPen / diameter << "%");
   }
   collisionPropertiesEvaluator->resetMaximumPenetration();

   //initial setup if no initial simulation has been carried out earlier
   if( !initialSimulationToAdjustFluidForcing )
   {
      // init solid volume fraction field
      if( dist == Distribution::DNearestNeighbor )
      {
         pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::NearestNeighborDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
         for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );
      } else {
         pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::KernelDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
         for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );
      }

      (*svfCommunicationScheme)();

      pe_coupling::discrete_particle_methods::ForceFieldResetter forceFieldResetter( forceFieldID );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  forceFieldResetter( &(*blockIt) );

      pe_coupling::discrete_particle_methods::GNSExternalForceToForceFieldAdder gnsExternalForceAdder( forceFieldID, svfFieldID, extForce );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  gnsExternalForceAdder( &(*blockIt) );
   }


   ///////// begin of GNS timeloop ///////////////////////

   if( addedMassCorr != AddedMassCorrelation::NoAM )
   {
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::FieldDataSwapper<Vec3Field_T>( velocityFieldID, swappedOldVelocityFieldID ), "Velocity Field Swap" );
   }

   if( addedMassCorr != AddedMassCorrelation::NoAM )
   {
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( oldVelocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Old Velocity Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(oldVelocityCommunicationScheme), "Old Velocity Field Communication" );
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::VelocityGradientFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( velocityGradientFieldID, oldVelocityFieldID, boundaryHandlingID ), "Velocity Gradient Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityGradientCommunicationScheme), "Velocity Gradient Field Communication" );
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::VelocityTotalTimeDerivativeFieldEvaluator( timeDerivativeVelocityFieldID, oldVelocityFieldID, swappedOldVelocityFieldID, velocityGradientFieldID, dt ), "Velocity Time Derivative Field Evaluation" );
   }

   // subcycling loop begin
   for( uint_t subcycle = 1; subcycle <= interactionSubCycles; ++subcycle )
   {
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" );
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSExternalForceToForceFieldAdder( forceFieldID, svfFieldID, extForce ), "Force Field Add" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (liftForceFieldToForceFieldAdder ), "Lift Force Field To Force Field Adder" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (amForceFieldToForceFieldAdder ), "AM Force Field To Force Field Adder" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Velocity Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCommunicationScheme), "Velocity Field Communication" );

      if( liftCorr != LiftCorrelation::NoLift )
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::VelocityCurlFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( velocityCurlFieldID, velocityFieldID, boundaryHandlingID ), "Velocity Curl Field Evaluation" )
                        << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCurlCommunicationScheme), "Velocity Curl Field Communication" );
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
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" );
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSExternalForceToForceFieldAdder( forceFieldID, svfFieldID, extForce ), "Force Field Add" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (liftForceFieldToForceFieldAdder ), "Lift Force Field To Force Field Adder" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (amForceFieldToForceFieldAdder ), "AM Force Field To Force Field Adder" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Velocity Field Evaluation" )
                        << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCommunicationScheme), "Velocity Field Communication" );
      }

      // evaluate Fdrag
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( dragForceFieldID ), "Drag Force Field Reset" )
                     << AfterFunction( dragAndPressureForceEvaluationFunction, "Fluid-Particle Interaction Force Evaluation" )
                     << AfterFunction( dragForceComm, "Drag Force Field Communication" );

      // ext forces on bodies
      timeloop.add() << Sweep( DummySweep(), "Dummy Sweep ")
                     << AfterFunction( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, Vector3<real_t>(0,0,- gravity * densityRatio * diameter * diameter * diameter * math::pi / real_t(6) )  ), "Gravitational Force Add" )
                     << AfterFunction( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, Vector3<real_t>(0,0,gravity * real_t(1) * diameter * diameter * diameter * math::pi / real_t(6) ) ), "Buoyancy Force (due to gravity) Add" )
                     << AfterFunction( pe_coupling::ForceOnBodiesAdder( blocks, bodyStorageID, Vector3<real_t>(0,0,extForce[2] * real_t(1) * diameter * diameter * diameter * math::pi / real_t(6) ) ), "Buoyancy Force (due to external fluid force) Add" )
                     << AfterFunction( pe_coupling::TimeStep( blocks, bodyStorageID, *cr, syncCall, dtInteractionSubCycle, peSubSteps, lubricationEvaluationFunction ), "Pe Time Step" );

      timeloop.add() << Sweep( DummySweep(), "Dummy Sweep ")
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

   timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" );
   timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSExternalForceToForceFieldAdder( forceFieldID, svfFieldID, extForce ), "Force Field Add" )
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

   // configure vtk output
   if( vtkIO )
   {
       shared_ptr<vtk::VTKOutput> pdfFieldVTKWriter = createFluidFieldVTKWriter( blocks, pdfFieldID, flagFieldID, vtkWriteFrequency, baseFolder );
       timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK (fluid field data)" );

       auto bodyVtkOutput   = make_shared<pe::SphereVtkOutput>( bodyStorageID, blocks->getBlockStorage() );
       auto bodyVTK   = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", vtkWriteFrequency, baseFolder );
       timeloop.addFuncAfterTimeStep( vtk::writeFiles( bodyVTK ), "VTK (sphere data)" );

       timeloop.addFuncAfterTimeStep( field::createVTKOutput<Vec3Field_T, float>( forceFieldID, *blocks, "force_field", vtkWriteFrequency, uint_t(0), false, baseFolder ), "VTK (force field)" );
       timeloop.addFuncAfterTimeStep( field::createVTKOutput<ScalarField_T, float>( svfFieldID, *blocks, "svf_field", vtkWriteFrequency, uint_t(0), false, baseFolder ), "VTK (svf field)" );
       timeloop.addFuncAfterTimeStep( field::createVTKOutput<ScalarField_T, float>( omegaFieldID, *blocks, "omega_field", vtkWriteFrequency, uint_t(0), false, baseFolder ), "VTK (omega field)" );

   }

   // execute simulation
   WcTimingPool timeloopTiming;

   real_t oldSettlingVel( real_t(0) );
   real_t curSettlingVel( real_t(0) );
   for( uint_t t = 0; t < timesteps; ++t )
   {
      timeloop.singleStep( timeloopTiming );
      curSettlingVel = quantityEvaluator->getSettlingVelocity();
      if( t > tStokes )
      {
         if( std::fabs((curSettlingVel - oldSettlingVel ) / curSettlingVel) < relativeVelDiffLimit )
         {
            WALBERLA_LOG_INFO_ON_ROOT("simulation terminated with relative difference of settling velocities of " << relativeVelDiffLimit << " after " << t << " time steps.");
            break;
         }
         oldSettlingVel = curSettlingVel;
      }
   }

   //log to file
   if( fileIO ) logSingleResultToFile( baseFolder+"/logSettlingVel.txt", actualSolidVolumeFraction, curSettlingVel, expectedVelocity, velUnhindered, ReynoldsNumber, tau );

   timeloopTiming.logResultOnRoot();

   real_t meanParticleVel = quantityEvaluator->getMeanParticleVelocity()[2];
   real_t meanFluidVel = quantityEvaluator->getMeanFluidVelocity()[2];

   WALBERLA_LOG_INFO_ON_ROOT("Result for solid volume fraction = " << actualSolidVolumeFraction << " with " << numSpheres << " spheres.");
   WALBERLA_LOG_INFO_ON_ROOT("diameter = " << diameter << ", tau = " << tau);
   WALBERLA_LOG_INFO_ON_ROOT(" - simulated settling velocity = " << curSettlingVel << ", us/uT = " << curSettlingVel / velUnhindered);
   WALBERLA_LOG_INFO_ON_ROOT(" - expected settling velocity = " << expectedVelocity << ", ue/uT = " << expectedVelocity / velUnhindered);
   WALBERLA_LOG_INFO_ON_ROOT("detailed overview:");
   WALBERLA_LOG_INFO_ON_ROOT(" - mean particle velocity  = " << meanParticleVel << " = " << meanParticleVel * dx_SI/dt_SI << " m/s ( " << std::fabs(meanParticleVel/curSettlingVel)*real_t(100) << "% of settling vel)");
   WALBERLA_LOG_INFO_ON_ROOT(" - mean fluid velocity     = " << meanFluidVel << " = " << meanFluidVel * dx_SI/dt_SI << " m/s ( " << std::fabs(meanFluidVel/curSettlingVel)*real_t(100) << "% of settling vel)");
   WALBERLA_LOG_INFO_ON_ROOT(" - mean relative velocity  = " << curSettlingVel << " = " << curSettlingVel * dx_SI/dt_SI << " m/s");
   WALBERLA_LOG_INFO_ON_ROOT(" - expected velocity (R&Z) = " << expectedVelocity << " = " << expectedVelocity * dx_SI/dt_SI << " m/s");
   
   real_t maxPen = collisionPropertiesEvaluator->getMaximumPenetrationInSimulation();
   WALBERLA_LOG_INFO_ON_ROOT(" - maximum penetration (maxPen) = " <<  maxPen << ", maxPen / D = " << real_t(100) * maxPen / diameter << "%");


   if ( fileIO ) {
      WALBERLA_ROOT_SECTION() {
         std::string fileName = baseFolder + "/logSettlingBehavior.txt";
         std::ofstream file;
         file.open(fileName.c_str(), std::ofstream::app);
         file.precision(8);

         file << actualSolidVolumeFraction << " \t" << numSpheres << " \t"
              << diameter << " \t" << tau << " \t" << extForceZ << " \t" << relativeVelDiffLimit << " \t"
              << interpol << " \t" << dist << " \t" << addedMassCorr << " \t" << liftCorr << " \t" << effVisc << " \t"
              << ((useTurbulenceModel) ? 1 : 0) << " \t"
              << interactionSubCycles << " \t" << peSubSteps << " \t" << lubricationCutOffDistance << " \t"
              << curSettlingVel << " \t" << expectedVelocity << " \t" << velUnhindered << " \t" << maxPen << " \n";
      }
   }

   return 0;
}

} //namespace hindered_settling_dynamics_dpm

int main( int argc, char **argv ){
   hindered_settling_dynamics_dpm::main(argc, argv);
}

