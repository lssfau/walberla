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
//! \file BiDisperseFluidizedBedDPM.cpp
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
#include "pe/utility/DestroyBody.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/VTKOutput.h"

#include <vector>
#include <iomanip>
#include <iostream>
#include <random>

namespace bidisperse_fluidized_bed_dpm
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

//#define OutletBC

// PDF field, flag field & body field
using TensorField_T = GhostLayerField<Matrix3<real_t>, 1>;
using Vec3Field_T = GhostLayerField<Vector3<real_t>, 1>;
using ScalarField_T = GhostLayerField<real_t, 1>;
using ForceModel_T = lbm::force_model::GuoField<Vec3Field_T>;

using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRTField<ScalarField_T>, false, ForceModel_T>;
using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

// boundary handling
using NoSlip_T = lbm::NoSlip<LatticeModel_T, flag_t>;
using Inflow_T = lbm::SimpleUBB<LatticeModel_T, flag_t>;

#ifdef OutletBC
typedef lbm::Outlet< LatticeModel_T, FlagField_T >                     Outflow_T;
#else
using Outflow_T = lbm::SimplePressure<LatticeModel_T, flag_t>;
#endif

using BoundaryHandling_T = BoundaryHandling<FlagField_T, Stencil_T, NoSlip_T, Inflow_T, Outflow_T>;

using BodyTypeTuple = std::tuple<pe::Plane, pe::Sphere> ;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID NoSlip_Flag  ( "no slip flag" );
const FlagUID Inflow_Flag  ( "inflow flag" );
const FlagUID Outflow_Flag ( "outflow flag" );

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

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID, real_t uInflow ) :
      flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), uInflow_( uInflow ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   real_t uInflow_;

}; // class MyBoundaryHandling


BoundaryHandling_T * MyBoundaryHandling::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   WALBERLA_ASSERT_NOT_NULLPTR( block );

   FlagField_T * flagField = block->getData< FlagField_T >( flagFieldID_ );
   PdfField_T *  pdfField  = block->getData< PdfField_T > ( pdfFieldID_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

#ifdef OutletBC
   BoundaryHandling_T * handling = new BoundaryHandling_T( "Boundary handling", flagField, fluid,
                                    NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                    Inflow_T( "Inflow", Inflow_Flag, pdfField, Vector3<real_t>(real_t(0),real_t(0),uInflow_) ),
                                    Outflow_T( "Outflow", Outflow_Flag, pdfField, flagField, fluid ) );
#else
   BoundaryHandling_T * handling = new BoundaryHandling_T( "Boundary handling", flagField, fluid,
                                    NoSlip_T( "NoSlip", NoSlip_Flag, pdfField ),
                                    Inflow_T( "Inflow", Inflow_Flag, pdfField, Vector3<real_t>(real_t(0),real_t(0),uInflow_) ),
                                    Outflow_T( "Outflow", Outflow_Flag, pdfField, real_t(1) ) );
#endif

   const auto noslip  = flagField->getFlag( NoSlip_Flag );
   const auto inflow  = flagField->getFlag( Inflow_Flag );
   const auto outflow = flagField->getFlag( Outflow_Flag );

   CellInterval domainBB = storage->getDomainCellBB();

   domainBB.zMin() -= cell_idx_c( FieldGhostLayers ); // 0
   domainBB.zMax() += cell_idx_c( FieldGhostLayers ); // cellsZ+1

   // BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   storage->transformGlobalToBlockLocalCellInterval( bottom, *block );
   handling->forceBoundary( inflow, bottom );

   // TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   storage->transformGlobalToBlockLocalCellInterval( top, *block );
   handling->forceBoundary( outflow, top );

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
                              real_t diameter1, real_t diameter2, real_t diameterAvg,
                              real_t solidVolumeFraction, pe::MaterialID & material )
{
   real_t domainVolume = generationDomain.volume();
   real_t totalSphereVolume = domainVolume * solidVolumeFraction;

   real_t percentageOfSpecies1 = real_t(1);
   if( diameter1 < diameter2 || diameter1 > diameter2 )
   {
      real_t effectiveMass1 = diameter1 * diameter1 * diameter1;
      real_t effectiveMass2 = diameter2 * diameter2 * diameter2;
      real_t effectiveMassAvg = diameterAvg * diameterAvg * diameterAvg;

      percentageOfSpecies1 = (effectiveMassAvg - effectiveMass2) / (effectiveMass1 - effectiveMass2);
   }

   WALBERLA_LOG_INFO_ON_ROOT("Creating "<< percentageOfSpecies1 * real_t(100) << "% of sphere type 1 with diameter = " << diameter1 <<
                             " and " << (real_t(1) - percentageOfSpecies1) * real_t(100) << "% of sphere type 2 with diameter = " << diameter2);

   real_t xParticle = real_t(0);
   real_t yParticle = real_t(0);
   real_t zParticle = real_t(0);
   real_t creationDiameter = real_t(0);

   real_t currentSphereVolume = real_t(0);
   uint_t numberOfSpheres = uint_t(0);

   while( currentSphereVolume < totalSphereVolume )
   {

      WALBERLA_ROOT_SECTION()
      {
         xParticle = math::realRandom<real_t>(generationDomain.xMin(), generationDomain.xMax());
         yParticle = math::realRandom<real_t>(generationDomain.yMin(), generationDomain.yMax());
         zParticle = math::realRandom<real_t>(generationDomain.zMin(), generationDomain.zMax());

         real_t speciesDecider = math::realRandom<real_t>(real_t(0), real_t(1));
         // if decider is in [0,...,percentageOfSpecies1], then species 1 is created, else species 2
         if( percentageOfSpecies1 > speciesDecider )
         {
            creationDiameter = diameter1;
         }
         else
         {
            creationDiameter = diameter2;
         }
      }

      WALBERLA_MPI_SECTION()
      {
         mpi::broadcastObject( xParticle );
         mpi::broadcastObject( yParticle );
         mpi::broadcastObject( zParticle );
         mpi::broadcastObject( creationDiameter );
      }

      pe::createSphere( globalBodyStorage, forest.getBlockStorage(), bodyStorageID, 0, Vector3<real_t>( xParticle, yParticle, zParticle ), creationDiameter * real_t(0.5), material );

      currentSphereVolume += math::pi / real_t(6) * creationDiameter * creationDiameter * creationDiameter;

      ++numberOfSpheres;
   }

   return numberOfSpheres;
}

class ForceOnBodiesAdder
{
public:

   ForceOnBodiesAdder( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & bodyStorageID,
                       const Vector3<real_t> & forcePerVolume )
         : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), forcePerVolume_( forcePerVolume )
   { }

   // set a force on all (only local, to avoid force duplication) bodies
   void operator()()
   {
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            real_t volume = bodyIt->getVolume();
            bodyIt->addForce ( forcePerVolume_ * volume );
         }
      }
   }

private:

   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   Vector3<real_t> forcePerVolume_;
};

class BodiesQuantityEvaluator
{
public:

   BodiesQuantityEvaluator( SweepTimeloop* timeloop, const shared_ptr<StructuredBlockStorage> & blockStorage,
                            const BlockDataID & bodyStorageID, const std::string & fileName, uint_t numSpheres, real_t totalMass ) :
         timeloop_( timeloop ), blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), fileName_( fileName ),
         numSpheres_( numSpheres ), totalMass_( totalMass )
   {
      std::ofstream file;
      file.open( fileName_.c_str() );
      file << "#t \tfX \tfY \tfZ \tvelX \tvelY \tvelZ \tCoMX \tCoMY \tCoMZ\n";
      file.close();
   }

   void operator()()
   {
      // get mean particle velocity
      Vector3<real_t> particleVelocity = getMeanParticleVelocity();
      // get force on particles
      Vector3<real_t> particleForce = getParticleForce();
      // get center of mass (average position weighted by mass proportion)
      Vector3<real_t> centerOfMass = getCenterOfMass();

      WALBERLA_ROOT_SECTION()
      {
         writeToFile( particleForce,  particleVelocity, centerOfMass );
      }
   }

private:

   Vector3<real_t> getMeanParticleVelocity()
   {
      Vector3<real_t> velocity(0.);
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            auto bodyVelocity = bodyIt->getLinearVel();
            velocity[0] += std::fabs(bodyVelocity[0]);
            velocity[1] += std::fabs(bodyVelocity[1]);
            velocity[2] += std::fabs(bodyVelocity[2]);
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( velocity[0], mpi::SUM );
         mpi::allReduceInplace( velocity[1], mpi::SUM );
         mpi::allReduceInplace( velocity[2], mpi::SUM );
      }
      return velocity / real_c(numSpheres_);
   }

   Vector3<real_t> getParticleForce()
   {
      Vector3<real_t> force(0.);
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
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

   Vector3<real_t> getCenterOfMass()
   {
      Vector3<real_t> centerOfMass(0.);
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin< pe::Sphere >( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end<pe::Sphere>(); ++bodyIt )
         {
            auto bodyMass = bodyIt->getMass();
            auto bodyPosition = bodyIt->getPosition();
            centerOfMass += ( bodyMass / totalMass_ ) * bodyPosition;
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( centerOfMass[0], mpi::SUM );
         mpi::allReduceInplace( centerOfMass[1], mpi::SUM );
         mpi::allReduceInplace( centerOfMass[2], mpi::SUM );
      }
      return centerOfMass;
   }

   void writeToFile( const Vector3<real_t> & particleForce, const Vector3<real_t> & meanParticleVel, const Vector3<real_t> & centerOfMass )
   {
      std::ofstream file;
      file.open( fileName_.c_str(), std::ofstream::app );
      file.precision(8);

      file << timeloop_->getCurrentTimeStep() << "\t "
           << particleForce[0] << "\t " << particleForce[1] << "\t " << particleForce[2] << "\t "
           << meanParticleVel[0] << "\t " << meanParticleVel[1] << "\t " << meanParticleVel[2] << "\t"
           << centerOfMass[0] << "\t " << centerOfMass[1] << "\t " << centerOfMass[2] << "\n";
      file.close();
   }

   SweepTimeloop* timeloop_;
   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   std::string fileName_;
   const uint_t numSpheres_;
   const real_t totalMass_;
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
   void resetMaximumPenetration()
   {
      maximumPenetration_ = real_t(0);
   }
private:
   pe::cr::ICR & collisionResponse_;
   real_t maximumPenetration_;
};

class SphereCounter
{
public:
   SphereCounter( const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & bodyStorageID, uint_t numSpheres)
         : blockStorage_( blockStorage ), bodyStorageID_( bodyStorageID ), numSpheres_( numSpheres )
   {}
   void operator()()
   {
      uint_t curNumSpheres = getNumSpheres();
      if( curNumSpheres != numSpheres_)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Warning: Number of spheres has changed: " << numSpheres_ << " -> " << curNumSpheres );
         numSpheres_ = curNumSpheres;
      }
   }

   void updateNumSpheres()
   {
      numSpheres_ = getNumSpheres();
   }

   uint_t getNumSpheres()
   {
      uint_t curNumSpheres = uint_t(0);
      for( auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt )
      {
         for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID_); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
         {
            ++curNumSpheres;
         }
      }
      WALBERLA_MPI_SECTION()
      {
         mpi::allReduceInplace( curNumSpheres, mpi::SUM );
      }
      return curNumSpheres;
   }

private:
   shared_ptr<StructuredBlockStorage> blockStorage_;
   const BlockDataID bodyStorageID_;
   uint_t numSpheres_;
};

void setBodyVelocities(const shared_ptr<StructuredBlockStorage> & blockStorage, const BlockDataID & bodyStorageID, const Vector3<real_t>& velocity )
{
   for( auto blockIt = blockStorage->begin(); blockIt != blockStorage->end(); ++blockIt )
   {
      for( auto bodyIt = pe::LocalBodyIterator::begin( *blockIt, bodyStorageID); bodyIt != pe::LocalBodyIterator::end(); ++bodyIt )
      {
         bodyIt->setLinearVel(velocity);
      }
   }

}

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
/*!\brief Simulation of a bidisperse fluidized bed with the discrete particle method.
 *
 * The simulation features a fluidized bed with spherical particles inside a rectangular column.
 * The domain size is [32 x 16 x 256] * d_avg
 * Since a bidisperse setup is used, i.e. the spheres have either diameter d1 or d2, one can specify d1, d2, and d_avg.
 * The volume-average diameter d_avg determines the fractional distribution of spheres with d1 and d2.
 * A monodisperse setting is obtained if d1 = d2 = d_avg.
 *
 * The fluidization is driven by a constant and uniform inflow from the bottom plane in positive z-direction.
 * A plane, not visible to the fluid, is placed above this inlet to prevent particles from directly interfering
 * with the inlet.
 * The Reynolds number and the Galileo number characterize the setup.
 *
 * Command line customization allows to switch on/off different parts of the DPS (discrete particle simulation)
 * algorithm.
 *
 * In a first step, the spherical particles are created with either d1 or d2 at random positions in the lower part of
 * the domain. Possible overlaps are resolved, artificially bounded by a horizontal plane at the top of the generation
 * domain, by carrying out PE-only steps.
 *
 * Afterwards, an initial, fluid-only, simulation is possible to equilibrate the fluid- and solid phase.
 * Since this is not needed for the default case, it is not used by default (numberOfInitialTimeSteps=0).
 *
 * Finally, the fully coupled simulation is carried out with the option to write vtk files.
 *
 * The algorithm, as well as the setup and the outcome are described in detail in
 * Rettinger, Ruede - "A Coupled Lattice Boltzmann Method and Discrete Element Method for Discrete Particle Simulations
 * of Particulate Flows" to be published
 *
 */
//*******************************************************************************************************************
int main( int argc, char **argv ) {
   debug::enterTestMode();

   mpi::Environment env(argc, argv);


   /////////////////////////////////////////////
   //                                         //
   //   Command Line Argument Customization   //
   //                                         //
   /////////////////////////////////////////////

   bool funcTest = false;
   bool vtkIO = false;
   uint_t vtkWriteFrequency = uint_t(0);
   uint_t numberOfTimeSteps = uint_t(50000);
   uint_t numberOfInitialTimeSteps = uint_t(0);
   std::string vtkBaseFolder = "vtk_out_BiDisperseFluidizedBedDPM";

   real_t densityRatio = real_t(2500) / real_t(1000);
   real_t ReynoldsNumber = real_t(10); // = rho_f * diameterAvg * uInflow / visc_f
   real_t GalileoNumber = real_t(28); // = sqrt((densityRatio - 1) * g * diameterAvg ) * diameterAvg / visc_f

   real_t diameter1 = real_t(0.7); // diameter of species 1
   real_t diameter2 = real_t(0.8); // diameter of species 2
   real_t diameterAvg = real_t(0.75); // (mass) average diameter of all particles

   uint_t interactionSubCycles = uint_t(2); // number of subcycles that involve evaluation of the interaction force
   uint_t peSubSteps = uint_t(10); // number of pe only calls in each subcycle
   real_t solidVolumeFraction = real_t(0.2); // in the generation domain (fraction of the whole domain)

   DPMethod dpm = DPMethod::GNS;
   Interpolation interpol = Interpolation::IKernel;
   Distribution dist = Distribution::DKernel;
   DragCorrelation dragCorr = DragCorrelation::Tenneti;
   LiftCorrelation liftCorr = LiftCorrelation::Saffman;
   AddedMassCorrelation addedMassCorr = AddedMassCorrelation::Finn;
   EffectiveViscosity effVisc = EffectiveViscosity::Eilers;

   bool useTurbulenceModel = true;
   const real_t smagorinskyConstant = real_t(0.1); //for turbulence model

   real_t lubricationCutOffDistance = diameterAvg; //0 switches it off

   for (int i = 1; i < argc; ++i) {
      if (std::strcmp(argv[i], "--funcTest") == 0) funcTest = true;
      else if (std::strcmp(argv[i], "--vtkIOFreq") == 0) vtkWriteFrequency = uint_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--numberOfTimeSteps") == 0) numberOfTimeSteps = uint_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--numberOfInitialTimeSteps") == 0) numberOfInitialTimeSteps = uint_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--densityRatio") == 0) densityRatio = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--Ga") == 0) GalileoNumber = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--Re") == 0) ReynoldsNumber = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--diameter1") == 0) diameter1 = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--diameter2") == 0) diameter2 = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--diameterAvg") == 0) diameterAvg = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--solidVolumeFraction") == 0) solidVolumeFraction = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--interactionSubCycles") == 0) interactionSubCycles = uint_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--peSubSteps") == 0) peSubSteps = uint_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--DPMethod") == 0) dpm = to_DPMethod(argv[++i]);
      else if (std::strcmp(argv[i], "--interpolation") == 0) interpol = to_interpolation(argv[++i]);
      else if (std::strcmp(argv[i], "--distribution") == 0) dist = to_distribution(argv[++i]);
      else if (std::strcmp(argv[i], "--dragCorrelation") == 0) dragCorr = to_dragCorrelation(argv[++i]);
      else if (std::strcmp(argv[i], "--liftCorrelation") == 0) liftCorr = to_liftCorrelation(argv[++i]);
      else if (std::strcmp(argv[i], "--addedMassCorrelation") == 0) addedMassCorr = to_addedMassCorrelation(argv[++i]);
      else if (std::strcmp(argv[i], "--effectiveViscosity") == 0) effVisc = to_effvisc(argv[++i]);
      else if (std::strcmp(argv[i], "--disableTurbulenceModel") == 0) useTurbulenceModel = false;
      else if (std::strcmp(argv[i], "--lubricationCutOff") == 0) lubricationCutOffDistance = real_c(std::atof(argv[++i]));
      else if (std::strcmp(argv[i], "--vtkBaseFolder") == 0) vtkBaseFolder = argv[++i];
      else WALBERLA_ABORT("Found invalid command line argument \"" << argv[i] << "\" - aborting...");
   }

   if (vtkWriteFrequency > 0) vtkIO = true;

   WALBERLA_CHECK(diameter1 <= real_t(1), "Diameter is not allowed to be > 1!");
   WALBERLA_CHECK(diameter2 <= real_t(1), "Diameter is not allowed to be > 1!");
   WALBERLA_CHECK((diameter1 <= diameterAvg && diameterAvg <= diameter2) ||
                  (diameter2 <= diameterAvg && diameterAvg <= diameter1),
                  "Average diameter has to be between diameter 1 and 2!");
   WALBERLA_CHECK(solidVolumeFraction > real_t(0), "Solid volume fraction has to be > 0!");
   WALBERLA_CHECK(solidVolumeFraction <= real_t(0.65), "Solid volume fraction is not allowed to be > 0.65!");
   WALBERLA_CHECK(interactionSubCycles > uint_t(0), "Number of interaction sub cycles has to be at least 1!");
   WALBERLA_CHECK(peSubSteps > uint_t(0), "Number of pe sub steps has to be at least 1!");
   WALBERLA_CHECK(lubricationCutOffDistance >= real_t(0), "Lubrication cut off distance has to be non-negative!");

   if (funcTest) {
      walberla::logging::Logging::instance()->setLogLevel(logging::Logging::LogLevel::WARNING);
   }
   if( !funcTest )
   {
      // create base directory if it does not yet exist
      filesystem::path tpath( vtkBaseFolder );
      if( !filesystem::exists( tpath ) )
         filesystem::create_directory( tpath );
   }


   ///////////////////////////////
   //                           //
   //   SIMULATION PROPERTIES   //
   //                           //
   ///////////////////////////////

   const uint_t xlength = uint_c(real_t(32) * diameterAvg);
   const uint_t ylength = uint_c(real_t(16) * diameterAvg);
   const uint_t zlength = uint_c(real_t(256) * diameterAvg);

   const real_t uInflow = real_t(0.01);
   const real_t viscosity = diameterAvg * uInflow / ReynoldsNumber;

   const real_t ug = GalileoNumber * viscosity / diameterAvg;
   const real_t gravitationalAcc = ug * ug / ((densityRatio - real_t(1)) * diameterAvg);

   const real_t omega = lbm::collision_model::omegaFromViscosity(viscosity);
   const real_t tau = real_t(1) / omega;

   const real_t dx = real_t(1);

   const real_t dt_DEM = real_t(1) / real_c(interactionSubCycles * peSubSteps);

   const real_t dt = real_t(1);
   const real_t dtInteractionSubCycle = dt / real_c(interactionSubCycles);
   const real_t dtBodyVelocityTimeDerivativeEvaluation = dtInteractionSubCycle;


   const uint_t timesteps = (funcTest) ? uint_t(3)
                                       : numberOfTimeSteps; // total number of time steps for the whole simulation
   const uint_t initialTimesteps = (funcTest) ? uint_t(0)
                                              : numberOfInitialTimeSteps; // total number of time steps for the initial simulation

   const Vector3<real_t> initialFluidVelocity(real_t(0), real_t(0), uInflow);


   if (!funcTest) {
      WALBERLA_LOG_INFO_ON_ROOT("Lx x Ly x Lz = " << xlength << " x " << ylength << " x " << zlength);
      WALBERLA_LOG_INFO_ON_ROOT("diameter1 = " << diameter1 << ", diameter2 = " << diameter2 << ", diameterAvg = " << diameterAvg);
      WALBERLA_LOG_INFO_ON_ROOT("Re = " << ReynoldsNumber << ", Ga = " << GalileoNumber);
      WALBERLA_LOG_INFO_ON_ROOT("tau = " << tau << ", omega = " << omega);
      WALBERLA_LOG_INFO_ON_ROOT("viscosity = " << viscosity);
      WALBERLA_LOG_INFO_ON_ROOT("gravity = " << gravitationalAcc);
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
      WALBERLA_LOG_INFO_ON_ROOT("dt_DEM = " << dt_DEM);
   }

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = uint_t(4);
   const uint_t YBlocks = uint_t(1);
   const uint_t ZBlocks = uint_t(8);

   const uint_t XCells = xlength / XBlocks;
   const uint_t YCells = ylength / YBlocks;
   const uint_t ZCells = zlength / ZBlocks;

   if (XBlocks * XCells != xlength ||
       YBlocks * YCells != ylength ||
       ZBlocks * ZCells != zlength) WALBERLA_ABORT("Domain decomposition failed!");

   auto blocks = blockforest::createUniformBlockGrid(XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells,
                                                     dx, 0, false, false,
                                                     false, false, false,
                                                     false);


   //write domain decomposition to file
   if( vtkIO )
   {
      vtk::writeDomainDecomposition( blocks, "domain_decomposition", vtkBaseFolder );
   }

   ////////
   // PE //
   ////////

   shared_ptr<pe::BodyStorage> globalBodyStorage = make_shared<pe::BodyStorage>();
   pe::SetBodyTypeIDs<BodyTypeTuple>::execute();
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling<BodyTypeTuple>(), "pe Body Storage");
   auto ccdID = blocks->addBlockData(pe::ccd::createHashGridsDataHandling(globalBodyStorage, bodyStorageID), "CCD");
   auto fcdID = blocks->addBlockData(
         pe::fcd::createGenericFCDDataHandling<BodyTypeTuple, pe::fcd::AnalyticCollideFunctor>(), "FCD");

   pe::cr::ICR *cr;
   pe::cr::DEM cr_dem(globalBodyStorage, blocks->getBlockStoragePointer(), bodyStorageID, ccdID, fcdID);
   cr = &cr_dem;

   const real_t restitutionCoeff = real_t(0.88);
   const real_t frictionCoeff = real_t(0.25);

   real_t sphereVolume = diameterAvg * diameterAvg * diameterAvg * math::pi / real_t(6); // based on avg. diameter
   const real_t particleMass = densityRatio * sphereVolume;
   const real_t Mij = particleMass * particleMass / (real_t(2) * particleMass);
   const real_t lnDryResCoeff = std::log(restitutionCoeff);
   const real_t collisionTime = real_t(0.5);
   const real_t stiffnessCoeff = math::pi * math::pi * Mij / (collisionTime * collisionTime *
                                 (real_t(1) - lnDryResCoeff * lnDryResCoeff / (math::pi * math::pi + lnDryResCoeff * lnDryResCoeff)));
   const real_t dampingCoeff = -real_t(2) * std::sqrt(Mij * stiffnessCoeff) *
                               (std::log(restitutionCoeff) / std::sqrt(math::pi * math::pi + (std::log(restitutionCoeff) * std::log(restitutionCoeff))));
   const real_t contactDuration = real_t(2) * math::pi * Mij / (std::sqrt(real_t(4) * Mij * stiffnessCoeff - dampingCoeff * dampingCoeff)); //formula from Uhlman

   WALBERLA_LOG_INFO_ON_ROOT("Created particle material with:\n"
                                   << " - coefficient of restitution = " << restitutionCoeff << "\n"
                                   << " - coefficient of friction = " << frictionCoeff << "\n"
                                   << " - stiffness coefficient kn = " << stiffnessCoeff << "\n"
                                   << " - damping coefficient cdn = " << dampingCoeff << "\n"
                                   << " - contact time Tc = " << contactDuration);

   auto peMaterial = pe::createMaterial("particleMat", densityRatio, restitutionCoeff, frictionCoeff, frictionCoeff,
                                        real_t(0), real_t(200), stiffnessCoeff, dampingCoeff, dampingCoeff);

   // create bounding planes
   pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(1, 0, 0), Vector3<real_t>(0, 0, 0), peMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(-1, 0, 0), Vector3<real_t>(real_c(xlength), 0, 0), peMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 1, 0), Vector3<real_t>(0, 0, 0), peMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, -1, 0), Vector3<real_t>(0, real_c(ylength), 0), peMaterial);

   // set the planes in z-direction slightly inwards the simulation domain to avoid direct interaction of the spheres with the in- and outflow BC
   real_t zOffset = real_t(4);
   pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 0, 1), Vector3<real_t>(0, 0, zOffset), peMaterial);
   pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 0, -1), Vector3<real_t>(0, 0, real_c(zlength) - zOffset), peMaterial);


   /////////////////
   // PE COUPLING //
   /////////////////

   // connect to pe
   const real_t overlap = real_t(1.5) * dx;
   auto syncCall = std::bind(pe::syncNextNeighbors<BodyTypeTuple>, std::ref(blocks->getBlockForest()),
                               bodyStorageID, static_cast<WcTimingTree *>(nullptr), overlap, false);
   shared_ptr<CollisionPropertiesEvaluator> collisionPropertiesEvaluator = walberla::make_shared<CollisionPropertiesEvaluator>(*cr);

   // create the spheres
   uint_t numSpheres = 0;
   {
      real_t radiusMax = real_t(0.5) * std::max(diameter1, diameter2);
      AABB generationDomain(radiusMax, radiusMax, radiusMax + zOffset,
                            real_c(xlength) - radiusMax, real_c(ylength) - radiusMax, real_t(64) * diameterAvg + zOffset);
      numSpheres = createSpheresRandomly(*blocks, *globalBodyStorage, bodyStorageID, generationDomain,
                                         diameter1, diameter2, diameterAvg, solidVolumeFraction, peMaterial);
      syncCall();

      const uint_t initialPeSteps = uint_t(50000);
      const real_t dt_DEM_init = collisionTime / real_t(uint_t(10) * peSubSteps);
      const real_t overlapLimit = real_t(0.05) * diameterAvg;

      WALBERLA_LOG_INFO_ON_ROOT("Sphere creation done -- " << numSpheres << " spheres created");
      WALBERLA_LOG_INFO_ON_ROOT(
            "resolving overlaps with goal all < " << overlapLimit / diameterAvg * real_t(100) << "%");

      auto boundingPlaneForInit = pe::createPlane(*globalBodyStorage, 0, Vector3<real_t>(0, 0, -1),
                                                  Vector3<real_t>(0, 0, generationDomain.zMax() + radiusMax),
                                                  peMaterial);

      for (uint_t pet = uint_t(1); pet <= initialPeSteps; ++pet) {
         cr->timestep(dt_DEM_init);
         syncCall();
         (*collisionPropertiesEvaluator)();
         real_t maxPen = collisionPropertiesEvaluator->getMaximumPenetrationInSimulation();
         if (maxPen < overlapLimit) {
            WALBERLA_LOG_INFO_ON_ROOT("Carried out " << pet << " DEM-only time steps to resolve initial overlaps");
            break;
         } else {
            if (pet % uint_t(200) == uint_t(0)) {
               WALBERLA_LOG_INFO_ON_ROOT(
                     pet << " - current max overlap = " << maxPen / diameterAvg * real_t(100) << "%");
            }
         }
         collisionPropertiesEvaluator->resetMaximumPenetration();
      }
      pe::destroyBodyBySID(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID,
                           boundingPlaneForInit->getSystemID());

      // set body velocities to zero
      setBodyVelocities(blocks, bodyStorageID, Vector3<real_t>(real_t(0)));

      // some spheres might have gone lost due to unfortunate random generation
      SphereCounter sphereNumberChecker(blocks, bodyStorageID, numSpheres);
      sphereNumberChecker.updateNumSpheres();
      numSpheres = sphereNumberChecker.getNumSpheres();

   }



   // log setup to file
   if (!funcTest && vtkIO) {
      WALBERLA_ROOT_SECTION() {
         std::string fileName = vtkBaseFolder + "/setup.txt";
         std::ofstream file;
         file.open(fileName.c_str());
         file.precision(8);

         file << "Lx x Ly x Lz = " << xlength << " x " << ylength << " x " << zlength << "\n"
              << "diameter1 = " << diameter1 << ", diameter2 = " << diameter2 << ", diameterAvg = " << diameterAvg << "\n"
              << "Re = " << ReynoldsNumber << ", Ga = " << GalileoNumber << "\n"
              << "tau = " << tau << ", omega = " << omega << "\n"
              << "viscosity = " << viscosity << "\n"
              << "gravity = " << gravitationalAcc << "\n"
              << "input solid volume fraction = " << solidVolumeFraction << "\n"
              << "discrete particle method = " << dpm_to_string(dpm) << "\n"
              << "interpolator = " << interpol << "\n"
              << "distribution = " << dist << "\n"
              << "dragCorrelation = " << dragCorr << "\n"
              << "addedMassCorrelation = " << addedMassCorr << "\n"
              << "liftCorrelation = " << liftCorr << "\n"
              << "effective viscosity = " << effVisc << "\n"
              << "turbulence model = " << (useTurbulenceModel ? "yes" : "no") << "\n"
              << "interaction sub cycles = " << interactionSubCycles << "\n"
              << "pe sub steps = " << peSubSteps << "\n"
              << "lubrication cut off distance = " << lubricationCutOffDistance << "\n"
              << "dt_DEM = " << dt_DEM << "\n"
              << "Material:\n"
              << " - coefficient of restitution = " << restitutionCoeff << "\n"
              << " - coefficient of friction = " << frictionCoeff << "\n"
              << " - stiffness coefficient kn = " << stiffnessCoeff << "\n"
              << " - damping coefficient cdn = " << dampingCoeff << "\n"
              << " - contact time Tc = " << contactDuration << "\n"
              << "number of spheres = " << numSpheres << "\n";
      }
   }


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
   BlockDataID omegaFieldID = field::addToStorage< ScalarField_T >( blocks, "omega field", omega, field::zyxf, FieldGhostLayers );

   // create the lattice model
   LatticeModel_T latticeModel = LatticeModel_T( omegaFieldID, ForceModel_T( forceFieldID ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage( blocks, "pdf field (zyxf)", latticeModel, initialFluidVelocity, real_t(1), FieldGhostLayers, field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   // add boundary handling & initialize outer domain boundaries
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >( MyBoundaryHandling( flagFieldID, pdfFieldID, uInflow ), "boundary handling" );

   // field to store fluid velolcity
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
            using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor>;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         } else if (dist == Distribution::DKernel) {
            using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor>;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         }
      } else if (interpol == Interpolation::IKernel) {
         if (dist == Distribution::DNearestNeighbor) {
            using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor>;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         } else if (dist == Distribution::DKernel) {
            using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor>;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         }
      } else if (interpol == Interpolation::ITrilinear) {
         if (dist == Distribution::DNearestNeighbor) {
            using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::NearestNeighborDistributor>;
            shared_ptr<IFE_T> forceEvaluatorPtr = make_shared<IFE_T>(blocks, dragForceFieldID, bodyStorageID,
                                                                     flagFieldID, Fluid_Flag,
                                                                     velocityFieldID, svfFieldID,
                                                                     pressureGradientFieldID, dragCorrelationFunction,
                                                                     viscosity);
            dragAndPressureForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
         } else if (dist == Distribution::DKernel) {
            using IFE_T = pe_coupling::discrete_particle_methods::InteractionForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::KernelDistributor>;
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
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor>;
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
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor>;
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
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, liftForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        velocityFieldID, velocityCurlFieldID, liftCorrelationFunction,
                                                                        viscosity );
         liftForceEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         using IFE_T = pe_coupling::discrete_particle_methods::LiftForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::KernelDistributor>;
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
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::NearestNeighborFieldInterpolator, field::KernelDistributor>;
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
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::KernelFieldInterpolator, field::KernelDistributor>;
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
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::NearestNeighborDistributor>;
         shared_ptr< IFE_T > forceEvaluatorPtr = make_shared< IFE_T > ( blocks, amForceFieldID, bodyStorageID, flagFieldID, Fluid_Flag,
                                                                        timeDerivativeVelocityFieldID, addedMassCorrelationFunction,
                                                                        bodyVelocityTimeDerivativeEvaluator );
         addedMassEvaluationFunction = std::bind(&IFE_T::operator(), forceEvaluatorPtr);
      }
      else if( dist == Distribution::DKernel )
      {
         using IFE_T = pe_coupling::discrete_particle_methods::AddedMassForceEvaluator<FlagField_T, field::TrilinearFieldInterpolator, field::KernelDistributor>;
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
      using LE_T = pe_coupling::discrete_particle_methods::LubricationForceEvaluator;
      shared_ptr<LE_T> lubEval = make_shared<LE_T>( blocks, globalBodyStorage, bodyStorageID, viscosity, lubricationCutOffDistance );
      lubricationEvaluationFunction = std::bind(&LE_T::operator(), lubEval);
   }
   else
   {
      lubricationEvaluationFunction = emptyFunction;
   }


   //////////////////////////////
   //                          //
   //    INITIAL SIMULATION    //
   //                          //
   //////////////////////////////

   {
      // init solid volume fraction field
      pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::KernelDistributor> svfEvaluator( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  svfEvaluator( &(*blockIt) );

      (*svfCommunicationScheme)();

      pe_coupling::discrete_particle_methods::ForceFieldResetter forceFieldResetter( forceFieldID );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  forceFieldResetter( &(*blockIt) );
   }

   if( initialTimesteps > uint_t(0) )
   {

      SweepTimeloop timeloop( blocks->getBlockStorage(), initialTimesteps );

      ///////// begin of GNS timeloop ///////////////////////

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSPressureFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( pressureFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Pressure Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(pressureCommunicationScheme), "Pressure Field Communication" );

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::PressureGradientFieldEvaluator<LatticeModel_T, BoundaryHandling_T>( pressureGradientFieldID, pressureFieldID, boundaryHandlingID ), "Pressure Gradient Field Evaluation" )
                     << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(pressureGradientCommunicationScheme), "Pressure Gradient Field Communication" );


      // subcycling loop begin
      for( uint_t subcycle = 1; subcycle <= interactionSubCycles; ++subcycle )
      {
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" )
                        << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> ( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID ), "Velocity Field Evaluation" )
                        << AfterFunction ( SharedFunctor<blockforest::communication::UniformBufferedScheme<stencil::D3Q27> >(velocityCommunicationScheme), "Velocity Field Communication" );

         // evaluate Fdrag
         timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( dragForceFieldID ), "Drag Force Field Reset" )
                        << AfterFunction( dragAndPressureForceEvaluationFunction, "Fluid-Particle Interaction Force Evaluation" )
                        << AfterFunction( dragForceComm, "Drag Force Field Communication" );

         // ext forces on bodies
         timeloop.add() << Sweep( DummySweep(), "Dummy Sweep ")
                        << AfterFunction( pe_coupling::ForceTorqueOnBodiesResetter( blocks, bodyStorageID  ), "Force on Bodies Reset" );

      }

      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" )
                     << AfterFunction( SharedFunctor<pe_coupling::discrete_particle_methods::AveragedInteractionForceFieldToForceFieldAdder> (dragForceFieldToForceFieldAdder ), "Drag Force Field To Force Field Adder" );

      timeloop.add() << Sweep( makeSharedSweep<pe_coupling::discrete_particle_methods::EffectiveViscosityFieldEvaluator>( effectiveViscosityEvaluator ), "Effective Viscosity Evaluation");
      if( useTurbulenceModel ) timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::GNSSmagorinskyLESField<LatticeModel_T>(blocks, omegaFieldID, pdfFieldID, svfFieldID, smagorinskyConstant), "Turbulence Model" );

      // execute GNS-LBM sweep, boundary handling, PDF communication
      auto sweep = pe_coupling::discrete_particle_methods::makeGNSSweep< LatticeModel_T, FlagField_T >( pdfFieldID, svfFieldID, flagFieldID, Fluid_Flag );

      timeloop.add() << Sweep( lbm::makeCollideSweep( sweep ), "GNS-LBM sweep (collide)" );

      timeloop.add() << BeforeFunction( pdfScheme, "LBM Communication" )
                     << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

      timeloop.add() << Sweep( lbm::makeStreamSweep( sweep ), "GNS-LBM sweep (stream)" );


      timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

      WALBERLA_LOG_INFO_ON_ROOT("Starting initial (GNS-LBM only) simulation with " << initialTimesteps << " time steps.")
      timeloop.run();

   }


   /////////////////////////////
   //                         //
   //    ACTUAL SIMULATION    //
   //                         //
   /////////////////////////////

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );

   {
      pe_coupling::discrete_particle_methods::ForceFieldResetter forceFieldResetter( forceFieldID );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  forceFieldResetter( &(*blockIt) );

      // evaluate fluid velocity field
      pe_coupling::discrete_particle_methods::GNSVelocityFieldEvaluator<LatticeModel_T, BoundaryHandling_T> velocityEvaluator( velocityFieldID, pdfFieldID, svfFieldID, boundaryHandlingID );
      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )  velocityEvaluator( &(*blockIt) );

   }

   // configure vtk output
   if( vtkIO )
   {
      shared_ptr<vtk::VTKOutput> pdfFieldVTKWriter = createFluidFieldVTKWriter( blocks, pdfFieldID, flagFieldID, vtkWriteFrequency, vtkBaseFolder );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK (fluid field data)" );

      auto bodyVtkOutput = make_shared<pe::SphereVtkOutput>( bodyStorageID, blocks->getBlockStorage() );
      auto bodyVTK = vtk::createVTKOutput_PointData( bodyVtkOutput, "bodies", vtkWriteFrequency, vtkBaseFolder );
      timeloop.addFuncBeforeTimeStep( vtk::writeFiles( bodyVTK ), "VTK (sphere data)" );

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
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::ForceFieldResetter( forceFieldID ), "Force Field Reset" )
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
                     << AfterFunction( dragForceComm, "Drag Force Field Communication" );

      // ext forces on bodies and PE step(s)
      timeloop.add() << Sweep( DummySweep(), "Dummy Sweep ")
                     << AfterFunction( BodiesQuantityEvaluator(&timeloop, blocks, bodyStorageID, vtkBaseFolder+"/quantityLogging.txt", numSpheres, real_t(numSpheres) * particleMass ), "Body Quantity Evaluator")
                     << AfterFunction( ForceOnBodiesAdder( blocks, bodyStorageID, Vector3<real_t>(0,0,- gravitationalAcc * ( densityRatio - real_t(1) ) )  ), "Gravitational and Buoyancy Force Add" )
                     << AfterFunction( pe_coupling::TimeStep( blocks, bodyStorageID, *cr, syncCall, dtInteractionSubCycle, peSubSteps, lubricationEvaluationFunction ), "Pe Time Step" );

      // update solid volume fraction field
      timeloop.add() << Sweep( pe_coupling::discrete_particle_methods::SolidVolumeFractionFieldEvaluator<FlagField_T, field::KernelDistributor> ( blocks, svfFieldID, bodyStorageID, flagFieldID, Fluid_Flag ), "Solid Volume Fraction Field Evaluation" )
                     << AfterFunction( SharedFunctor< pe_coupling::discrete_particle_methods::CombinedReductionFieldCommunication<ScalarField_T> >(svfCommunicationScheme), "Solid Volume Fraction Field Communication" );

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

   timeloop.addFuncAfterTimeStep( SphereCounter(blocks, bodyStorageID, numSpheres ), "Sphere Counter" );

   timeloop.addFuncAfterTimeStep( RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );


   WALBERLA_LOG_INFO_ON_ROOT("Starting actual (fully coupled) simulation with " << timesteps << " time steps.")

   // execute simulation
   WcTimingPool timeloopTiming;

   for( uint_t t = 0; t < timesteps; ++t )
   {
      timeloop.singleStep( timeloopTiming );
   }

   timeloopTiming.logResultOnRoot();

   return 0;
}

} //namespace bidisperse_fluidized_bed_dpm

int main( int argc, char **argv ){
   bidisperse_fluidized_bed_dpm::main(argc, argv);
}

