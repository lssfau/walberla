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
//! \file ForceBetweenTwoStationaryObjects.cpp
//! \ingroup lbm_mesapd_coupling
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

#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/SimpleBB.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/utility/ResetHydrodynamicForceTorqueKernel.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/DataTypes.h"

#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/ParticleSelector.h"

#include "stencil/D3Q27.h"

#include "timeloop/SweepTimeloop.h"

#include <iomanip>
#include <iostream>

namespace force_between_two_stationary_objects
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

using LatticeModel_T = lbm::D3Q19< lbm::collision_model::TRT>;
using LatticeModelComp_T = lbm::D3Q19< lbm::collision_model::TRT, true>;


using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag   ( "fluid" );
const FlagUID BB_Flag   ( "obstacle BB" );
const FlagUID CLI_Flag  ( "obstacle CLI" );

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////

template <typename LM_T, typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using PdfField_T = lbm::PdfField<LM_T>;
   using SBB_T = lbm_mesapd_coupling::SimpleBB< LM_T, FlagField_T, ParticleAccessor_T >;
   using CLI_T = lbm_mesapd_coupling::CurvedLinear< LM_T, FlagField_T, ParticleAccessor_T >;
   using Type = BoundaryHandling< FlagField_T, typename LM_T::Stencil, SBB_T, CLI_T >;

   MyBoundaryHandling( const BlockDataID & flagFieldID, const BlockDataID & pdfFieldID,
                       const BlockDataID & particleFieldID, const shared_ptr<ParticleAccessor_T>& ac) :
         flagFieldID_( flagFieldID ), pdfFieldID_( pdfFieldID ), particleFieldID_( particleFieldID ), ac_( ac ) {}

   Type * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR( block );
      WALBERLA_ASSERT_NOT_NULLPTR( storage );

      auto * flagField     = block->getData< FlagField_T >( flagFieldID_ );
      auto *  pdfField     = block->getData< PdfField_T > ( pdfFieldID_ );
      auto * particleField = block->getData< lbm_mesapd_coupling::ParticleField_T > ( particleFieldID_ );

      const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

      Type * handling = new Type( "moving obstacle boundary handling", flagField, fluid,
                                  SBB_T("SBB_BB", BB_Flag,  pdfField, flagField, particleField, ac_, fluid, *storage, *block ),
                                  CLI_T("CLI_BB", CLI_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ) );

      handling->fillWithDomain( FieldGhostLayers );

      return handling;
   }

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;
};

void createPlaneSetup(const shared_ptr<mesa_pd::data::ParticleStorage> & ps, const shared_ptr<mesa_pd::data::ShapeStorage> & ss,
                      const math::AABB & simulationDomain, Vector3<real_t> velocity)
{
   // create bounding planes
   mesa_pd::data::Particle p0 = *ps->create(true);
   p0.setPosition(simulationDomain.minCorner());
   p0.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,1) ));
   p0.setOwner(mpi::MPIManager::instance()->rank());
   p0.setType(0);
   p0.setLinearVelocity(velocity);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p1 = *ps->create(true);
   p1.setPosition(simulationDomain.maxCorner());
   p1.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,0,-1) ));
   p1.setOwner(mpi::MPIManager::instance()->rank());
   p1.setType(0);
   p1.setLinearVelocity(velocity);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p1.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p2 = *ps->create(true);
   p2.setPosition(simulationDomain.minCorner());
   p2.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(1,0,0) ));
   p2.setOwner(mpi::MPIManager::instance()->rank());
   p2.setType(0);
   p2.setLinearVelocity(velocity);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p2.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p3 = *ps->create(true);
   p3.setPosition(simulationDomain.maxCorner());
   p3.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(-1,0,0) ));
   p3.setOwner(mpi::MPIManager::instance()->rank());
   p3.setType(0);
   p3.setLinearVelocity(velocity);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p3.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p4 = *ps->create(true);
   p4.setPosition(simulationDomain.minCorner());
   p4.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,1,0) ));
   p4.setOwner(mpi::MPIManager::instance()->rank());
   p4.setType(0);
   p4.setLinearVelocity(velocity);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p4.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

   mesa_pd::data::Particle p5 = *ps->create(true);
   p5.setPosition(simulationDomain.maxCorner());
   p5.setShapeID(ss->create<mesa_pd::data::HalfSpace>( Vector3<real_t>(0,-1,0) ));
   p5.setOwner(mpi::MPIManager::instance()->rank());
   p5.setType(0);
   p5.setLinearVelocity(velocity);
   mesa_pd::data::particle_flags::set(p5.getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p5.getFlagsRef(), mesa_pd::data::particle_flags::FIXED);

}

template< typename LM_T, typename ParticleAccessor_T>
void createSimulationSetup( shared_ptr< StructuredBlockForest > blocks, shared_ptr<ParticleAccessor_T> accessor,
                            bool useSBB, shared_ptr<mesa_pd::data::ParticleStorage> ps, Vector3<real_t> velocity,
                            SweepTimeloop & timeloop )
{
   real_t omega = real_t(1);

   // create the lattice model
   LM_T latticeModel = LM_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );

   // add PDF field
   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LM_T >( blocks, "pdf field (zyxf)", latticeModel, velocity, real_t(1), uint_t(1), field::zyxf );

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   // add particle field
   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::zyxf, FieldGhostLayers );

   // add boundary handling
   using BoundaryHandling_T = typename MyBoundaryHandling<LM_T,ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<LM_T,ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );

   // boundary handling sweep (does the hydro force calculations and the no-slip treatment)
   timeloop.add() << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingID ), "Boundary Handling" );

   // LBM
   auto sweep = lbm::makeCellwiseSweep< LM_T, FlagField_T >( pdfFieldID, flagFieldID, Fluid_Flag );
   timeloop.add() << Sweep( makeSharedSweep(sweep), "cell-wise LB sweep" );

   // map particles
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);
   if( useSBB )
   {
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, movingParticleMappingKernel, *accessor, BB_Flag);
   }else{
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, movingParticleMappingKernel, *accessor, CLI_Flag);
   }
}


//*******************************************************************************************************************
/*!\brief Testcase that checks the force between two stationary objects (sphere-sphere and sphere-wall) in a quiescent fluid.
 *
 * Since the fluid has zero velocity, as well as the objects, the force has to be zero!
 *
 * The objects have only a small distance.
 * Thus, not all cells around the sphere contain fluid and are thus skipped by the boundary handling and force computation.
 * As a result, some force components are missing that balance the total force such that it adds up to zero.
 * As such, one would need to reconstruct the PDF values inside the missing (from the perspective of one of the obstacles)
 * fluid cells to then obtain the balance again.
 * This is e.g. done in
 *  - Ernst, Dietzel, Sommerfeld - A lattice Boltzmann method for simulating transport and agglomeration of resolved particles
 *    doi: 10.1007/s00707-013-0923-1 (2013), Section 2.2.4
 *  - in the work of Simon Bogner
 *
 * This is avoided here via the following:
 * ( also see the docu in the beginning of lbm/PdfField.h)
 * When using an incompressible LBM model in waLBerla, all stored PDF values are simply the offset from the default weights.
 * I.e. in a quiescent fluid with density 1, all of them are zero.
 * Consequently, the force contributions computed by the momentum exchange method are also all zero in that case.
 * It thus does not matter if some of them are missing due to a nearby obstacle because only zero values would be missing.
 *
 * This changes when using a compressible model, as then the stored PDF values are the regular ones in waLBerla.
 * As such, one would see the missing force contributions in this test.
 * To compensate this, the PDF values are explicitly normalized like for the incompressible case in the MEM boundary conditions.
 * Then, the same arguments as for the incompressible flow hold.
 *
 * This can also be seen as an implicit reconstruction of the PDFs in the missing cells with feq(0, vec(0)).
 * As such, it can/will lead to errors if the velocity of the obstacle covering the missing cells is different from the
 * surrounding flow. In that case, an explicit reconstruction like in Ernst et al. should be applied, which however will
 * change the algorithm, as the computation of the force components can probbaly no longer be fused with the boundary handling.
 *
 * Tested variants:
 *  - boundary condition: SimpleBB, CLI
 *  - compressible and incompressible lattice model
 *  - different surface distances
 *  - sphere-sphere and sphere-wall
 *  - different system velocity, i.e. velocity of ALL things is constant but can be non-zero -> zero force expected
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

   bool useCompressible = false;
   bool useSBB = false;
   bool useSphereWallSetup = false;
   real_t surfaceDistance = real_t(0.1);
   real_t systemVelocity = real_t(0);
   uint_t timesteps = uint_t(10);
   real_t radius = real_t(5);

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--useCompressible"  ) == 0 ) { useCompressible  = true; continue;}
      if( std::strcmp( argv[i], "--useSBB"  ) == 0 ) { useSBB  = true; continue;}
      if( std::strcmp( argv[i], "--useSphereWallSetup"  ) == 0 ) { useSphereWallSetup  = true; continue;}
      if( std::strcmp( argv[i], "--surfaceDistance"    ) == 0 ) { surfaceDistance = real_c(std::atof( argv[++i])); continue;}
      if( std::strcmp( argv[i], "--systemVelocity"    ) == 0 ) { systemVelocity = real_c(std::atof( argv[++i])); continue;}
      if( std::strcmp( argv[i], "--radius"    ) == 0 ) { radius = real_c(std::atof( argv[++i])); continue;}
      if( std::strcmp( argv[i], "--timesteps"    ) == 0 ) { timesteps = uint_c(std::atof( argv[++i])); continue;}
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }


   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   const uint_t length = uint_t(real_t(4) * radius);
   const Vector3<real_t> velocity(systemVelocity, real_t(0), real_t(0));

   ///////////////////////////
   // BLOCK STRUCTURE SETUP //
   ///////////////////////////

   const uint_t XBlocks = uint_t( 1 );
   const uint_t YBlocks = uint_t( 1 );
   const uint_t ZBlocks = uint_t( 1 );
   const uint_t XCells = length / XBlocks;
   const uint_t YCells = length / YBlocks;
   const uint_t ZCells = length / ZBlocks;

   auto blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells, uint_t(1), true,
                                                      false, false, false );


   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());

   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = make_shared<ParticleAccessor_T >(ps, ss);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );

   createPlaneSetup(ps, ss, walberla::math::AABB(real_t(0), real_t(0), real_t(0), real_c(length), real_c(length), real_c(length)), velocity);

   walberla::id_t sphereID;
   if(useSphereWallSetup)
   {
      // create sphere close to right wall
      Vector3<real_t> position1 (real_c(length) - radius - surfaceDistance,
                                 real_c(length) * real_c(0.5),
                                 real_c(length) * real_c(0.5));
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(position1);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         p.setLinearVelocity(velocity);
         sphereID = p.getUid();
      }
   } else {
      // create two spheres
      Vector3<real_t> position1 (real_c(length) * real_c(0.5) - radius - surfaceDistance*real_t(0.5),
                                 real_c(length) * real_c(0.5),
                                 real_c(length) * real_c(0.5));
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(position1);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         p.setLinearVelocity(velocity);
         sphereID = p.getUid();
      }


      Vector3<real_t> position2 (real_c(length) * real_c(0.5) + radius + surfaceDistance*real_t(0.5),
                                 real_c(length) * real_c(0.5),
                                 real_c(length) * real_c(0.5));
      {
         mesa_pd::data::Particle&& p = *ps->create();
         p.setPosition(position2);
         p.setInteractionRadius(radius);
         p.setOwner(mpi::MPIManager::instance()->rank());
         p.setShapeID(sphereShape);
         p.setLinearVelocity(velocity);
      }
   }

   // create the timeloop
   SweepTimeloop timeloop( blocks->getBlockStorage(), timesteps );
   if(useCompressible) createSimulationSetup<LatticeModelComp_T>(blocks, accessor, useSBB, ps, velocity , timeloop);
   else                createSimulationSetup<LatticeModel_T>(blocks, accessor, useSBB, ps, velocity, timeloop);


   for(uint_t t = 0; t < timeloop.getNrOfTimeSteps(); ++t)
   {
      // some timesteps
      timeloop.singleStep();

      // check force
      size_t idx = accessor->uidToIdx(sphereID);
      auto hydrodynamicForce = accessor->getHydrodynamicForce(idx);

      //WALBERLA_LOG_INFO(hydrodynamicForce);

      WALBERLA_CHECK_FLOAT_EQUAL(hydrodynamicForce, Vector3<real_t>(real_t(0)), "Found non-zero force");


      lbm_mesapd_coupling::ResetHydrodynamicForceTorqueKernel resetHydrodynamicForceTorque;
      ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, resetHydrodynamicForceTorque, *accessor );

   }

   return 0;

}

} //namespace force_between_two_stationary_objects

int main( int argc, char **argv ){
   force_between_two_stationary_objects::main(argc, argv);
}
