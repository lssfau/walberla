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
//! \file InspectionProbe.cpp
//! \ingroup lbm_mesapd_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "blockforest/Initialization.h"

#include "boundary/all.h"

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/all.h"
#include "core/timing/RemainingTimeLogger.h"
#include "core/logging/all.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"

#include "lbm/boundary/all.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D3Q19.h"

#include "lbm_mesapd_coupling/mapping/ParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/MovingParticleMapping.h"
#include "lbm_mesapd_coupling/momentum_exchange_method/boundary/CurvedLinear.h"
#include "lbm_mesapd_coupling/utility/ParticleSelector.h"
#include "lbm_mesapd_coupling/utility/InspectionProbe.h"
#include "lbm_mesapd_coupling/DataTypes.h"

#include "mesa_pd/data/ParticleAccessorWithShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/ShapeStorage.h"
#include "mesa_pd/data/DataTypes.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/kernel/DoubleCast.h"
#include "mesa_pd/kernel/ParticleSelector.h"

#include "timeloop/SweepTimeloop.h"

namespace inspection_probe
{

///////////
// USING //
///////////

using namespace walberla;
using walberla::uint_t;

using LatticeModel_T = lbm::D3Q19< lbm::collision_model::TRT>;

using Stencil_T = LatticeModel_T::Stencil;
using PdfField_T = lbm::PdfField<LatticeModel_T>;

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

const uint_t FieldGhostLayers = 1;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag( "fluid" );
const FlagUID MO_Flag( "moving obstacle" );

/////////////////////////////////////
// BOUNDARY HANDLING CUSTOMIZATION //
/////////////////////////////////////
template <typename ParticleAccessor_T>
class MyBoundaryHandling
{
public:

   using MO_T = lbm_mesapd_coupling::CurvedLinear< LatticeModel_T, FlagField_T, ParticleAccessor_T >;
   using Type = BoundaryHandling< FlagField_T, Stencil_T, MO_T >;

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
                                  MO_T( "MO", MO_Flag, pdfField, flagField, particleField, ac_, fluid, *storage, *block ) );

      handling->fillWithDomain( FieldGhostLayers );

      return handling;
   }

private:

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID particleFieldID_;

   shared_ptr<ParticleAccessor_T> ac_;
};


//////////
// MAIN //
//////////
int main( int argc, char **argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   bool printToScreen = false;
   bool printSurroundingState = false;
   std::string outputFileName = "";

   for( int i = 1; i < argc; ++i )
   {
      if( std::strcmp( argv[i], "--printToScreen"  ) == 0 ) { printToScreen  = true; continue; }
      if( std::strcmp( argv[i], "--printSurroundingState"  ) == 0 ) { printSurroundingState  = true; continue; }
      if( std::strcmp( argv[i], "--outputFileName" ) == 0 ) { outputFileName = argv[++i]; continue; }
      WALBERLA_ABORT("Unrecognized command line argument found: " << argv[i]);
   }

   const uint_t length = uint_t(32);
   const real_t radius = real_t(5);

   const uint_t XBlocks = uint_t( 1 );
   const uint_t YBlocks = uint_t( 1 );
   const uint_t ZBlocks = uint_t( 1 );
   const uint_t XCells = length / XBlocks;
   const uint_t YCells = length / YBlocks;
   const uint_t ZCells = length / ZBlocks;

   auto blocks = blockforest::createUniformBlockGrid( XBlocks, YBlocks, ZBlocks, XCells, YCells, ZCells, uint_t(1), true,
                                                      true, true, true );


   mesa_pd::domain::BlockForestDomain domain(blocks->getBlockForestPointer());

   auto ps = std::make_shared<mesa_pd::data::ParticleStorage>(1);
   auto ss = std::make_shared<mesa_pd::data::ShapeStorage>();
   using ParticleAccessor_T = mesa_pd::data::ParticleAccessorWithShape;
   auto accessor = make_shared<ParticleAccessor_T >(ps, ss);
   auto sphereShape = ss->create<mesa_pd::data::Sphere>( radius );

   Vector3<real_t> centerPosition (real_c(length) * real_c(0.5));
   {
      mesa_pd::data::Particle&& p = *ps->create();
      p.setPosition(centerPosition);
      p.setInteractionRadius(radius);
      p.setOwner(mpi::MPIManager::instance()->rank());
      p.setShapeID(sphereShape);

   }

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( real_t(1) ) );

   BlockDataID pdfFieldID = lbm::addPdfFieldToStorage< LatticeModel_T >( blocks, "pdf field (zyxf)", latticeModel,
                                                                         Vector3< real_t >( real_t(0) ), real_t(1),
                                                                         uint_t(1), field::zyxf );

   BlockDataID flagFieldID = field::addFlagFieldToStorage<FlagField_T>( blocks, "flag field" );

   BlockDataID particleFieldID = field::addToStorage<lbm_mesapd_coupling::ParticleField_T>( blocks, "particle field", accessor->getInvalidUid(), field::zyxf, FieldGhostLayers );

   using BoundaryHandling_T = MyBoundaryHandling<ParticleAccessor_T>::Type;
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(MyBoundaryHandling<ParticleAccessor_T>( flagFieldID, pdfFieldID, particleFieldID, accessor), "boundary handling" );

   // map sphere
   lbm_mesapd_coupling::MovingParticleMappingKernel<BoundaryHandling_T> movingParticleMappingKernel(blocks, boundaryHandlingID, particleFieldID);
   ps->forEachParticle(false, mesa_pd::kernel::SelectAll(), *accessor, movingParticleMappingKernel, *accessor, MO_Flag);

   Vector3<real_t> probeLocation = centerPosition - Vector3<real_t>(radius + real_t(1), real_t(0), real_t(0));
   lbm_mesapd_coupling::InspectionProbe<PdfField_T, BoundaryHandling_T, ParticleAccessor_T> probe(probeLocation, blocks, pdfFieldID, boundaryHandlingID, particleFieldID, accessor, printToScreen, printSurroundingState, outputFileName );

   real_t rhoAtLocation;
   Vector3<real_t> velocityAtLocation;
   probe(rhoAtLocation, velocityAtLocation);

   WALBERLA_ASSERT_FLOAT_EQUAL(rhoAtLocation, real_t(1));
   WALBERLA_ASSERT_FLOAT_EQUAL(velocityAtLocation[0], real_t(0));
   WALBERLA_ASSERT_FLOAT_EQUAL(velocityAtLocation[1], real_t(0));
   WALBERLA_ASSERT_FLOAT_EQUAL(velocityAtLocation[2], real_t(0));

   return EXIT_SUCCESS;
}

} // namespace inspection_probe

int main( int argc, char **argv ){
   inspection_probe::main(argc, argv);
}
