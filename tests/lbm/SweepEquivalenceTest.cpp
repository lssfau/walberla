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
//! \file SweepEquivalenceTest.cpp
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimpleUBB.h"
#include "lbm/communication/PdfFieldPackInfo.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/D2Q9.h"
#include "lbm/lattice_model/D3Q19.h"
#include "lbm/lattice_model/D3Q27.h"
#include "lbm/sweeps/CellwiseSweep.h"
#include "lbm/sweeps/SplitPureSweep.h"
#include "lbm/sweeps/SplitSweep.h"
#include "lbm/sweeps/SweepWrappers.h"

#include "blockforest/Initialization.h"
#include "blockforest/communication/UniformBufferedScheme.h"

#include "boundary/BoundaryHandling.h"

#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/SharedSweep.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"

#include "timeloop/SweepTimeloop.h"

#include <type_traits>


//#define TEST_USES_VTK_OUTPUT
#ifdef TEST_USES_VTK_OUTPUT
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"
#include "vtk/VTKOutput.h"
#endif

//#define DEVEL_OUTPUT
#ifdef DEVEL_OUTPUT
#include "core/math/Sample.h"
#endif





namespace walberla {

using flag_t = walberla::uint64_t;
using FlagField_T = FlagField<flag_t>;

const FlagUID  Fluid_Flag( "fluid" );
const FlagUID    UBB_Flag( "velocity bounce back" );
const FlagUID NoSlip_Flag( "no slip" );

const uint_t FieldSize        = uint_t(10);
const uint_t FieldGhostLayers = uint_t(1);
const real_t GlobalOmega      = real_t(1.4);
const real_t GlobalLambdaE    = real_t(1.8);
const real_t GlobalLambdaD    = real_t(1.7);



template< typename LatticeModel_T >
class MyBoundaryHandling
{
public:

   typedef lbm::NoSlip< LatticeModel_T, flag_t >    NoSlip_T;
   typedef lbm::SimpleUBB< LatticeModel_T, flag_t > UBB_T;

   typedef BoundaryHandling< FlagField_T, typename LatticeModel_T::Stencil, NoSlip_T, UBB_T > BoundaryHandling_T;

   MyBoundaryHandling( const std::string & id, const BlockDataID & flagField, const BlockDataID & pdfField, const real_t velocity ) :
      id_( id ), flagField_( flagField ), pdfField_( pdfField ), velocity_( velocity ) {}

   BoundaryHandling_T * operator()( IBlock* const block, const StructuredBlockStorage* const storage ) const;

private:

   const std::string id_;

   const BlockDataID flagField_;
   const BlockDataID  pdfField_;

   const real_t velocity_;

};

template< typename LatticeModel_T >
typename MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T *
MyBoundaryHandling<LatticeModel_T>::operator()( IBlock * const block, const StructuredBlockStorage * const storage ) const
{
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   WALBERLA_ASSERT_NOT_NULLPTR( block );
   WALBERLA_ASSERT_NOT_NULLPTR( storage );

   FlagField_T * flagField = block->getData< FlagField_T >( flagField_ );
   PdfField_T *   pdfField = block->getData< PdfField_T > (  pdfField_ );

   const auto fluid = flagField->flagExists( Fluid_Flag ) ? flagField->getFlag( Fluid_Flag ) : flagField->registerFlag( Fluid_Flag );

   BoundaryHandling_T * handling = new BoundaryHandling_T( std::string("boundary handling ")+id_, flagField, fluid,
                                                           NoSlip_T( "no slip", NoSlip_Flag, pdfField ),
                                                           UBB_T( "velocity bounce back", UBB_Flag, pdfField, velocity_, real_c(0), real_c(0) ) );

   // Couette flow -> periodic in x- and y-direction!

   CellInterval domainBB = storage->getDomainCellBB();
   storage->transformGlobalToBlockLocalCellInterval( domainBB, *block );

   domainBB.xMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.xMax() += cell_idx_c( FieldGhostLayers );
   domainBB.yMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.yMax() += cell_idx_c( FieldGhostLayers );
   domainBB.zMin() -= cell_idx_c( FieldGhostLayers );
   domainBB.zMax() += cell_idx_c( FieldGhostLayers );

   // no slip BOTTOM
   CellInterval bottom( domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(), domainBB.zMin() );
   handling->forceBoundary( NoSlip_Flag, bottom );

   // velocity bounce back TOP
   CellInterval top( domainBB.xMin(), domainBB.yMin(), domainBB.zMax(), domainBB.xMax(), domainBB.yMax(), domainBB.zMax() );
   handling->forceBoundary( UBB_Flag, top );

   handling->fillWithDomain( domainBB );

   return handling;
}



template< typename LatticeModel_T  >
void addTest( shared_ptr< StructuredBlockForest > & blocks, SweepTimeloop & timeloop, std::vector< BlockDataID > & fieldIds,
              const LatticeModel_T & latticeModel, field::Layout layout, const BlockDataID & flagFieldId, const real_t velocity,
              const char * fieldName )
{
   fieldIds.push_back( lbm::addPdfFieldToStorage( blocks, std::string("pdf field ") + std::string(fieldName),
                                                  latticeModel, Vector3<real_t>( velocity, velocity / real_t(2), velocity / real_t(4) ), real_t(1),
                                                  FieldGhostLayers, layout ) );

   BlockDataID boundaryHandlingId = blocks->addStructuredBlockData< typename MyBoundaryHandling< LatticeModel_T >::BoundaryHandling_T >(
            MyBoundaryHandling< LatticeModel_T >( fieldName, flagFieldId, fieldIds.back(), velocity ), std::string("boundary handling ") + std::string(fieldName) );

   blockforest::communication::UniformBufferedScheme< typename LatticeModel_T::CommunicationStencil > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LatticeModel_T > >( fieldIds.back() ) );

   timeloop.add()
      << BeforeFunction( communication, std::string("communication ") + std::string(fieldName) )
      << Sweep( MyBoundaryHandling<LatticeModel_T>::BoundaryHandling_T::getBlockSweep( boundaryHandlingId ), std::string("boundary sweep ") + std::string(fieldName) );
}



template< typename LatticeModel_T, class Enable = void >
struct AddTest;

template< typename LatticeModel_T  >
struct AddTest< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                             lbm::collision_model::SRT_tag >::value >::type >
{
   static void add( shared_ptr< StructuredBlockForest > & blocks, SweepTimeloop & timeloop, std::vector< BlockDataID > & fieldIds,
                    field::Layout layout, const BlockDataID & flagFieldId, const real_t velocity, const char * fieldName )
   {
      LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( GlobalOmega ) );
      addTest( blocks, timeloop, fieldIds, latticeModel, layout, flagFieldId, velocity, fieldName );
   }
};

template< typename LatticeModel_T  >
struct AddTest< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                             lbm::collision_model::TRT_tag >::value >::type >
{
   static void add( shared_ptr< StructuredBlockForest > & blocks, SweepTimeloop & timeloop, std::vector< BlockDataID > & fieldIds,
                    field::Layout layout, const BlockDataID & flagFieldId, const real_t velocity, const char * fieldName )
   {
      LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT( GlobalOmega, GlobalOmega ) );
      addTest( blocks, timeloop, fieldIds, latticeModel, layout, flagFieldId, velocity, fieldName );
   }
};

template< typename LatticeModel_T  >
struct AddTest< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel,
                                                                             lbm::collision_model::D3Q19MRT >::value >::type >
{
   static void add( shared_ptr< StructuredBlockForest > & blocks, SweepTimeloop & timeloop, std::vector< BlockDataID > & fieldIds,
                    field::Layout layout, const BlockDataID & flagFieldId, const real_t velocity, const char * fieldName )
   {
      LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT( GlobalOmega, GlobalOmega, GlobalOmega,
                                                                                    GlobalOmega, GlobalOmega, GlobalOmega ) );
      addTest( blocks, timeloop, fieldIds, latticeModel, layout, flagFieldId, velocity, fieldName );
   }
};



template< typename LatticeModel_T, class Enable = void >
struct AddTRTTest;

template< typename LatticeModel_T  >
struct AddTRTTest< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel::tag,
                                                                                lbm::collision_model::TRT_tag >::value >::type >
{
   static void add( shared_ptr< StructuredBlockForest > & blocks, SweepTimeloop & timeloop, std::vector< BlockDataID > & fieldIds,
                    field::Layout layout, const BlockDataID & flagFieldId, const real_t velocity, const char * fieldName )
   {
      LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT( GlobalLambdaE, GlobalLambdaD ) );
      addTest( blocks, timeloop, fieldIds, latticeModel, layout, flagFieldId, velocity, fieldName );
   }
};

template< typename LatticeModel_T  >
struct AddTRTTest< LatticeModel_T, typename std::enable_if< std::is_same< typename LatticeModel_T::CollisionModel,
                                                                                lbm::collision_model::D3Q19MRT >::value >::type >
{
   static void add( shared_ptr< StructuredBlockForest > & blocks, SweepTimeloop & timeloop, std::vector< BlockDataID > & fieldIds,
                    field::Layout layout, const BlockDataID & flagFieldId, const real_t velocity, const char * fieldName )
   {
      LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT::constructTRT( GlobalLambdaE, GlobalLambdaD ) );
      addTest( blocks, timeloop, fieldIds, latticeModel, layout, flagFieldId, velocity, fieldName );
   }
};



template< typename LatticeModel1_T, typename LatticeModel2_T >
void check( const shared_ptr< StructuredBlockForest > & blocks, const BlockDataID & fieldId1, const BlockDataID & fieldId2 )
{
   using PdfField1_T = lbm::PdfField< LatticeModel1_T >;
   using PdfField2_T = lbm::PdfField< LatticeModel2_T >;

   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      PdfField1_T * referenceField = block->template getData< PdfField1_T >( fieldId1 );
      PdfField2_T *          field = block->template getData< PdfField2_T >( fieldId2 );

      const auto & id1 = blocks->getBlockDataIdentifier( fieldId1 );
      const auto & id2 = blocks->getBlockDataIdentifier( fieldId2 );

      const std::string msg = std::string("Check failed for fields with block data ID \"") + id1 + std::string("\" and \"") + id2 + std::string("\"");

      #ifdef DEVEL_OUTPUT
      math::Sample samples[4];
      #endif

      for( cell_idx_t z = cell_idx_t(0); z < cell_idx_c(FieldSize); ++z )
         for( cell_idx_t y = cell_idx_t(0); y < cell_idx_c(FieldSize); ++y )
            for( cell_idx_t x = cell_idx_t(0); x < cell_idx_c(FieldSize); ++x )
            {
               Vector3< real_t > velocityReference;
               Vector3< real_t > velocity;

               real_t rhoReference = referenceField->getDensityAndVelocity( velocityReference, x, y, z );
               real_t rho          =          field->getDensityAndVelocity( velocity         , x, y, z );

               #ifdef DEVEL_OUTPUT
               samples[0].insert( std::fabs( velocityReference[0] - velocity[0] ) );
               samples[1].insert( std::fabs( velocityReference[1] - velocity[1] ) );
               samples[2].insert( std::fabs( velocityReference[2] - velocity[2] ) );
               samples[3].insert( std::fabs( rhoReference - rho ) );
               #else
               #ifdef __INTEL_COMPILER // std::math::float_distance causes a segmentation fault with Intel compiler on our test machine...
               WALBERLA_CHECK( floatIsEqual( velocityReference[0], velocity[0] ), msg );
               WALBERLA_CHECK( floatIsEqual( velocityReference[1], velocity[1] ), msg );
               WALBERLA_CHECK( floatIsEqual( velocityReference[2], velocity[2] ), msg );
               WALBERLA_CHECK( floatIsEqual( rhoReference, rho ), msg );
               #else
               WALBERLA_CHECK_FLOAT_EQUAL( velocityReference[0], velocity[0], msg );
               WALBERLA_CHECK_FLOAT_EQUAL( velocityReference[1], velocity[1], msg );
               WALBERLA_CHECK_FLOAT_EQUAL( velocityReference[2], velocity[2], msg );
               WALBERLA_CHECK_FLOAT_EQUAL( rhoReference, rho, msg );
               #endif
               #endif
            }

      #ifdef DEVEL_OUTPUT
      #ifdef __INTEL_COMPILER // std::math::float_distance causes a segmentation fault with Intel compiler on our test machine...
      WALBERLA_LOG_DEVEL( "Velocity (x): " << samples[0].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK( floatIsEqual( samples[0].range(), 0.0 ) );
      WALBERLA_LOG_DEVEL( "Velocity (y): " << samples[1].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK( floatIsEqual( samples[1].range(), 0.0 ) );
      WALBERLA_LOG_DEVEL( "Velocity (z): " << samples[2].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK( floatIsEqual( samples[2].range(), 0.0 ) );
      WALBERLA_LOG_DEVEL( "Density: " << samples[3].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK( floatIsEqual( samples[3].range(), 0.0 ) );
      #else
      WALBERLA_LOG_DEVEL( "Velocity (x): " << samples[0].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK_FLOAT_EQUAL( samples[0].range(), 0.0 );
      WALBERLA_LOG_DEVEL( "Velocity (y): " << samples[1].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK_FLOAT_EQUAL( samples[1].range(), 0.0 );
      WALBERLA_LOG_DEVEL( "Velocity (z): " << samples[2].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK_FLOAT_EQUAL( samples[2].range(), 0.0 );
      WALBERLA_LOG_DEVEL( "Density: " << samples[3].format( "[%min, %max], %mean, %med" ) );
      WALBERLA_CHECK_FLOAT_EQUAL( samples[3].range(), 0.0 );
      #endif
      #endif
   }
}



int main( int argc, char ** argv )
{
   debug::enterTestMode();

   mpi::Environment env( argc, argv );

   if( MPIManager::instance()->numProcesses() != 1 )
      WALBERLA_ABORT( "The number of processes must be equal to 1!" );

   auto blocks = blockforest::createUniformBlockGrid( uint_t(1), uint_t(1), uint_t(1),
                                                      FieldSize, FieldSize, FieldSize,
                                                      real_c(1.0), true,
                                                      true, true, false ); // periodicty

   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   #ifdef TEST_USES_VTK_OUTPUT
   SweepTimeloop timeloop( blocks->getBlockStorage(), uint_t(101) );
   #else
   SweepTimeloop timeloop( blocks->getBlockStorage(), uint_t(10) );
   #endif

   std::vector< std::vector< BlockDataID > > fieldIds;

   const real_t velocity = real_t(0.05);

   #ifdef TEST_USES_VTK_OUTPUT
   auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, "pdf_field", uint_t(10) );
   #endif

   ///////////////////////////
   // D3Q19, incompressible //
   ///////////////////////////

   fieldIds.emplace_back( );

   // SRT

   typedef lbm::D3Q19< lbm::collision_model::SRT, false > D3Q19_SRT_INCOMP;

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 SRT incomp zyxf cell-wise)" );

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q19_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q19 SRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q19 SRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_SRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D3Q19 SRT incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_SRT_INCOMP, float > >( fieldIds.back().back(), "density (D3Q19 SRT incomp zyxf cell-wise)" ) );
   #endif

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 SRT incomp fzyx cell-wise)" );

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q19_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q19 SRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q19 SRT incomp fzyx cell-wise - separate stream+collide)" );
   }

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT incomp zyxf split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                              "LB stream & collide (D3Q19 SRT incomp zyxf split)" );

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT incomp fzyx split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                              "LB stream & collide (D3Q19 SRT incomp fzyx split)" );

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT incomp zyxf split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_SRT_INCOMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 SRT incomp zyxf split pure)" );

   AddTest< D3Q19_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT incomp fzyx split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_SRT_INCOMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 SRT incomp fzyx split pure)" );

   // TRT

   typedef lbm::D3Q19< lbm::collision_model::TRT, false > D3Q19_TRT_INCOMP;

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 TRT incomp zyxf cell-wise)" );

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q19 TRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q19 TRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_TRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D3Q19 TRT incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_TRT_INCOMP, float > >( fieldIds.back().back(), "density (D3Q19 TRT incomp zyxf cell-wise)" ) );
   #endif

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 TRT incomp fzyx cell-wise)" );

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q19 TRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q19 TRT incomp fzyx cell-wise - separate stream+collide)" );
   }

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT incomp zyxf split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                              "LB stream & collide (D3Q19 TRT incomp zyxf split)" );

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT incomp fzyx split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                              "LB stream & collide (D3Q19 TRT incomp fzyx split)" );

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT incomp zyxf split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_TRT_INCOMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 TRT incomp zyxf split pure)" );

   AddTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT incomp fzyx split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_TRT_INCOMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 TRT incomp fzyx split pure)" );

   // MRT

   typedef lbm::D3Q19< lbm::collision_model::D3Q19MRT, false > D3Q19_MRT_INCOMP;

   AddTest< D3Q19_MRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 MRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_MRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 MRT incomp zyxf cell-wise)" );

   AddTest< D3Q19_MRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 MRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q19_MRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q19 MRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q19 MRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_MRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D3Q19 MRT incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_MRT_INCOMP, float > >( fieldIds.back().back(), "density (D3Q19 MRT incomp zyxf cell-wise)" ) );
   #endif

   AddTest< D3Q19_MRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 MRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_MRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 MRT incomp fzyx cell-wise)" );

   AddTest< D3Q19_MRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 MRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q19_MRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q19 MRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q19 MRT incomp fzyx cell-wise - separate stream+collide)" );
   }
   
   /////////////////////////
   // D3Q19, compressible //
   /////////////////////////

   fieldIds.emplace_back( );

   // SRT

   typedef lbm::D3Q19< lbm::collision_model::SRT, true > D3Q19_SRT_COMP;

   AddTest< D3Q19_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT comp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_SRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q19 SRT comp zyxf cell-wise)" );

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_SRT_COMP, float > >( fieldIds.back().back(), "velocity (D3Q19 SRT comp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_SRT_COMP, float > >( fieldIds.back().back(), "density (D3Q19 SRT comp zyxf cell-wise)" ) );
   #endif

   AddTest< D3Q19_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT comp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_SRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q19 SRT comp fzyx cell-wise)" );

   AddTest< D3Q19_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT comp zyxf split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_SRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                            "LB stream & collide (D3Q19 SRT comp zyxf split)" );

   AddTest< D3Q19_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT comp fzyx split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_SRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                            "LB stream & collide (D3Q19 SRT comp fzyx split)" );

   AddTest< D3Q19_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 SRT comp zyxf split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_SRT_COMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 SRT comp zyxf split pure)" );

   AddTest< D3Q19_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 SRT comp fzyx split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_SRT_COMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 SRT comp fzyx split pure)" );

   // TRT

   typedef lbm::D3Q19< lbm::collision_model::TRT, true > D3Q19_TRT_COMP;

   AddTest< D3Q19_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT comp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q19 TRT comp zyxf cell-wise)" );

   AddTest< D3Q19_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT comp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q19 TRT comp fzyx cell-wise)" );

   AddTest< D3Q19_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT comp zyxf split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_TRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                            "LB stream & collide (D3Q19 TRT comp zyxf split)" );

   AddTest< D3Q19_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT comp fzyx split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_TRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                            "LB stream & collide (D3Q19 TRT comp fzyx split)" );

   AddTest< D3Q19_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 TRT comp zyxf split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_TRT_COMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 TRT comp zyxf split pure)" );

   AddTest< D3Q19_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 TRT comp fzyx split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_TRT_COMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 TRT comp fzyx split pure)" );

   ///////////////////////////
   // D3Q27, incompressible //
   ///////////////////////////

   fieldIds.emplace_back( );

   // SRT

   typedef lbm::D3Q27< lbm::collision_model::SRT, false > D3Q27_SRT_INCOMP;

   AddTest< D3Q27_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q27 SRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q27 SRT incomp zyxf cell-wise)" );

   AddTest< D3Q27_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q27 SRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q27_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q27 SRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q27 SRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q27_SRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D3Q27 SRT incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q27_SRT_INCOMP, float > >( fieldIds.back().back(), "density (D3Q27 SRT incomp zyxf cell-wise)" ) );
   #endif

   AddTest< D3Q27_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q27 SRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q27 SRT incomp fzyx cell-wise)" );

   AddTest< D3Q27_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q27 SRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q27_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q27 SRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q27 SRT incomp fzyx cell-wise - separate stream+collide)" );
   }

   // TRT

   typedef lbm::D3Q27< lbm::collision_model::TRT, false > D3Q27_TRT_INCOMP;

   AddTest< D3Q27_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q27 TRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q27 TRT incomp zyxf cell-wise)" );

   AddTest< D3Q27_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q27 TRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q27_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q27 TRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q27 TRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   AddTest< D3Q27_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q27 TRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q27 TRT incomp fzyx cell-wise)" );

   AddTest< D3Q27_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q27 TRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D3Q27_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D3Q27 TRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D3Q27 TRT incomp fzyx cell-wise - separate stream+collide)" );
   }

   /////////////////////////
   // D3Q27, compressible //
   /////////////////////////

   fieldIds.emplace_back( );

   // SRT

   typedef lbm::D3Q27< lbm::collision_model::SRT, true > D3Q27_SRT_COMP;

   AddTest< D3Q27_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q27 SRT comp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_SRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q27 SRT comp zyxf cell-wise)" );

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q27_SRT_COMP, float > >( fieldIds.back().back(), "velocity (D3Q27 SRT comp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q27_SRT_COMP, float > >( fieldIds.back().back(), "density (D3Q27 SRT comp zyxf cell-wise)" ) );
   #endif

   AddTest< D3Q27_SRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q27 SRT comp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_SRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q27 SRT comp fzyx cell-wise)" );

   // TRT

   typedef lbm::D3Q27< lbm::collision_model::TRT, true > D3Q27_TRT_COMP;

   AddTest< D3Q27_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q27 TRT comp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_TRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q27 TRT comp zyxf cell-wise)" );

   AddTest< D3Q27_TRT_COMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q27 TRT comp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q27_TRT_COMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                    "LB stream & collide (D3Q27 TRT comp fzyx cell-wise)" );

   ////////////////////////////
   // TRT <-> MRT COMPARISON //
   ////////////////////////////

   fieldIds.emplace_back( );

   // TRT

   AddTRTTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 (real) TRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 (real) TRT incomp zyxf cell-wise)" );

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_TRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D3Q19 (real) TRT incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_TRT_INCOMP, float > >( fieldIds.back().back(), "density (D3Q19 (real) TRT incomp zyxf cell-wise)" ) );
   #endif

   AddTRTTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 (real) TRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 (real) TRT incomp fzyx cell-wise)" );

   AddTRTTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 (real) TRT incomp zyxf split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                              "LB stream & collide (D3Q19 (real) TRT incomp zyxf split)" );

   AddTRTTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 (real) TRT incomp fzyx split)" );
   timeloop.add() << Sweep( lbm::SplitSweep< D3Q19_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ),
                                                                              "LB stream & collide (D3Q19 (real) TRT incomp fzyx split)" );

   AddTRTTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 (real) TRT incomp zyxf split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_TRT_INCOMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 (real) TRT incomp zyxf split pure)" );

   AddTRTTest< D3Q19_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 (real) TRT incomp fzyx split pure)" );
   timeloop.add() << Sweep( lbm::SplitPureSweep< D3Q19_TRT_INCOMP >( fieldIds.back().back() ), "LB stream & collide (D3Q19 (real) TRT incomp fzyx split pure)" );

   // MRT

   typedef lbm::D3Q19< lbm::collision_model::D3Q19MRT, false > D3Q19_MRT_INCOMP;

   AddTRTTest< D3Q19_MRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D3Q19 MRT (TRT) incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_MRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 MRT (TRT) incomp zyxf cell-wise)" );

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D3Q19_MRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D3Q19 MRT (TRT) incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D3Q19_MRT_INCOMP, float > >( fieldIds.back().back(), "density (D3Q19 MRT (TRT) incomp zyxf cell-wise)" ) );
   #endif

   AddTRTTest< D3Q19_MRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D3Q19 MRT (TRT) incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D3Q19_MRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                      "LB stream & collide (D3Q19 MRT (TRT) incomp fzyx cell-wise)" );

   //////////////////////////
   // D2Q9, incompressible //
   //////////////////////////

   fieldIds.emplace_back( );

   // SRT

   typedef lbm::D2Q9< lbm::collision_model::SRT, false > D2Q9_SRT_INCOMP;

   AddTest< D2Q9_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D2Q9 SRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D2Q9_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                     "LB stream & collide (D2Q9 SRT incomp zyxf cell-wise)" );

   AddTest< D2Q9_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D2Q9 SRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D2Q9_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D2Q9 SRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D2Q9 SRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   #ifdef TEST_USES_VTK_OUTPUT
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::VelocityVTKWriter< D2Q9_SRT_INCOMP, float > >( fieldIds.back().back(), "velocity (D2Q9 SRT incomp zyxf cell-wise)" ) );
   pdfFieldVTKWriter->addCellDataWriter( make_shared< lbm::DensityVTKWriter < D2Q9_SRT_INCOMP, float > >( fieldIds.back().back(), "density (D2Q9 SRT incomp zyxf cell-wise)" ) );
   #endif

   AddTest< D2Q9_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D2Q9 SRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D2Q9_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                     "LB stream & collide (D2Q9 SRT incomp fzyx cell-wise)" );

   AddTest< D2Q9_SRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D2Q9 SRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D2Q9_SRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D2Q9 SRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D2Q9 SRT incomp fzyx cell-wise - separate stream+collide)" );
   }

   // TRT

   typedef lbm::D2Q9< lbm::collision_model::TRT, false > D2Q9_TRT_INCOMP;

   AddTest< D2Q9_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D2Q9 TRT incomp zyxf cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D2Q9_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                     "LB stream & collide (D2Q9 TRT incomp zyxf cell-wise)" );

   AddTest< D2Q9_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::zyxf, flagFieldId, velocity, "(D2Q9 TRT incomp zyxf cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D2Q9_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D2Q9 TRT incomp zyxf cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D2Q9 TRT incomp zyxf cell-wise - separate stream+collide)" );
   }

   AddTest< D2Q9_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D2Q9 TRT incomp fzyx cell-wise)" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< D2Q9_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag ) ),
                                                                                                     "LB stream & collide (D2Q9 TRT incomp fzyx cell-wise)" );

   AddTest< D2Q9_TRT_INCOMP >::add( blocks, timeloop, fieldIds.back(), field::fzyx, flagFieldId, velocity, "(D2Q9 TRT incomp fzyx cell-wise - separate stream+collide)" );
   {
      auto sweep = lbm::makeCellwiseSweep< D2Q9_TRT_INCOMP, FlagField_T >( fieldIds.back().back(), flagFieldId, Fluid_Flag );
      timeloop.add() << Sweep( makeStreamSweep( sweep ), "LB stream (D2Q9 TRT incomp fzyx cell-wise - separate stream+collide)" );
      timeloop.add() << Sweep( makeCollideSweep( sweep ), "LB collide (D2Q9 TRT incomp fzyx cell-wise - separate stream+collide)" );
   }

   #ifdef TEST_USES_VTK_OUTPUT
   timeloop.addFuncAfterTimeStep( vtk::writeFiles( pdfFieldVTKWriter ), "VTK" );
   #endif

   WcTimingPool timeloopTiming;
   timeloop.run( timeloopTiming );
   timeloopTiming.logResultOnRoot();

   ///////////////////////////
   // D3Q19, incompressible //
   ///////////////////////////

   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][1] );
   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][2] );
   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][3] );
   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][4] );
   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][5] );
   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][6] );
   check< D3Q19_SRT_INCOMP, D3Q19_SRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][7] );

   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][8]  );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][9]  );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][10] );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][11] );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][12] );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][13] );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][14] );
   check< D3Q19_SRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][15] );

   check< D3Q19_SRT_INCOMP, D3Q19_MRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][16] );
   check< D3Q19_SRT_INCOMP, D3Q19_MRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][17] );
   check< D3Q19_SRT_INCOMP, D3Q19_MRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][18] );
   check< D3Q19_SRT_INCOMP, D3Q19_MRT_INCOMP >( blocks, fieldIds[0][0], fieldIds[0][19] );

   /////////////////////////
   // D3Q19, compressible //
   /////////////////////////

   check< D3Q19_SRT_COMP, D3Q19_SRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][1] );
   check< D3Q19_SRT_COMP, D3Q19_SRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][2] );
   check< D3Q19_SRT_COMP, D3Q19_SRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][3] );
   check< D3Q19_SRT_COMP, D3Q19_SRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][4] );
   check< D3Q19_SRT_COMP, D3Q19_SRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][5] );

   check< D3Q19_SRT_COMP, D3Q19_TRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][6]  );
   check< D3Q19_SRT_COMP, D3Q19_TRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][7]  );
   check< D3Q19_SRT_COMP, D3Q19_TRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][8]  );
   check< D3Q19_SRT_COMP, D3Q19_TRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][9]  );
   check< D3Q19_SRT_COMP, D3Q19_TRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][10] );
   check< D3Q19_SRT_COMP, D3Q19_TRT_COMP >( blocks, fieldIds[1][0], fieldIds[1][11] );

   ///////////////////////////
   // D3Q27, incompressible //
   ///////////////////////////

   check< D3Q27_SRT_INCOMP, D3Q27_SRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][1] );
   check< D3Q27_SRT_INCOMP, D3Q27_SRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][2] );
   check< D3Q27_SRT_INCOMP, D3Q27_SRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][3] );

   check< D3Q27_SRT_INCOMP, D3Q27_TRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][4] );
   check< D3Q27_SRT_INCOMP, D3Q27_TRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][5] );
   check< D3Q27_SRT_INCOMP, D3Q27_TRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][6] );
   check< D3Q27_SRT_INCOMP, D3Q27_TRT_INCOMP >( blocks, fieldIds[2][0], fieldIds[2][7] );

   /////////////////////////
   // D3Q27, compressible //
   /////////////////////////

   check< D3Q27_SRT_COMP, D3Q27_SRT_COMP >( blocks, fieldIds[3][0], fieldIds[3][1] );

   check< D3Q27_SRT_COMP, D3Q27_TRT_COMP >( blocks, fieldIds[3][0], fieldIds[3][2] );
   check< D3Q27_SRT_COMP, D3Q27_TRT_COMP >( blocks, fieldIds[3][0], fieldIds[3][3] );

   ////////////////////////////
   // TRT <-> MRT COMPARISON //
   ////////////////////////////

   check< D3Q19_TRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][1] );
   check< D3Q19_TRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][2] );
   check< D3Q19_TRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][3] );
   check< D3Q19_TRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][4] );
   check< D3Q19_TRT_INCOMP, D3Q19_TRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][5] );

   check< D3Q19_TRT_INCOMP, D3Q19_MRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][6] );
   check< D3Q19_TRT_INCOMP, D3Q19_MRT_INCOMP >( blocks, fieldIds[4][0], fieldIds[4][7] );

   //////////////////////////
   // D2Q9, incompressible //
   //////////////////////////

   check< D2Q9_SRT_INCOMP, D2Q9_SRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][1] );
   check< D2Q9_SRT_INCOMP, D2Q9_SRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][2] );
   check< D2Q9_SRT_INCOMP, D2Q9_SRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][3] );

   check< D2Q9_SRT_INCOMP, D2Q9_TRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][4] );
   check< D2Q9_SRT_INCOMP, D2Q9_TRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][5] );
   check< D2Q9_SRT_INCOMP, D2Q9_TRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][6] );
   check< D2Q9_SRT_INCOMP, D2Q9_TRT_INCOMP >( blocks, fieldIds[5][0], fieldIds[5][7] );

   return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
