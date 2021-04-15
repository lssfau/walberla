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
//! \file PermeabilityTest.cpp
//! \author Tobias Schruff <tobias.schruff@gmail.com>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "geometry/all.h"
#include "lbm/all.h"
#include "timeloop/all.h"


namespace walberla {

using flag_t = walberla::uint8_t;
using FlagField_T = FlagField<flag_t>;

enum Scenario { BCC };

// constants for analytical BCC permeability calculation, values taken from Sangani & Acrivos (1982)
static const double qs[31] = { 0.1000000E+01, 0.1575834E+01, 0.2483254E+01, 0.3233022E+01, 0.4022864E+01, 0.4650320E+01,
                               0.5281412E+01, 0.5826374E+01, 0.6258376E+01, 0.6544504E+01, 0.6878396E+01, 0.7190839E+01,
                               0.7268068E+01, 0.7304025E+01, 0.7301217E+01, 0.7236410E+01, 0.7298014E+01, 0.7369849E+01,
                               0.7109497E+01, 0.6228418E+01, 0.5235796E+01, 0.4476874E+01, 0.3541982E+01, 0.2939353E+01,
                               0.3935484E+01, 0.5179097E+01, 0.3959872E+01, 0.2227627E+01, 0.3393390E+01, 0.4491369E+01,
                               0.2200686E+01 };


struct Setup
{
   Setup( int argc, char ** argv )
   {
      epsilon = real_c( 0.01 );
      omega = real_t( 0.3 );
      rhoDiff = real_t( 0.00001 );
      kappa = real_t( 0.85 );
      length = uint_t( 32 );
      timesteps = uint_t( 10000 );
      checkFrequency = uint_t( 100 );
      scenario = BCC;
      collisionModel = std::string( "TRT" );

      for( int i = 1; i < argc; ++i )
      {
              if( std::strcmp( argv[i], "--omega" )          == 0 ) omega          = real_c( std::atof( argv[++i] ) );
         else if( std::strcmp( argv[i], "--epsilon" )        == 0 ) epsilon        = real_c( std::atof( argv[++i] ) );
         else if( std::strcmp( argv[i], "--dRho" )           == 0 ) rhoDiff        = real_c( std::atof( argv[++i] ) );
         else if( std::strcmp( argv[i], "--kappa" )          == 0 ) kappa          = real_c( std::atof( argv[++i] ) );
         else if( std::strcmp( argv[i], "--length" )         == 0 ) length         = uint_c( std::atoi( argv[++i] ) );
         else if( std::strcmp( argv[i], "--timesteps" )      == 0 ) timesteps      = uint_c( std::atoi( argv[++i] ) );
         else if( std::strcmp( argv[i], "--checkFrequency" ) == 0 ) checkFrequency = uint_c( std::atoi( argv[++i] ) );
      // else if( std::strcmp( argv[i], "--scenario" )       == 0 ) scenario       = std::string( argv[++i] );
         else if( std::strcmp( argv[i], "--collisionModel" ) == 0 ) collisionModel = std::string( argv[++i] );
         else WALBERLA_ABORT( "Unexpected argument " << argv[i] << "! Aborting ..." );
      }
   }

   // allowed deviation between simulation and analytical value
   real_t epsilon;
   // LBM omega value
   real_t omega;
   // density (pressure) difference applied at domain borders
   real_t rhoDiff;
   // control parameter for BCC scenario
   real_t kappa;
   // side length of domain
   uint_t length;
   // number of simulation time steps
   uint_t timesteps;
   // permeability check frequency
   uint_t checkFrequency;
   // simulation scenario (currently only BCC is implemented)
   Scenario scenario;
   // block distribution
   Vector3<uint_t> blocks;
   // LBM collision model (one of SRT, TRT, MRT, or CUM)
   std::string collisionModel;
};


// calculates (semi-)analytical permeability value for given setup
real_t permeability( Setup setup )
{
   // BCC implementation
   const real_t L = real_c(setup.length);
   const real_t r = real_c(std::sqrt(real_t(3.0))) / real_t(4) * L * setup.kappa;

   real_t drag( 0.0 );
   for( uint_t i = 0; i < 31; i++ )
      drag += real_c(qs[i]) * real_c(std::pow( setup.kappa, real_c(i) ));

   return ( L * L * L ) / ( real_t(6) * math::pi * r * real_t(2) * drag );
}


template< typename CollisionModel_T >
BlockDataID initPdfField( const shared_ptr<StructuredBlockForest> & blocks, real_t omega );

template< >
BlockDataID initPdfField< lbm::collision_model::SRT >( const shared_ptr<StructuredBlockForest> & blocks, real_t omega )
{
   using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::SRT( omega ) );
   return lbm::addPdfFieldToStorage( blocks, "PDF Field (SRT)", latticeModel, Vector3<real_t>(), real_t(1) );
}

template< >
BlockDataID initPdfField< lbm::collision_model::TRT >( const shared_ptr<StructuredBlockForest> & blocks, real_t omega )
{
   using LatticeModel_T = lbm::D3Q19<lbm::collision_model::TRT>;

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::TRT::constructWithMagicNumber( omega ) );
   return lbm::addPdfFieldToStorage( blocks, "PDF Field (TRT)", latticeModel, Vector3<real_t>(), real_t(1) );
}

template< >
BlockDataID initPdfField< lbm::collision_model::D3Q19MRT >( const shared_ptr<StructuredBlockForest> & blocks, real_t omega )
{
   using LatticeModel_T = lbm::D3Q19<lbm::collision_model::D3Q19MRT>;

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q19MRT::constructPanWithMagicNumber( omega ) );
   return lbm::addPdfFieldToStorage( blocks, "PDF Field (MRT)", latticeModel, Vector3<real_t>(), real_t(1) );
}

template< >
BlockDataID initPdfField< lbm::collision_model::D3Q27Cumulant >( const shared_ptr<StructuredBlockForest> & blocks, real_t omega )
{
   typedef lbm::D3Q27< lbm::collision_model::D3Q27Cumulant, true > LatticeModel_T;

   LatticeModel_T latticeModel = LatticeModel_T( lbm::collision_model::D3Q27Cumulant( omega ) );
   return lbm::addPdfFieldToStorage( blocks, "PDF Field (Cumulant)", latticeModel, Vector3<real_t>(), real_t(1) );
}


template< typename LatticeModel_T >
BlockDataID initBoundaryHandling( shared_ptr<StructuredBlockForest> & blocks, const BlockDataID & pdfFieldId, const BlockDataID & flagFieldId, const FlagUID & fluid, Setup setup )
{
   typedef lbm::DefaultBoundaryHandlingFactory< LatticeModel_T, FlagField_T > BHFactory_T;
   using BoundaryHandling_T = typename BHFactory_T::BoundaryHandling;

   BlockDataID boundaryHandlingId = BHFactory_T::addBoundaryHandlingToStorage( blocks, "boundary handling", flagFieldId, pdfFieldId, fluid,
                                                                               Vector3<real_t>(),
                                                                               Vector3<real_t>(),
                                                                               real_c(1),
                                                                               real_c(1) + setup.rhoDiff );
   std::vector< geometry::Sphere > spheres;

   switch( setup.scenario )
   {
   case BCC:

      const real_t L = real_c( setup.length );
      const real_t r = real_c(std::sqrt(real_c(3))) / real_c(4) * L * setup.kappa;

      // spheres in all eight corners of the domain
      spheres.emplace_back( Vector3<real_t>( 0, 0, 0 ), r );
      spheres.emplace_back( Vector3<real_t>( L, 0, 0 ), r );
      spheres.emplace_back( Vector3<real_t>( 0, L, 0 ), r );
      spheres.emplace_back( Vector3<real_t>( 0, 0, L ), r );
      spheres.emplace_back( Vector3<real_t>( L, L, 0 ), r );
      spheres.emplace_back( Vector3<real_t>( L, 0, L ), r );
      spheres.emplace_back( Vector3<real_t>( 0, L, L ), r );
      spheres.emplace_back( Vector3<real_t>( L, L, L ), r );
      // and one sphere in the middle
      spheres.emplace_back( Vector3<real_t>( L / real_c(2), L / real_c(2), L / real_c(2) ), r );

      break;
   }

   // set no-slip bcs inside the spheres
   geometry::initializer::BoundaryFromBody< BoundaryHandling_T > bodyInitializer( *blocks, boundaryHandlingId );
   for( size_t i = 0; i < spheres.size(); ++i )
      bodyInitializer.template init< geometry::Sphere >( spheres.at( i ), BHFactory_T::getNoSlip() );

   // set top and bottom pressure bcs
   geometry::initializer::BoundaryFromDomainBorder< BoundaryHandling_T > borderInitializer( *blocks, boundaryHandlingId );
   borderInitializer.init( BHFactory_T::getPressure0(), stencil::T, cell_idx_t(-1) );
   borderInitializer.init( BHFactory_T::getPressure1(), stencil::B, cell_idx_t(-1) );

   // initialize remaining cells with fluid
   geometry::setNonBoundaryCellsToDomain< BoundaryHandling_T >( *blocks, boundaryHandlingId );

   return boundaryHandlingId;
}


// determines the permeability for a given setup and compares to respective analytical solution
template< typename LM_T >
int setupAndExecute( Setup setup )
{
   using PdfField_T = lbm::PdfField< LM_T >;
   typedef lbm::DefaultBoundaryHandlingFactory< LM_T, FlagField_T > BHFactory_T;
   using BoundaryHandling_T = typename BHFactory_T::BoundaryHandling;

   Vector3<uint_t> blockLength;
   blockLength[0] = setup.length / setup.blocks[0];
   blockLength[1] = setup.length / setup.blocks[1];
   blockLength[2] = setup.length / setup.blocks[2];

   auto blocks = blockforest::createUniformBlockGrid( setup.blocks[0], setup.blocks[1], setup.blocks[2],
                                                      blockLength[0] , blockLength[1] , blockLength[2],
                                                      real_c(1),
                                                      true,
                                                      true, true, false );

   // create fields
   BlockDataID pdfFieldId  = initPdfField< typename LM_T::CollisionModel >( blocks, setup.omega );
   BlockDataID flagFieldId = field::addFlagFieldToStorage< FlagField_T >( blocks, "flag field" );

   const FlagUID fluid( "Fluid" );

   // create and initialize boundary handling
   BlockDataID boundaryHandlingId = initBoundaryHandling< LM_T >( blocks, pdfFieldId, flagFieldId, fluid, setup );

   // create communication for PdfField
   blockforest::communication::UniformBufferedScheme< typename LM_T::CommunicationStencil > communication( blocks );
   communication.addPackInfo( make_shared< lbm::PdfFieldPackInfo< LM_T > >( pdfFieldId ) );

   // create time loop
   SweepTimeloop timeloop( blocks->getBlockStorage(), setup.timesteps );

   // add LBM sweep and communication to time loop
   timeloop.add() << BeforeFunction( communication,                                                                           "LBM Communication" )
                  << Sweep( BoundaryHandling_T::getBlockSweep( boundaryHandlingId ),                                          "LBM Boundary Handling" );
   timeloop.add() << Sweep( makeSharedSweep( lbm::makeCellwiseSweep< LM_T, FlagField_T >( pdfFieldId, flagFieldId, fluid ) ), "LBM Stream & Collide" );

   const real_t tau( real_c(1) / setup.omega );
   const real_t nu( ( tau - real_c(0.5) ) / real_c(3) );

   typedef lbm::evaluations::Permeability< PdfField_T, BoundaryHandling_T > Permeability_T;
   shared_ptr< Permeability_T > permEval = make_shared< Permeability_T >( nu, pdfFieldId, boundaryHandlingId, fluid, blocks );
   permEval->init( blocks->getDomainCellBB(), uint_t(2), setup.checkFrequency );

   // add permeability evaluation to timeloop
   timeloop.addFuncAfterTimeStep( SharedFunctor< Permeability_T >( permEval ), "Permeability Evaluation" );
   // log remaining time
   timeloop.addFuncAfterTimeStep( timing::RemainingTimeLogger( timeloop.getNrOfTimeSteps() ), "Remaining Time Logger" );

   timeloop.run();

   // calculate semi-analytical permeability value
   const real_t K = permeability( setup );

   WALBERLA_CHECK_LESS( std::abs( permEval->currentValue() - K ) / K, setup.epsilon );

   return EXIT_SUCCESS;
}


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   mpi::Environment env( argc, argv );

   const int np = mpi::MPIManager::instance()->numProcesses();

   Setup setup( argc, argv );

   if( np == 1 )
      setup.blocks = Vector3<uint_t>( 1 );
   else if( np == 4 )
      setup.blocks = Vector3<uint_t>( 2, 2, 1 );
   else if( np == 8 )
      setup.blocks = Vector3<uint_t>( 2, 2, 2 );
   else
      WALBERLA_ABORT( "This app can only be executed with 1, 4, or 8 processes! Aborting ..." );

   std::string collisionModel( "TRT" );

   try {
      if( std::strcmp( setup.collisionModel.c_str(), "SRT" ) == 0 )
         return setupAndExecute< lbm::D3Q19< lbm::collision_model::SRT          , false > >( setup );
      else if( std::strcmp( setup.collisionModel.c_str(), "TRT" ) == 0 )
         return setupAndExecute< lbm::D3Q19< lbm::collision_model::TRT          , false > >( setup );
      else if( std::strcmp( setup.collisionModel.c_str(), "MRT" ) == 0 )
         return setupAndExecute< lbm::D3Q19< lbm::collision_model::D3Q19MRT     , false > >( setup );
      else if( std::strcmp( setup.collisionModel.c_str(), "CUM" ) == 0 )
         return setupAndExecute< lbm::D3Q27< lbm::collision_model::D3Q27Cumulant, true  > >( setup );
      else
         WALBERLA_ABORT( "Unexpected LBM collision model: " << setup.collisionModel << "! Aborting ..." )
   } catch( const std::exception & e ) {
      WALBERLA_ABORT( "Unexpected error: " << e.what() << "! Aborting ..." );
   }
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
