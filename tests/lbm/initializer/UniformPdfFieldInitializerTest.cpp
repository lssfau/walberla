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
//! \file UniformPdfFieldInitializerTest.cpp
//! \author Tobias Schruff <tobias.schruff@gmail.com>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"
#include "field/all.h"
#include "lbm/all.h"

#include <vector>
#include <string>
#include <sstream>


namespace walberla {


using LatticeModel_T = lbm::D3Q19<lbm::collision_model::SRT>;
using PdfField_T = lbm::PdfField<LatticeModel_T>;
using PdfFieldInitializer_T = lbm::initializer::PdfFieldInitializer<LatticeModel_T>;


struct DensityInit
{
   real_t operator()( const Cell & globalCell ) const
   {
      return real_c( globalCell.x() + globalCell.y() + globalCell.z() );
   }
};


struct VelocityInit
{
   Vector3<real_t> operator()( const Cell & globalCell ) const
   {
      return Vector3<real_t>( real_c( globalCell.x() ), real_c( globalCell.y() ), real_c( globalCell.z() ) );
   }
};

struct DensityAndVelocityInit
{
   std::vector<real_t> operator()( const Cell & globalCell ) const
   {
      DensityInit densityInit;
      VelocityInit velocityInit;

      real_t density = densityInit( globalCell );
      Vector3<real_t> velocity = velocityInit( globalCell );

      std::vector<real_t> result;
      result.push_back( density );
      result.push_back( velocity[0] );
      result.push_back( velocity[1] );
      result.push_back( velocity[2] );

      return result;
   }
};


void writeVTK( const std::string & identifier, const shared_ptr<StructuredBlockForest> & blocks, const BlockDataID & pdfFieldId )
{
    auto pdfFieldVTKWriter = vtk::createVTKOutput_BlockData( blocks, identifier, 1u, uint_t(1) );

    blockforest::communication::UniformBufferedScheme<stencil::D3Q27> pdfGhostLayerSync( blocks );
    pdfGhostLayerSync.addPackInfo( make_shared< field::communication::PackInfo<PdfField_T> >( pdfFieldId ) );
    pdfFieldVTKWriter->addBeforeFunction( pdfGhostLayerSync );

    auto velocityWriter = make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >( pdfFieldId, "Velocity (Lattice)" );
    auto  densityWriter = make_shared< lbm::DensityVTKWriter < LatticeModel_T, float > >( pdfFieldId, "Density (Lattice)" );

    pdfFieldVTKWriter->addCellDataWriter( velocityWriter );
    pdfFieldVTKWriter->addCellDataWriter( densityWriter );

    pdfFieldVTKWriter->write();
}


void testDensityInit( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks, bool writeVtk = false )
{
   lbm::initializer::PdfFieldInitializer< LatticeModel_T > initializer( pdfFieldId, blocks );

   DensityInit densityInit;
   initializer.initDensity( densityInit );

   if( writeVtk ) writeVTK( "density_init", blocks, pdfFieldId );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto pdfField = blockIt->getData<PdfField_T>( pdfFieldId );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField )

      for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Cell globalCell( cellIt.cell() );
         blocks->transformBlockLocalToGlobalCell( globalCell, *blockIt );

         WALBERLA_CHECK_FLOAT_EQUAL( pdfField->getDensity( cellIt.cell() ), densityInit( globalCell ) )
      }
   }
}


void testVelocityInit( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks, bool writeVtk = false )
{
   lbm::initializer::PdfFieldInitializer< LatticeModel_T > initializer( pdfFieldId, blocks );

   VelocityInit velocityInit;
   initializer.initVelocity( velocityInit );

   if( writeVtk ) writeVTK( "velocity_init", blocks, pdfFieldId );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto pdfField = blockIt->getData<PdfField_T>( pdfFieldId );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField )

      for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Cell globalCell( cellIt.cell() );
         blocks->transformBlockLocalToGlobalCell( globalCell, *blockIt );

         const Vector3<real_t> initV  = velocityInit( globalCell );
         const Vector3<real_t> fieldV = pdfField->getVelocity( cellIt.cell() );

         WALBERLA_CHECK_FLOAT_EQUAL( initV[0], fieldV[0] )
         WALBERLA_CHECK_FLOAT_EQUAL( initV[1], fieldV[1] )
         WALBERLA_CHECK_FLOAT_EQUAL( initV[2], fieldV[2] )
      }
   }
}


void testDensityAndVelocityInit( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks, bool writeVtk = false )
{
   lbm::initializer::PdfFieldInitializer< LatticeModel_T > initializer( pdfFieldId, blocks );

   DensityAndVelocityInit densityAndVelocityInit;
   initializer.initDensityAndVelocity( densityAndVelocityInit );

   if( writeVtk ) writeVTK( "density_and_velocity_init", blocks, pdfFieldId );

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      auto pdfField = blockIt->getData<PdfField_T>( pdfFieldId );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField )

      for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
      {
         Cell globalCell( cellIt.cell() );
         blocks->transformBlockLocalToGlobalCell( globalCell, *blockIt );

         Vector3<real_t> velocity;
         const real_t density = pdfField->getDensityAndVelocity( velocity, cellIt.cell() );
         const std::vector<real_t> control = densityAndVelocityInit( globalCell );

         WALBERLA_CHECK_FLOAT_EQUAL( control[0], density     )
         WALBERLA_CHECK_FLOAT_EQUAL( control[1], velocity[0] )
         WALBERLA_CHECK_FLOAT_EQUAL( control[2], velocity[1] )
         WALBERLA_CHECK_FLOAT_EQUAL( control[3], velocity[2] )
      }
   }
}


void testDensityAndVelocityInitFromConfig( const BlockDataID & pdfFieldId, const shared_ptr<StructuredBlockForest> & blocks, const shared_ptr<Config> & config, bool writeVtk = false )
{
   PdfFieldInitializer_T initializer( pdfFieldId, blocks );

   Config::Blocks initBlocks;
   config->getBlocks( "Initializer", initBlocks );

   uint_t configCounter( 1 );

   for( auto configIt = initBlocks.begin(); configIt != initBlocks.end(); ++configIt )
   {
      initializer.initFromConfig( *configIt );

      if( writeVtk )
      {
         std::ostringstream ss;
         ss << "config_init_" << configCounter;
         writeVTK( ss.str(), blocks, pdfFieldId );
      }

      lbm::initializer::ExprSystemInitFunction initControl( blocks );
      initControl.parse( *configIt );

      for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
      {
         auto pdfField = blockIt->getData<PdfField_T>( pdfFieldId );
         WALBERLA_ASSERT_NOT_NULLPTR( pdfField )

         for( auto cellIt = pdfField->beginWithGhostLayerXYZ(); cellIt != pdfField->end(); ++cellIt )
         {
            Cell globalCell( cellIt.cell() );
            blocks->transformBlockLocalToGlobalCell( globalCell, *blockIt );

            Vector3<real_t> velocity;
            const real_t density = pdfField->getDensityAndVelocity( velocity, cellIt.cell() );
            const std::vector<real_t> control = initControl( globalCell );

            WALBERLA_CHECK_FLOAT_EQUAL( control[0], density     )
            WALBERLA_CHECK_FLOAT_EQUAL( control[1], velocity[0] )
            WALBERLA_CHECK_FLOAT_EQUAL( control[2], velocity[1] )
            WALBERLA_CHECK_FLOAT_EQUAL( control[3], velocity[2] )
         }
      }

      configCounter++;
   }
}


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   Environment env( argc, argv );
   if( !env.config() )
      WALBERLA_ABORT("You have to specify a configuration file!")

   auto blocks = blockforest::createUniformBlockGrid( 2u, 2u, 1u, 5u, 5u, 10u, real_c(1) );

   LatticeModel_T lm( real_c(1) );

   BlockDataID pdfFieldId = lbm::addPdfFieldToStorage( blocks, "pdf_field", lm );

   testDensityInit                     ( pdfFieldId, blocks              , false );
   testVelocityInit                    ( pdfFieldId, blocks              , false );
   testDensityAndVelocityInit          ( pdfFieldId, blocks              , false );
   testDensityAndVelocityInitFromConfig( pdfFieldId, blocks, env.config(), false );

   return 0;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}