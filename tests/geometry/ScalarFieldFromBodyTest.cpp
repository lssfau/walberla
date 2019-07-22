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
//! \file ScalarFieldFromBodyTest.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Tests setup of scalar field from configuration file
//
//======================================================================================================================

#include "geometry/bodies/Ellipsoid.h"
#include "geometry/bodies/Sphere.h"
#include "geometry/initializer/InitializationManager.h"
#include "geometry/initializer/OverlapFieldFromBody.h"

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/math/Vector3.h"

#include "field/AddToStorage.h"

#include "gui/Gui.h"

#include "timeloop/SweepTimeloop.h"

#include <fstream>
#include <memory>


namespace walberla {
using namespace geometry;

const uint_t confBlockCount []      = { 1, 1, 1 };
const uint_t confCells []           = { 30, 30, 30 };

const bool   useGui = false;


typedef GhostLayerField<real_t,1> ScalarField;

//======================================================================================================================
//
//  Helper Functions
//
//======================================================================================================================


void resetField( StructuredBlockStorage & blocks, BlockDataID & fieldID, real_t value = 1.0 )
{
   for( auto b = blocks.begin(); b != blocks.end(); ++b )
   {
      GhostLayerField<real_t,1> * field = b->getData<GhostLayerField<real_t,1> >( fieldID );
      field->setWithGhostLayer( value );
   }
}


real_t getVolume( const StructuredBlockStorage & storage, BlockDataID fillLevelField )
{
   real_t volume = 0;

   for ( auto block = storage.begin(); block != storage.end(); ++block )
   {
      auto fl = block->getData<GhostLayerField<real_t,1> >( fillLevelField );

      for ( auto cell = fl->begin(); cell != fl->end(); ++cell )
         volume += ( 1 - *cell );
   }

   const real_t cellVolume = storage.dx() * storage.dy() * storage.dz();

   return volume * cellVolume;
}



//======================================================================================================================
//
//  Tests
//
//======================================================================================================================


void boxTest( StructuredBlockStorage & storage,
              BlockDataID fieldID,
              initializer::InitializationManager & manager,
              initializer::OverlapFieldFromBody & initializer )
{
   const Vector3<real_t> min ( real_c(5.2),  real_c(11.0), real_c(3.2) );
   const Vector3<real_t> max ( real_c(10.5), real_c(14.2), real_c(17.2) );
   const Vector3<real_t> diff = max - min;
   const real_t expectedVolume = diff[0] * diff[1] * diff[2];

   resetField( storage, fieldID );

   // Write Sample configuration file - and initialize ellipsoid using the file
   {
      std::ofstream configFile ( "sampleBox.dat" );
      configFile << "Geometry  {   FreeSurface   {  \n";

      configFile << "Box { bubble; \n "
                 << "shape " << "Box;"
                 << "min " << min << "; "
                 << "max " << max << "; "
                 << "} ";

      configFile << "}  } \n";
   }

   Config cfg;
   cfg.readParameterFile( "sampleBox.dat" );

   manager.init( cfg.getBlock("Geometry") );

   real_t volumeFromFileInit = getVolume( storage, fieldID );

   resetField( storage, fieldID );

   // Manual initialization
   initializer.init( AABB(min[0], min[1], min[2], max[0],max[1], max[2] ), false );
   real_t volumeFromManualInit = getVolume( storage, fieldID );


   WALBERLA_LOG_RESULT("Box initialization test (Manual,File,Expected): ("
                << volumeFromManualInit << ","
                << volumeFromFileInit << ','
                << expectedVolume << ")" );

   WALBERLA_CHECK_FLOAT_EQUAL( volumeFromManualInit, volumeFromFileInit );
   WALBERLA_CHECK_FLOAT_EQUAL( volumeFromManualInit, expectedVolume );

}

void ellipsoidTest( StructuredBlockStorage & storage,
                    BlockDataID fieldID,
                    initializer::InitializationManager & manager,
                    initializer::OverlapFieldFromBody & initializer )
{
   const Vector3<real_t> midpoint ( real_t(15), real_t(15), real_t(15) );
   const Vector3<real_t> radii (  real_t(5), real_t(3), real_t(4) );
   const Vector3<real_t> axis1 (  real_t(1), real_t(1), real_t(0) );
   const Vector3<real_t> axis2 ( -real_t(1), real_t(1), real_t(0) );

   const real_t expectedVolume = real_t(4) / real_t(3) * radii[0] * radii[1] * radii[2] * math::pi;

   resetField( storage, fieldID );

   // Write Sample configuration file - and initialize ellipsoid using the file
   {
      std::ofstream configFile ( "sampleEllipsoid.dat" );
      configFile << "Geometry  {   FreeSurface   {  \n";

      configFile << "ellipsoid { bubble; midpoint " << midpoint << ";\n "
                 << "shape"  << " ellipsoid" << "; "
                 << "radii " << radii << "; "
                 << "axis1 " << axis1 << "; "
                 << "axis2 " << axis2 << "; "
                 << "} ";

      configFile << "}  } \n";
   }

   Config cfg;
   cfg.readParameterFile( "sampleEllipsoid.dat" );

   manager.init( cfg.getBlock("Geometry") );

   real_t volumeFromFileInit = getVolume( storage, fieldID );

   resetField( storage, fieldID );

   // Manual initialization

   initializer.init( Ellipsoid( midpoint, axis1, axis2, radii ), false );
   real_t volumeFromManualInit = getVolume( storage, fieldID );


   WALBERLA_LOG_RESULT("Ellipsoid initialization test (Manual,File,Expected): ("
                << volumeFromManualInit << ","
                << volumeFromFileInit << ','
                << expectedVolume << ")" );

   WALBERLA_CHECK_FLOAT_EQUAL( volumeFromManualInit, volumeFromFileInit );

   WALBERLA_CHECK_LESS ( std::abs( expectedVolume - volumeFromManualInit ),  0.2 );
}



void sphereTest( StructuredBlockStorage & storage,
                 BlockDataID fieldID,
                 initializer::InitializationManager & manager,
                 initializer::OverlapFieldFromBody & initializer )
{
   const Vector3<real_t> midpoint ( real_t(15), real_t(15), real_t(15) );
   const real_t radius = real_t(5);
   const real_t expectedVolume = real_t(4) / real_t(3) * radius * radius * radius * math::pi;

   resetField( storage, fieldID );

   // Write Sample configuration file - and initialize sphere using the file
   {
      std::ofstream configFile ( "sampleSphere.dat" );
      configFile << "Geometry  {   FreeSurface   {  \n";
      configFile << "Sphere { shape sphere; bubble; midpoint " << midpoint << "; radius " << radius << "; } ";
      configFile << "}  } \n";
   }
   Config cfg;
   cfg.readParameterFile( "sampleSphere.dat" );

   manager.init( cfg.getBlock("Geometry") );

   real_t volumeFromFileInit = getVolume( storage, fieldID );

   resetField( storage, fieldID );


   // Manual initialization

   initializer.init( Sphere( midpoint, radius), false );
   real_t volumeFromManualInit = getVolume( storage, fieldID );


   WALBERLA_LOG_RESULT("Sphere initialization test (Manual,File,Expected): ("
                << volumeFromManualInit << ","
                << volumeFromFileInit << ','
                << expectedVolume << ")" );

   WALBERLA_CHECK_FLOAT_EQUAL( volumeFromManualInit, volumeFromFileInit );

   WALBERLA_CHECK_LESS ( std::abs( expectedVolume - volumeFromManualInit ),  0.1 );
}



int main( int argc, char ** argv )
{
   debug::enterTestMode();

   auto mpiManager = MPIManager::instance();
   mpiManager->initializeMPI( &argc, &argv );

   using blockforest::createUniformBlockGrid;
   auto blocks = createUniformBlockGrid(confBlockCount[0],confBlockCount[1],confBlockCount[2], //blocks
                                        confCells[0],confCells[1],confCells[2],    //cells
                                        real_c(0.7),                                       //dx
                                        true,                                      //one block per process
                                        false,false,false );                       //periodicity


   BlockDataID scalarFieldID = field::addToStorage<ScalarField> ( blocks, "ScalarField" );

   // Geometry Initialization from config file


   using namespace geometry::initializer;

   auto geometryInitializationManager = std::make_shared<InitializationManager> ( blocks->getBlockStorage() );
   auto freeSurfaceInitializer        = std::make_shared<OverlapFieldFromBody> ( *blocks, scalarFieldID, "drop", "bubble" );

   geometryInitializationManager->registerInitializer( "FreeSurface", freeSurfaceInitializer );

   sphereTest   ( *blocks, scalarFieldID, *geometryInitializationManager, *freeSurfaceInitializer );
   ellipsoidTest( *blocks, scalarFieldID, *geometryInitializationManager, *freeSurfaceInitializer );
   boxTest      ( *blocks, scalarFieldID, *geometryInitializationManager, *freeSurfaceInitializer );


   if ( useGui )
   {
      SweepTimeloop timeloop ( blocks, 100 );
      GUI gui ( timeloop, blocks, argc, argv );
      gui.run();
   }

  return EXIT_SUCCESS;
}
} // namespace walberla

int main( int argc, char* argv[] )
{
  return walberla::main( argc, argv );
}
