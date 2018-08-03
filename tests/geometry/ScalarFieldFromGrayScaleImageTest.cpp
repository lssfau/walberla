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
//! \file ScalarFieldFromGrayScaleImageTest.cpp
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "geometry/initializer/InitializationManager.h"
#include "geometry/initializer/ScalarFieldFromGrayScaleImage.h"

#include "blockforest/Initialization.h"

#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
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


int main( int argc, char ** argv )
{
   debug::enterTestMode();

   auto mpiManager = MPIManager::instance();
   mpiManager->initializeMPI( &argc, &argv );

   using blockforest::createUniformBlockGrid;
   auto blocks = createUniformBlockGrid(confBlockCount[0],confBlockCount[1],confBlockCount[2], //blocks
                                        confCells[0],confCells[1],confCells[2],    //cells
                                        real_c(0.7),                                       //dx
                                        false,                                      //one block per process
                                        false,false,false );                       //periodicity


   BlockDataID scalarFieldID = field::addToStorage<GhostLayerField<real_t,1> > ( blocks, "ScalarField" );

   // Geometry Initialization from config file
   using namespace geometry::initializer;

   auto geometryInitializationManager = std::make_shared<InitializationManager> ( blocks->getBlockStorage() );
   auto freeSurfaceInitializer        = std::make_shared<ScalarFieldFromGrayScaleImage> ( *blocks, scalarFieldID );

   geometryInitializationManager->registerInitializer( "FreeSurfaceImage", freeSurfaceInitializer );

   {
      std::ofstream configFile ( "sampleImage.dat" );
      configFile << "Geometry  {   FreeSurfaceImage   {  \n";

      configFile << "file test.png; \n";
      configFile << "extrusionCoordinate 0;";
      configFile << "lowerExtrusionLimit 5;";
      configFile << "upperExtrusionLimit 10;";
      configFile << "xOffset 5; ";
      configFile << "yOffset -5; ";

      configFile << "}  } \n";
   }
   Config cfg;
   cfg.readParameterFile( "sampleImage.dat" );
   geometryInitializationManager->init( cfg.getBlock("Geometry") );


   if ( useGui )
   {
      SweepTimeloop timeloop ( blocks, 100 );
      GUI gui ( timeloop, blocks, argc, argv );
      gui.run();
   }

   return EXIT_SUCCESS;
}
}

int main( int argc, char ** argv ){
   return walberla::main(argc, argv);
}
