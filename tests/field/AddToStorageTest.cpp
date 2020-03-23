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
//! \file CollectTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/Environment.h"

#include "blockforest/Initialization.h"

#include "field/AddToStorage.h"


namespace walberla {


int main( int argc, char ** argv )
{
   debug::enterTestMode();
   walberla::Environment walberlaEnv( argc, argv );

   auto  blocks = blockforest::createUniformBlockGrid( 2, 2, 2,       // blocks in x,y,z
                                                       4, 4, 4,       // nr of cells per block
                                                       1.0,           // dx
                                                       false
                                                       );
   typedef GhostLayerField<Vector3<uint_t>,1> VectorField;
   typedef GhostLayerField<uint_t, 3> FlattenedField;
   BlockDataID fieldID = field::addToStorage<VectorField>( blocks, "Field" );
   BlockDataID flattenedID = field::addFlattenedShallowCopyToStorage<VectorField>( blocks, fieldID, "flattened Field");

   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      VectorField * field = blockIt->getData<VectorField>( fieldID );

      for( auto cellIt = field->beginWithGhostLayerXYZ(); cellIt != field->end(); ++cellIt )
      {
         for( uint_t f = 0; f < 3; ++f)
         {
            uint_t val = uint_t(&(*cellIt)[f]);
            (*cellIt)[f] = val;
         }
      }
   }
   
   BlockDataID copyID = field::addCloneToStorage<VectorField>( blocks, fieldID, "copied Field");
   
   for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
   {
      VectorField * field = blockIt->getData<VectorField>( fieldID );
      FlattenedField * flattened = blockIt->getData<FlattenedField>( flattenedID );
      VectorField * copy = blockIt->getData<VectorField>( copyID );

      for( auto cellIt = field->beginWithGhostLayerXYZ(); cellIt != field->end(); ++cellIt )
      {
         for( uint_t f = 0; f < 3; ++f)
         {
            WALBERLA_CHECK_EQUAL((*cellIt)[f], copy->get(cellIt.x(), cellIt.y(), cellIt.z())[f]);
            WALBERLA_CHECK_UNEQUAL(&(*cellIt)[f], &(copy->get(cellIt.x(), cellIt.y(), cellIt.z())[f]));
            WALBERLA_CHECK_EQUAL(&(*cellIt)[f], &(flattened->get(cellIt.x(), cellIt.y(), cellIt.z(), f)));
         }
      }
   }
   
   return 0;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc,argv);
}
