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
//! \file BodyVtkOutput.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/RigidBody.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/uid/SUID.h"
#include "domain_decomposition/BlockStorage.h"
#include "vtk/Base64Writer.h"
#include "vtk/PointDataSource.h"
#include "vtk/UtilityFunctions.h"

#include <vector>

namespace walberla {
namespace pe {

class DefaultBodyVTKOutput : public vtk::PointDataSource
{
public:
   DefaultBodyVTKOutput( ConstBlockDataID storageID, const BlockStorage & blockStorage)
      : storageID_( storageID )
      , blockStorage_( blockStorage ) { }

   std::vector< Attributes > getAttributes() const;

   void configure();

   std::vector< Vector3< real_t > > getPoints();

   inline void push( std::ostream& os , const uint_t /*data*/, const uint_t point, const uint_t component );
   inline void push( vtk::Base64Writer& b64, const uint_t /*data*/, const uint_t point, const uint_t component );

private:

   ConstBlockDataID storageID_;
   const BlockStorage & blockStorage_;
   std::vector< ConstBodyID > bodies_;
};


void DefaultBodyVTKOutput::push( std::ostream& os, const uint_t data, const uint_t point, const uint_t component )
{
   WALBERLA_ASSERT_LESS( point, bodies_.size() );
   WALBERLA_ASSERT_LESS( component, 3u );

   switch( data )
   {
   case 0:
      vtk::toStream( os, numeric_cast< float >( bodies_.at( point )->getLinearVel()[component] ) );
      break;
   case 1:
      vtk::toStream( os, numeric_cast< int >( bodies_.at( point )->MPITrait.getOwner().rank_ ) );
      break;
   case 2:
      vtk::toStream( os, 0 );
      break;
   }

}


void DefaultBodyVTKOutput::push( vtk::Base64Writer& b64, const uint_t data, const uint_t point, const uint_t component )
{
   WALBERLA_ASSERT_LESS( point, bodies_.size() );
   WALBERLA_ASSERT_LESS( component, 3u );

   switch( data )
   {
   case 0:
      b64 << numeric_cast< float >( bodies_.at( point )->getLinearVel()[component] );
      break;
   case 1:
      b64 << numeric_cast< int >( bodies_.at( point )->MPITrait.getOwner().rank_ );
      break;
   case 2:
      b64 << 0;
      break;
   }
}


} // namespace pe
} // namespace walberla

