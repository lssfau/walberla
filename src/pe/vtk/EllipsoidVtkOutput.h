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
//! \file EllipsoidVtkOutput.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "pe/rigidbody/RigidBody.h"
#include "pe/rigidbody/Ellipsoid.h"

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/BlockStorage.h"

#include "vtk/Base64Writer.h"
#include "vtk/PointDataSource.h"
#include "vtk/UtilityFunctions.h"

#include <vector>
#include <array>

namespace walberla {
namespace pe {

class EllipsoidVtkOutput : public vtk::PointDataSource
{
public:
   EllipsoidVtkOutput( ConstBlockDataID storageID, const BlockStorage & blockStorage)
      : storageID_( storageID )
      , blockStorage_( blockStorage ) { }

   std::vector< Attributes > getAttributes() const override;

   void configure() override;

   std::vector< Vector3< real_t > > getPoints() override;

   inline void push( std::ostream& os , const uint_t /*data*/, const uint_t point, const uint_t component ) override;
   inline void push( vtk::Base64Writer& b64, const uint_t /*data*/, const uint_t point, const uint_t component ) override;

private:

   ConstBlockDataID storageID_;
   const BlockStorage & blockStorage_;
   std::vector< Ellipsoid const * > bodies_;
   std::vector< std::array<real_t,6> > tensorGlyphs_;
};


void EllipsoidVtkOutput::push( std::ostream& os, const uint_t data, const uint_t point, const uint_t component )
{
   WALBERLA_ASSERT_LESS( point, bodies_.size() );
   WALBERLA_ASSERT_LESS( component, 6u );

   switch( data )
   {
   case 0:
      if( bodies_.at( point )->hasInfiniteMass() )
      {
         vtk::toStream( os, std::numeric_limits< float >::infinity() );
      }
      else
      {
         vtk::toStream( os, numeric_cast< float >(bodies_.at( point )->getMass()) );
      }
      break;
   case 1:
      vtk::toStream( os, numeric_cast< float >(tensorGlyphs_.at(point)[component]) );
      break;
   case 2:
      vtk::toStream( os, numeric_cast< float >(bodies_.at( point )->getLinearVel()[component]) );
      break;
   case 3:
      vtk::toStream( os, numeric_cast< int >(bodies_.at( point )->MPITrait.getOwner().rank_) );
      break;
   case 4:
      vtk::toStream( os, numeric_cast< id_t >(bodies_.at( point )->getTopSuperBody()->getSystemID()) );
      break;
   case 5:
      vtk::toStream( os, numeric_cast< id_t >(bodies_.at( point )->getTopSuperBody()->getID()) );
      break;
   }

}



void EllipsoidVtkOutput::push( vtk::Base64Writer& b64, const uint_t data, const uint_t point, const uint_t component )
{
   WALBERLA_ASSERT_LESS( point, bodies_.size() );
   WALBERLA_ASSERT_LESS( component, 6u );

   switch( data )
   {
   case 0:
      if( bodies_.at( point )->hasInfiniteMass() )
      {
         b64 << std::numeric_limits< float >::infinity();
      }
      else
      {
         b64 << numeric_cast< float >(bodies_.at( point )->getMass());
      }
      break;
   case 1:
      b64 << numeric_cast< float >(tensorGlyphs_.at(point)[component]);
      break;
   case 2:
      b64 << numeric_cast< float >(bodies_.at( point )->getLinearVel()[component]);
      break;
   case 3:
      b64 << numeric_cast< int >(bodies_.at( point )->MPITrait.getOwner().rank_);
      break;
   case 4:
      b64 << numeric_cast< id_t >(bodies_.at( point )->getTopSuperBody()->getSystemID());
      break;
   case 5:
      b64 << numeric_cast< id_t >(bodies_.at( point )->getTopSuperBody()->getID());
      break;
   }
}


} // namespace pe
} // namespace walberla

