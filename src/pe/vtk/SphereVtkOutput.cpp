

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
//! \file SphereVtkOutput.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "SphereVtkOutput.h"
#include "pe/rigidbody/BodyStorage.h"
#include "pe/rigidbody/Squirmer.h"

#include "core/selectable/IsSetSelected.h"
#include "core/uid/GlobalState.h"


namespace walberla {
namespace pe {

std::vector< SphereVtkOutput::Attributes > SphereVtkOutput::getAttributes() const
{
   std::vector< Attributes > attributes;
   attributes.emplace_back( vtk::typeToString< float >(), "mass", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< float >(), "radius", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< float >(), "velocity", uint_c(3) );
   attributes.emplace_back( vtk::typeToString< float >(), "orientation", uint_c(3) );
   attributes.emplace_back( vtk::typeToString< int >(),   "rank", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< id_t >(),   "id", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< id_t >(),  "uid", uint_c(1) );

   return attributes;
}

void SphereVtkOutput::configure()
{
   bodies_.clear();
   for( auto& block : blockStorage_ )
   {

      const BodyStorage& localStorage = (*(block.getData<const Storage>( storageID_ )))[0];

      for( auto& body : localStorage )
      {
         if (body.getTypeID() == Sphere::getStaticTypeID() || body.getTypeID() == Squirmer::getStaticTypeID())
            bodies_.push_back( static_cast<Sphere const *> (&body) );
         if (body.getTypeID() == Union<Sphere>::getStaticTypeID())
         {
            auto un = static_cast<Union<Sphere> const * > (&body);
            for( auto it2 = un->begin(); it2 != un->end(); ++it2 )
            {
               if (it2->getTypeID() == Sphere::getStaticTypeID())
                  bodies_.push_back( static_cast<ConstSphereID> (it2.getBodyID()) );
            }
         }
         if (body.getTypeID() == Union<Squirmer>::getStaticTypeID())
         {
            auto un = static_cast<Union<Squirmer> const * > (&body);
            for( auto it2 = un->begin(); it2 != un->end(); ++it2 )
            {
               if (it2->getTypeID() == Squirmer::getStaticTypeID())
                  bodies_.push_back( static_cast<ConstSphereID> (it2.getBodyID()) );
            }
         }
         if (body.getTypeID() == Union<Sphere,Squirmer>::getStaticTypeID())
         {
            auto un = static_cast<Union<Sphere,Squirmer> const * > (&body);
            for( auto it2 = un->begin(); it2 != un->end(); ++it2 )
            {
               if (it2->getTypeID() == Sphere::getStaticTypeID() || it2->getTypeID() == Squirmer::getStaticTypeID())
                  bodies_.push_back( static_cast<ConstSphereID> (it2.getBodyID()) );
            }
         }
         if (body.getTypeID() == Union<Squirmer,Sphere>::getStaticTypeID())
         {
            auto un = static_cast<Union<Squirmer,Sphere> const * > (&body);
            for( auto it2 = un->begin(); it2 != un->end(); ++it2 )
            {
               if (it2->getTypeID() == Sphere::getStaticTypeID() || it2->getTypeID() == Squirmer::getStaticTypeID())
                  bodies_.push_back( static_cast<ConstSphereID> (it2.getBodyID()) );
            }
         }
      }
   }
}

std::vector< math::Vector3< real_t > > SphereVtkOutput::getPoints()
{
   std::vector< math::Vector3< real_t > > result( bodies_.size() );
   auto resultIt = result.begin();
   for( auto it = bodies_.begin(); it != bodies_.end(); ++it, ++resultIt )
   {
      *resultIt = (*it)->getPosition();
   }
   return result;
}

} // namespace pe
} // namespace walberla

