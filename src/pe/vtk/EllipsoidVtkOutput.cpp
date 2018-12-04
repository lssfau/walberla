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
//! \file EllipsoidVtkOutput.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "EllipsoidVtkOutput.h"
#include "pe/rigidbody/BodyStorage.h"

#include "core/selectable/IsSetSelected.h"
#include "core/uid/GlobalState.h"


namespace walberla {
namespace pe {

std::vector< EllipsoidVtkOutput::Attributes > EllipsoidVtkOutput::getAttributes() const
{
   std::vector< Attributes > attributes;
   attributes.emplace_back( vtk::typeToString< float >(), "mass", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< float >(), "tensorGlyph", uint_c(6) );
   attributes.emplace_back( vtk::typeToString< float >(), "velocity", uint_c(3) );
   attributes.emplace_back( vtk::typeToString< int >(),   "rank", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< id_t >(),  "id", uint_c(1) );
   attributes.emplace_back( vtk::typeToString< id_t >(),  "uid", uint_c(1) );

   return attributes;
}

void EllipsoidVtkOutput::configure()
{
   bodies_.clear();
   tensorGlyphs_.clear();

   for( auto& block : blockStorage_ )
   {

      const BodyStorage& localStorage = (*(block.getData<const Storage>( storageID_ )))[0];

      for( auto& body : localStorage )
      {
         if (body.getTypeID() == Ellipsoid::getStaticTypeID())
         {
            auto ellipsoid = static_cast<ConstEllipsoidID> (&body);
            bodies_.push_back(ellipsoid);

            // compute tensor glyph for visualization with ParaView (tensorGlyph)
            Mat3 rotMat = ellipsoid->getRotation();
            Vector3<real_t> directionVectorX(rotMat[0], rotMat[3], rotMat[6]);
            Vector3<real_t> directionVectorY(rotMat[1], rotMat[4], rotMat[7]);
            Vector3<real_t> directionVectorZ(rotMat[2], rotMat[5], rotMat[8]);
            Vector3<real_t> semiAxes = ellipsoid->getSemiAxes();
            Mat3 axa = math::dyadicProduct(directionVectorX, directionVectorX);
            Mat3 bxb = math::dyadicProduct(directionVectorY, directionVectorY);
            Mat3 cxc = math::dyadicProduct(directionVectorZ, directionVectorZ);
            Mat3 tensor = axa * semiAxes[0] + bxb * semiAxes[1] + cxc * semiAxes[2];
            // use symmetry to only write 6 of the 9 elements: XX YY ZZ XY YZ XZ
            tensorGlyphs_.push_back({{tensor(0,0), tensor(1,1), tensor(2,2),
                                      tensor(0,1), tensor(1,2), tensor(0,2)}});
         }
      }
   }
}

std::vector< math::Vector3< real_t > > EllipsoidVtkOutput::getPoints()
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

