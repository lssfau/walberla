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
//! \file ShapeTypes.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/Sphere.h>

namespace walberla {
namespace mesa_pd {
namespace data {

const int Box::SHAPE_TYPE                ;
const int CylindricalBoundary::SHAPE_TYPE;
const int HalfSpace::SHAPE_TYPE          ;
const int Ellipsoid::SHAPE_TYPE          ;
const int Sphere::SHAPE_TYPE             ;

} //namespace data
} //namespace mesa_pd
} //namespace walberla
