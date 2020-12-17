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
//! \file
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

//======================================================================================================================
//
//  THIS FILE IS GENERATED - PLEASE CHANGE THE TEMPLATE !!!
//
//======================================================================================================================

#pragma once

#include <mesa_pd/data/DataTypes.h>
#include <mesa_pd/data/ParticleStorage.h>
#include <mesa_pd/data/shape/BaseShape.h>
{%- for shape in particle.shapes %}
#include <mesa_pd/data/shape/{{shape}}.h>
{%- endfor %}

#include <core/Abort.h>
#include <core/debug/Debug.h>
#include <core/math/AABB.h>

#include <memory>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace data {

struct ShapeStorage
{
   std::vector<std::unique_ptr<BaseShape>> shapes {};

   template <typename ShapeT, typename... Args>
   size_t create(Args&&... args);

   template <typename ReturnType, typename func>
   ReturnType singleDispatch( ParticleStorage& ps, size_t idx, func& f );

   template <typename ReturnType, typename func>
   ReturnType doubleDispatch( ParticleStorage& ps, size_t idx, size_t idy, func& f );
};
//Make sure that no two different shapes have the same unique identifier!
{%- for shape1 in particle.shapes %}
{%- for shape2 in particle.shapes[loop.index:] %}
static_assert( {{shape1}}::SHAPE_TYPE != {{shape2}}::SHAPE_TYPE, "Shape types have to be different!" );
{%- endfor %}
{%- endfor %}

template <typename ShapeT, typename... Args>
size_t ShapeStorage::create(Args&&... args)
{
   shapes.push_back(std::make_unique<ShapeT>(std::forward<Args>(args)...));
   return shapes.size() - 1;
}

template <typename ReturnType, typename func>
ReturnType ShapeStorage::singleDispatch( ParticleStorage& ps, size_t idx, func& f )
{
   WALBERLA_ASSERT_LESS( idx, ps.size() );

   switch (shapes[ps.getShapeID(idx)]->getShapeType())
   {
      {%- for shape in particle.shapes %}
      case {{shape}}::SHAPE_TYPE : return f(ps, idx, *static_cast<{{shape}}*>(shapes[ps.getShapeID(idx)].get()));
      {%- endfor %}
      default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idx)]->getShapeType() << ") could not be determined!");
   }
}

template <typename ReturnType, typename func>
ReturnType ShapeStorage::doubleDispatch( ParticleStorage& ps, size_t idx, size_t idy, func& f )
{
   WALBERLA_ASSERT_LESS( idx, ps.size() );
   WALBERLA_ASSERT_LESS( idy, ps.size() );

   switch (shapes[ps.getShapeID(idx)]->getShapeType())
   {
      {%- for shape1 in particle.shapes %}
      case {{shape1}}::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            {%- for shape2 in particle.shapes %}
            case {{shape2}}::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<{{shape1}}*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<{{shape2}}*>(shapes[ps.getShapeID(idy)].get()));
            {%- endfor %}
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      {%- endfor %}
      default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idx)]->getShapeType() << ") could not be determined!");
   }
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla
