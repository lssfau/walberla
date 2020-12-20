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
#include <mesa_pd/data/shape/Sphere.h>
#include <mesa_pd/data/shape/HalfSpace.h>
#include <mesa_pd/data/shape/CylindricalBoundary.h>
#include <mesa_pd/data/shape/Box.h>
#include <mesa_pd/data/shape/Ellipsoid.h>
#include <mesa_pd/data/shape/ConvexPolyhedron.h>

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
static_assert( Sphere::SHAPE_TYPE != HalfSpace::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Sphere::SHAPE_TYPE != CylindricalBoundary::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Sphere::SHAPE_TYPE != Box::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Sphere::SHAPE_TYPE != Ellipsoid::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Sphere::SHAPE_TYPE != ConvexPolyhedron::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( HalfSpace::SHAPE_TYPE != CylindricalBoundary::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( HalfSpace::SHAPE_TYPE != Box::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( HalfSpace::SHAPE_TYPE != Ellipsoid::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( HalfSpace::SHAPE_TYPE != ConvexPolyhedron::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( CylindricalBoundary::SHAPE_TYPE != Box::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( CylindricalBoundary::SHAPE_TYPE != Ellipsoid::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( CylindricalBoundary::SHAPE_TYPE != ConvexPolyhedron::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Box::SHAPE_TYPE != Ellipsoid::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Box::SHAPE_TYPE != ConvexPolyhedron::SHAPE_TYPE, "Shape types have to be different!" );
static_assert( Ellipsoid::SHAPE_TYPE != ConvexPolyhedron::SHAPE_TYPE, "Shape types have to be different!" );

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
      case Sphere::SHAPE_TYPE : return f(ps, idx, *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()));
      case HalfSpace::SHAPE_TYPE : return f(ps, idx, *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()));
      case CylindricalBoundary::SHAPE_TYPE : return f(ps, idx, *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()));
      case Box::SHAPE_TYPE : return f(ps, idx, *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()));
      case Ellipsoid::SHAPE_TYPE : return f(ps, idx, *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()));
      case ConvexPolyhedron::SHAPE_TYPE : return f(ps, idx, *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()));
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
      case Sphere::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            case Sphere::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idy)].get()));
            case HalfSpace::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idy)].get()));
            case CylindricalBoundary::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idy)].get()));
            case Box::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idy)].get()));
            case Ellipsoid::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idy)].get()));
            case ConvexPolyhedron::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idy)].get()));
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      case HalfSpace::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            case Sphere::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idy)].get()));
            case HalfSpace::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idy)].get()));
            case CylindricalBoundary::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idy)].get()));
            case Box::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idy)].get()));
            case Ellipsoid::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idy)].get()));
            case ConvexPolyhedron::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idy)].get()));
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      case CylindricalBoundary::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            case Sphere::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idy)].get()));
            case HalfSpace::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idy)].get()));
            case CylindricalBoundary::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idy)].get()));
            case Box::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idy)].get()));
            case Ellipsoid::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idy)].get()));
            case ConvexPolyhedron::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idy)].get()));
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      case Box::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            case Sphere::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idy)].get()));
            case HalfSpace::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idy)].get()));
            case CylindricalBoundary::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idy)].get()));
            case Box::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idy)].get()));
            case Ellipsoid::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idy)].get()));
            case ConvexPolyhedron::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idy)].get()));
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      case Ellipsoid::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            case Sphere::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idy)].get()));
            case HalfSpace::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idy)].get()));
            case CylindricalBoundary::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idy)].get()));
            case Box::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idy)].get()));
            case Ellipsoid::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idy)].get()));
            case ConvexPolyhedron::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idy)].get()));
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      case ConvexPolyhedron::SHAPE_TYPE :
         switch (shapes[ps.getShapeID(idy)]->getShapeType())
         {
            case Sphere::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Sphere*>(shapes[ps.getShapeID(idy)].get()));
            case HalfSpace::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<HalfSpace*>(shapes[ps.getShapeID(idy)].get()));
            case CylindricalBoundary::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<CylindricalBoundary*>(shapes[ps.getShapeID(idy)].get()));
            case Box::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Box*>(shapes[ps.getShapeID(idy)].get()));
            case Ellipsoid::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<Ellipsoid*>(shapes[ps.getShapeID(idy)].get()));
            case ConvexPolyhedron::SHAPE_TYPE : return f(ps,
                                                   idx,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idx)].get()),
                                                   idy,
                                                   *static_cast<ConvexPolyhedron*>(shapes[ps.getShapeID(idy)].get()));
            default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idy)]->getShapeType() << ") could not be determined!");
         }
      default : WALBERLA_ABORT("Shape type (" << shapes[ps.getShapeID(idx)]->getShapeType() << ") could not be determined!");
   }
}

} //namespace data
} //namespace mesa_pd
} //namespace walberla