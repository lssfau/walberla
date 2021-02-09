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
//! \file Sphere.h
//! \ingroup geometry
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "BodyOverlapFunctions.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

namespace walberla {
namespace geometry {

class AbstractBody {
public:
   virtual ~AbstractBody() = default;
   virtual bool contains (const Vector3<real_t> & point ) const = 0;
   virtual FastOverlapResult fastOverlapCheck ( const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx ) const = 0;
   virtual FastOverlapResult fastOverlapCheck ( const AABB & box ) const = 0;
};

// this is an adaptor from static to dynamic polymorphism
template<typename Body>
class DynamicBody : public AbstractBody 
{
public:
   DynamicBody( const Body & b )
      : body_(b)
   {}
   
   bool contains (const Vector3<real_t> & point ) const override
   {
        return geometry::contains( body_, point );
   }
   FastOverlapResult fastOverlapCheck ( const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx ) const override
   {
      return geometry::fastOverlapCheck( body_, cellMidpoint, dx );
   }
   FastOverlapResult fastOverlapCheck ( const AABB & box ) const override
   {
      return geometry::fastOverlapCheck( body_, box);
   }
private:
   const Body body_;
};

template<typename Body>
shared_ptr<DynamicBody<Body> > make_DynamicBody(Body body)
{
   return make_shared<DynamicBody<Body> >( DynamicBody<Body>(body) );
}

//===================================================================================================================
//
//  Body concept implementation
//
//===================================================================================================================


template<>
inline FastOverlapResult fastOverlapCheck ( const AbstractBody & body, const AABB & box )
{
   return body.fastOverlapCheck( box );
}

template<>
inline FastOverlapResult fastOverlapCheck ( const AbstractBody & body, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx )
{
   return body.fastOverlapCheck( cellMidpoint, dx );
}

template<>
inline bool contains ( const AbstractBody & body, const Vector3<real_t> & point )
{
   return body.contains( point );
}



} // namespace geometry
} // namespace walberla
