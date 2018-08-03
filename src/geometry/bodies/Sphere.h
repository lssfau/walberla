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
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BodyOverlapFunctions.h"
#include "core/config/Config.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"


namespace walberla {
namespace geometry {


   //*******************************************************************************************************************
   /*!
   * Class representing a Sphere
   *
   * Implements body concept, i.e. a fastOverlapCheck and a contains() function, which can be used
   * to initialize fields.
   * For details about the body concept see BodyOverlapFunctions.h
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   class Sphere
   {
   public:

      explicit Sphere( const Vector3<real_t> & midp, real_t rad );
      Sphere( const Sphere & o ) = default;

      void setMidpoint( const Vector3<real_t> & point ) { midpoint_ = point;         updateBoxes(); }
      void setMidpoint( real_t newVal, uint_t coord )   { midpoint_[coord] = newVal; updateBoxes(); }
      void setRadius  ( real_t newRadius)               { radius_ = newRadius;       updateBoxes(); }

      const Vector3<real_t> & midpoint() const { return midpoint_; }
                     real_t   radius()   const { return radius_;   }


      const AABB & boundingBox() const { return boundingBox_; }
      const AABB & innerBox()    const { return innerBox_;    }


   private:
      /// Recalculates boundingBox_ and innerBox_
      void updateBoxes();

      Vector3<real_t> midpoint_;
      real_t          radius_;

      AABB boundingBox_;
      AABB innerBox_;
   };


   // Body concept
   template<> FastOverlapResult fastOverlapCheck ( const Sphere & sphere, const AABB & box );
   template<> FastOverlapResult fastOverlapCheck ( const Sphere & sphere, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx );
   template<> bool contains ( const Sphere & sphere, const Vector3<real_t> & point );



} // namespace geometry
} // namespace walberla



