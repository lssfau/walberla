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
//! \file Ellipsoid.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BodyOverlapFunctions.h"
#include "core/config/Config.h"
#include "core/math/AABB.h"
#include "core/math/Matrix3.h"
#include "core/math/Vector3.h"


namespace walberla {
namespace geometry {


   //*******************************************************************************************************************
   /*! Class representing an Ellipsoid in 3D
   *
   * Implements body concept, i.e. a fastOverlapCheck and a contains() function, which can be used
   * to initialize fields.
   * For details about the body concept see BodyOverlapFunctions.h
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   class Ellipsoid
   {
   public:

      /*************************************************************************************************************//**
      *  Ellipsoid constructor
      *
      *  \param midpoint    midpoint of the ellipsoid
      *  \param axis1       the first semi-axis of the ellipsoid ( associated with radius1 )
      *                     not required to be normalized, is normalized internally
      *  \param axis2       second semi-axis (must not be parallel to axis1 ) should ideally be perpendicular to axis1
      *                     if not perpendicular, axis2 is made perpendicular to axis1 using one Gram-Schmitt step.
      *                     The third semi-axis is computed as cross-product of axis1 and axis2
      *  \param radii       the length of the semi axes
      *****************************************************************************************************************/
      explicit Ellipsoid( const Vector3<real_t> & midpoint,
               Vector3<real_t> axis1, Vector3<real_t> axis2, const Vector3<real_t>& radii );


      const AABB & boundingBox() const { return boundingBox_; }

      bool contains( const Vector3<real_t> & point ) const;

      const Vector3<real_t> midpoint() const { return midpoint_; }

      real_t maxRadius() const { return maxRadius_; }
      real_t minRadius() const { return minRadius_; }

   private:


      void updateBoundingBox();

      Vector3<real_t> midpoint_;
      Matrix3<real_t> mat_;             //< (x-midpoint_)' * mat_ * (x-midpoint) < 1   defines ellipsoid interior
      Matrix3<real_t> rotationMatrix_;  //< required to update bounding box


      Vector3<real_t> radii_;          //< stored explicitly for updating bounding box
      real_t maxRadius_;
      real_t minRadius_;

      AABB boundingBox_;
   };


   // Body concept
   template<> FastOverlapResult fastOverlapCheck ( const Ellipsoid & e, const AABB & box );
   template<> FastOverlapResult fastOverlapCheck ( const Ellipsoid & e, const Vector3<real_t> & cellMidpoint, const Vector3<real_t> & dx );
   template<> bool contains ( const Ellipsoid & ellipsoid, const Vector3<real_t> & point );


} // namespace geometry
} // namespace walberla



