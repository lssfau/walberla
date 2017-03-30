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
//! \file Torus.h
//! \ingroup geometry
//! \author Matthias Markl <matthias.markl@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "BodyOverlapFunctions.h"

namespace walberla {
namespace geometry {


   //*******************************************************************************************************************
   /*!
   * Class representing a Torus
   *
   * Implements body concept, i.e. a fastOverlapCheck and a contains() function, which can be used
   * to initialize fields.
   * For details about the body concept see BodyOverlapFunctions.h
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   class Torus
   {
   public:

      explicit Torus( const Vector3<real_t> & _midPnt, const Vector3<real_t> & _normal, real_t _radius, real_t _distance );

      void setMidPoint( const Vector3<real_t> & _midPnt ) { midPnt_ = _midPnt;       update(); }
      void setMidPoint( real_t newVal, uint_t coord     ) { midPnt_[coord] = newVal; update(); }

      void setNormal( const Vector3<real_t> & _normal ) { normal_ = _normal;       update(); }
      void setNormal( real_t newVal, uint_t coord     ) { normal_[coord] = newVal; update(); }

      void setRadius  ( const real_t value ){ radius_ = value; update(); }
      void setDistance( const real_t value ){ distnc_ = value; update(); }

      const Vector3<real_t> & normal() const { return normal_; }
      const Vector3<real_t> & defTan() const { return defTan_; }
      const Vector3<real_t> & comTan() const { return comTan_; }
      const Vector3<real_t> & midPnt() const { return midPnt_; }

      real_t radius()   const { return radius_; }
      real_t distance() const { return distnc_; }
      
      const AABB & boundingBox() const { return boundingBox_; }
      
   private:
      /// Recalculates boundingBox_ and dir_
      void update();

      Vector3<real_t> midPnt_;
      Vector3<real_t> normal_;
      real_t          radius_;
      real_t          distnc_;

      Vector3<real_t> defTan_;
      Vector3<real_t> comTan_;
      
      AABB boundingBox_;
   };
   
   // Body concept
   template<> FastOverlapResult fastOverlapCheck ( const Torus & torus, const AABB & box );
   //template<> FastOverlapResult fastOverlapCheck ( const Torus & torus, const Vector3<real_t> & cellMidpoint, real_t dx );
   template<> bool contains ( const Torus & torus, const Vector3<real_t> & point );
   
} // namespace geometry
} // namespace walberla
