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
//! \file Cylinder.h
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
   * Class representing a Cylinder
   *
   * Implements body concept, i.e. a fastOverlapCheck and a contains() function, which can be used
   * to initialize fields.
   * For details about the body concept see BodyOverlapFunctions.h
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   class Cylinder
   {
   public:

      explicit Cylinder( const Vector3<real_t> & _start, const Vector3<real_t> & _end, real_t _rad );

      void setStart( const Vector3<real_t> & point ) { start_ = point;          update(); }
      void setStart( real_t newVal, uint_t coord   ) { start_[coord] = newVal;  update(); }

      void setEnd( const Vector3<real_t> & point ) { end_ = point;          update(); }
      void setEnd( real_t newVal, uint_t coord   ) { end_[coord] = newVal;  update(); }

      void setRadius( const real_t value ){ rad_ = value; update(); }

      const Vector3<real_t> & start()  const { return start_; }
      const Vector3<real_t> & end()    const { return end_;   }
      const Vector3<real_t> & dir()    const { return dir_;   }
          
      real_t radius() const { return rad_; }
      real_t length() const { return len_; }
      
      const AABB & boundingBox() const { return boundingBox_; }
      
   private:
      /// Recalculates boundingBox_ and dir_
      void update();

      Vector3<real_t> start_;
      Vector3<real_t> end_;
      real_t          rad_;

      Vector3<real_t> dir_;
      real_t          len_;

      AABB boundingBox_;
   };
   
   // Body concept
   template<> FastOverlapResult fastOverlapCheck ( const Cylinder & cyl, const AABB & box );
   //template<> FastOverlapResult fastOverlapCheck ( const Cylinder & cyl, const Vector3<real_t> & cellMidpoint, real_t dx );
   template<> bool contains ( const Cylinder & cyl, const Vector3<real_t> & point );
   
} // namespace geometry
} // namespace walberla
