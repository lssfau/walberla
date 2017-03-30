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

   
   //*******************************************************************************************************************
   /*!
   * Class representing a union, difference, etc. between two bodies
   *
   * \ingroup geometry
   *
   */
   //*******************************************************************************************************************
   
   template<typename A, typename B>
   class BodyLogicalOperationBinary
   {
   public:

      explicit BodyLogicalOperationBinary( const shared_ptr<A> &a, const shared_ptr<B> &b )
         : a_(a), b_(b)
      {}
      
      const A& getA() const
      {
         return *a_;
      }
      
      const B& getB() const
      {
         return *b_;
      }

   private:
      const shared_ptr<A> a_;
      const shared_ptr<B> b_;
   };
   
   template<typename A, typename B>
   class BodyLogicalAND : public BodyLogicalOperationBinary<A,B>
   {
   public:
      explicit BodyLogicalAND( const shared_ptr<A> &a, const shared_ptr<B> &b ): BodyLogicalOperationBinary<A,B>(a,b) {};
   };
   template<typename A, typename B>
   class BodyLogicalOR : public BodyLogicalOperationBinary<A,B>
   {
   public:
      explicit BodyLogicalOR ( const shared_ptr<A> &a, const shared_ptr<B> &b ): BodyLogicalOperationBinary<A,B>(a,b) {};
   };
   template<typename A, typename B>
   class BodyLogicalXOR : public BodyLogicalOperationBinary<A,B>
   {
   public:
      explicit BodyLogicalXOR( const shared_ptr<A> &a, const shared_ptr<B> &b ): BodyLogicalOperationBinary<A,B>(a,b) {};
   };
   
   template<typename A>
   class BodyLogicalNOT
   {
   public:

      explicit BodyLogicalNOT( const shared_ptr<A> &a )
         : a_(a)
      {}
      
      const A& getA() const
      {
         return *a_;
      }

   private:
      const shared_ptr<A> a_;
   };


   // Body concept
   template<typename A>
   bool contains ( const BodyLogicalNOT<A> & body, const Vector3<real_t> & point )
   {
      return !contains( body.getA(), point);
   }
   template<typename A, typename B>
   bool contains ( const BodyLogicalAND<A,B> & body, const Vector3<real_t> & point )
   {
      return contains(body.getA(), point) && contains(body.getB(), point);
   }
   template<typename A, typename B>
   bool contains ( const BodyLogicalOR<A,B> & body, const Vector3<real_t> & point )
   {
      return contains(body.getA(), point) || contains(body.getB(), point);
   }
   template<typename A, typename B>
   bool contains ( const BodyLogicalXOR<A,B> & body, const Vector3<real_t> & point )
   {
      return contains(body.getA(), point) != contains(body.getB(), point);
   }



} // namespace geometry
} // namespace walberla
