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
//! \file BodyCast.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <pe/rigidbody/RigidBody.h>

#include <tuple>
#include <type_traits>

namespace walberla {
namespace pe {

template < typename TypeList, typename Functor, typename ReturnType, int N = std::tuple_size<TypeList>::value - 1 >
class SingleCast
{
public:
   static ReturnType execute(const id_t typeID, Functor& func)
   {
      using CastBodyType = typename std::tuple_element<N, TypeList>::type;
      static_assert(std::is_base_of<RigidBody, CastBodyType>::value, "only downcasting allowed!");
      if (CastBodyType::getStaticTypeID() == typeID)
      {
         CastBodyType* bd = NULL;
         return func( static_cast<CastBodyType *>( bd ) );
      } else
      {
         return SingleCast<TypeList, Functor, ReturnType, N - 1>::execute(typeID, func);
      }
   }

   static ReturnType execute(RigidBody* bd, Functor& func){
      using CastBodyType = typename std::tuple_element<N, TypeList>::type;
      static_assert(std::is_base_of<RigidBody, CastBodyType>::value, "only downcasting allowed!");
      if (CastBodyType::getStaticTypeID() == bd->getTypeID())
      {
         return func( static_cast<CastBodyType *>(bd) );
      } else
      {
         return SingleCast<TypeList, Functor, ReturnType, N - 1>::execute(bd, func);
      }
   }

   static ReturnType execute(const RigidBody* bd, Functor& func){
      using CastBodyType = typename std::tuple_element<N, TypeList>::type;
      static_assert(std::is_base_of<RigidBody, CastBodyType>::value, "only downcasting allowed!");
      if (CastBodyType::getStaticTypeID() == bd->getTypeID())
      {
         return func( static_cast<const CastBodyType *>(bd) );
      } else
      {
         return SingleCast<TypeList, Functor, ReturnType, N - 1>::execute(bd, func);
      }
   }
};

template < typename TypeList, typename Functor, typename ReturnType >
struct SingleCast< TypeList, Functor, ReturnType, -1 >{
   static ReturnType execute(const id_t typeID, Functor& /*func*/)
   {
      WALBERLA_ABORT("SingleCast: BodyType could not be determined (" << typeID << ")");
      return ReturnType();
   }

   static ReturnType execute(RigidBody* bd, Functor& /*func*/)
   {
      WALBERLA_ABORT("SingleCast: Type of body " << bd->getSystemID() << " could not be determined (" << bd->getTypeID() << ")");
      return ReturnType();
   }
   static ReturnType execute(const RigidBody* bd, Functor& /*func*/)
   {
      WALBERLA_ABORT("SingleCast: Type of body " << bd->getSystemID() << " could not be determined (" << bd->getTypeID() << ")");
      return ReturnType();
   }
};

template < typename TypeListA,
           typename TypeListB,
           typename Functor,
           typename ReturnType,
           int N = std::tuple_size<TypeListA>::value - 1 >
class DoubleCast
{
private:
   template< typename BodyType1 >
   struct SingleCastFunctor
   {
      BodyType1* a_;
      Functor& func_;

      SingleCastFunctor(BodyType1* a, Functor& func) : a_(a), func_(func) {}

      template< typename BodyType2 >
      ReturnType operator()( BodyType2* bd) { return func_( a_, bd); }
   };
   template< typename BodyType1 >
   struct SingleCastConstFunctor
   {
      BodyType1 const * a_;
      Functor& func_;

      SingleCastConstFunctor(BodyType1 const * a, Functor& func) : a_(a), func_(func) {}

      template< typename BodyType2 >
      ReturnType operator()( BodyType2 const * bd) { return func_( a_, bd); }
   };
public:
   static ReturnType execute(RigidBody* a, RigidBody* b, Functor& func){
      using CastBodyType = typename std::tuple_element<N, TypeListA>::type;
      static_assert(std::is_base_of<RigidBody, CastBodyType>::value, "only downcasting allowed!");
      if (CastBodyType::getStaticTypeID() == a->getTypeID())
      {
         SingleCastFunctor<CastBodyType> singleFunc( static_cast<CastBodyType *>(a), func);
         return SingleCast<TypeListB, SingleCastFunctor<CastBodyType>, ReturnType>::execute(b, singleFunc );
      } else
      {
         return DoubleCast<TypeListA, TypeListB, Functor, ReturnType, N - 1>::execute(a, b, func);
      }
   }

   static ReturnType execute(const RigidBody* a, const RigidBody* b, Functor& func){
      using CastBodyType = typename std::tuple_element<N, TypeListA>::type;
      static_assert(std::is_base_of<RigidBody, CastBodyType>::value, "only downcasting allowed!");
      if (CastBodyType::getStaticTypeID() == a->getTypeID())
      {
         SingleCastConstFunctor<CastBodyType> singleFunc( static_cast<CastBodyType *>(a), func);
         return SingleCast<TypeListB, SingleCastConstFunctor<CastBodyType>, ReturnType>::execute(b, singleFunc );
      } else
      {
         return DoubleCast<TypeListA, TypeListB, Functor, ReturnType, N - 1>::execute(a, b, func);
      }
   }
};

template < typename TypeListA, typename TypeListB, typename Functor, typename ReturnType >
struct DoubleCast< TypeListA, TypeListB, Functor, ReturnType, -1 >{
   static ReturnType execute(RigidBody* a, RigidBody* /*b*/, Functor& /*func*/)
   {
      WALBERLA_ABORT("DoubleCast: Type of body " << a->getSystemID() << " could not be determined (" << a->getTypeID() << ")");
      return ReturnType();
   }
   static ReturnType execute(const RigidBody* a, const RigidBody* /*b*/, Functor& /*func*/)
   {
      WALBERLA_ABORT("DoubleCast: Type of body " << a->getSystemID() << " could not be determined (" << a->getTypeID() << ")");
      return ReturnType();
   }
};

}  // namespace pe
}  // namespace walberla
