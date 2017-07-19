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
//! \file BodyDispatcher.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <pe/rigidbody/RigidBody.h>

#include <boost/tuple/tuple.hpp>

namespace walberla {
namespace pe {

template < typename TypeList, typename Functor, typename ReturnType >
class SingleDispatch
{
public:
   static ReturnType execute(const id_t typeID, Functor& func){
      static_assert(boost::is_base_of<RigidBody, typename TypeList::head_type>::value, "only downcasting allowed!");
      if (TypeList::head_type::getStaticTypeID() == typeID)
      {
         typedef typename TypeList::head_type BodyType;
         BodyType* bd = NULL;
         return func( static_cast<BodyType *>( bd ) );
      } else
      {
         return SingleDispatch<typename TypeList::tail_type, Functor, ReturnType>::execute(typeID, func);
      }
   }

   static ReturnType execute(RigidBody* bd, Functor& func){
      static_assert(boost::is_base_of<RigidBody, typename TypeList::head_type>::value, "only downcasting allowed!");
      if (TypeList::head_type::getStaticTypeID() == bd->getTypeID())
      {
         typedef typename TypeList::head_type BodyType;
         return func( static_cast<BodyType *>(bd) );
      } else
      {
         return SingleDispatch<typename TypeList::tail_type, Functor, ReturnType>::execute(bd, func);
      }
   }

   static ReturnType execute(const RigidBody* bd, Functor& func){
      static_assert(boost::is_base_of<RigidBody, typename TypeList::head_type>::value, "only downcasting allowed!");
      if (TypeList::head_type::getStaticTypeID() == bd->getTypeID())
      {
         typedef typename TypeList::head_type BodyType;
         return func( static_cast<const BodyType *>(bd) );
      } else
      {
         return SingleDispatch<typename TypeList::tail_type, Functor, ReturnType>::execute(bd, func);
      }
   }
};

template < typename Functor, typename ReturnType >
struct SingleDispatch< boost::tuples::null_type, Functor, ReturnType >{
   static ReturnType execute(const id_t typeID, Functor& /*func*/)
   {
      WALBERLA_ABORT("SingleDispatch: BodyType could not be determined (" << typeID << ")");
      return ReturnType();
   }

   static ReturnType execute(RigidBody* bd, Functor& /*func*/)
   {
      WALBERLA_ABORT("SingleDispatch: Type of body " << bd->getSystemID() << " could not be determined (" << bd->getTypeID() << ")");
      return ReturnType();
   }
   static ReturnType execute(const RigidBody* bd, Functor& /*func*/)
   {
      WALBERLA_ABORT("SingleDispatch: Type of body " << bd->getSystemID() << " could not be determined (" << bd->getTypeID() << ")");
      return ReturnType();
   }
};

}  // namespace pe
}  // namespace walberla
