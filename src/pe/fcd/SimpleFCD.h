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
//! \file IFCD.h
//! \ingroup FineCollisionDetection
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "IFCD.h"

#include "pe/collision/Collide.h"

#include <boost/type_traits/is_base_of.hpp>
#include <boost/tuple/tuple.hpp>

namespace walberla{
namespace pe{
namespace fcd {

template < typename BodyA, typename BodyTypeTuple, typename BaseT >
struct SingleDispatch{
   static void execute(RigidBody* a, RigidBody* b, Contacts& contacts){
      static_assert(boost::is_base_of<RigidBody, typename BodyTypeTuple::head_type>::value, "only downcasting allowed!");
      if (BodyTypeTuple::head_type::getStaticTypeID() == b->getTypeID())
      {
         typedef typename BodyTypeTuple::head_type BodyB;

         auto bd1 = static_cast<BodyA *>(a);
         auto bd2 = static_cast<BodyB *>(b);

         collide(bd1, bd2, contacts);
      } else
      {
         SingleDispatch<BodyA, typename BodyTypeTuple::tail_type, BaseT>::execute(a, b, contacts);
      }
   }
};

template < typename BodyA, typename BaseT>
struct SingleDispatch< BodyA, boost::tuples::null_type, BaseT>{
   static void execute(RigidBody* /*a*/, RigidBody* b, Contacts& /*contacts*/){
      WALBERLA_ABORT("SingleDispatch: Type of body " << b->getSystemID() << " could not be determined (" << b->getTypeID() << ")");
   }
};

template < typename BodyTypeTuple, typename BaseT = BodyTypeTuple>
struct DoubleDispatch{
   static void execute(RigidBody* a, RigidBody* b, Contacts& contacts){
      // Force a defined order of collision detection across processes
      if( b->getSystemID() < a->getSystemID() )
         std::swap( a, b );
      static_assert(boost::is_base_of<RigidBody, typename BodyTypeTuple::head_type>::value, "only downcasting allowed!");
      if (BodyTypeTuple::head_type::getStaticTypeID() == a->getTypeID()) {
         SingleDispatch<typename BodyTypeTuple::head_type, BaseT, BaseT>::execute(a, b, contacts);
      } else {
         DoubleDispatch<typename BodyTypeTuple::tail_type, BaseT>::execute(a, b, contacts);
      }
   }
};

template < typename BaseT>
struct DoubleDispatch< boost::tuples::null_type, BaseT>{
   static void execute(RigidBody* /*a*/, RigidBody* /*b*/, Contacts& /*contacts*/){
      WALBERLA_ABORT("DoubleDispatch: Type could not be determined");
   }
};

template <typename BodyTypeTuple>
class SimpleFCD : public IFCD{
public:
   virtual Contacts& generateContacts(PossibleContacts& possibleContacts)
   {
      contacts_.clear();
      for (auto it = possibleContacts.begin(); it != possibleContacts.end(); ++it)
      {
         DoubleDispatch<BodyTypeTuple>::execute(it->first, it->second, contacts_);
      }
      return contacts_;
   }
};

}
}
}
