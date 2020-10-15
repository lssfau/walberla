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
//! \file GJKEPACollideFunctor.h
//! \author Tobias Leemann
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once


#include "pe/Types.h"
#include "pe/collision/EPA.h"
#include "pe/collision/GJK.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Union.h"
#include "pe/Thresholds.h"
#include "pe/utility/BodyCast.h"

namespace walberla{
namespace pe{
namespace fcd {
namespace gjkepa{

   //function for all single rigid bodies.
   template<typename Container>
   inline bool generateContacts(GeomPrimitive *a, GeomPrimitive *b, Container& contacts_);

   //Planes
   template<typename Container>
   inline bool generateContacts(Plane *a, GeomPrimitive *b, Container& contacts_);

   template<typename Container>
   inline bool generateContacts(GeomPrimitive *a, Plane *b, Container& contacts_);

   template< typename Container>
   inline bool generateContacts(Plane *a, Plane *b, Container& contacts_);

   //Unions
   template<typename... BodyTypesA, typename BodyB, typename Container>
   inline bool generateContacts(Union<BodyTypesA...> *a, BodyB *b, Container& contacts_);

   template<typename BodyA, typename... BodyTypesB, typename Container>
   inline bool generateContacts(BodyA *a, Union<BodyTypesB...> *b, Container& contacts_);

   template<typename... BodyTypesA, typename... BodyTypesB, typename Container>
   inline bool generateContacts(Union<BodyTypesA...> *a, Union<BodyTypesB...>  *b, Container& contacts_);

   //Union and Plane
   template<typename... BodyTypesA, typename Container>
   inline bool generateContacts(Union<BodyTypesA...> *a, Plane *b, Container& contacts_);

   template<typename... BodyTypesB, typename Container>
   inline bool generateContacts(Plane *a, Union<BodyTypesB...> *b, Container& contacts_);
}

/* Iterative Collide Functor for contact Generation with iterative collision detection (GJK and EPA algorithms).
 * Usage: fcd::GenericFCD<BodyTypes..., fcd::GJKEPACollideFunctor> testFCD;
 * testFCD.generateContacts(...);
 */
template <typename Container>
struct GJKEPACollideFunctor
{
   Container& contacts_;

   GJKEPACollideFunctor(Container& contacts) : contacts_(contacts) {}

   template< typename BodyType1, typename BodyType2 >
   bool operator()( BodyType1* bd1, BodyType2* bd2) {
      using namespace gjkepa;
      return generateContacts(bd1, bd2, contacts_);
   }
};

template <typename BodyType1, typename Container>
struct GJKEPASingleCollideFunctor
{
   BodyType1* bd1_;
   Container& contacts_;

   GJKEPASingleCollideFunctor(BodyType1* bd1, Container& contacts) : bd1_(bd1), contacts_(contacts) {}

   template< typename BodyType2 >
   bool operator()( BodyType2* bd2) {
      using namespace gjkepa;
      return generateContacts( bd1_, bd2, contacts_);
   }
};


namespace gjkepa{

   //function for all single rigid bodies.
   template<typename Container>
   inline bool generateContacts(GeomPrimitive *a, GeomPrimitive *b, Container& contacts_){
      Vec3 normal;
      Vec3 contactPoint;
      real_t penetrationDepth;

      real_t margin = real_t(1e-4);
      GJK gjk;
      if(gjk.doGJKmargin(*a, *b, margin)){
         //2. If collision is possible perform EPA.
         EPA epa;
         epa.useSphereOptimization(true);
         if(epa.doEPAmargin(*a, *b, gjk, normal, contactPoint, penetrationDepth, margin)){
            contacts_.push_back( Contact(a, b, contactPoint, normal, penetrationDepth) );
            return true;
         }else{
            return false;
         }

      }else{
         return false;
      }
   }

   //Planes
   template<typename Container>
   inline bool generateContacts(Plane *a, GeomPrimitive *b, Container& contacts_){
      Vec3 normal;
      Vec3 contactPoint;
      real_t penetrationDepth;

      Vec3 support_dir = -a->getNormal();
      // We now have a direction facing to the "wall".
      // Compute support point of body b in this direction. This will be the furthest point overlapping.
      Vec3 contactp = b->support(support_dir);
      real_t pdepth = contactp * a->getNormal() - a->getDisplacement();
      if(pdepth < contactThreshold){ //We have a collision
         normal = support_dir;
         penetrationDepth = pdepth;
         contactPoint = contactp + real_t(0.5) * penetrationDepth * normal;
         contacts_.push_back( Contact(a, b, contactPoint, normal, penetrationDepth) );
         return true;
      }else{ //No collision
         return false;
      }
   }

   template<typename Container>
   inline bool generateContacts(GeomPrimitive *a, Plane *b, Container& contacts_){
      return generateContacts(b, a, contacts_);
   }

   //Planes cannot collide with each other
   template< typename Container>
   inline bool generateContacts(Plane*, Plane*, Container&){
      return false;
   }

   //Unions
   template<typename... BodyTypesA, typename BodyB, typename Container>
   inline bool generateContacts(Union<BodyTypesA...> *a, BodyB *b, Container& contacts_){
      GJKEPASingleCollideFunctor<BodyB, Container> func(b, contacts_);
      bool collision = false;
      for( auto it=a->begin(); it!=a->end(); ++it )
      {
         collision |= SingleCast<std::tuple<BodyTypesA...>, GJKEPASingleCollideFunctor<BodyB, Container>, bool>::execute(it.getBodyID(), func);
      }
      return collision;
   }

   template<typename BodyA, typename... BodyTypesB, typename Container>
   inline bool generateContacts(BodyA *a, Union<BodyTypesB...> *b, Container& contacts_){
      return generateContacts(b, a, contacts_);
   }

   template<typename... BodyTypesA, typename... BodyTypesB, typename Container>
   inline bool generateContacts(Union<BodyTypesA...> *a, Union<BodyTypesB...>  *b, Container& contacts_){
      GJKEPACollideFunctor<Container> func(contacts_);
      bool collision = false;
      for( auto it1=a->begin(); it1!=a->end(); ++it1 )
      {
         for( auto it2=b->begin(); it2!=b->end(); ++it2 )
         {
            collision |= DoubleCast<std::tuple<BodyTypesA...>, std::tuple<BodyTypesB...>, GJKEPACollideFunctor<Container>, bool>::execute(it1.getBodyID(), it2.getBodyID(), func);
         }
      }
      return collision;
   }

   //Union and Plane (these calls are ambigous if not implemented seperatly)
   template<typename... BodyTypesA, typename Container>
   inline bool generateContacts(Union<BodyTypesA...> *a, Plane *b, Container& contacts_){
      GJKEPASingleCollideFunctor<Plane, Container> func(b, contacts_);
      bool collision = false;
      for( auto it=a->begin(); it!=a->end(); ++it )
      {
         collision |= SingleCast<std::tuple<BodyTypesA...>, GJKEPASingleCollideFunctor<Plane, Container>, bool>::execute(it.getBodyID(), func);
      }
      return collision;
   }

   template<typename... BodyTypesB, typename Container>
   inline bool generateContacts(Plane *a, Union<BodyTypesB...> *b, Container& contacts_){
      return generateContacts(b, a, contacts_);
   }


} //namespace gjkepa


} //fcd
} //pe
} //walberla
