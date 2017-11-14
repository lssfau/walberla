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
//! \file IterativeFCD.h
//! \author Tobias Leemann <tobias.leemann@fau.de>
//
//======================================================================================================================


#pragma once

#include "IFCD.h"
#include "pe/Types.h"
#include "pe/collision/EPA.h"
#include "pe/collision/GJK.h"
#include "pe/rigidbody/Plane.h"
#include "pe/rigidbody/Plane.h"
#include <pe/Thresholds.h>
#include <boost/tuple/tuple.hpp>

namespace walberla{
namespace pe{
namespace fcd {

/**
 * \ingroup pe
 * \brief Implementation of the IFCD Interface with an iterative detection.
 *
 * This Class implements a collision detection compatible with the
 * IFCD interface. It detects collisions using an iterative technique which employs
 * the GJK and EPA algorithms.
 */
template <typename BodyTypeTuple>
class IterativeFCD : public IFCD {
   typedef Union< boost::tuple<BodyTypeTuple>> UnionGenericType;

private:
   /*!\brief This computes a collision between two rigid bodies.
   *
   * \param a The first body (must not be a Plane or Union).
   * \param b The second body (must not be a Plane or Union).
   * \param normal Returns the normal vector of the collision.
   * \param contactPoint Returns a point of contact btw. the bodies.
   * \param penetrationDepth Returns the depth.
   * \return True, if the bodies collide.
   *
   * The Return parameters remain untouched if there is no collision.
   * This function performes collision detection with bodies enlarged by a margin for numerical reasons.
   * See "Collision detection in interactive 3D environments" by Gino van den Bergen
   * for a detailed explanation of the algorithm used.
   */
   bool performIterativeDetection(GeomPrimitive &a, GeomPrimitive &b, Vec3& normal, Vec3& contactPoint, real_t& penetrationDepth){
      real_t margin = real_t(1e-4);
      GJK gjk;
      if(gjk.doGJKcontactThreshold(a, b, margin)){
         //2. If collision is possible perform EPA.
         EPA epa;
         return epa.doEPAmargin(a, b, gjk, normal, contactPoint, penetrationDepth, margin);
      }else{
         return false;
      }
   }

   /*!\brief This computes a collision between a Plane and another body.
   *
   * \param pl The Plane.
   * \param b The second body (must not be a Plane or Union).
   * \param normal Returns the normal vector of the collision.
   * \param contactPoint Returns a point of contact btw. the bodies.
   * \param penetrationDepth Returns the depth.
   * \return True, if the Plane and the body collide.
   *
   * The Return parameters remain untouched if there is no collision.
   */
   bool performPlaneDetection(Plane &pl, GeomPrimitive &b, Vec3& normal, Vec3& contactPoint, real_t& penetrationDepth){
      Vec3 support_dir = -pl.getNormal();

      // We now have a direction facing to the "wall".
      // Compute support point of body b in this direction. This will be the furthest point overlapping.

      Vec3 contactp = b.support(support_dir);
      //std::cerr << contactp << std::endl;
      //std::cerr << pl.getDisplacement() <<std::endl;
      real_t pdepth = contactp * pl.getNormal() - pl.getDisplacement();
      if(pdepth < contactThreshold){ //We have a collision
         normal = support_dir;
         penetrationDepth = pdepth;
         contactPoint = contactp + real_t(0.5) * penetrationDepth * normal;
         return true;
      }else{ //No collision
         return false;
      }
   }

   /*!\brief This function adds contacts for single bodies (not unions)
   *
   * \param bd1 The first body.
   * \param bd2 The second body.
   *
   * This function is called by generateContactsPair() on all normal (single)
   * body pairs and on each single bodys of the union.
   * It checks if one of the bodies is a plane and calls the corresponding
   * detection function. If collision is detected a Contact is generated and
   * added to the vector.
   */
   void generateContactsSingle(BodyID bd1, BodyID bd2){
      //Temporary Storage for data
      Vec3 normal;
      Vec3 contactPoint;
      real_t penetrationDepth;

      //Check if one primitive is a plane
      if(bd1->getTypeID() == Plane::getStaticTypeID()){
         if(!(bd2->getTypeID() == Plane::getStaticTypeID())){
            //First object is a plane, second is not
            Plane* pl = static_cast<Plane*>(bd1);
            GeomPrimitive *gb = static_cast<GeomPrimitive*>(bd2);
            if (performPlaneDetection(*pl, *gb, normal, contactPoint, penetrationDepth)){
               contacts_.push_back( Contact(pl, gb, contactPoint, normal, penetrationDepth) );
            }
         }else{
            return; //Plane plane collisions cannot be handled here.
         }
      }else if(bd2->getTypeID() == Plane::getStaticTypeID()){
         if(!(bd1->getTypeID() == Plane::getStaticTypeID())){
            //Second object is a plane, first is not
            Plane* pl = static_cast<Plane*>(bd2);
            GeomPrimitive *ga = static_cast<GeomPrimitive*>(bd1);
            if (performPlaneDetection(*pl, *ga, normal, contactPoint, penetrationDepth)){
               contacts_.push_back( Contact(ga, pl, contactPoint, -normal, penetrationDepth) );
            }
         }
      }else{ // Both objects are geometries.
         GeomPrimitive *ga = static_cast<GeomPrimitive*>(bd1);
         GeomPrimitive *gb = static_cast<GeomPrimitive*>(bd2);

         if (performIterativeDetection(*ga, *gb, normal, contactPoint, penetrationDepth)){
            contacts_.push_back( Contact(ga, gb, contactPoint, normal, penetrationDepth) );
         }
      }
   }

   /* Resolve pairs with unions and call generateContacts single for each pair of non-united
    * particles.
    */
   void generateContactsPair(BodyID bd1, BodyID bd2){
      bool hasSubbodies1 = bd1->hasSubBodies();
      if( hasSubbodies1 ){
         UnionGenericType* u1 = nullptr;
         u1= static_cast<UnionGenericType*>(bd1);
         //Loop over possible subbodies of first rigid body
         auto it1 = u1->begin();
         while(it1 != u1->end()){
            generateContactsPair(*it1, bd2);
            it1++;
         }
      }else{
         //bd1 is a single body
         bool hasSubbodies2 = bd2->hasSubBodies();
         if( hasSubbodies2 ){
            UnionGenericType* u2 = nullptr;
            u2= static_cast<UnionGenericType*>(bd2);
            auto it2 = u2->begin();
            while(it2 != u2->end()){
               generateContactsPair(bd1, *it2);
               it2++;
            }
         }else{
            //bd1 and bd2 are single bodies
            generateContactsSingle(bd1, bd2);
         }
      }
   }

public:
   /*!\brief Checks each body pair in the given vector and generates
   * a vector of contacts.
   *
   * \param possibleContacts Vector of body pairs thay may collide.
   * \return Vector of contacts with information of all contacts generated.
   */
   virtual Contacts& generateContacts(PossibleContacts& possibleContacts)
   {
      contacts_.clear();
      for (auto it = possibleContacts.begin(); it != possibleContacts.end(); ++it)
      {
         generateContactsPair(it->first, it->second);
      }
      return contacts_;
   }

};


}
}
}
