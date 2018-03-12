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
//! \file Shading.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <pe/raytracing/Color.h>
#include <pe/raytracing/ShadingParameters.h>


namespace walberla {
namespace pe {
namespace raytracing {
inline ShadingParameters defaultBodyTypeDependentShadingParams (const BodyID body);
inline ShadingParameters defaultShadingParams (const BodyID body);
inline ShadingParameters blackShadingParams (const BodyID body);
inline ShadingParameters whiteShadingParams (const BodyID body);
inline ShadingParameters lightGreyShadingParams (const BodyID body);
inline ShadingParameters greyShadingParams (const BodyID body);
inline ShadingParameters darkGreyShadingParams (const BodyID body);
inline ShadingParameters redShadingParams (const BodyID body);
inline ShadingParameters blueShadingParams (const BodyID body);
inline ShadingParameters violetShadingParams (const BodyID body);

inline ShadingParameters defaultBodyTypeDependentShadingParams (const BodyID body) {
   auto bodyTypeID = body->getTypeID();
   
   if (bodyTypeID == Plane::getStaticTypeID()) {
      return lightGreyShadingParams(body).makeMatte();
   } else if (bodyTypeID == Sphere::getStaticTypeID()) {
      return redShadingParams(body).makeGlossy();
   } else if (bodyTypeID == Capsule::getStaticTypeID()) {
      return blueShadingParams(body).makeGlossy();
   } else if (bodyTypeID == Box::getStaticTypeID()) {
      return violetShadingParams(body);
   } else {
      return defaultShadingParams(body);
   }
}

inline ShadingParameters defaultShadingParams (const BodyID body) {
   return greyShadingParams(body);
}
   
inline ShadingParameters whiteShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(1, 1, 1),
                       Color(0.9, 0.9, 0.9),
                       Color(0, 0, 0),
                       0);
   return s;
}
   
inline ShadingParameters blackShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0, 0, 0),
                       Color(0, 0, 0),
                       Color(0.1, 0.1, 0.1),
                       0);
   return s;
}

inline ShadingParameters lightGreyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.82, 0.82, 0.82),
                       Color(0.5, 0.5, 0.5),
                       Color(0, 0, 0),
                       0);
   return s;
}

inline ShadingParameters greyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.5, 0.5, 0.5),
                       Color(0.4, 0.4, 0.4),
                       Color(0.1, 0.1, 0.1),
                       0);
   return s;
}

inline ShadingParameters darkGreyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.2, 0.2, 0.2),
                       Color(0.06, 0.06, 0.06),
                       Color(0.1, 0.1, 0.1),
                       0);
   return s;
}
   
inline ShadingParameters redShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(1, 0, 0),
             Color(0.5, 0, 0),
             Color(0.1, 0.1, 0.1),
             0);
   return s;
}

inline ShadingParameters greenShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0, 0.72, 0),
                       Color(0, 0.41, 0),
                       Color(0.1, 0.1, 0.1),
                       0);
   return s;
}

inline ShadingParameters blueShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.15, 0.44, 0.91),
             Color(0, 0, 0.4),
             Color(0.1, 0.1, 0.1),
             0);
   return s;
}
   
inline ShadingParameters yellowShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(1, 0.96, 0),
                       Color(0.5, 0.48, 0),
                       Color(0, 0, 0),
                       0);
   return s;
}

inline ShadingParameters violetShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(0.6, 0, 0.9),
                       Color(0.5, 0, 0.8),
                       Color(0, 0, 0),
                       0);
   return s;
}
}
}
}
