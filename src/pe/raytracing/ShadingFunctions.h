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
//! \file ShadingFunctions.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <pe/raytracing/Color.h>
#include <pe/raytracing/ShadingParameters.h>
#include <core/mpi/MPIWrapper.h>
#include <core/mpi/MPIManager.h>

namespace walberla {
namespace pe {
namespace raytracing {

inline ShadingParameters defaultBodyTypeDependentShadingParams (const BodyID body);
inline ShadingParameters processRankDependentShadingParams (const BodyID body);
inline ShadingParameters defaultShadingParams (const BodyID body);
inline ShadingParameters blackShadingParams (const BodyID body);
inline ShadingParameters whiteShadingParams (const BodyID body);
inline ShadingParameters lightGreyShadingParams (const BodyID body);
inline ShadingParameters greyShadingParams (const BodyID body);
inline ShadingParameters darkGreyShadingParams (const BodyID body);
inline ShadingParameters redShadingParams (const BodyID body);
inline ShadingParameters greenShadingParams (const BodyID body);
inline ShadingParameters blueShadingParams (const BodyID body);
inline ShadingParameters violetShadingParams (const BodyID body);
inline ShadingParameters yellowShadingParams (const BodyID body);

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
   } else if (bodyTypeID == Ellipsoid::getStaticTypeID()) {
      return yellowShadingParams(body).makeGlossy(60);
   } else {
      return defaultShadingParams(body);
   }
}

inline ShadingParameters processRankDependentShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   int numProcesses = mpi::MPIManager::instance()->numProcesses();
   int rank = mpi::MPIManager::instance()->rank();
   
   real_t hue = real_t(360) * real_t(rank)/real_t(numProcesses);
   Color color = Color::colorFromHSV(hue, real_t(1), real_t(0.9));
   
   return ShadingParameters(color,
                            color*real_t(0.5),
                            Color(0,0,0),
                            real_t(0));
}
   
inline ShadingParameters defaultShadingParams (const BodyID body) {
   return greyShadingParams(body);
}
   
inline ShadingParameters whiteShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(1), real_t(1), real_t(1)),
                       Color(real_t(0.9), real_t(0.9), real_t(0.9)),
                       Color(real_t(0), real_t(0), real_t(0)),
                       real_t(0));
   return s;
}
   
inline ShadingParameters blackShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0), real_t(0), real_t(0)),
                       Color(real_t(0), real_t(0), real_t(0)),
                       Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                       real_t(0));
   return s;
}

inline ShadingParameters lightGreyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0.82), real_t(0.82), real_t(0.82)),
                       Color(real_t(0.5), real_t(0.5), real_t(0.5)),
                       Color(real_t(0), real_t(0), real_t(0)),
                       real_t(0));
   return s;
}

inline ShadingParameters greyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0.5), real_t(0.5), real_t(0.5)),
                       Color(real_t(0.4), real_t(0.4), real_t(0.4)),
                       Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                       real_t(0));
   return s;
}

inline ShadingParameters darkGreyShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0.2), real_t(0.2), real_t(0.2)),
                       Color(real_t(0.06), real_t(0.06), real_t(0.06)),
                       Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                       real_t(0));
   return s;
}
   
inline ShadingParameters redShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(1), real_t(0), real_t(0)),
                       Color(real_t(0.5), real_t(0), real_t(0)),
                       Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                       real_t(0));
   return s;
}

inline ShadingParameters greenShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0), real_t(0.72), real_t(0)),
                       Color(real_t(0), real_t(0.41), real_t(0)),
                       Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                       real_t(0));
   return s;
}

inline ShadingParameters blueShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0.15), real_t(0.44), real_t(0.91)),
                       Color(real_t(0), real_t(0), real_t(0.4)),
                       Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                       real_t(0));
   return s;
}
   
inline ShadingParameters yellowShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(1), real_t(0.96), real_t(0)),
                       Color(real_t(0.5), real_t(0.48), real_t(0)),
                       Color(real_t(0), real_t(0), real_t(0)),
                       real_t(0));
   return s;
}

inline ShadingParameters violetShadingParams (const BodyID body) {
   WALBERLA_UNUSED(body);
   ShadingParameters s(Color(real_t(0.6), real_t(0), real_t(0.9)),
                       Color(real_t(0.5), real_t(0), real_t(0.8)),
                       Color(real_t(0), real_t(0), real_t(0)),
                       real_t(0));
   return s;
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
