#pragma once

#include "core/math/Vector3.h"

namespace walberla {
namespace pe {
namespace raytracing {
   class Ray {
   public:
      Vec3 origin;
      Vec3 direction; // e.g. (1,0,0) for ray in x direction
      Vec3 inv_direction;
      int sign[3];
      
      Ray (Vec3 origin_, Vec3 direction_) : origin(origin_), direction(direction_) {
         inv_direction = Vec3(1/direction[0], 1/direction[1], 1/direction[2]);
         sign[0] = (inv_direction[0] < 0);
         sign[1] = (inv_direction[1] < 0);
         sign[2] = (inv_direction[2] < 0);
      }
   };
}
}
}
