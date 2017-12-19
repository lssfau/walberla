#pragma once

#include <pe/basic.h>

namespace walberla {
namespace pe {
namespace raytracing {
   class Ray {
   public:
      Vec3 origin; // coordinates of screen pixel
      Vec3 direction; // e.g. (1,0,0) for ray in x direction
      
      Ray (Vec3 origin_, Vec3 direction_) : origin(origin_), direction(direction_) {
         
      }
   };
}
}
}
