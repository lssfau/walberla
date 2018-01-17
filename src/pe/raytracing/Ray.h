#pragma once

#include "core/math/Vector3.h"

namespace walberla {
namespace pe {
namespace raytracing {
   class Ray {
   private:
      Vec3 direction_;
      Vec3 inv_direction_;
      int sign[3];
   public:
      Vec3 origin;

      Ray () {
         Ray (Vec3(0,0,0), Vec3(1,0,0));
      }
      
      Ray (Vec3 origin_, Vec3 direction) : direction_(direction), origin(origin_) {
         WALBERLA_ASSERT(realIsEqual(direction_.length(), real_t(1)));
         calcInvDirection();
      }
      
      void setDirection (Vec3 direction) {
         direction_ = direction;
         calcInvDirection();
      }
      
      Vec3& getDirection () {
         return direction_;
      }
      
      Vec3& getInvDirection () {
         return inv_direction_;
      }
      
      int* getInvDirectionSigns () {
         return sign;
      }

      void calcInvDirection () {
         inv_direction_ = Vec3(1/direction_[0], 1/direction_[1], 1/direction_[2]);
         sign[0] = (inv_direction_[0] < 0);
         sign[1] = (inv_direction_[1] < 0);
         sign[2] = (inv_direction_[2] < 0);
      }
   };
}
}
}
