#pragma once
#include "core/DataTypes.h"
#include "core/math/Vector3.h"

namespace walberla {
namespace pe {
namespace raytracing {
   class Ray {
   private:
      Vec3 direction_;
      Vec3 inv_direction_;
      Vector3<int8_t> sign_;
      Vec3 origin_;
      
   public:
      Ray () {
         Ray (Vec3(0,0,0), Vec3(1,0,0));
      }
      
      Ray (Vec3 origin, Vec3 direction) {
         setDirection(direction);
         setOrigin(origin);
      }
      
      inline void setOrigin (const Vec3& origin) {
         origin_ = origin;
      }
      
      inline const Vec3& getOrigin () const {
         return origin_;
      }
      
      inline void setDirection (const Vec3& direction) {
         // im kommentar verweis auf normalisierung
         WALBERLA_CHECK_FLOAT_EQUAL(direction.length(), real_t(1));
         direction_ = direction;
         calcInvDirection();
      }
      
      inline const Vec3& getDirection () const {
         return direction_;
      }
      
      inline const Vec3& getInvDirection () const {
         return inv_direction_;
      }
      
      inline const Vector3<int8_t>& getInvDirectionSigns () const {
         return sign_;
      }

      inline void calcInvDirection () {
         inv_direction_ = Vec3(1/direction_[0], 1/direction_[1], 1/direction_[2]);
         sign_[0] = int8_c(inv_direction_[0] < 0);
         sign_[1] = int8_c(inv_direction_[1] < 0);
         sign_[2] = int8_c(inv_direction_[2] < 0);
      }
   };
}
}
}
