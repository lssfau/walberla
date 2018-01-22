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
//! \file Raytracer.cpp
//! \author Lukas Werner
//
//======================================================================================================================

#include "Raytracer.h"
#include "Ray.h"
#include "Intersects.h"

using namespace walberla;

struct BodyIntersectionInfo {
   size_t imageX; // viewing plane pixel coordinates to ...
   size_t imageY; // ... identify ray by pixel it intersected
   walberla::id_t bodySystemID; // body which was hit
   real_t t; // distance from camera to intersection point on body
};

struct Coordinates {
   size_t x;
   size_t y;
};

struct CoordinatesComparator {
   bool operator() (const Coordinates& lhs, const Coordinates& rhs) const {
      // standard lexicographical ordering
      return (lhs.x < rhs.x) || (lhs.x == rhs.x && lhs.y < rhs.y);
   }
};

real_t deg2rad(real_t deg) {
   return deg * math::M_PI / real_t(180.0);
}

namespace walberla {
namespace pe {
namespace raytracing {
/*!\brief Instantiation constructor for the Box class.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Storage ID of the block data storage the raytracer operates on.
 * \param pixelsHorizontal Horizontal amount of pixels of the generated image.
 * \param pixelsVertical Vertical amount of pixels of the generated image.
 * \param fov_vertical Vertical field-of-view of the camera.
 * \param cameraPosition Position of the camera in the global world frame.
 * \param lookAtPoint Point the camera looks at in the global world frame.
 * \param upVector Vector indicating the upwards direction of the camera.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage>& forest, BlockDataID storageID,
                     uint8_t pixelsHorizontal, uint8_t pixelsVertical,
                     real_t fov_vertical,
                     const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector)
   : forest_(forest), storageID_(storageID),
   pixelsHorizontal_(pixelsHorizontal), pixelsVertical_(pixelsVertical),
   fov_vertical_(fov_vertical),
   cameraPosition_(cameraPosition), lookAtPoint_(lookAtPoint), upVector_(upVector)
{
   // eye coordinate system setup
   n = (cameraPosition_ - lookAtPoint_).getNormalized();
   u = (upVector_ % n).getNormalized();
   v = n % u;
   
   // viewing plane setup
   d = (cameraPosition_ - lookAtPoint_).length();
   aspectRatio = real_t(pixelsHorizontal_) / real_t(pixelsVertical_);
   viewingPlaneHeight = tan(deg2rad(fov_vertical_)/real_t(2.)) * real_t(2.) * d;
   viewingPlaneWidth = viewingPlaneHeight * aspectRatio;
   viewingPlaneOrigin = lookAtPoint_ - u*viewingPlaneWidth/real_t(2.) - v*viewingPlaneHeight/real_t(2.);
   
   pixelWidth = viewingPlaneWidth / real_c(pixelsHorizontal_);
   pixelHeight = viewingPlaneHeight / real_c(pixelsVertical);
}

/*!\brief Does one raytracing step.
 *
 * \param timestep The timestep after which the raytracing starts.
 */
template <typename BodyTypeTuple>
void Raytracer::rayTrace(const size_t timestep) const {
   std::map<Coordinates, real_t, CoordinatesComparator> tBuffer; // t values for each pixel
   std::map<Coordinates, walberla::id_t, CoordinatesComparator> idBuffer; // ids of the intersected body for each pixel
   std::vector<BodyIntersectionInfo> intersections; // contains for each pixel information about an intersection, if existent
   
   std::map<Coordinates, BodyIntersectionInfo, CoordinatesComparator> localPixelIntersectionMap;
   
   real_t t, t_closest;
   walberla::id_t id_closest;
   RigidBody* body_closest = NULL;
   Ray ray(cameraPosition_, Vec3(1,0,0));
   IntersectsFunctor func(ray, t);
   for (size_t x = 0; x < pixelsHorizontal_; x++) {
      for (size_t y = 0; y < pixelsVertical_; y++) {
         //WALBERLA_LOG_INFO(x << "/" << y);
         Vec3 pixelLocation = viewingPlaneOrigin + u*(real_c(x)+real_t(0.5))*pixelWidth + v*(real_c(y)+real_t(0.5))*pixelHeight;
         Vec3 direction = (pixelLocation - cameraPosition_).getNormalized();
         ray.setDirection(direction);
         
         t_closest = INFINITY;
         id_closest = 0;
         body_closest = NULL;
         for (auto blockIt = forest_->begin(); blockIt != forest_->end(); ++blockIt) {
            // blockIt->getAABB();
            /*const AABB& blockAabb = blockIt->getAABB();
             if (!intersects(blockAabb, ray, t)) {
             continue;
             }*/
            for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID_); bodyIt != LocalBodyIterator::end(); ++bodyIt) {
               bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(*bodyIt, func);
               
               if (intersects && t < t_closest) {
                  // body was shot by ray and currently closest to camera
                  t_closest = t;
                  id_closest = bodyIt->getID();
                  body_closest = *bodyIt;
               }
            }
         }
         
         //std::cout << (t_closest != INFINITY ? size_t(t_closest) : 0) << " ";
         
         Coordinates c = {
            x,
            y
         };
         
         if (!realIsIdentical(t_closest, INFINITY) && body_closest != NULL) {
            BodyIntersectionInfo intersectionInfo = {
               x,
               y,
               body_closest->getSystemID(),
               t_closest
            };
            intersections.push_back(intersectionInfo);
            localPixelIntersectionMap[c] = intersectionInfo;
         }
         
         tBuffer[c] = t_closest;
         idBuffer[c] = id_closest;
      }
      //std::cout << std::endl;
   }
   
   // intersections synchronisieren
   mpi::SendBuffer sendBuffer;
   for (auto& info: intersections) {
      sendBuffer << info.imageX << info.imageY << info.bodySystemID << info.t;
   }
   
   std::vector<BodyIntersectionInfo> gatheredIntersections;
   
   std::map<walberla::id_t, bool> visibleBodyIDs;
   
   //std::map<Coordinates, BodyIntersectionInfo, CoordinatesComparator> pixelIntersectionMap;
   
   mpi::RecvBuffer recvBuffer;
   mpi::allGathervBuffer(sendBuffer, recvBuffer);
   while (!recvBuffer.isEmpty()) {
      BodyIntersectionInfo info;
      
      recvBuffer >> info.imageX;
      recvBuffer >> info.imageY;
      recvBuffer >> info.bodySystemID;
      recvBuffer >> info.t;
      
      Coordinates c = {
         info.imageX,
         info.imageY
      };
      
      /*if (pixelIntersectionMap.find(c) == pixelIntersectionMap.end()) {
       // map didnt contain coordinates
       pixelIntersectionMap.insert(std::make_pair(c, info));
       } else {
       // map already contains info at coordinates, check if current info is closer
       BodyIntersectionInfo& existingInfo = pixelIntersectionMap.at(c);
       if (existingInfo.t < info.t) {
       pixelIntersectionMap[c] = info;
       }
       }*/
      auto it = localPixelIntersectionMap.find(c);
      if (it != localPixelIntersectionMap.end()) {
         // there was a local hit at coordinate c
         BodyIntersectionInfo& localInfo = localPixelIntersectionMap.at(c);
         if (localInfo.t < info.t) {
            localPixelIntersectionMap.erase(it);
         }
      }
      
      //gatheredIntersections.push_back(info);
   }
   
   for (auto& info: localPixelIntersectionMap) {
      visibleBodyIDs[info.second.bodySystemID] = true;
   }
   
   WALBERLA_LOG_INFO("#particles visible: " << visibleBodyIDs.size());
   //WALBERLA_LOG_INFO_ON_ROOT("#gatheredIntersections: " << gatheredIntersections.size());
   
#ifdef __APPLE__
   //mpi::MPIRank rank = mpi::MPIManager::instance()->rank();
   //writeTBufferToFile(tBuffer, idBuffer, pixelsHorizontal_, pixelsVertical_, "/Users/ng/Desktop/walberla/tbuffer_" + std::to_string(rank) + ".ppm");
#endif
}

}
}
}
