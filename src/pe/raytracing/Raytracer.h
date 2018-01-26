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
//! \file Raytracer.h
//! \author Lukas Werner
//
//======================================================================================================================

#pragma once

#include <pe/basic.h>
#include <pe/Types.h>
#include <core/math/Vector3.h>
#include <core/mpi/all.h>
#include <core/config/Config.h>
#include <boost/filesystem.hpp>
#include <core/timing/TimingTree.h>
#include "Ray.h"
#include "Intersects.h"

using namespace walberla;
using namespace walberla::pe;
using namespace walberla::timing;

namespace walberla {
namespace pe {
namespace raytracing {

/*!\brief Contains information about a ray-body intersection.
 */
struct BodyIntersectionInfo {
   size_t imageX;                //!< viewing plane x pixel coordinate the ray belongs to.
   size_t imageY;                //!< viewing plane y pixel coordinate the ray belongs to.
   walberla::id_t bodySystemID;  //!< system ID of body which was hit.
   real_t t;                     //!< distance from camera to intersection point on body.
};

class Raytracer {
public:
   /*!\name Constructors */
   //@{
   explicit Raytracer(const shared_ptr<BlockStorage> forest, BlockDataID storageID,
                      const shared_ptr<BodyStorage> globalBodyStorage,
                      size_t pixelsHorizontal, size_t pixelsVertical,
                      real_t fov_vertical,
                      const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                      real_t blockAABBIntersectionPadding = real_t(0.0));
   explicit Raytracer(const shared_ptr<BlockStorage> forest, BlockDataID storageID,
                      const shared_ptr<BodyStorage> globalBodyStorage,
                      const Config::BlockHandle& config);
   //@}

private:
   /*!\name Member variables */
   //@{
   const shared_ptr<BlockStorage> forest_; //!< The BlockForest the raytracer operates on.
   BlockDataID storageID_;    /*!< The storage ID of the block data storage the raytracer operates
                               on.*/
   const shared_ptr<BodyStorage> globalBodyStorage_; //!< The global body storage the raytracer operates on.
   
   size_t pixelsHorizontal_;  //!< The horizontal amount of pixels of the generated image.
   size_t pixelsVertical_;    //!< The vertical amount of pixels of the generated image.
   real_t fov_vertical_;      //!< The vertical field-of-view of the camera.
   Vec3 cameraPosition_;      //!< The position of the camera in the global world frame.
   Vec3 lookAtPoint_;         /*!< The point the camera looks at in the global world frame,
                               marks the center of the view plane.*/
   Vec3 upVector_;            //!< The vector indicating the upwards direction of the camera.
   real_t blockAABBIntersectionPadding_; /*!< The padding applied in block AABB intersection pretesting, as
                                          some objects within a block might protrude from the block's AABB.*/
   
   bool tBufferOutputEnabled_; //!< Enable / disable dumping the tbuffer to a file
   std::string tBufferOutputDirectory_; //!< Path to the tbuffer output directory
   //@}
   
   Vec3 n_;                   //!< The normal vector of viewing plane.
   Vec3 u_;                   //!< The vector spanning the viewing plane in the "right direction".
   Vec3 v_;                   //!< The vector spanning the viewing plane in the "up direction".
   real_t d_;                 //!< The the distance from camera to viewing plane.
   real_t aspectRatio_;       //!< The aspect ratio of the generated image and viewing plane.
   real_t viewingPlaneHeight_; //!< The viewing plane height in the world frame.
   real_t viewingPlaneWidth_; //!< The viewing plane width in the world frame.
   Vec3 viewingPlaneOrigin_;  //!< The origin of the viewing plane.
   real_t pixelWidth_;        //!< The width of a pixel of the generated image in the viewing plane.
   real_t pixelHeight_;       //!< The height of a pixel of the generated image in the viewing plane.
   WcTimingPool tp_;

public:
   /*!\name Get functions */
   //@{
   inline size_t getPixelsHorizontal() const;
   inline size_t getPixelsVertical() const;
   inline real_t getFOVVertical() const;
   inline const Vec3& getCameraPosition() const;
   inline const Vec3& getLookAtPoint() const;
   inline const Vec3& getUpVector() const;
   inline bool getTBufferOutputEnabled() const;
   inline const std::string& getTBufferOutputDirectory() const;
   //@}

   /*!\name Set functions */
   //@{
   inline void setTBufferOutputEnabled(const bool enabled);
   inline void setTBufferOutputDirectory(const std::string& path);
   //@}
   
   /*!\name Functions */
   //@{
   void setupView_();
   template <typename BodyTypeTuple>
   void rayTrace(const size_t timestep);
   
private:
   void writeTBufferToFile(const std::vector<real_t>& tBuffer, size_t timestep, bool isGlobalImage = false) const;
   void writeTBufferToFile(const std::vector<real_t>& tBuffer, const std::string& fileName) const;
   bool isPlaneVisible(const PlaneID plane, const Ray& ray) const;
   size_t coordinateToArrayIndex(size_t x, size_t y) const;
   //@}
};
   
/*!\brief Returns the horizontal amount of pixels of the generated image.
 *
 * \return The horizontal amount of pixels of the generated image.
 */
inline size_t Raytracer::getPixelsHorizontal() const {
   return pixelsHorizontal_;
}

/*!\brief Returns the vertical amount of pixels of the generated image.
 *
 * \return The vertical amount of pixels of the generated image.
 */
inline size_t Raytracer::getPixelsVertical() const {
   return pixelsVertical_;
}

/*!\brief Returns the vertical field-of-view of the camera.
 *
 * \return The vertical field-of-view of the camera.
 */
inline real_t Raytracer::getFOVVertical() const {
   return fov_vertical_;
}

/*!\brief Returns the position of the camera in the global world frame.
 *
 * \return The position of the camera.
 *
 * Returns the position of the camera in the global world frame.
 */
inline const Vec3& Raytracer::getCameraPosition() const {
   return cameraPosition_;
}

/*!\brief Returns the point the camera looks at in the global world frame.
 *
 * \return The looked at point.
 *
 * Returns the point the camera looks at in the global world frame, which also marks the center of
 * the view plane.
 */
inline const Vec3& Raytracer::getLookAtPoint() const {
   return lookAtPoint_;
}

/*!\brief Returns the vector indicating the upwards direction of the camera.
 *
 * \return The upwards vector of the camera.
 *
 * Returns the vector indicating the upwards direction of the camera.
 */
inline const Vec3& Raytracer::getUpVector() const {
   return upVector_;
}
   
/*!\brief Returns true if tbuffer output to a file is enabled.
 *
 * \return True if tbuffer output enabled, otherwise false.
 */
inline bool Raytracer::getTBufferOutputEnabled() const {
   return tBufferOutputEnabled_;
}

/*!\brief Returns the directory where the tbuffer output file will be saved in.
 *
 * \return Path to the tbuffer output directory.
 */
inline const std::string& Raytracer::getTBufferOutputDirectory() const {
   return tBufferOutputDirectory_;
}
   
/*!\brief Enabled / disable outputting the tBuffer to a file in the specified directory.
 * \param enabled Set to true / false to enable / disable tbuffer output.
 */
inline void Raytracer::setTBufferOutputEnabled(const bool enabled) {
   tBufferOutputEnabled_ = enabled;
}

/*!\brief Enable / disable outputting the tBuffer to a file in the specified directory.
 * \param enabled Set to true / false to enable / disable tbuffer output.
 */
inline void Raytracer::setTBufferOutputDirectory(const std::string& path) {
   namespace fs = boost::filesystem;
   
   fs::path dir (path);
   WALBERLA_CHECK(fs::exists(dir) && fs::is_directory(dir), "Tbuffer output directory " << path << " is invalid.");
   
   tBufferOutputDirectory_ = path;
}

/*!\brief Does one raytracing step.
 *
 * \param timestep The timestep after which the raytracing starts.
 *
 * \attention Planes will not get rendered if their normal and the rays direction point in the approximately
 * same direction. See Raytracer::isPlaneVisible() for further information.
 */
template <typename BodyTypeTuple>
void Raytracer::rayTrace(const size_t timestep) {
   real_t inf = std::numeric_limits<real_t>::max();

   std::vector<real_t> tBuffer(pixelsVertical_ * pixelsHorizontal_, inf);
   std::vector<BodyIntersectionInfo> intersections; // contains for each pixel information about an intersection, if existent
   
   real_t t, t_closest;
   walberla::id_t id_closest;
   RigidBody* body_closest = NULL;
   Ray ray(cameraPosition_, Vec3(1,0,0));
   IntersectsFunctor func(ray, t);
   tp_["Raytracing"].start();
   for (size_t x = 0; x < pixelsHorizontal_; x++) {
      for (size_t y = 0; y < pixelsVertical_; y++) {
         Vec3 pixelLocation = viewingPlaneOrigin_ + u_*(real_c(x)+real_t(0.5))*pixelWidth_ + v_*(real_c(y)+real_t(0.5))*pixelHeight_;
         Vec3 direction = (pixelLocation - cameraPosition_).getNormalized();
         ray.setDirection(direction);
         
         t_closest = inf;
         id_closest = 0;
         body_closest = NULL;
         for (auto blockIt = forest_->begin(); blockIt != forest_->end(); ++blockIt) {
#ifndef DISABLE_BLOCK_AABB_INTERSECTION_PRECHECK
            const AABB& blockAabb = blockIt->getAABB();
            if (!intersects(blockAabb, ray, t, blockAABBIntersectionPadding_)) {
               continue;
            }
#endif
            for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID_); bodyIt != LocalBodyIterator::end(); ++bodyIt) {
               if (bodyIt->getTypeID() == Plane::getStaticTypeID()) {
                  PlaneID plane = (Plane*)(*bodyIt);
                  if (!isPlaneVisible(plane, ray)) {
                     continue;
                  }
               }
               
               bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(*bodyIt, func);
               
               if (intersects && t < t_closest) {
                  // body was shot by ray and currently closest to camera
                  t_closest = t;
                  id_closest = bodyIt->getID();
                  body_closest = *bodyIt;
               }
            }
         }
         
         // only iterate over global body storage in one process.
         // optimization required, e.g. split up global bodies over all processes.
         WALBERLA_ROOT_SECTION() {
            for( auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt ) {
               if (bodyIt->getTypeID() == Plane::getStaticTypeID()) {
                  PlaneID plane = (Plane*)(*bodyIt);
                  if (!isPlaneVisible(plane, ray)) {
                     continue;
                  }
               }
               
               bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(*bodyIt, func);
               
               if (intersects && t < t_closest) {
                  // body was shot by ray and currently closest to camera
                  t_closest = t;
                  id_closest = bodyIt->getID();
                  body_closest = *bodyIt;
               }
            }
         }
         
         if (!realIsIdentical(t_closest, inf) && body_closest != NULL) {
            BodyIntersectionInfo intersectionInfo = {
               x,
               y,
               body_closest->getSystemID(),
               t_closest
            };
            intersections.push_back(intersectionInfo);
         }
         
         tBuffer[coordinateToArrayIndex(x, y)] = t_closest;
      }
   }
   tp_["Raytracing"].end();
   
   tp_["Reduction"].start();
   
   // intersections synchronisieren
   mpi::SendBuffer sendBuffer;
   for (auto& info: intersections) {
      sendBuffer << info.imageX << info.imageY << info.bodySystemID << info.t;
   }
   
   int gatheredIntersectionCount = 0;
   
   std::vector<real_t> fullImage(pixelsHorizontal_ * pixelsVertical_, inf);
   
   mpi::RecvBuffer recvBuffer;
   mpi::allGathervBuffer(sendBuffer, recvBuffer);
   BodyIntersectionInfo info;
   while (!recvBuffer.isEmpty()) {
      recvBuffer >> info.imageX;
      recvBuffer >> info.imageY;
      recvBuffer >> info.bodySystemID;
      recvBuffer >> info.t;
      
      size_t i = coordinateToArrayIndex(info.imageX, info.imageY);
      real_t currentFullImageT = fullImage[i];
      if (currentFullImageT > info.t) {
         fullImage[i] = info.t;
      }
      
      gatheredIntersectionCount++;
   }
   
   tp_["Reduction"].end();
   
   //WALBERLA_LOG_INFO("#particles visible: " << visibleBodyIDs.size());
   WALBERLA_LOG_INFO_ON_ROOT("#gathered intersections: " << gatheredIntersectionCount);
   
   auto tpReduced = tp_.getReduced();
   WALBERLA_ROOT_SECTION() {
      WALBERLA_LOG_INFO("Raytracing timing:");
      tpReduced->print(std::cout);
   }
   
   if (getTBufferOutputEnabled()) {
      writeTBufferToFile(tBuffer, timestep);
      WALBERLA_ROOT_SECTION() {
         writeTBufferToFile(fullImage, timestep, true);
      }
   }
}

}
}
}
