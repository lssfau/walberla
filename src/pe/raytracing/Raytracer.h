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
#include "Lighting.h"

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
   Vec3 color;                   //!< color computed for the pixel.
};

class Raytracer {
public:
   /*!\name Constructors */
   //@{
   explicit Raytracer(const shared_ptr<BlockStorage> forest, BlockDataID storageID,
                      const shared_ptr<BodyStorage> globalBodyStorage,
                      uint16_t pixelsHorizontal, uint16_t pixelsVertical,
                      real_t fov_vertical,
                      const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                      const Lighting& lighting,
                      const Color& backgroundColor = Color(0.1, 0.1, 0.1),
                      real_t blockAABBIntersectionPadding = real_t(0.0));
   explicit Raytracer(const shared_ptr<BlockStorage> forest, BlockDataID storageID,
                      const shared_ptr<BodyStorage> globalBodyStorage,
                      const Config::BlockHandle& config);
   //@}

private:
   /*!\name Member variables */
   //@{
   const shared_ptr<BlockStorage> forest_; //!< The BlockForest the raytracer operates on.
   BlockDataID storageID_;    //!< The storage ID of the block data storage the raytracer operates on.
   const shared_ptr<BodyStorage> globalBodyStorage_; //!< The global body storage the raytracer operates on.
   
   uint16_t pixelsHorizontal_;  //!< The horizontal amount of pixels of the generated image.
   uint16_t pixelsVertical_;    //!< The vertical amount of pixels of the generated image.
   real_t fov_vertical_;      //!< The vertical field-of-view of the camera.
   Vec3 cameraPosition_;      //!< The position of the camera in the global world frame.
   Vec3 lookAtPoint_;         /*!< The point the camera looks at in the global world frame,
                               marks the center of the view plane.*/
   Vec3 upVector_;            //!< The vector indicating the upwards direction of the camera.
   Lighting lighting_;        //!< The lighting of the scene.
   Color backgroundColor_;    //!< Background color of the scene.
   real_t blockAABBIntersectionPadding_; /*!< The padding applied in block AABB intersection pretesting, as
                                          some objects within a block might protrude from the block's AABB.*/
   
   bool tBufferOutputEnabled_; //!< Enable / disable dumping the tbuffer to file.
   std::string tBufferOutputDirectory_; //!< Path to the tbuffer output directory.
   
   bool imageOutputEnabled_;  //!< Enable / disable writing images to file.
   bool localImageOutputEnabled_; //!< Enable / disable writing images of the local process to file.
   std::string imageOutputDirectory_; //!< Path to the image output directory.
   
   int8_t filenameTimestepWidth_; /*!< Width of the timestep number in output filenames.
                                  * Use e.g. 5 for ranges from 1 to 99 999: Will result in
                                  * filenames like image_00001.png up to image_99999.png. */
   int8_t filenameRankWidth_;  //!< Width of the mpi rank part in a filename.
   //@}
   
   /*!\name Member variables for raytracing geometry */
   //@{
   Vec3 n_;                   //!< The normal vector of the viewing plane.
   Vec3 u_;                   //!< The vector spanning the viewing plane in the "right direction".
   Vec3 v_;                   //!< The vector spanning the viewing plane in the "up direction".
   real_t d_;                 //!< The the distance from camera to viewing plane.
   real_t aspectRatio_;       //!< The aspect ratio of the generated image and viewing plane.
   real_t viewingPlaneHeight_; //!< The viewing plane height in the world frame.
   real_t viewingPlaneWidth_; //!< The viewing plane width in the world frame.
   Vec3 viewingPlaneOrigin_;  //!< The origin of the viewing plane.
   real_t pixelWidth_;        //!< The width of a pixel of the generated image in the viewing plane.
   real_t pixelHeight_;       //!< The height of a pixel of the generated image in the viewing plane.
   //@}
   
   WcTimingPool tp_;

public:
   /*!\name Get functions */
   //@{
   inline uint16_t getPixelsHorizontal() const;
   inline uint16_t getPixelsVertical() const;
   inline real_t getFOVVertical() const;
   inline const Vec3& getCameraPosition() const;
   inline const Vec3& getLookAtPoint() const;
   inline const Vec3& getUpVector() const;
   inline const Color& getBackgroundColor() const;
   inline bool getTBufferOutputEnabled() const;
   inline const std::string& getTBufferOutputDirectory() const;
   inline bool getImageOutputEnabled() const;
   inline bool getLocalImageOutputEnabled() const;
   inline const std::string& getImageOutputDirectory() const;
   inline int8_t getFilenameTimestepWidth() const;
   //@}

   /*!\name Set functions */
   //@{
   inline void setBackgroundColor(const Color& color);
   inline void setTBufferOutputEnabled(const bool enabled);
   inline void setTBufferOutputDirectory(const std::string& path);
   inline void setImageOutputEnabled(const bool enabled);
   inline void setLocalImageOutputEnabled(const bool enabled);
   inline void setImageOutputDirectory(const std::string& path);
   inline void setFilenameTimestepWidth(int8_t width);
   //@}
   
   /*!\name Functions */
   //@{
   void setupView_();
   void setupFilenameRankWidth_();
   template <typename BodyTypeTuple>
   void rayTrace(const size_t timestep);
   
private:
   std::string getOutputFilename(const std::string& base, size_t timestep, bool isGlobalImage) const;
   void writeTBufferToFile(const std::vector<real_t>& tBuffer, size_t timestep, bool isGlobalImage = false) const;
   void writeTBufferToFile(const std::vector<real_t>& tBuffer, const std::string& fileName) const;
   void writeImageBufferToFile(const std::vector<Color>& imageBuffer, size_t timestep, bool isGlobalImage = false) const;
   void writeImageBufferToFile(const std::vector<Color>& imageBuffer, const std::string& fileName) const;
   
   inline bool isPlaneVisible(const PlaneID plane, const Ray& ray) const;
   inline size_t coordinateToArrayIndex(size_t x, size_t y) const;
   
   inline Color getColor(const BodyID body, const Ray& ray, real_t t, const Vec3& n) const;
   //@}
};
   
/*!\brief Returns the horizontal amount of pixels of the generated image.
 *
 * \return The horizontal amount of pixels of the generated image.
 */
inline uint16_t Raytracer::getPixelsHorizontal() const {
   return pixelsHorizontal_;
}

/*!\brief Returns the vertical amount of pixels of the generated image.
 *
 * \return The vertical amount of pixels of the generated image.
 */
inline uint16_t Raytracer::getPixelsVertical() const {
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

/*!\brief Returns the background color of the scene.
 *
 * \return Color.
 *
 * Returns the background color of the scene.
 */
inline const Color& Raytracer::getBackgroundColor() const {
   return backgroundColor_;
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
   
/*!\brief Returns true if image output to a file is enabled.
 *
 * \return True if image output enabled, otherwise false.
 */
inline bool Raytracer::getImageOutputEnabled() const {
   return imageOutputEnabled_;
}

/*!\brief Returns true if local image output to a file is enabled.
 *
 * \return True if local image output enabled, otherwise false.
 */
inline bool Raytracer::getLocalImageOutputEnabled() const {
   return localImageOutputEnabled_;
}
   
/*!\brief Returns the directory where the images will be saved in.
 *
 * \return Path to the image output directory.
 */
inline const std::string& Raytracer::getImageOutputDirectory() const {
   return imageOutputDirectory_;
}

/*!\brief Returns width of the timestep number in output filenames.
 * \return Width of the timestep part in filenames.
 */
inline int8_t Raytracer::getFilenameTimestepWidth() const {
   return filenameTimestepWidth_;
}

/*!\brief Set the background color of the scene.
 *
 * \param color New background color.
 */
inline void Raytracer::setBackgroundColor(const Color& color) {
   backgroundColor_ = color;
}
   
/*!\brief Enable / disable outputting the tBuffer to a file in the specified directory.
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
   
/*!\brief Enable / disable outputting images in the specified directory.
 * \param enabled Set to true / false to enable / disable image output.
 */
inline void Raytracer::setImageOutputEnabled(const bool enabled) {
   imageOutputEnabled_ = enabled;
}

/*!\brief Enable / disable outputting local images in the specified directory.
 * \param enabled Set to true / false to enable / disable image output.
 */
inline void Raytracer::setLocalImageOutputEnabled(const bool enabled) {
   localImageOutputEnabled_ = enabled;
}
   
/*!\brief Enable / disable outputting images in the specified directory.
 * \param enabled Set to true / false to enable / disable image output.
 */
inline void Raytracer::setImageOutputDirectory(const std::string& path) {
   namespace fs = boost::filesystem;
   
   fs::path dir (path);
   WALBERLA_CHECK(fs::exists(dir) && fs::is_directory(dir), "Image output directory " << path << " is invalid.");
   
   imageOutputDirectory_ = path;
}

/*!\brief Set width of timestep number in output filenames.
 * \param width Width of timestep part in a filename.
 */
inline void Raytracer::setFilenameTimestepWidth(int8_t width) {
   filenameTimestepWidth_ = width;
}
   
/*!\brief Checks if a plane should get rendered.
 * \param plane Plane to check for visibility.
 * \param ray Ray which is intersected with plane.
 *
 * Checks if a plane should get rendered by comparing the planes normal and the ray direction.
 * If the rays direction vectors projection on the planes normal is positive, the plane is considered invisible.
 */
inline bool Raytracer::isPlaneVisible(const PlaneID plane, const Ray& ray) const {
   return plane->getNormal() * ray.getDirection() < 0;
}

/*!\brief Converts a coordinate to an array index.
 * \param x X component of the coordinate.
 * \param y Y component of the coordinate.
 * \return Array index.
 */
inline size_t Raytracer::coordinateToArrayIndex(size_t x, size_t y) const {
   return y*pixelsHorizontal_ + x;
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

   int numProcesses = mpi::MPIManager::instance()->numProcesses();
   int rank = mpi::MPIManager::instance()->rank();
   
   std::vector<real_t> tBuffer(pixelsVertical_ * pixelsHorizontal_, inf);
   std::vector<Color> imageBuffer(pixelsVertical_ * pixelsHorizontal_);
   std::vector<BodyIntersectionInfo> intersections; // contains for each pixel information about an intersection, if existent
   
   real_t t, t_closest;
   Vec3 n;
   Vec3 n_closest;
   BodyID body_closest = NULL;
   Ray ray(cameraPosition_, Vec3(1,0,0));
   IntersectsFunctor func(ray, t, n);
   tp_["Raytracing"].start();
   for (size_t x = 0; x < pixelsHorizontal_; x++) {
      for (size_t y = 0; y < pixelsVertical_; y++) {
         Vec3 pixelLocation = viewingPlaneOrigin_ + u_*(real_c(x)+real_t(0.5))*pixelWidth_ + v_*(real_c(y)+real_t(0.5))*pixelHeight_;
         Vec3 direction = (pixelLocation - cameraPosition_).getNormalized();
         ray.setDirection(direction);
         
         n.reset();
         t_closest = inf;
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
                  body_closest = *bodyIt;
                  n_closest = n;
               }
            }
         }
         
         int i = 0;
         for(auto bodyIt: *globalBodyStorage_) {
            i++;
            // distribute global objects more or less evenly on all processes
            if (((i-1)%numProcesses) != rank) {
               continue;
            }
            
            if (bodyIt->getTypeID() == Plane::getStaticTypeID()) {
               PlaneID plane = (PlaneID)bodyIt;
               if (!isPlaneVisible(plane, ray)) {
                  continue;
               }
            }
            
            bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(bodyIt, func);
            
            if (intersects && t < t_closest) {
               // body was shot by ray and currently closest to camera
               t_closest = t;
               body_closest = bodyIt;
               n_closest = n;
            }
         }
         
         if (!realIsIdentical(t_closest, inf) && body_closest != NULL) {
            Color color = getColor(body_closest, ray, t_closest, n_closest);
            imageBuffer[coordinateToArrayIndex(x, y)] = color;
            BodyIntersectionInfo intersectionInfo = {
               x,
               y,
               body_closest->getSystemID(),
               t_closest,
               color
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
      sendBuffer << info.imageX << info.imageY
      << info.bodySystemID << info.t
      << info.color[0] << info.color[1] << info.color[2];
   }
   int gatheredIntersectionCount = 0;
   std::vector<real_t> fullTBuffer(pixelsHorizontal_ * pixelsVertical_, inf);
   std::vector<Color> fullImageBuffer(pixelsHorizontal_ * pixelsVertical_, backgroundColor_);
   mpi::RecvBuffer recvBuffer;
   
   mpi::gathervBuffer(sendBuffer, recvBuffer, 0);
   //mpi::allGathervBuffer(sendBuffer, recvBuffer);
   
   WALBERLA_ROOT_SECTION() {
      BodyIntersectionInfo info;
      while (!recvBuffer.isEmpty()) {
         recvBuffer >> info.imageX;
         recvBuffer >> info.imageY;
         recvBuffer >> info.bodySystemID;
         recvBuffer >> info.t;
         recvBuffer >> info.color[0];
         recvBuffer >> info.color[1];
         recvBuffer >> info.color[2];
         
         size_t i = coordinateToArrayIndex(info.imageX, info.imageY);
         real_t currentFullTBufferT = fullTBuffer[i];
         if (currentFullTBufferT > info.t) {
            fullTBuffer[i] = info.t;
            fullImageBuffer[i] = info.color;
         }
         
         gatheredIntersectionCount++;
      }
   }
   
   tp_["Reduction"].end();
   
   //WALBERLA_LOG_INFO("#particles visible: " << visibleBodyIDs.size());
   WALBERLA_LOG_INFO_ON_ROOT("#gathered intersections: " << gatheredIntersectionCount);
   
   auto tpReduced = tp_.getReduced();
   WALBERLA_ROOT_SECTION() {
      WALBERLA_LOG_INFO("Raytracing timing:");
      tpReduced->print(std::cout);
   }
   
   if (getImageOutputEnabled()) {
      if (getLocalImageOutputEnabled()) {
         writeImageBufferToFile(imageBuffer, timestep);
      }
      WALBERLA_ROOT_SECTION() {
         writeImageBufferToFile(fullImageBuffer, timestep, true);
      }
   }
   
   if (getTBufferOutputEnabled()) {
      writeTBufferToFile(tBuffer, timestep);
      WALBERLA_ROOT_SECTION() {
         writeTBufferToFile(fullTBuffer, timestep, true);
      }
   }
}

/*!\brief Computes the color for a certain intersection.
 *
 * \param body Intersected body.
 * \param Ray Ray which intersected the body.
 * \param t Distance from eye to intersection point.
 * \param n Intersection normal at the intersection point.
 *
 * \return Computed color.
 */
inline Color Raytracer::getColor(const BodyID body, const Ray& ray, real_t t, const Vec3& n) const {
   //----
   Color diffuseColor(0.6, 0, 0.9);
   Color specularColor(0.8, 0.8, 0.8);
   Color ambientColor(0.5, 0, 0.8);
   real_t shininess = 100;

   if (body->getTypeID() == Plane::getStaticTypeID()) {
      diffuseColor = Color(real_t(0.7), real_t(0.7), real_t(0.7));
      ambientColor.set(real_t(0.5), real_t(0.5), real_t(0.5));
      specularColor.set(real_t(0), real_t(0), real_t(0));
      shininess = real_t(0);
   }
   if (body->getTypeID() == Sphere::getStaticTypeID()) {
      diffuseColor = Color(real_t(0.98), real_t(0.1), real_t(0.1));
      ambientColor.set(real_t(0.6), real_t(0.05), real_t(0.05));
      specularColor.set(real_t(1), real_t(1), real_t(1));
      shininess = real_t(30);
   }
   if (body->getTypeID() == Capsule::getStaticTypeID()) {
      diffuseColor = Color(real_t(0.15), real_t(0.44), real_t(0.91));
      ambientColor.set(real_t(0), real_t(0), real_t(0.3));
      specularColor.set(real_t(1), real_t(1), real_t(1));
      shininess = real_t(20);
   }
   //----
   
   const Vec3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
   Vec3 lightDirection = lighting_.pointLightOrigin - intersectionPoint;
   lightDirection = lightDirection.getNormalized();
   
   real_t lambertian = std::max(real_t(0), lightDirection * n);
   
   real_t specular = real_t(0);
   
   if (lambertian > 0) {
      // Blinn-Phong
      Vec3 viewDirection = -ray.getDirection();
      Vec3 halfDirection = (lightDirection + viewDirection).getNormalized();
      real_t specularAngle = std::max(halfDirection * n, real_t(0));
      specular = real_c(pow(specularAngle, shininess));
   }
   
   Color color = lighting_.ambientColor.mulComponentWise(ambientColor)
      + lighting_.diffuseColor.mulComponentWise(diffuseColor)*lambertian
      + lighting_.specularColor.mulComponentWise(specularColor)*specular;
   
   // Capping of color channels to 1.
   // Capping instead of scaling will make specular highlights stronger.
   color.clamp();

   return color;
}
}
}
}
