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

#include <core/config/Config.h>
#include <core/Filesystem.h>
#include <core/math/Vector3.h>
#include <core/mpi/Datatype.h>
#include <core/timing/TimingTree.h>

#include <pe/ccd/ICCD.h>
#include <pe/ccd/HashGrids.h>
#include <pe/raytracing/Ray.h>
#include <pe/raytracing/Intersects.h>
#include <pe/raytracing/Lighting.h>
#include <pe/raytracing/ShadingFunctions.h>
#include <pe/Types.h>

#include <cstddef>
#include <functional>

namespace walberla {
namespace pe {
namespace raytracing {

/*!\brief Contains information about a ray-body intersection.
 */
struct BodyIntersectionInfo {
   unsigned int imageX;          //!< Viewing plane x pixel coordinate the ray belongs to.   -> MPI_UNSIGNED
   unsigned int imageY;          //!< Viewing plane y pixel coordinate the ray belongs to.   -> MPI_UNSIGNED
   walberla::id_t bodySystemID;  //!< System ID of body which was hit.                       -> MPI_UNSIGNED_LONG_LONG
   double t;                     //!< Distance from camera to intersection point on body.    -> MPI_DOUBLE
   double r;                     //!< Red value for the pixel.                               -> MPI_DOUBLE
   double g;                     //!< Green value for the pixel.                             -> MPI_DOUBLE
   double b;                     //!< Blue value for the pixel.                              -> MPI_DOUBLE
};

inline bool defaultIsBodyVisible(const BodyID body) {
   WALBERLA_UNUSED(body);
   return true;
}

class Raytracer {
public:
   /*!\brief Which method to use when reducing the process-local image to a global one.
    */
   enum ReductionMethod {
      MPI_REDUCE,    //!< Reduce info from all processes onto root (assembling happens during reduction).
      MPI_GATHER     //!< Gather info from all processes onto root process and assemble global image there.
   };
   /*!\brief Which algorithm to use when doing ray-object intersection finding.
    */
   enum Algorithm {
      RAYTRACE_HASHGRIDS,              //!< Use hashgrids to find ray-body intersections.
      RAYTRACE_NAIVE,                  //!< Use the brute force approach of checking all objects for intersection testing.
      RAYTRACE_COMPARE_BOTH,           //!< Compare both methods and check for pixel errors.
      RAYTRACE_COMPARE_BOTH_STRICTLY   //!< Same as RAYTRACE_COMPARE_BOTH but abort if errors found.
   };
   
   /*!\name Constructors */
   //@{
   explicit Raytracer(const shared_ptr<BlockStorage>& forest, const BlockDataID storageID,
                      const shared_ptr<BodyStorage>& globalBodyStorage,
                      const BlockDataID ccdID,
                      uint16_t pixelsHorizontal, uint16_t pixelsVertical,
                      real_t fov_vertical, uint16_t antiAliasFactor,
                      const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                      const Lighting& lighting,
                      const Color& backgroundColor = Color(real_t(0.1), real_t(0.1), real_t(0.1)),
                      const std::function<ShadingParameters (const BodyID)>& bodyToShadingParamsFunc = defaultBodyTypeDependentShadingParams,
                      const std::function<bool (const BodyID)>& isBodyVisibleFunc = defaultIsBodyVisible);

   explicit Raytracer(const shared_ptr<BlockStorage>& forest, const BlockDataID storageID,
                      const shared_ptr<BodyStorage>& globalBodyStorage,
                      const BlockDataID ccdID,
                      const Config::BlockHandle& config,
                      const std::function<ShadingParameters (const BodyID)>& bodyToShadingParamsFunction = defaultBodyTypeDependentShadingParams,
                      const std::function<bool (const BodyID)>& isBodyVisibleFunc = defaultIsBodyVisible);
   //@}

private:
   /*!\name Member variables */
   //@{
   const shared_ptr<BlockStorage> forest_; //!< The BlockForest the raytracer operates on.
   const BlockDataID storageID_;    //!< The storage ID of the block data storage the raytracer operates on.
   const shared_ptr<BodyStorage> globalBodyStorage_; //!< The global body storage the raytracer operates on.
   const BlockDataID ccdID_;  //!< The ID of the hash grids block data.
   
   uint16_t pixelsHorizontal_;  //!< The horizontal amount of pixels of the generated image.
   uint16_t pixelsVertical_;    //!< The vertical amount of pixels of the generated image.
   real_t fov_vertical_;      //!< The vertical field-of-view of the camera.
   uint16_t antiAliasFactor_; /*!< Factor used for oversampling. Should be between 1 (fast, but jagged edges)
                               * and 4 (16 times slower, very smooth edges).*/
   Vec3 cameraPosition_;      //!< The position of the camera in the global world frame.
   Vec3 lookAtPoint_;         /*!< The point the camera looks at in the global world frame,
                               marks the center of the view plane.*/
   Vec3 upVector_;            //!< The vector indicating the upwards direction of the camera.
   Lighting lighting_;        //!< The lighting of the scene.
   Color backgroundColor_;    //!< Background color of the scene.
   
   bool imageOutputEnabled_;  //!< Enable / disable writing images to file.
   bool localImageOutputEnabled_; //!< Enable / disable writing images of the local process to file.
   std::string imageOutputDirectory_; //!< Path to the image output directory.
   
   uint8_t filenameTimestepWidth_; /*!< Width of the timestep number in output filenames.
                                   * Use e.g. 5 for ranges from 1 to 99 999: Will result in
                                   * filenames like image_00001.png up to image_99999.png. */
   uint8_t filenameRankWidth_;  //!< Width of the mpi rank part in a filename.
   bool confinePlanesToDomain_; //!< Enable to render only the parts of planes within the simulation domain.
   std::function<ShadingParameters (const BodyID)> bodyToShadingParamsFunc_; /*!< Function which returns a
                                                                              * ShadingParameters struct for the
                                                                              * specified body. */
   std::function<bool (const BodyID)> isBodyVisibleFunc_; /*!< Function which returns a boolean indicating if
                                                           * a given body should be visible in the final image. */
   Algorithm raytracingAlgorithm_;  //!< Algorithm to use while intersection testing.
   ReductionMethod reductionMethod_; //!< Reduction method used for assembling the image from all processes.
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
   
   MPI_Op bodyIntersectionInfo_reduction_op;
   MPI_Datatype bodyIntersectionInfo_mpi_type;
   
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
   inline bool getImageOutputEnabled() const;
   inline bool getLocalImageOutputEnabled() const;
   inline const std::string& getImageOutputDirectory() const;
   inline uint8_t getFilenameTimestepWidth() const;
   inline bool getConfinePlanesToDomain() const;
   //@}

   /*!\name Set functions */
   //@{
   inline void setBackgroundColor(const Color& color);
   inline void setImageOutputEnabled(const bool enabled);
   inline void setLocalImageOutputEnabled(const bool enabled);
   inline void setImageOutputDirectory(const std::string& path);
   inline void setFilenameTimestepWidth(uint8_t width);
   inline void setRaytracingAlgorithm(Algorithm algorithm);
   inline void setReductionMethod(ReductionMethod reductionMethod);
   inline void setConfinePlanesToDomain(bool confinePlanesToOrigin);
   //@}
   
   /*!\name Functions */
   //@{
   template <typename BodyTypeTuple>
   void generateImage(const size_t timestep, WcTimingTree* tt = nullptr );
   
   void setupView_();
   void setupFilenameRankWidth_();
   void setupMPI_();
   
private:
   void localOutput(const std::vector<BodyIntersectionInfo>& intersectionsBuffer, size_t timestep,
                    WcTimingTree* tt = nullptr);
   void output(const std::vector<BodyIntersectionInfo>& intersectionsBuffer, size_t timestep,
               WcTimingTree* tt = nullptr);

   std::string getOutputFilename(const std::string& base, size_t timestep, bool isGlobalImage) const;
   void writeImageToFile(const std::vector<BodyIntersectionInfo>& intersectionsBuffer,
                         size_t timestep, bool isGlobalImage = false) const;
   void writeImageToFile(const std::vector<BodyIntersectionInfo>& intersectionsBuffer,
                         const std::string& fileName) const;
   
   void syncImageUsingMPIReduce(std::vector<BodyIntersectionInfo>& intersectionsBuffer, WcTimingTree* tt = nullptr);
   void syncImageUsingMPIGather(std::vector<BodyIntersectionInfo>& intersections,
                                std::vector<BodyIntersectionInfo>& intersectionsBuffer, WcTimingTree* tt = nullptr);
   
   inline bool isPlaneVisible(const PlaneID plane, const Ray& ray) const;
   inline size_t coordinateToArrayIndex(size_t x, size_t y) const;
   
   template <typename BodyTypeTuple>
   inline void traceRayInGlobalBodyStorage(const Ray& ray, BodyID& body_closest, real_t& t_closest, Vec3& n_closest) const;
   template <typename BodyTypeTuple>
   inline void traceRayNaively(const Ray& ray, BodyID& body_closest, real_t& t_closest, Vec3& n_closest) const;
   template <typename BodyTypeTuple>
   inline void traceRayInHashGrids(const Ray& ray, BodyID& body_closest, real_t& t_closest, Vec3& n_closest) const;

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
inline uint8_t Raytracer::getFilenameTimestepWidth() const {
   return filenameTimestepWidth_;
}

/*!\brief Returns whether the rendering of planes gets confined in the simulation domain.
 * \return True if the rendering of planes gets confined in the simulation domain.
 */
inline bool Raytracer::getConfinePlanesToDomain() const {
   return confinePlanesToDomain_;
}

/*!\brief Set the background color of the scene.
 *
 * \param color New background color.
 */
inline void Raytracer::setBackgroundColor(const Color& color) {
   backgroundColor_ = color;
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
   filesystem::path dir (path);
   WALBERLA_CHECK(filesystem::exists(dir) && filesystem::is_directory(dir), "Image output directory " << path << " is invalid.");
   
   imageOutputDirectory_ = path;
}

/*!\brief Set width of timestep number in output filenames.
 * \param width Width of timestep part in a filename.
 */
inline void Raytracer::setFilenameTimestepWidth(uint8_t width) {
   filenameTimestepWidth_ = width;
}

/*!\brief Set the algorithm to use while ray tracing.
 * \param algorithm One of RAYTRACE_HASHGRIDS, RAYTRACE_NAIVE, RAYTRACE_COMPARE_BOTH, RAYTRACE_COMPARE_BOTH_STRICTLY (abort on errors).
 */
inline void Raytracer::setRaytracingAlgorithm(Algorithm algorithm) {
   raytracingAlgorithm_ = algorithm;
}

/*!\brief Set the algorithm to use while reducing.
 * \param reductionMethod One of MPI_GATHER or MPI_REDUCE (latter one only works on MPI builds).
 */
inline void Raytracer::setReductionMethod(ReductionMethod reductionMethod) {
   reductionMethod_ = reductionMethod;
}

/*!\brief Set if the rendering of planes should get confined to the simulation domain.
 * \param confinePlanesToOrigin True if the rendering of planes should get confined to the simulation domain.
 */
inline void Raytracer::setConfinePlanesToDomain(bool confinePlanesToDomain) {
   confinePlanesToDomain_ = confinePlanesToDomain;
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
   return y*pixelsHorizontal_*antiAliasFactor_ + x;
}

/*!\brief Traces a ray in the global body storage and finds the closest ray-body intersection.
 * \param ray Ray which is shot.
 * \param body_closest Reference where the closest body will be stored in.
 * \param t_closest Reference where the distance of the currently closest body is stored in,
                    will get updated if closer intersection found.
 * \param n_closest Reference where the intersection normal will be stored in.
 */
template <typename BodyTypeTuple>
inline void Raytracer::traceRayInGlobalBodyStorage(const Ray& ray, BodyID& body_closest, real_t& t_closest, Vec3& n_closest) const {
   WALBERLA_ROOT_SECTION(){
      real_t t = std::numeric_limits<real_t>::max();
      Vec3 n;
      
      IntersectsFunctor func(ray, t, n);
      
      for(auto bodyIt = globalBodyStorage_->begin(); bodyIt != globalBodyStorage_->end(); ++bodyIt)
      {
         if (!isBodyVisibleFunc_(bodyIt.getBodyID()))
         {
            continue;
         }
         
         bool isPlane = (bodyIt->getTypeID() == Plane::getStaticTypeID());
         if (isPlane)
         {
            PlaneID plane = (PlaneID)bodyIt.getBodyID();
            if (!isPlaneVisible(plane, ray))
            {
               continue;
            }
         }
         
         bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(bodyIt.getBodyID(), func);
         if (intersects && t < t_closest)
         {
            if (isPlane && confinePlanesToDomain_)
            {
               Vec3 intersectionPoint = ray.getOrigin()+ray.getDirection()*t;
               if (!forest_->getDomain().contains(intersectionPoint, real_t(1e-8)))
               {
                  continue;
               }
            }
            // body was shot by ray and is currently closest to camera
            t_closest = t;
            body_closest = bodyIt.getBodyID();
            n_closest = n;
         }
      }
   }
}

/*!\brief Traces a ray naively and finds the closest ray-body intersection.
 * \param ray Ray which is shot.
 * \param body_closest Reference where the closest body will be stored in.
 * \param t_closest Reference where the distance of the currently closest body is stored in,
                    will get updated if closer intersection found.
 * \param n_closest Reference where the intersection normal will be stored in.
 */
template <typename BodyTypeTuple>
inline void Raytracer::traceRayNaively(const Ray& ray, BodyID& body_closest, real_t& t_closest, Vec3& n_closest) const {
   real_t t = std::numeric_limits<real_t>::max();
   Vec3 n;
   
   IntersectsFunctor func(ray, t, n);
   
   for (auto blockIt = forest_->begin(); blockIt != forest_->end(); ++blockIt) {
      for (auto bodyIt = LocalBodyIterator::begin(*blockIt, storageID_); bodyIt != LocalBodyIterator::end(); ++bodyIt) {
         if (!isBodyVisibleFunc_(bodyIt.getBodyID())) {
            continue;
         }
         
         bool intersects = SingleCast<BodyTypeTuple, IntersectsFunctor, bool>::execute(bodyIt.getBodyID(), func);
         if (intersects && t < t_closest) {
            // body was shot by ray and is currently closest to camera
            t_closest = t;
            body_closest = bodyIt.getBodyID();
            n_closest = n;
         }
      }
   }
}

/*!\brief Traces a ray in the global body storage and finds the closest ray-body intersection.
 * \param ray Ray which is shot.
 * \param body_closest Reference where the closest body will be stored in.
 * \param t_closest Reference where the distance of the currently closest body is stored in,
                    will get updated if closer intersection found.
 * \param n_closest Reference where the intersection normal will be stored in.
 */
template <typename BodyTypeTuple>
inline void Raytracer::traceRayInHashGrids(const Ray& ray, BodyID& body_closest, real_t& t_closest, Vec3& n_closest) const {
   real_t t = std::numeric_limits<real_t>::max();
   Vec3 n;
   
   for (auto blockIt = forest_->begin(); blockIt != forest_->end(); ++blockIt) {
      const AABB& blockAABB = blockIt->getAABB();
      const ccd::HashGrids* hashgrids = blockIt->uncheckedFastGetData<ccd::HashGrids>(ccdID_);
      BodyID body = hashgrids->getClosestBodyIntersectingWithRay<BodyTypeTuple>(ray, blockAABB, t, n,
                                                                                isBodyVisibleFunc_);
      if (body != nullptr && t < t_closest) {
         t_closest = t;
         body_closest = body;
         n_closest = n;
      }
   }
}
   
/*!\brief Does one raytracing step.
 *
 * \param timestep The timestep after which the raytracing starts.
 *
 * \attention Planes will not get rendered if their normal and the rays direction point in the approximately
 * same direction. See Raytracer::isPlaneVisible() for further information.
 */
template <typename BodyTypeTuple>
void Raytracer::generateImage(const size_t timestep, WcTimingTree* tt) {
   if (tt != nullptr) tt->start("Raytracing");
   const real_t realMax = std::numeric_limits<real_t>::max();
   
   std::vector<BodyIntersectionInfo> intersections;
   // contains for each pixel information about an intersection:
   size_t bufferSize = (pixelsVertical_*antiAliasFactor_)*(pixelsHorizontal_*antiAliasFactor_);
   std::vector<BodyIntersectionInfo> intersectionsBuffer(bufferSize);

   if (raytracingAlgorithm_ == RAYTRACE_HASHGRIDS || raytracingAlgorithm_ == RAYTRACE_COMPARE_BOTH
      || raytracingAlgorithm_ == RAYTRACE_COMPARE_BOTH_STRICTLY) {
      if (tt != nullptr) tt->start("HashGrids Update");
      for (auto blockIt = forest_->begin(); blockIt != forest_->end(); ++blockIt) {
         ccd::HashGrids* hashgrids = blockIt->getData<ccd::HashGrids>(ccdID_);
         hashgrids->update();
      }
      if (tt != nullptr) tt->stop("HashGrids Update");
   }
   
   real_t t, t_closest;
   Vec3 n;
   Vec3 n_closest;
   BodyID body_closest = nullptr;
   Ray ray(cameraPosition_, Vec3(1,0,0));
   IntersectsFunctor func(ray, t, n);
   bool isErrorneousPixel = false;
   uint_t pixelErrors = 0;
   std::map<BodyID, std::unordered_set<BodyID>> correctToIncorrectBodyIDsMap;
   
   if (tt != nullptr) tt->start("Intersection Testing");
   for (size_t x = 0; x < pixelsHorizontal_*antiAliasFactor_; x++) {
      for (size_t y = 0; y < pixelsVertical_*antiAliasFactor_; y++) {
         Vec3 pixelLocation = viewingPlaneOrigin_ + u_*(real_c(x)+real_t(0.5))*pixelWidth_ + v_*(real_c(y)+real_t(0.5))*pixelHeight_;
         Vec3 direction = (pixelLocation - cameraPosition_).getNormalized();
         ray.setDirection(direction);
         
         n.reset();
         t_closest = realMax;
         body_closest = nullptr;
         
         if (raytracingAlgorithm_ == RAYTRACE_HASHGRIDS) {
            traceRayInHashGrids<BodyTypeTuple>(ray, body_closest, t_closest, n_closest);
         } else if (raytracingAlgorithm_ == RAYTRACE_NAIVE) {
            traceRayNaively<BodyTypeTuple>(ray, body_closest, t_closest, n_closest);
         } else {
            traceRayInHashGrids<BodyTypeTuple>(ray, body_closest, t_closest, n_closest);
            BodyID hashgrids_body_closest = body_closest;
            
            t_closest = realMax;
            body_closest = nullptr;
            traceRayNaively<BodyTypeTuple>(ray, body_closest, t_closest, n_closest);
            
            if (body_closest != hashgrids_body_closest) {
               correctToIncorrectBodyIDsMap[body_closest].insert(hashgrids_body_closest);
               isErrorneousPixel = true;
               ++pixelErrors;
            }
         }
         
         traceRayInGlobalBodyStorage<BodyTypeTuple>(ray, body_closest, t_closest, n_closest);
         
         BodyIntersectionInfo& intersectionInfo = intersectionsBuffer[coordinateToArrayIndex(x, y)];
         intersectionInfo.imageX = uint32_t(x);
         intersectionInfo.imageY = uint32_t(y);
         
         if (!realIsIdentical(t_closest, realMax) && body_closest != nullptr) {
            Color color = getColor(body_closest, ray, t_closest, n_closest);
            if (isErrorneousPixel) {
               color = Color(1,0,0);
               isErrorneousPixel = false;
            }
            
            intersectionInfo.bodySystemID = body_closest->getSystemID();
            intersectionInfo.t = t_closest;
            intersectionInfo.r = color[0];
            intersectionInfo.g = color[1];
            intersectionInfo.b = color[2];
            
            intersections.push_back(intersectionInfo);
         } else {
            intersectionInfo.bodySystemID = 0;
            intersectionInfo.t = realMax;
            intersectionInfo.r = backgroundColor_[0];
            intersectionInfo.g = backgroundColor_[1];
            intersectionInfo.b = backgroundColor_[2];
         }
      }
   }
   if (tt != nullptr) tt->stop("Intersection Testing");

   if (raytracingAlgorithm_ == RAYTRACE_COMPARE_BOTH || raytracingAlgorithm_ == RAYTRACE_COMPARE_BOTH_STRICTLY) {
      if (pixelErrors > 0) {
         WALBERLA_LOG_WARNING(pixelErrors << " pixel errors found!");
         
         std::stringstream ss;
         for (auto it: correctToIncorrectBodyIDsMap) {
            const BodyID correctBody = it.first;
            if (it.first != nullptr) {
               ss << " correct body: " << correctBody->getID() << "(" << correctBody->getHash() << ")";
            } else {
               ss << " no body naively found";
            }
            ss << ", hashgrids found:";
            for (auto incorrectBody: it.second) {
               ss << " ";
               if (incorrectBody != nullptr) {
                  ss << incorrectBody->getID() << "(" << incorrectBody->getHash();
               } else {
                  ss << "NULL";
               }
            }
            ss << std::endl;
            ss << "  minCorner: " << correctBody->getAABB().minCorner() << ", maxCorner: " << correctBody->getAABB().maxCorner();
            ss << std::endl;
         }
         WALBERLA_LOG_WARNING("Problematic bodies: " << std::endl << ss.str());
         
         if (raytracingAlgorithm_ == RAYTRACE_COMPARE_BOTH_STRICTLY) {
            WALBERLA_ABORT("Pixel errors found, aborting due to strict comparison option.");
         }
      } else {
         WALBERLA_LOG_INFO("No pixel errors found.");
      }
   }
   
   localOutput(intersectionsBuffer, timestep, tt);
   
   // Reduction with different methods only makes sense if actually using MPI.
   // Besides that, the MPI reduce routine does not compile without MPI.
   WALBERLA_MPI_SECTION() {
      switch(reductionMethod_) {
         case MPI_REDUCE:
            syncImageUsingMPIReduce(intersectionsBuffer, tt);
            break;
         case MPI_GATHER:
            syncImageUsingMPIGather(intersections, intersectionsBuffer, tt);
            break;
      }
   } else {
      syncImageUsingMPIGather(intersections, intersectionsBuffer, tt);
   }
   
   output(intersectionsBuffer, timestep, tt);
   
   if (tt != nullptr) tt->stop("Raytracing");
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
   const ShadingParameters shadingParams = bodyToShadingParamsFunc_(body);
   
   const Vec3 intersectionPoint = ray.getOrigin() + ray.getDirection() * t;
   Vec3 lightDirection = lighting_.pointLightOrigin - intersectionPoint;
   lightDirection = lightDirection.getNormalized();
   
   real_t lambertian = std::max(real_t(0), lightDirection * n);
   
   real_t specular = real_t(0);
   
   if (lambertian > 0 || realIsEqual(lambertian, real_t(0))) {
      // Blinn-Phong
      Vec3 viewDirection = -ray.getDirection();
      Vec3 halfDirection = (lightDirection + viewDirection).getNormalized();
      real_t specularAngle = std::max(halfDirection * n, real_t(0));
      specular = real_c(pow(specularAngle, shadingParams.shininess));
   }
   
   Color color = lighting_.ambientColor.mulComponentWise(shadingParams.ambientColor)
      + lighting_.diffuseColor.mulComponentWise(shadingParams.diffuseColor)*lambertian
      + lighting_.specularColor.mulComponentWise(shadingParams.specularColor)*specular;
   
   // Capping of color channels to 1.
   // Capping instead of scaling will make specular highlights stronger.
   color.clamp();

   return color;
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
