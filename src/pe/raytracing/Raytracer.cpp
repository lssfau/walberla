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

using namespace walberla;

real_t deg2rad(real_t deg) {
   return deg * math::M_PI / real_t(180.0);
}

namespace walberla {
namespace pe {
namespace raytracing {
/*!\brief Instantiation constructor for the Raytracer class.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Storage ID of the block data storage the raytracer operates on.
 * \param pixelsHorizontal Horizontal amount of pixels of the generated image.
 * \param pixelsVertical Vertical amount of pixels of the generated image.
 * \param fov_vertical Vertical field-of-view of the camera.
 * \param cameraPosition Position of the camera in the global world frame.
 * \param lookAtPoint Point the camera looks at in the global world frame.
 * \param upVector Vector indicating the upwards direction of the camera.
 * \param blockAABBIntersectionPadding The padding applied in block AABB intersection pretesting. Usually not required.
 *                                     Set it to the value of the farthest distance a object might protrude from
 *                                     its containing block.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage> forest, BlockDataID storageID,
                     const shared_ptr<BodyStorage> globalBodyStorage,
                     size_t pixelsHorizontal, size_t pixelsVertical,
                     real_t fov_vertical,
                     const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                     real_t blockAABBIntersectionPadding)
   : forest_(forest), storageID_(storageID), globalBodyStorage_(globalBodyStorage),
   pixelsHorizontal_(pixelsHorizontal), pixelsVertical_(pixelsVertical),
   fov_vertical_(fov_vertical),
   cameraPosition_(cameraPosition), lookAtPoint_(lookAtPoint), upVector_(upVector),
   blockAABBIntersectionPadding_(blockAABBIntersectionPadding),
   tBufferOutputEnabled_(false)
{
   setupView_();
}

/*!\brief Instantiation constructor for the Raytracer class using a config object for view setup.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Storage ID of the block data storage the raytracer operates on.
 * \param config Config block for the raytracer.
 *
 * The config block has to contain image_x (int), image_y (int), fov_vertical (real, in degrees)
 * and tbuffer_output_directory (string) parameters. Additionally a vector of reals
 * for each of cameraPosition, lookAt and the upVector. Optional is blockAABBIntersectionPadding (real).
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage> forest, BlockDataID storageID,
                     const shared_ptr<BodyStorage> globalBodyStorage,
                     const Config::BlockHandle& config)
   : forest_(forest), storageID_(storageID), globalBodyStorage_(globalBodyStorage) {
   WALBERLA_CHECK(config.isValid(), "No valid config passed to raytracer");
   
   pixelsHorizontal_ = config.getParameter<size_t>("image_x");
   pixelsVertical_ = config.getParameter<size_t>("image_y");
   fov_vertical_ = config.getParameter<real_t>("fov_vertical");
   
   if (config.isDefined("tbuffer_output_directory")) {
      setTBufferOutputEnabled(true);
      setTBufferOutputDirectory(config.getParameter<std::string>("tbuffer_output_directory"));
      WALBERLA_LOG_INFO_ON_ROOT("t buffers will be written to " << getTBufferOutputDirectory() << ".");
   }
   
   cameraPosition_ = config.getParameter<Vec3>("cameraPosition");
   lookAtPoint_ = config.getParameter<Vec3>("lookAt");
   upVector_ = config.getParameter<Vec3>("upVector");
   
   blockAABBIntersectionPadding_ = config.getParameter<real_t>("blockAABBIntersectionPadding", real_t(0.0));

   setupView_();
}

/*!\brief Utility function for setting up the view plane and calculating required variables.
 */
void Raytracer::setupView_() {
   // eye coordinate system setup
   n_ = (cameraPosition_ - lookAtPoint_).getNormalized();
   u_ = (upVector_ % n_).getNormalized();
   v_ = n_ % u_;
   
   // viewing plane setup
   d_ = (cameraPosition_ - lookAtPoint_).length();
   aspectRatio_ = real_t(pixelsHorizontal_) / real_t(pixelsVertical_);
   viewingPlaneHeight_ = tan(deg2rad(fov_vertical_)/real_t(2.)) * real_t(2.) * d_;
   viewingPlaneWidth_ = viewingPlaneHeight_ * aspectRatio_;
   viewingPlaneOrigin_ = lookAtPoint_ - u_*viewingPlaneWidth_/real_t(2.) - v_*viewingPlaneHeight_/real_t(2.);
   
   pixelWidth_ = viewingPlaneWidth_ / real_c(pixelsHorizontal_);
   pixelHeight_ = viewingPlaneHeight_ / real_c(pixelsVertical_);
}

/*!\brief Writes the tBuffer to a file in the tBuffer output directory.
 * \param tBuffer Buffer with t values as generated in rayTrace(...).
 */
void Raytracer::writeTBufferToFile(const std::map<Coordinates, real_t, CoordinatesComparator>& tBuffer, const size_t timestep) const {
   mpi::MPIRank rank = mpi::MPIManager::instance()->rank();
   std::string fileName = "tbuffer_" + std::to_string(timestep) + "+" + std::to_string(rank) + ".ppm";
   writeTBufferToFile(tBuffer, fileName);
}

/*!\brief Writes the tBuffer to a file in the tBuffer output directory.
 * \param tBuffer Buffer with t values as generated in rayTrace(...).
 * \param fileName Name of the output file.
 */
void Raytracer::writeTBufferToFile(const std::map<Coordinates, real_t, CoordinatesComparator>& tBuffer, const std::string& fileName) const {
   namespace fs = boost::filesystem;
   
   real_t inf = std::numeric_limits<real_t>::max();
   
   real_t t_max = 1;
   real_t t_min = inf;
   for (size_t x = 0; x < pixelsHorizontal_; x++) {
      for (size_t y = 0; y < pixelsVertical_; y++) {
         Coordinates c = {x, y};
         real_t t = tBuffer.at(c);
         if (t > t_max && !realIsIdentical(t, inf)) {
            t_max = t;
         }
         if (t < t_min) {
            t_min = t;
         }
      }
   }
   if (realIsIdentical(t_min, inf)) t_min = 0;
   
   fs::path dir (getTBufferOutputDirectory());
   fs::path file (fileName);
   fs::path fullPath = dir / file;
   
   std::ofstream ofs(fullPath.string<std::string>(), std::ios::out | std::ios::binary);
   ofs << "P6\n" << pixelsHorizontal_ << " " << pixelsVertical_ << "\n255\n";
   for (size_t y = pixelsVertical_-1; y > 0; y--) {
      for (size_t x = 0; x < pixelsHorizontal_; x++) {
         Coordinates c = {x, y};
         char r = 0, g = 0, b = 0;
         real_t t = tBuffer.at(c);
         if (realIsIdentical(t, inf)) {
            r = g = b = (char)0;
         } else {
            r = g = b = (char)(240 * (1-(t-t_min)/(t_max-t_min)));
         }
         ofs << r << g << b;
      }
   }
   
   ofs.close();
}

/*!\brief Checks if a plane should get rendered.
 * \param plane Plane to check for visibility.
 * \param ray Ray which is intersected with plane.
 *
 * Checks if a plane should get rendered by comparing the planes normal and the ray direction.
 * If the rays direction vectors projection on the planes normal is positive, the plane is considered invisible.
 */
bool Raytracer::isPlaneVisible(const PlaneID plane, const Ray& ray) const {
   return plane->getNormal() * ray.getDirection() < 0;
}

}
}
}
