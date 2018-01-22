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
   setupView_();
}

/*!\brief Instantiation constructor for the Raytracer class using a config object for view setup.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Storage ID of the block data storage the raytracer operates on.
 * \param config Config block for the raytracer.
 *
 * The config block has to contain image_x (int), image_y (int) and fov_vertical (real, in degrees)
 * parameters, additionally one block with x, y and z values (real) for each of camera,
 * lookAt and the upVector.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage> & forest, BlockDataID storageID,
                     const Config::BlockHandle& config) : forest_(forest), storageID_(storageID) {
   pixelsHorizontal_ = config.getParameter<uint8_t>("image_x");
   pixelsVertical_ = config.getParameter<uint8_t>("image_y");
   fov_vertical_ = config.getParameter<real_t>("fov_vertical");
   
   WALBERLA_CHECK(config.isDefined("camera"), "No camera block found in config");
   const Config::BlockHandle cameraConf = config.getBlock("camera");
   cameraPosition_ = Vec3(cameraConf.getParameter<real_t>("x"),
                          cameraConf.getParameter<real_t>("y"),
                          cameraConf.getParameter<real_t>("z"));
   
   WALBERLA_CHECK(config.isDefined("lookAt"), "No lookAt block found in config");
   const Config::BlockHandle lookAtConf = config.getBlock("lookAt");
   lookAtPoint_ = Vec3(lookAtConf.getParameter<real_t>("x"),
                       lookAtConf.getParameter<real_t>("y"),
                       lookAtConf.getParameter<real_t>("z"));
   
   WALBERLA_CHECK(config.isDefined("upVector"), "No upVector block found in config");
   const Config::BlockHandle upVectorConf = config.getBlock("upVector");
   upVector_ = Vec3(upVectorConf.getParameter<real_t>("x"),
                    upVectorConf.getParameter<real_t>("y"),
                    upVectorConf.getParameter<real_t>("z"));
   
   setupView_();
}
   
void Raytracer::setupView_() {
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
   pixelHeight = viewingPlaneHeight / real_c(pixelsVertical_);
}

}
}
}
