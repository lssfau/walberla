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

namespace walberla {
namespace pe {
namespace raytracing {

class Raytracer {
public:
   /*!\name Constructors */
   //@{
   explicit Raytracer(const shared_ptr<BlockStorage> & forest, BlockDataID storageID,
                      uint8_t pixelsHorizontal, uint8_t pixelsVertical,
                      real_t fov_vertical,
                      const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector);
   //@}

private:
   /*!\name Member variables */
   //@{
   const shared_ptr<BlockStorage>& forest_; //!< The BlockForest the raytracer operates on.
   BlockDataID storageID_;    /*!< The storage ID of the block data storage the raytracer operates
                               on.*/
   uint8_t pixelsHorizontal_; //!< The horizontal amount of pixels of the generated image.
   uint8_t pixelsVertical_;   //!< The vertical amount of pixels of the generated image.
   real_t fov_vertical_;      //!< The vertical field-of-view of the camera.
   Vec3 cameraPosition_;      //!< The position of the camera in the global world frame.
   Vec3 lookAtPoint_;         /*!< The point the camera looks at in the global world frame,
                               marks the center of the view plane.*/
   Vec3 upVector_;            //!< The vector indicating the upwards direction of the camera.
   //@}
   
   Vec3 n; // normal vector of viewing plane
   Vec3 u; // u and ...
   Vec3 v; // ... v span the viewing plane
   real_t d; // distance from camera to viewing plane
   real_t aspectRatio; // aspect ratio of the generated image and viewing plane
   real_t viewingPlaneHeight; // viewing plane height
   real_t viewingPlaneWidth; // viewing plane width
   Vec3 viewingPlaneOrigin; // origin of the viewing plane
   real_t pixelWidth; // width of a pixel of the generated image in the viewing plane
   real_t pixelHeight; // height of a pixel of the generated image in the viewing plane

public:
   /*!\name Get functions */
   //@{
   inline uint8_t       getPixelsHorizontal()   const;
   inline uint8_t       getPixelsVertical()     const;
   inline real_t        getFOVVertical()        const;
   inline const Vec3&   getCameraPosition()     const;
   inline const Vec3&   getLookAtPoint()        const;
   inline const Vec3&   getUpVector()           const;
   //@}
   
   /*!\name Functions */
   //@{
   template <typename BodyTypeTuple>
   void rayTrace(const size_t timestep)   const;
   //@}
};
   
/*!\brief Returns the horizontal amount of pixels of the generated image.
 *
 * \return The horizontal amount of pixels of the generated image.
 */
inline uint8_t Raytracer::getPixelsHorizontal() const {
   return pixelsHorizontal_;
}

/*!\brief Returns the vertical amount of pixels of the generated image.
 *
 * \return The vertical amount of pixels of the generated image.
 */
inline uint8_t Raytracer::getPixelsVertical() const {
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
}
}
}
