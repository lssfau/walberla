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
#include "geometry/structured/extern/lodepng.h"

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
 * \param backgroundColor Background color of the scene.
 * \param blockAABBIntersectionPadding The padding applied in block AABB intersection pretesting. Usually not required.
 *                                     Set it to the value of the farthest distance a object might protrude from
 *                                     its containing block.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage> forest, const BlockDataID storageID,
                     const shared_ptr<BodyStorage> globalBodyStorage,
                     const BlockDataID ccdID,
                     uint16_t pixelsHorizontal, uint16_t pixelsVertical,
                     real_t fov_vertical,
                     const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                     const Lighting& lighting,
                     const Color& backgroundColor,
                     real_t blockAABBIntersectionPadding,
                     std::function<ShadingParameters (const BodyID)> bodyToShadingParamsFunction)
   : forest_(forest), storageID_(storageID), globalBodyStorage_(globalBodyStorage), ccdID_(ccdID),
   pixelsHorizontal_(pixelsHorizontal), pixelsVertical_(pixelsVertical),
   fov_vertical_(fov_vertical),
   cameraPosition_(cameraPosition), lookAtPoint_(lookAtPoint), upVector_(upVector),
   lighting_(lighting),
   backgroundColor_(backgroundColor),
   blockAABBIntersectionPadding_(blockAABBIntersectionPadding),
   tBufferOutputEnabled_(false),
   imageOutputEnabled_(false),
   localImageOutputEnabled_(false),
   filenameTimestepWidth_(5),
   bodyToShadingParamsFunction_(bodyToShadingParamsFunction) {
   
   setupView_();
   setupFilenameRankWidth_();
}

/*!\brief Instantiation constructor for the Raytracer class using a config object for view setup.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Storage ID of the block data storage the raytracer operates on.
 * \param config Config block for the raytracer.
 *
 * The config block has to contain image_x (int), image_y (int), fov_vertical (real, in degrees)
 * and tbuffer_output_directory (string) parameters. Additionally a vector of reals
 * for each of cameraPosition, lookAt and the upVector. Optional is blockAABBIntersectionPadding (real) and backgroundColor (Vec3).
 * To output both process local and global tbuffers after raytracing, set tbuffer_output_directory (string).
 * For image output after raytracing, set image_output_directory (string); for local image output additionally set
 * local_image_output_enabled (bool) to true. outputFilenameTimestepZeroPadding (int) sets zero padding for timesteps of output filenames.
 * For the lighting a config block named Lighting has to be defined, information about its contents is in Lighting.h.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage> forest, const BlockDataID storageID,
                     const shared_ptr<BodyStorage> globalBodyStorage,
                     const BlockDataID ccdID,
                     const Config::BlockHandle& config,
                     std::function<ShadingParameters (const BodyID)> bodyToShadingParamsFunction)
   : forest_(forest), storageID_(storageID), globalBodyStorage_(globalBodyStorage), ccdID_(ccdID),
   bodyToShadingParamsFunction_(bodyToShadingParamsFunction) {
   WALBERLA_CHECK(config.isValid(), "No valid config passed to raytracer");
   
   pixelsHorizontal_ = config.getParameter<uint16_t>("image_x");
   pixelsVertical_ = config.getParameter<uint16_t>("image_y");
   fov_vertical_ = config.getParameter<real_t>("fov_vertical");
   
   if (config.isDefined("tbuffer_output_directory")) {
      setTBufferOutputEnabled(true);
      setTBufferOutputDirectory(config.getParameter<std::string>("tbuffer_output_directory"));
      WALBERLA_LOG_INFO_ON_ROOT("t buffers will be written to " << getTBufferOutputDirectory() << ".");
   }
   
   setLocalImageOutputEnabled(config.getParameter<bool>("local_image_output_enabled", false));
      
   if (config.isDefined("image_output_directory")) {
      setImageOutputEnabled(true);
      setImageOutputDirectory(config.getParameter<std::string>("image_output_directory"));
      WALBERLA_LOG_INFO_ON_ROOT("Images will be written to " << getImageOutputDirectory() << ".");
   } else if (getLocalImageOutputEnabled()) {
      WALBERLA_ABORT("Cannot enable local image output without image_output_directory parameter being set.");
   }
      
   filenameTimestepWidth_ = config.getParameter<uint8_t>("filenameTimestepWidth", uint8_t(5));
   
   cameraPosition_ = config.getParameter<Vec3>("cameraPosition");
   lookAtPoint_ = config.getParameter<Vec3>("lookAt");
   upVector_ = config.getParameter<Vec3>("upVector");
   lighting_ = Lighting(config.getBlock("Lighting"));
   backgroundColor_ = config.getParameter<Color>("backgroundColor", Vec3(0.1, 0.1, 0.1));

   blockAABBIntersectionPadding_ = config.getParameter<real_t>("blockAABBIntersectionPadding", real_t(0.0));
      
   setupView_();
   setupFilenameRankWidth_();
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

/*!\brief Utility function for initializing the attribute filenameRankWidth.
 */
void Raytracer::setupFilenameRankWidth_() {
   int numProcesses = mpi::MPIManager::instance()->numProcesses();
   filenameRankWidth_ = uint8_c(log10(numProcesses))+1;
}

/*!\brief Generates the filename for output files.
 * \param base String that precedes the timestap and rank info.
 * \param timestep Timestep this image is from.
 * \param isGlobalImage Whether this image is the fully stitched together one.
 */
std::string Raytracer::getOutputFilename(const std::string& base, size_t timestep, bool isGlobalImage) const {
   uint_t maxTimestep = uint_c(pow(10, filenameTimestepWidth_));
   WALBERLA_CHECK(timestep < maxTimestep, "Raytracer only supports outputting " << (maxTimestep-1) << " timesteps for the configured filename timestep width.");
   mpi::MPIRank rank = mpi::MPIManager::instance()->rank();
   std::stringstream fileNameStream;
   fileNameStream << base << "_";
   fileNameStream << std::setfill('0') << std::setw(int_c(filenameTimestepWidth_)) << timestep; // add timestep
   fileNameStream << "+";
   if (isGlobalImage) {
      fileNameStream << "global";
   } else {
      fileNameStream << std::setfill('0') << std::setw(int_c(filenameRankWidth_)) << std::to_string(rank); // add rank
   }
   fileNameStream << ".png"; // add extension
   return fileNameStream.str();
}

/*!\brief Writes the tBuffer to a file in the tBuffer output directory.
 * \param tBuffer Buffer with t values as generated in rayTrace(...).
 * \param timestep Timestep this image is from.
 * \param isGlobalImage Whether this image is the fully stitched together one.
 */
void Raytracer::writeTBufferToFile(const std::vector<real_t>& tBuffer, size_t timestep, bool isGlobalImage) const {
   writeTBufferToFile(tBuffer, getOutputFilename("tbuffer", timestep, isGlobalImage));
}

/*!\brief Writes the tBuffer to a file in the tBuffer output directory.
 * \param tBuffer Buffer with t values as generated in rayTrace(...).
 * \param fileName Name of the output file.
 */
void Raytracer::writeTBufferToFile(const std::vector<real_t>& tBuffer, const std::string& fileName) const {
   namespace fs = boost::filesystem;
   
   real_t inf = std::numeric_limits<real_t>::max();
   
   real_t t_max = 1;
   real_t t_min = inf;
   for (size_t x = 0; x < pixelsHorizontal_; x++) {
      for (size_t y = 0; y < pixelsVertical_; y++) {
         size_t i = coordinateToArrayIndex(x, y);
         real_t t = tBuffer[i];
         //if (t > t_max && !realIsIdentical(t, inf)) {
         //   t_max = t;
         //}
         if (t < t_min) {
            t_min = t;
         }
      }
   }
   if (realIsIdentical(t_min, inf)) t_min = 0;
   
   t_max = forest_->getDomain().maxDistance(cameraPosition_);
   
   fs::path dir (getTBufferOutputDirectory());
   fs::path file (fileName);
   fs::path fullPath = dir / file;
   
   std::vector<u_char> lodeTBuffer(pixelsHorizontal_*pixelsVertical_);
   
   uint32_t l = 0;
   for (size_t y = pixelsVertical_-1; y > 0; y--) {
      for (size_t x = 0; x < pixelsHorizontal_; x++) {
         size_t i = coordinateToArrayIndex(x, y);
         u_char g = 0;
         real_t t = tBuffer[i];
         if (realIsIdentical(t, inf)) {
            g = (u_char)0;
         } else {
            real_t t_scaled = (1-(t-t_min)/(t_max-t_min));
            g = (u_char)(255 * std::max(std::min(t_scaled, real_t(1)), real_t(0)));
         }
         lodeTBuffer[l] = g;
         l++;
      }
   }
   
   uint32_t error = lodepng::encode(fullPath.string(), lodeTBuffer, getPixelsHorizontal(), getPixelsVertical(), LCT_GREY);
   if(error) {
      WALBERLA_LOG_WARNING("lodePNG error " << error << " when trying to save tbuffer file to " << fullPath.string() << ": " << lodepng_error_text(error));
   }
}

/*!\brief Writes a image buffer to a file in the image output directory.
 * \param imageBuffer Buffer with color vectors.
 * \param timestep Timestep this image is from.
 * \param isGlobalImage Whether this image is the fully stitched together one.
 */
void Raytracer::writeImageBufferToFile(const std::vector<Color>& imageBuffer, size_t timestep, bool isGlobalImage) const {
   writeImageBufferToFile(imageBuffer, getOutputFilename("image", timestep, isGlobalImage));
}

/*!\brief Writes the image buffer to a file in the image output directory.
 * \param imageBuffer Buffer with color vectors.
 * \param fileName Name of the output file.
 */
void Raytracer::writeImageBufferToFile(const std::vector<Color>& imageBuffer, const std::string& fileName) const {
   namespace fs = boost::filesystem;

   fs::path dir (getImageOutputDirectory());
   fs::path file (fileName);
   fs::path fullPath = dir / file;
   
   std::vector<u_char> lodeImageBuffer(pixelsHorizontal_*pixelsVertical_*3);
   
   uint32_t l = 0;
   for (size_t y = pixelsVertical_-1; y > 0; y--) {
      for (size_t x = 0; x < pixelsHorizontal_; x++) {
         size_t i = coordinateToArrayIndex(x, y);
         const Color& color = imageBuffer[i];
         u_char r = (u_char)(255 * color[0]);
         u_char g = (u_char)(255 * color[1]);
         u_char b = (u_char)(255 * color[2]);
         lodeImageBuffer[l] = r;
         lodeImageBuffer[l+1] = g;
         lodeImageBuffer[l+2] = b;
         l+=3;
      }
   }
   
   uint32_t error = lodepng::encode(fullPath.string(), lodeImageBuffer, getPixelsHorizontal(), getPixelsVertical(), LCT_RGB);
   if(error) {
      WALBERLA_LOG_WARNING("lodePNG error " << error << " when trying to save image file to " << fullPath.string() << ": " << lodepng_error_text(error));
   }
}
}
}
}
