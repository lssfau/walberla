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

#include <core/mpi/Gatherv.h>
#include <core/mpi/MPIManager.h>
#include <core/mpi/RecvBuffer.h>
#include <core/mpi/SendBuffer.h>

#include "lodepng.h"

namespace walberla {
namespace pe {
namespace raytracing {
   
void BodyIntersectionInfo_Comparator_MPI_OP( BodyIntersectionInfo *in, BodyIntersectionInfo *inout, const int *len, MPI_Datatype *dptr) {
   WALBERLA_UNUSED(dptr);
   for (int i = 0; i < *len; ++i) {
      if (in->bodySystemID != 0 && inout->bodySystemID != 0) {
         WALBERLA_ASSERT(in->imageX == inout->imageX && in->imageY == inout->imageY, "coordinates of infos do not match: " << in->imageX << "/" << in->imageY << " and " << inout->imageX << "/" << inout->imageY);
      }
      
      if ((in->t < inout->t && in->bodySystemID != 0) || (inout->bodySystemID == 0 && in->bodySystemID != 0)) {
         // info in "in" is closer than the one in "inout" -> update inout to values of in
         inout->imageX = in->imageX;
         inout->imageY = in->imageY;
         inout->bodySystemID = in->bodySystemID;
         inout->t = in->t;
         inout->r = in->r;
         inout->g = in->g;
         inout->b = in->b;
      }
      
      in++;
      inout++;
   }
}

/*!\brief Instantiation constructor for the Raytracer class.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Block data ID for the storage the raytracer operates on.
 * \param globalBodyStorage Pointer to the global body storage.
 * \param ccdID Block data ID for HashGrids.
 * \param pixelsHorizontal Horizontal amount of pixels of the generated image.
 * \param pixelsVertical Vertical amount of pixels of the generated image.
 * \param fov_vertical Vertical field-of-view of the camera.
 * \param cameraPosition Position of the camera in the global world frame.
 * \param lookAtPoint Point the camera looks at in the global world frame.
 * \param upVector Vector indicating the upwards direction of the camera.
 * \param backgroundColor Background color of the scene.
 * \param bodyToShadingParamsFunc A function mapping a BodyID to ShadingParameters for this body.
 *                                This can be used to customize the color and shading of bodies.
 * \param isBodyVisibleFunc A function which returns a boolean indicating if a given body should be visible
 *                          in the final image.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage>& forest, const BlockDataID storageID,
                     const shared_ptr<BodyStorage>& globalBodyStorage,
                     const BlockDataID ccdID,
                     uint16_t pixelsHorizontal, uint16_t pixelsVertical,
                     real_t fov_vertical, uint16_t antiAliasFactor,
                     const Vec3& cameraPosition, const Vec3& lookAtPoint, const Vec3& upVector,
                     const Lighting& lighting,
                     const Color& backgroundColor,
                     const std::function<ShadingParameters (const BodyID)>& bodyToShadingParamsFunc,
                     const std::function<bool (const BodyID)>& isBodyVisibleFunc)
   : forest_(forest), storageID_(storageID), globalBodyStorage_(globalBodyStorage), ccdID_(ccdID),
   pixelsHorizontal_(pixelsHorizontal), pixelsVertical_(pixelsVertical),
   fov_vertical_(fov_vertical), antiAliasFactor_(antiAliasFactor),
   cameraPosition_(cameraPosition), lookAtPoint_(lookAtPoint), upVector_(upVector),
   lighting_(lighting),
   backgroundColor_(backgroundColor),
   imageOutputEnabled_(true),
   localImageOutputEnabled_(false),
   imageOutputDirectory_("."),
   filenameTimestepWidth_(5),
   confinePlanesToDomain_(true),
   bodyToShadingParamsFunc_(bodyToShadingParamsFunc),
   isBodyVisibleFunc_(isBodyVisibleFunc),
   raytracingAlgorithm_(RAYTRACE_HASHGRIDS),
   reductionMethod_(MPI_REDUCE) {
   
   setupView_();
   setupFilenameRankWidth_();
   WALBERLA_MPI_SECTION() {
      setupMPI_();
   }
}

/*!\brief Instantiation constructor for the Raytracer class using a config object for view setup.
 *
 * \param forest BlockForest the raytracer operates on.
 * \param storageID Storage ID of the block data storage the raytracer operates on.
 * \param globalBodyStorage Pointer to the global body storage.
 * \param ccdID Block data ID for HashGrids.
 * \param config Config block for the raytracer.
 * \param bodyToShadingParamsFunc A function mapping a BodyID to ShadingParameters for this body.
 *                                This can be used to customize the color and shading of bodies.
 * \param isBodyVisibleFunc A function which returns a boolean indicating if a given body should be visible
 *                          in the final image.
 *
 * The config block has to contain image_x (int), image_y (int) and fov_vertical (real, in degrees).
 * Additionally a vector of reals for each of cameraPosition, lookAt and the upVector for the view setup are required.
 * Optional is antiAliasFactor (uint, usually between 1 and 4) for supersampling and backgroundColor (Vec3).
 * For image output after raytracing, set image_output_directory (string); for local image output additionally set
 * local_image_output_enabled (bool) to true. outputFilenameTimestepZeroPadding (int) sets the zero padding
 * for timesteps in output filenames.
 * For the lighting a config block within the Raytracer config block named Lighting has to be defined,
 * information about its contents is in the Lighting class.
 */
Raytracer::Raytracer(const shared_ptr<BlockStorage>& forest, const BlockDataID storageID,
                     const shared_ptr<BodyStorage>& globalBodyStorage,
                     const BlockDataID ccdID,
                     const Config::BlockHandle& config,
                     const std::function<ShadingParameters (const BodyID)>& bodyToShadingParamsFunc,
                     const std::function<bool (const BodyID)>& isBodyVisibleFunc)
   : forest_(forest), storageID_(storageID), globalBodyStorage_(globalBodyStorage), ccdID_(ccdID),
   bodyToShadingParamsFunc_(bodyToShadingParamsFunc),
   isBodyVisibleFunc_(isBodyVisibleFunc),
   raytracingAlgorithm_(RAYTRACE_HASHGRIDS),
   reductionMethod_(MPI_REDUCE) {
   WALBERLA_CHECK(config.isValid(), "No valid config passed to raytracer");
   
   pixelsHorizontal_ = config.getParameter<uint16_t>("image_x");
   pixelsVertical_ = config.getParameter<uint16_t>("image_y");
   fov_vertical_ = config.getParameter<real_t>("fov_vertical");
   antiAliasFactor_ = config.getParameter<uint16_t>("antiAliasFactor", 1);
   
   setLocalImageOutputEnabled(config.getParameter<bool>("local_image_output_enabled", false));
      
   if (config.isDefined("image_output_directory")) {
      setImageOutputEnabled(true);
      setImageOutputDirectory(config.getParameter<std::string>("image_output_directory", "."));
      WALBERLA_LOG_INFO_ON_ROOT("Images will be written to " << getImageOutputDirectory() << ".");
   } else if (getLocalImageOutputEnabled()) {
      WALBERLA_ABORT("Cannot enable local image output without image_output_directory parameter being set.");
   }
      
   filenameTimestepWidth_ = config.getParameter<uint8_t>("filenameTimestepWidth", uint8_t(5));
   confinePlanesToDomain_ = config.getParameter<bool>("confinePlanesToDomain", true);
   
   cameraPosition_ = config.getParameter<Vec3>("cameraPosition");
   lookAtPoint_ = config.getParameter<Vec3>("lookAt");
   upVector_ = config.getParameter<Vec3>("upVector");
   lighting_ = Lighting(config.getBlock("Lighting"));
   backgroundColor_ = config.getParameter<Color>("backgroundColor", Vec3(real_t(0.1), real_t(0.1), real_t(0.1)));

   std::string raytracingAlgorithm = config.getParameter<std::string>("raytracingAlgorithm", "RAYTRACE_HASHGRIDS");
   if (raytracingAlgorithm == "RAYTRACE_HASHGRIDS") {
      setRaytracingAlgorithm(RAYTRACE_HASHGRIDS);
   } else if (raytracingAlgorithm == "RAYTRACE_NAIVE") {
      setRaytracingAlgorithm(RAYTRACE_NAIVE);
   } else if (raytracingAlgorithm == "RAYTRACE_COMPARE_BOTH") {
      setRaytracingAlgorithm(RAYTRACE_COMPARE_BOTH);
   }
      
   std::string reductionMethod = config.getParameter<std::string>("reductionMethod", "MPI_REDUCE");
   if (reductionMethod == "MPI_REDUCE") {
      setReductionMethod(MPI_REDUCE);
   } else if (reductionMethod == "MPI_GATHER") {
      setReductionMethod(MPI_GATHER);
   }
      
   setupView_();
   setupFilenameRankWidth_();
   WALBERLA_MPI_SECTION() {
      setupMPI_();
   }
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
   real_t fov_vertical_rad = fov_vertical_ * math::pi / real_t(180.0);
   viewingPlaneHeight_ = real_c(tan(fov_vertical_rad/real_t(2.))) * real_t(2.) * d_;
   viewingPlaneWidth_ = viewingPlaneHeight_ * aspectRatio_;
   viewingPlaneOrigin_ = lookAtPoint_ - u_*viewingPlaneWidth_/real_t(2.) - v_*viewingPlaneHeight_/real_t(2.);
   
   pixelWidth_ = viewingPlaneWidth_ / real_c(pixelsHorizontal_*antiAliasFactor_);
   pixelHeight_ = viewingPlaneHeight_ / real_c(pixelsVertical_*antiAliasFactor_);
}

/*!\brief Utility function for initializing the attribute filenameRankWidth.
 */
void Raytracer::setupFilenameRankWidth_() {
   int numProcesses = mpi::MPIManager::instance()->numProcesses();
   filenameRankWidth_ = uint8_c(log10(numProcesses)+1);
}

/*!\brief Utility function for setting up the MPI datatype and operation.
 */
void Raytracer::setupMPI_() {
   MPI_Op_create((MPI_User_function *)BodyIntersectionInfo_Comparator_MPI_OP, true, &bodyIntersectionInfo_reduction_op);
   
   const int nblocks = 7;
   int blocklengths[nblocks] = {1,1,1,1,1,1,1};
   MPI_Datatype types[nblocks] = {
      MPI_UNSIGNED, // for coordinate
      MPI_UNSIGNED, // for coordinate
      MPI_UNSIGNED_LONG_LONG, // for id
      MPI_DOUBLE, // for distance
      MPI_DOUBLE, // for color
      MPI_DOUBLE, // for color
      MPI_DOUBLE // for color
   };
   MPI_Aint displacements[nblocks];
   displacements[0] = offsetof(BodyIntersectionInfo, imageX);
   displacements[1] = offsetof(BodyIntersectionInfo, imageY);
   displacements[2] = offsetof(BodyIntersectionInfo, bodySystemID);
   displacements[3] = offsetof(BodyIntersectionInfo, t);
   displacements[4] = offsetof(BodyIntersectionInfo, r);
   displacements[5] = offsetof(BodyIntersectionInfo, g);
   displacements[6] = offsetof(BodyIntersectionInfo, b);
   
   MPI_Datatype tmp_type;
   MPI_Type_create_struct(nblocks, blocklengths, displacements, types, &tmp_type);
   
   MPI_Aint lb;
   MPI_Aint extent;
   MPI_Type_get_extent( tmp_type, &lb, &extent );
   MPI_Type_create_resized( tmp_type, lb, extent, &bodyIntersectionInfo_mpi_type );
   
   MPI_Type_commit(&bodyIntersectionInfo_mpi_type);
}
   
/*!\brief Generates the filename for output files.
 * \param base String that precedes the timestep and rank info.
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
   WALBERLA_MPI_SECTION() {
      // Appending the rank to the filename only makes sense if actually using MPI.
      fileNameStream << "+";
      if (isGlobalImage) {
         fileNameStream << "global";
      } else {
         fileNameStream << std::setfill('0') << std::setw(int_c(filenameRankWidth_)) << std::to_string(rank); // add rank
      }
   }
   fileNameStream << ".png"; // add extension
   return fileNameStream.str();
}

/*!\brief Writes the image of the current intersection buffer to a file in the image output directory.
 * \param intersectionsBuffer Buffer with intersection info for each pixel.
 * \param timestep Timestep this image is from.
 * \param isGlobalImage Whether this image is the fully stitched together one.
 */
void Raytracer::writeImageToFile(const std::vector<BodyIntersectionInfo>& intersectionsBuffer, size_t timestep,
                                 bool isGlobalImage) const {
   writeImageToFile(intersectionsBuffer, getOutputFilename("image", timestep, isGlobalImage));
}

/*!\brief Writes the image of the current intersection buffer to a file in the image output directory.
 * \param intersectionsBuffer Buffer with intersection info for each pixel.
 * \param fileName Name of the output file.
 */
void Raytracer::writeImageToFile(const std::vector<BodyIntersectionInfo>& intersectionsBuffer,
                                 const std::string& fileName) const {
   filesystem::path dir = getImageOutputDirectory();
   filesystem::path file (fileName);
   filesystem::path fullPath = dir / file;
   
   std::vector<uint8_t> lodeImageBuffer(pixelsHorizontal_*pixelsVertical_*3);
   
   uint32_t l = 0;
   real_t patchSize = real_c(antiAliasFactor_*antiAliasFactor_);
   for (int y = pixelsVertical_-1; y >= 0; y--) {
      for (uint32_t x = 0; x < pixelsHorizontal_; x++) {
         real_t r_sum = 0;
         real_t g_sum = 0;
         real_t b_sum = 0;
         for (uint32_t ay = uint32_c(y)*antiAliasFactor_; ay < (uint32_c(y+1))*antiAliasFactor_; ay++) {
            for (uint32_t ax = x*antiAliasFactor_; ax < (x+1)*antiAliasFactor_; ax++) {
               size_t i = coordinateToArrayIndex(ax, ay);
               r_sum += real_c(intersectionsBuffer[i].r);
               g_sum += real_c(intersectionsBuffer[i].g);
               b_sum += real_c(intersectionsBuffer[i].b);
            }
         }
         uint8_t r = (uint8_t)(255 * (r_sum/patchSize));
         uint8_t g = (uint8_t)(255 * (g_sum/patchSize));
         uint8_t b = (uint8_t)(255 * (b_sum/patchSize));
         
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


/*!\brief Conflate the intersectionsBuffer of each process onto the root process using MPI_Reduce.
 * \param intersectionsBuffer Buffer containing all intersections for entire image (including non-hits).
 * \param tt Optional TimingTree.
 *
 * This function conflates the intersectionsBuffer of each process onto the root process using the MPI_Reduce
 * routine. It requires sending intersection info structs for the entire image instead of only the ones of the hits.
 *
 * \attention This function only works on MPI builds due to the explicit usage of MPI routines.
 */
void Raytracer::syncImageUsingMPIReduce(std::vector<BodyIntersectionInfo>& intersectionsBuffer, WcTimingTree* tt) {
   WALBERLA_NON_MPI_SECTION() {
      WALBERLA_UNUSED(intersectionsBuffer);
      WALBERLA_UNUSED(tt);
      WALBERLA_ABORT("Cannot call MPI reduce on a non-MPI build due to usage of MPI-specific code.");
   }
   
   WALBERLA_MPI_BARRIER();
   if (tt != nullptr) tt->start("Reduction");
   int rank = mpi::MPIManager::instance()->rank();

   const int recvRank = 0;
   if( rank == recvRank ) {
      MPI_Reduce(MPI_IN_PLACE,
                 &intersectionsBuffer[0], int_c(intersectionsBuffer.size()),
                 bodyIntersectionInfo_mpi_type, bodyIntersectionInfo_reduction_op,
                 recvRank, MPI_COMM_WORLD);
   } else {
      MPI_Reduce(&intersectionsBuffer[0], nullptr, int_c(intersectionsBuffer.size()),
                 bodyIntersectionInfo_mpi_type, bodyIntersectionInfo_reduction_op,
                 recvRank, MPI_COMM_WORLD);
   }
   
   WALBERLA_MPI_BARRIER();
   if (tt != nullptr) tt->stop("Reduction");
}
  
/*!\brief Conflate the intersectionsBuffer of each process onto the root process using MPI_Gather.
 * \param intersections Intersections to conflate.
 * \param intersectionsBuffer Buffer containing intersections.
 * \param tt Optional TimingTree.
 *
 * This function conflates the intersectionsBuffer of each process onto the root process using the MPI_Gather
 * routine. It only sends information for hits.
 */
void Raytracer::syncImageUsingMPIGather(std::vector<BodyIntersectionInfo>& intersections, std::vector<BodyIntersectionInfo>& intersectionsBuffer, WcTimingTree* tt) {
   WALBERLA_MPI_BARRIER();
   if (tt != nullptr) tt->start("Reduction");
   
   mpi::SendBuffer sendBuffer;
   for (auto& info: intersections) {
      sendBuffer << info.imageX << info.imageY
      << info.bodySystemID << info.t
      << info.r << info.g << info.b;
   }

   mpi::RecvBuffer recvBuffer;
   mpi::gathervBuffer(sendBuffer, recvBuffer, 0);
   
   WALBERLA_ROOT_SECTION() {
      BodyIntersectionInfo info;
      while (!recvBuffer.isEmpty()) {
         recvBuffer >> info.imageX;
         recvBuffer >> info.imageY;
         recvBuffer >> info.bodySystemID;
         recvBuffer >> info.t;
         recvBuffer >> info.r;
         recvBuffer >> info.g;
         recvBuffer >> info.b;
         
         size_t i = coordinateToArrayIndex(info.imageX, info.imageY);
         
         if (intersectionsBuffer[i].bodySystemID == 0 || info.t < intersectionsBuffer[i].t) {
            intersectionsBuffer[i] = info;
         }
      }
   }
   
   WALBERLA_MPI_BARRIER();
   if (tt != nullptr) tt->stop("Reduction");
}

void Raytracer::localOutput(const std::vector<BodyIntersectionInfo>& intersectionsBuffer, size_t timestep, WcTimingTree* tt) {
   if (getImageOutputEnabled()) {
      if (getLocalImageOutputEnabled()) {
         if (tt != nullptr) tt->start("Local Output");
         writeImageToFile(intersectionsBuffer, timestep);
         if (tt != nullptr) tt->stop("Local Output");
      }
   }
}

void Raytracer::output(const std::vector<BodyIntersectionInfo>& intersectionsBuffer, size_t timestep, WcTimingTree* tt) {
   if (tt != nullptr) tt->start("Output");
   WALBERLA_ROOT_SECTION() {
      if (getImageOutputEnabled()) {
         writeImageToFile(intersectionsBuffer, timestep, true);
      }
   }
   if (tt != nullptr) tt->stop("Output");
}

} //namespace raytracing
} //namespace pe
} //namespace walberla
