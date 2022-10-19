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
//! \file   ShapeGeneration.h
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/config/Config.h"
#include "core/Filesystem.h"

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"
#include "mesa_pd/data/shape/Ellipsoid.h"
#include "mesa_pd/data/shape/ConvexPolyhedron.h"

#include "mesh_common/TriangleMeshes.h"
#include "mesh_common/MeshIO.h"
#include "mesh_common/MeshOperations.h"
#include "mesh_common/PolyMeshes.h"

#include "Utility.h"

#include <random>

namespace walberla {
namespace mesa_pd {

Vec3 getSemiAxesFromMesh(mesh::TriangleMesh& mesh)
{
   auto volume = mesh::computeVolume(mesh);
   mesh::translate(mesh, mesh::toWalberla(-mesh::computeCentroid(mesh)));
   auto inertiaTensor = mesh::computeInertiaTensor(mesh);
   auto mass = volume * 1_r; // mesh has unit mass with density 1
   return semiAxesFromInertiaTensor(inertiaTensor, mass);
}

std::vector<Vec3> extractSemiAxesFromMeshFiles(const std::vector<std::string> & meshFiles)
{
   std::vector<Vec3> semiAxesVector;

   for(const auto& meshFileName : meshFiles)
   {
      mesh::TriangleMesh mesh;
      mesh::readAndBroadcast<mesh::TriangleMesh>(meshFileName, mesh);
      auto volume = mesh::computeVolume(mesh);

      auto semiAxes = getSemiAxesFromMesh(mesh);

      WALBERLA_LOG_INFO_ON_ROOT("Read mesh file: " << meshFileName << " (" << mesh.n_vertices() << " vertices, "
                                                   << mesh.n_faces()
                                                   << " faces, volume = " << volume << ", semi axes = " << semiAxes << ")");

      semiAxesVector.push_back(semiAxes);
   }

   return semiAxesVector;
}


std::vector<std::string> getMeshFilesFromPath(const std::string & meshPath)
{
   std::vector<std::string> meshNames;

   if(filesystem::is_directory(filesystem::path(meshPath)))
   {
      // assuming this folder contains the mesh files
      for(const auto& entry : filesystem::directory_iterator(meshPath)) {
         std::string meshFileName = entry.path();
         if(meshFileName.find(".obj") == std::string::npos)
         {
            // open mesh seemingly can only read .obj files reliably, so skip all others
            WALBERLA_LOG_WARNING_ON_ROOT("Ignoring mesh file " << meshFileName << ", since not of .obj file type.");
            continue;
         }
         meshNames.push_back(meshFileName);
      }
   } else {
      // assuming path is a single (mesh) file
      meshNames.push_back(meshPath);
   }

   return meshNames;
}

enum ScaleMode { sphereEquivalent, sieveLike };

ScaleMode str_to_scaleMode(const std::string & str)
{
   if(str == "sphereEquivalent") return ScaleMode::sphereEquivalent;
   else if(str == "sieveLike") return ScaleMode::sieveLike;
   else WALBERLA_ABORT("Unknown shape scale mode " << str);
}

using FormParameters = Vector3<real_t>; // contains L, I, S

class NormalizedFormGenerator
{
public:
   virtual FormParameters get() = 0;
   virtual FormParameters getNormalFormParameters() const = 0;
   virtual real_t getMaxDiameterScalingFactor() const = 0; // used to estimate maximum particle extent
   virtual real_t getNormalVolume() const = 0;  // volume of a typical particle with diameter 1
   virtual bool generatesSingleShape() const = 0; // generate a single shape or multiple shapes?
   virtual ~NormalizedFormGenerator() = default;
};

class ConstFormGenerator : public NormalizedFormGenerator
{
public:
   explicit ConstFormGenerator(){}
   FormParameters get() override { return FormParameters(1_r, 1_r, 1_r); }
   FormParameters getNormalFormParameters() const override { return FormParameters(1_r, 1_r, 1_r); }
   real_t getMaxDiameterScalingFactor() const override { return 1_r; }
   real_t getNormalVolume() const override { return 1_r; }
   bool generatesSingleShape() const override { return true; }
};

class SampleFormGenerator : public NormalizedFormGenerator
{
public:
   SampleFormGenerator(const std::vector<Vec3> & semiAxesVector,
                       ScaleMode scaleMode) :
         gen_(static_cast<unsigned long>(walberla::mpi::MPIManager::instance()->rank()))
   {

      maxDiameterScalingFactor_ = 1_r;
      normalVolume_ = 0_r;
      for(const auto& semiAxes : semiAxesVector)
      {
         FormParameters normalizedFormParameters;
         switch(scaleMode)
         {
            case ScaleMode::sphereEquivalent:
            {
               real_t ellipsoidScalingFactor = std::cbrt( 1_r / ( 8_r * semiAxes[0] * semiAxes[1] * semiAxes[2] ) ); // = V_sphere / V_ellipsoid
               normalizedFormParameters = 2_r * semiAxes * ellipsoidScalingFactor;
               break;
            }
            case ScaleMode::sieveLike:
            {
               auto scalingFactor = sizeFromSemiAxes(semiAxes);
               normalizedFormParameters = 2_r * semiAxes / scalingFactor;
               break;
            }
            default:
               WALBERLA_ABORT("Unknown shape scale mode");
         }

         auto volume = 1_r / 6_r * math::pi * normalizedFormParameters[0] * normalizedFormParameters[1] * normalizedFormParameters[2];

         maxDiameterScalingFactor_ = std::max(maxDiameterScalingFactor_, normalizedFormParameters.max());
         normalVolume_ += volume;

         normalizedFormParametersVector_.push_back(normalizedFormParameters);
      }

      WALBERLA_CHECK(!normalizedFormParametersVector_.empty(), "No shape samples provided");

      auto numFormParameters = normalizedFormParametersVector_.size();
      normalVolume_ /= real_c(numFormParameters); // use average normal volume here as estimation

      normalFormParameters_ = std::accumulate(normalizedFormParametersVector_.begin(),
                                              normalizedFormParametersVector_.end(),
                                              FormParameters(0_r)) / real_t(numFormParameters); // = average of form parameters

      dist_ = std::uniform_int_distribution<uint_t>(0, numFormParameters-1);
   }

   FormParameters get() override
   {
      return normalizedFormParametersVector_[dist_(gen_)];
   }

   FormParameters getNormalFormParameters() const override { return normalFormParameters_; }
   real_t getMaxDiameterScalingFactor() const override { return maxDiameterScalingFactor_; }
   real_t getNormalVolume() const override { return normalVolume_; }
   bool generatesSingleShape() const override { return normalizedFormParametersVector_.size() == 1; }

private:
   std::vector<FormParameters> normalizedFormParametersVector_;
   std::mt19937 gen_;
   std::uniform_int_distribution<uint_t> dist_;

   real_t maxDiameterScalingFactor_;
   real_t normalVolume_;
   FormParameters normalFormParameters_;

};

// assumes individual normal distribution of form factors elongation (= I/L) and flatness (= S/I)
// NOTE: since S <= I <= L, the combination elongation & equancy would not work trivially as it follows: equancy <= elongation,
// these two parameters are not be completely independent
// limits both values to interval (eps,1]
// -> best way would be to use a truncated Normal distribution https://en.wikipedia.org/wiki/Truncated_normal_distribution
// however, formula is more complicated and no library for sampling is available
class DistributionFormGenerator : public NormalizedFormGenerator
{
public:
   DistributionFormGenerator( real_t elongationMean, real_t elongationStdDev,
                              real_t flatnessMean, real_t flatnessStdDev,
                              ScaleMode scaleMode) :
         scaleMode_(scaleMode), gen_(static_cast<unsigned long>(walberla::mpi::MPIManager::instance()->rank()))
   {
      elongationDist_ = std::normal_distribution<real_t>(elongationMean, elongationStdDev);
      flatnessDist_ = std::normal_distribution<real_t>(flatnessMean, flatnessStdDev);

      normalFormParameters_ = getNormalizedFormParameters(elongationMean, flatnessMean);
      normalVolume_ = 1_r / 6_r * math::pi * normalFormParameters_[0] * normalFormParameters_[1] * normalFormParameters_[2];

      WALBERLA_LOG_INFO_ON_ROOT("Shape generation from distribution with mean size-normalized form parameters (L,I,S) = " << normalFormParameters_);

      real_t extremeElongation = std::max(eps_, elongationMean - 5_r * elongationStdDev); // covers almost all values
      real_t extremeFlatness = std::max(eps_, flatnessMean - 5_r * flatnessStdDev);
      auto extremeFormParams = getNormalizedFormParameters(extremeElongation, extremeFlatness);
      maxDiameterScalingFactor_ = extremeFormParams.max();

   }

   FormParameters get() override
   {
      // variant: re-roll until within bounds
      real_t elongation = elongationDist_(gen_);
      while(elongation < eps_ || elongation > 1_r) elongation = elongationDist_(gen_);
      real_t flatness = flatnessDist_(gen_);
      while(flatness < eps_ || flatness > 1_r) flatness = flatnessDist_(gen_);

      // variant: cap values to min/max
      //real_t elongation = std::min(std::max(eps_, elongationFromDist), 1_r);
      //real_t flatness = std::min(std::max(eps_, flatnessFromDist), 1_r);

      return getNormalizedFormParameters(elongation, flatness);
   }

   FormParameters getNormalFormParameters() const override { return normalFormParameters_; }
   real_t getNormalVolume() const override { return normalVolume_; }
   real_t getMaxDiameterScalingFactor() const override { return maxDiameterScalingFactor_; }
   bool generatesSingleShape() const override { return false; }

private:

   FormParameters getNormalizedFormParameters(real_t elongation, real_t flatness)
   {
      real_t I = 1_r;
      switch(scaleMode_)
      {
         case ScaleMode::sieveLike:
            I = std::sqrt(2_r) / std::sqrt(1_r + flatness * flatness); break;
         case ScaleMode::sphereEquivalent:
            I = std::cbrt(elongation / flatness); break;
      }
      return FormParameters( I / elongation, I, I * flatness);
   }

   ScaleMode scaleMode_;
   std::mt19937 gen_;
   std::normal_distribution<real_t> elongationDist_;
   std::normal_distribution<real_t> flatnessDist_;
   real_t normalVolume_;
   real_t maxDiameterScalingFactor_;
   FormParameters normalFormParameters_;

   const real_t eps_ = 0.1_r;
};


class ShapeGenerator
{
public:
   virtual void setShape(real_t diameter, real_t maximumAllowedInteractionRadius, data::ParticleStorage::baseShape_type& shape, real_t& interactionRadius) = 0;
   virtual real_t getMaxDiameterScalingFactor() = 0; // used to estimate maximum particle extent
   virtual real_t getNormalVolume() = 0;  // volume of a typical particle with diameter 1
   virtual FormParameters getNormalFormParameters() = 0; // form (L,I,S) of a typical particle with diameter 1
   virtual bool generatesSingleShape() = 0; // generate a single shape or multiple shapes?
   virtual ~ShapeGenerator() = default;
};

class SphereGenerator : public ShapeGenerator
{
public:
   explicit SphereGenerator()
   {
      maxDiameterScalingFactor_ = 1_r;
      normalVolume_ = math::pi / 6_r;
      WALBERLA_LOG_INFO_ON_ROOT("Will create spheres");
   }

   void setShape(real_t diameter, real_t maximumAllowedInteractionRadius, data::ParticleStorage::baseShape_type& shape, real_t& interactionRadius) override
   {
      real_t radius = diameter * 0.5_r;
      if(radius > maximumAllowedInteractionRadius) radius = maximumAllowedInteractionRadius;
      shape = std::make_shared<data::Sphere>(radius);
      interactionRadius = radius;
   }

   real_t getMaxDiameterScalingFactor() override { return maxDiameterScalingFactor_; }
   real_t getNormalVolume() override { return normalVolume_; }
   FormParameters getNormalFormParameters() override {return FormParameters(1_r);}
   bool generatesSingleShape() override { return true; }

private:
   real_t maxDiameterScalingFactor_;
   real_t normalVolume_;
};


class EllipsoidGenerator : public ShapeGenerator
{
public:
   EllipsoidGenerator(std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator) :
         normalizedFormGenerator_(std::move(normalizedFormGenerator)){}

   void setShape(real_t diameter, real_t maximumAllowedInteractionRadius, data::ParticleStorage::baseShape_type& shape, real_t& interactionRadius) override
   {
      Vector3<real_t> semiAxes = 0.5_r * normalizedFormGenerator_->get() * diameter;
      sortVector(semiAxes); // we want to have particles that are longest in z-direction

      real_t maxAxis = semiAxes.max();
      if(maxAxis > maximumAllowedInteractionRadius)
      {
         semiAxes *= (maximumAllowedInteractionRadius / maxAxis);
      }

      shape = std::make_shared<data::Ellipsoid>(semiAxes);
      interactionRadius = semiAxes.max();
   }

   real_t getMaxDiameterScalingFactor() override { return normalizedFormGenerator_->getMaxDiameterScalingFactor(); }
   real_t getNormalVolume() override { return normalizedFormGenerator_->getNormalVolume(); }
   FormParameters getNormalFormParameters() override {return normalizedFormGenerator_->getNormalFormParameters();}
   bool generatesSingleShape() override { return normalizedFormGenerator_->generatesSingleShape(); }

private:
   std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator_;
};


class MeshesGenerator : public ShapeGenerator
{
public:
   MeshesGenerator(const std::vector<std::string> & meshFiles, ScaleMode scaleMode, std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator) :
         normalizedFormGenerator_(std::move(normalizedFormGenerator)), gen_(static_cast<unsigned long>(walberla::mpi::MPIManager::instance()->rank()))
   {
      WALBERLA_CHECK(!meshFiles.empty());

      maxDiameterScalingFactor_ = 1_r;
      normalVolume_ = 0_r;
      for(const auto& meshFileName : meshFiles)
      {
         mesh::TriangleMesh mesh;
         mesh::readAndBroadcast<mesh::TriangleMesh>(meshFileName, mesh);
         mesh::translate(mesh, mesh::toWalberla(-mesh::computeCentroid(mesh)));

         if(normalizedFormGenerator_->generatesSingleShape())
         {
            WALBERLA_LOG_INFO_ON_ROOT("Read mesh file: " << meshFileName << " (" << mesh.n_vertices() << " vertices, "
                                                         << mesh.n_faces()
                                                         << " faces, volume = " << mesh::computeVolume(mesh) << ")");

            // all shape scaling is handled as given by the meshes and thus by this class internally
            switch( scaleMode )
            {
               case ScaleMode::sphereEquivalent:
               {
                  real_t meshScalingFactor = std::cbrt(math::pi / (6_r * mesh::computeVolume(mesh)));
                  mesh::scale(mesh, Vec3(meshScalingFactor));
                  break;
               }
               case ScaleMode::sieveLike:
               {
                  WALBERLA_LOG_INFO_ON_ROOT("Using sieve-like scaling! Assuming that the mesh is scaled such that a scaling with the mesh size results in the correct size!");
                  break;
               }
               default: WALBERLA_ABORT("Unknown shape scale mode");
            }

            maxDiameterScalingFactor_ = std::max( maxDiameterScalingFactor_, 2_r * mesh::computeBoundingSphereRadius(mesh, mesh::computeCentroid(mesh)) );
            normalVolume_ += mesh::computeVolume(mesh);

         } else
         {
            // meshes are scaled via form factors provided from outside
            // -> scale mesh such that the semi axes are 0.5, i.e. removing all shape features from it

            auto semiAxes = getSemiAxesFromMesh(mesh);
            auto meshScalingFactor = Vec3(1_r / (2_r * semiAxes[0]), 1_r / (2_r * semiAxes[1]), 1_r / (2_r * semiAxes[2]));
            mesh::scale(mesh, meshScalingFactor);

            maxDiameterScalingFactor_ = std::max( maxDiameterScalingFactor_,
                                                  2_r * mesh::computeBoundingSphereRadius(mesh, mesh::computeCentroid(mesh)) * normalizedFormGenerator_->getMaxDiameterScalingFactor() );
            auto meshCopy = mesh;
            auto normalFormParameters = normalizedFormGenerator_->getNormalFormParameters();
            sortVector(normalFormParameters);
            mesh::scale(meshCopy, normalFormParameters);
            normalVolume_ += mesh::computeVolume(meshCopy);
         }

         particleMeshes_.push_back(mesh);

      }

      WALBERLA_CHECK(!particleMeshes_.empty());

      normalVolume_ /= real_c(particleMeshes_.size()); // use average normal mesh volume here as estimation

      dist_ = std::uniform_int_distribution<uint_t>(0, particleMeshes_.size()-1);

      WALBERLA_LOG_INFO_ON_ROOT("Read and stored in total " << particleMeshes_.size() << " meshes, and will randomly pick one for each generated particle.")
   }

   void setShape(real_t diameter, real_t maximumAllowedInteractionRadius, data::ParticleStorage::baseShape_type& shape, real_t& interactionRadius) override
   {
      auto meshCopy = particleMeshes_[dist_(gen_)];
      auto normalizedFormParameters = normalizedFormGenerator_->get();
      sortVector(normalizedFormParameters); // we want to have particles that are longest in z-direction

      mesh::scale(meshCopy, normalizedFormParameters * diameter);
      auto convexPolyPtr = std::make_shared<data::ConvexPolyhedron>(meshCopy);
      convexPolyPtr->updateMeshQuantities();
      real_t boundingRadius = convexPolyPtr->getBoundingSphereRadius();

      if(boundingRadius > maximumAllowedInteractionRadius)
      {
         // scale to match limiting radius
         mesh::scale(meshCopy, Vec3(maximumAllowedInteractionRadius / boundingRadius) );
         convexPolyPtr = std::make_shared<data::ConvexPolyhedron>(meshCopy);
         convexPolyPtr->updateMeshQuantities();
      }
      shape = convexPolyPtr;
      interactionRadius = convexPolyPtr->getBoundingSphereRadius();
   }

   real_t getMaxDiameterScalingFactor() override { return maxDiameterScalingFactor_; }
   real_t getNormalVolume() override { return normalVolume_; }
   FormParameters getNormalFormParameters() override {return normalizedFormGenerator_->getNormalFormParameters();}
   bool generatesSingleShape() override { return particleMeshes_.size() == 1 && normalizedFormGenerator_->generatesSingleShape(); }

private:
   std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator_;

   std::vector<mesh::TriangleMesh> particleMeshes_;
   std::mt19937 gen_;
   std::uniform_int_distribution<uint_t> dist_;

   real_t maxDiameterScalingFactor_;
   real_t normalVolume_;
};

class UnscaledMeshesPerFractionGenerator : public ShapeGenerator
{
public:
   UnscaledMeshesPerFractionGenerator(const Config::BlockHandle & shapeConfig, const std::vector<real_t> & massFractions) :
      gen_(static_cast<unsigned long>(walberla::mpi::MPIManager::instance()->rank()))
   {

      auto meshesConfig = shapeConfig.getBlock("UnscaledMeshesPerFraction");
      std::string meshesTopFolder = meshesConfig.getParameter<std::string>("folder");

      std::vector<real_t> avgVolumesPerFraction(massFractions.size(), 0_r);

      maxDiameterScalingFactor_ = 1_r;
      generatesSingleShape_ = true;

      auto meshesTopFolderPath = filesystem::path(meshesTopFolder);
      for(uint_t fractionIdx = 0; fractionIdx < massFractions.size(); ++fractionIdx)
      {
         auto meshesFolder = meshesTopFolderPath / std::to_string(fractionIdx);

         if(!filesystem::exists(meshesFolder))
         {
            WALBERLA_ABORT("Path " << meshesFolder.string() << " expected but does not exist.");
         }

         std::vector<mesh::TriangleMesh> meshesVector;

         for (const auto &entry : filesystem::directory_iterator(meshesFolder)) {
            std::string meshFileName = entry.path();
            if (meshFileName.find(".obj") == std::string::npos) {
               // open mesh seemingly can only read .obj files reliably, so skip all others
               continue;
            }
            mesh::TriangleMesh mesh;
            mesh::readAndBroadcast<mesh::TriangleMesh>(meshFileName, mesh);

            WALBERLA_LOG_INFO_ON_ROOT("Read mesh file: " << meshFileName << " (" << mesh.n_vertices() << " vertices, "
                                                         << mesh.n_faces()
                                                         << " faces, volume = " << mesh::computeVolume(mesh) << ")");
            mesh::translate(mesh, mesh::toWalberla(-mesh::computeCentroid(mesh)));
            meshesVector.push_back(mesh);

            avgVolumesPerFraction[fractionIdx] += mesh::computeVolume(mesh);

            /*
            auto meshAABB = mesh::computeAABB(mesh);
            auto aabbHeight = meshAABB.zSize();
            auto aabbHorizontalDiagonal = std::sqrt(meshAABB.xSize()*meshAABB.xSize() + meshAABB.ySize()*meshAABB.ySize());
            maxDiameterScalingFactor_ = std::max(maxDiameterScalingFactor_, aabbHeight / aabbHorizontalDiagonal);
            */
            maxDiameterScalingFactor_ = std::max( maxDiameterScalingFactor_, 2_r * mesh::computeBoundingSphereRadius(mesh, mesh::computeCentroid(mesh)) );
         }

         if(meshesVector.size() != 1) generatesSingleShape_ = false;

         WALBERLA_CHECK(!meshesVector.empty(), "No meshes found in folder " << meshesFolder << ". Provide at least one per folder.");
         avgVolumesPerFraction[fractionIdx]  /= real_c(meshesVector.size()); // use average normal mesh volume here as estimation
         distsPerFraction_.emplace_back(std::uniform_int_distribution<uint_t>(0, meshesVector.size()-1));
         particleMeshPerFractionVector_.emplace_back(meshesVector);
      }

      auto particleNumbers = transferMassFractionsToParticleNumbersFromAvgVolumes(massFractions, avgVolumesPerFraction);
      auto totalParticles = std::accumulate(particleNumbers.begin(), particleNumbers.end(), real_t(0));
      std::string outString = "Particle probabilities per size fraction: ";
      for(const auto & p : particleNumbers) outString += std::to_string(p/totalParticles) + " | ";
      WALBERLA_LOG_INFO_ON_ROOT(outString);
      distForSizeFraction_ = std::discrete_distribution<uint_t>(particleNumbers.begin(), particleNumbers.end());

      WALBERLA_LOG_INFO_ON_ROOT("Read and stored in total " << particleMeshPerFractionVector_.size() << " size fractions and their meshes, and will randomly pick one for each to-be generated size fraction.")
   }

   void setShape(real_t /*diameter*/, real_t maximumAllowedInteractionRadius, data::ParticleStorage::baseShape_type& shape, real_t& interactionRadius) override
   {
      auto sizeFractionIdx = distForSizeFraction_(gen_);
      auto meshCopy = particleMeshPerFractionVector_[sizeFractionIdx][distsPerFraction_[sizeFractionIdx](gen_)];
      auto convexPolyPtr = std::make_shared<data::ConvexPolyhedron>(meshCopy);
      convexPolyPtr->updateMeshQuantities();
      shape = convexPolyPtr;
      interactionRadius = convexPolyPtr->getBoundingSphereRadius();
      WALBERLA_CHECK_GREATER_EQUAL(maximumAllowedInteractionRadius, interactionRadius, "Particle shape larger than allowed radius")
   }

   real_t getMaxDiameterScalingFactor() override { return maxDiameterScalingFactor_; }
   real_t getNormalVolume() override { return 1_r; } //dummy value
   FormParameters getNormalFormParameters() override {return FormParameters(1_r);} // dummy value
   bool generatesSingleShape() override { return generatesSingleShape_; }

private:
   std::vector<std::vector<mesh::TriangleMesh>> particleMeshPerFractionVector_;
   std::mt19937 gen_;
   std::discrete_distribution<uint_t> distForSizeFraction_;
   std::vector<std::uniform_int_distribution<uint_t>> distsPerFraction_;

   real_t maxDiameterScalingFactor_;
   bool generatesSingleShape_;
};


} // namespace msa_pd
} // namespace walberla
