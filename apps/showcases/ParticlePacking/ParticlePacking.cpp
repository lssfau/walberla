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
//! \file   ParticlePacking.cpp
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================


#include "blockforest/Initialization.h"
#include "blockforest/BlockForest.h"
#include "blockforest/loadbalancing/InfoCollection.h"
#include "blockforest/loadbalancing/DynamicCurve.h"
#include "blockforest/loadbalancing/weight_assignment/WeightAssignmentFunctor.h"
#include "core/Environment.h"
#include "core/grid_generator/SCIterator.h"
#include "core/grid_generator/HCPIterator.h"
#include "core/math/all.h"
#include "core/timing/TimingTree.h"
#include "vtk/VTKOutput.h"

#include "mesa_pd/collision_detection/AnalyticContactDetection.h"
#include "mesa_pd/collision_detection/GeneralContactDetection.h"

#include "mesa_pd/common/ParticleFunctions.h"

#include "mesa_pd/data/ContactStorage.h"
#include "mesa_pd/data/ContactAccessor.h"
#include "mesa_pd/data/ParticleAccessorWithBaseShape.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/HashGrids.h"

#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/domain/BlockForestDataHandling.h"
#include "mesa_pd/domain/InfoCollection.h"

#include "mesa_pd/kernel/AssocToBlock.h"
#include "mesa_pd/kernel/InsertParticleIntoLinkedCells.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/DetectAndStoreContacts.h"
#include "mesa_pd/kernel/InitContactsForHCSITS.h"
#include "mesa_pd/kernel/InitParticlesForHCSITS.h"
#include "mesa_pd/kernel/IntegrateParticlesHCSITS.h"
#include "mesa_pd/kernel/HCSITSRelaxationStep.h"
#include "mesa_pd/kernel/SemiImplicitEuler.h"
#include "mesa_pd/kernel/LinearSpringDashpot.h"

#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/ClearGhostOwnerSync.h"
#include "mesa_pd/mpi/ClearNextNeighborSync.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/ReduceContactHistory.h"
#include "mesa_pd/mpi/BroadcastProperty.h"
#include "mesa_pd/mpi/SyncNextNeighborsBlockForest.h"
#include "mesa_pd/mpi/SyncGhostOwners.h"
#include "mesa_pd/mpi/notifications/VelocityUpdateNotification.h"
#include "mesa_pd/mpi/notifications/VelocityCorrectionNotification.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/notifications/NumContactNotification.h"

#include "mesa_pd/sorting/LinearizedCompareFunctor.h"

#include "mesa_pd/vtk/ParticleVtkOutput.h"
#include "mesa_pd/vtk/TensorGlyph.h"
#include "mesa_pd/vtk/ConvexPolyhedron/MeshParticleVTKOutput.h"
#include "mesa_pd/vtk/ConvexPolyhedron/data_sources/SurfaceVelocityVertexDataSource.h"

#include "sqlite/SQLite.h"

#include <iostream>
#include <random>
#include <chrono>
#include <utility>
#include <vector>
#include <memory>

#include "Utility.h"
#include "Evaluation.h"
#include "DiameterDistribution.h"
#include "ShapeGeneration.h"

namespace walberla {
namespace mesa_pd {


kernel::HCSITSRelaxationStep::RelaxationModel relaxationModelFromString(const std::string & model)
{
   if(model == "InelasticFrictionlessContact") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::InelasticFrictionlessContact;
   if(model == "ApproximateInelasticCoulombContactByDecoupling") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::ApproximateInelasticCoulombContactByDecoupling;
   if(model == "ApproximateInelasticCoulombContactByOrthogonalProjections") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::ApproximateInelasticCoulombContactByOrthogonalProjections;
   if(model == "InelasticCoulombContactByDecoupling") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::InelasticCoulombContactByDecoupling;
   if(model == "InelasticCoulombContactByOrthogonalProjections") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::InelasticCoulombContactByOrthogonalProjections;
   if(model == "InelasticGeneralizedMaximumDissipationContact") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::InelasticGeneralizedMaximumDissipationContact;
   if(model == "InelasticProjectedGaussSeidel") 
      return kernel::HCSITSRelaxationStep::RelaxationModel::InelasticProjectedGaussSeidel;
   WALBERLA_ABORT("Unknown relaxation model " << model);
}

class GenerationPositionChecker
{
 public:
   using PositionDescription_T = std::pair<Vec3, real_t>;

   explicit GenerationPositionChecker(const std::string & domainSetup, real_t domainWidth, real_t largeParticleDiameterThreshold, bool alwaysGenerateLargeParticles)
      : domainSetup_(domainSetup), domainWidth_(domainWidth), largeParticleDiameterThreshold_(largeParticleDiameterThreshold),
        alwaysGenerateLargeParticles_(alwaysGenerateLargeParticles)
   {}

   void addParticleIfLarge(const Vec3& pos, real_t diameter)
   {
      if(largeParticleDiameterThreshold_ > 0_r && diameter > largeParticleDiameterThreshold_)
         positionsOfLargeParticles_.emplace_back(PositionDescription_T(pos, diameter));
   }

   bool isGenerationAllowed(const Vec3& pos, real_t diameter)
   {
      if( alwaysGenerateLargeParticles_ && diameter > largeParticleDiameterThreshold_) return true;

      for(const auto& descriptionPair : positionsOfLargeParticles_)
      {
         // check for containment
         auto largePos = descriptionPair.first;
         auto largeDiameter = descriptionPair.second;
         auto smallestDistanceSqr = (largePos - pos).sqrLength();


         if( diameter > largeParticleDiameterThreshold_ &&
             realIsEqual(largePos[0], pos[0]) &&
             realIsEqual(largePos[1], pos[1]) &&
             realIsEqual(largePos[2], pos[2])) continue; // large particle checks itself

         if(domainSetup_ == "periodic")
         {
            // need to check all periodic copies
            for(int deltaXFac = -1; deltaXFac <= 1; ++deltaXFac)
            {
               for(int deltaYFac = -1; deltaYFac <= 1; ++deltaYFac)
               {
                  auto posDelta = Vec3(domainWidth_ * real_c(deltaXFac), domainWidth_ * real_c(deltaYFac), 0_r);
                  smallestDistanceSqr = std::min(smallestDistanceSqr, (largePos + posDelta - pos).sqrLength());
               }
            }
         }

         auto distanceThreshold = (largeDiameter + diameter) * 0.5_r;
         if( smallestDistanceSqr < distanceThreshold * distanceThreshold) return false;

      }
      return true;
   }

 private:
   std::string domainSetup_;
   real_t domainWidth_;
   real_t largeParticleDiameterThreshold_;
   bool alwaysGenerateLargeParticles_;
   std::vector<PositionDescription_T> positionsOfLargeParticles_;
};


class ParticleCreator
{
public:
   ParticleCreator(const std::shared_ptr<data::ParticleStorage> & particleStorage, const std::shared_ptr<domain::IDomain> & particleDomain,
                   const AABB & simulationDomain, const std::string & domainSetup, real_t particleDensity, bool scaleGenerationSpacingWithForm ) :
                   particleStorage_(particleStorage), particleDomain_(particleDomain), simulationDomain_(simulationDomain),
                   domainSetup_(domainSetup), particleDensity_(particleDensity), scaleGenerationSpacingWithForm_(scaleGenerationSpacingWithForm),
                   gen_(static_cast<unsigned long>(walberla::mpi::MPIManager::instance()->rank()))
   {  }

   void createParticles( real_t zMin, real_t zMax, real_t spacing,
                         const std::unique_ptr<DiameterGenerator> & diameterGenerator, const std::unique_ptr<ShapeGenerator> & shapeGenerator,
                         real_t initialVelocity, real_t maximumAllowedInteractionRadius, real_t generationPositionCheckerLargeDiameterThreshold )
   {
      // this scaling is done to flexibly change the generation scaling in x,y, and z direction, based on the average form
      auto spacingScaling = (scaleGenerationSpacingWithForm_) ? shapeGenerator->getNormalFormParameters() / shapeGenerator->getNormalFormParameters()[1]  // divide by I for normalization
                                                              : Vec3(1_r); // no scaling (= equal spacing) in all directions
      sortVector(spacingScaling); // S, I, L
      Vec3 invScaling(1_r/spacingScaling[0], 1_r/spacingScaling[1], 1_r/spacingScaling[2]);

      AABB creationDomain(simulationDomain_.xMin()*invScaling[0], simulationDomain_.yMin()*invScaling[1], zMin*invScaling[2],
                          simulationDomain_.xMax()*invScaling[0], simulationDomain_.yMax()*invScaling[1], zMax*invScaling[2]);
      Vec3 pointOfReference(0,0,(zMax+zMin)*0.5_r*invScaling[2]);

      bool alwaysGenerateLargeParticles = true;
      GenerationPositionChecker generationPositionChecker(domainSetup_, simulationDomain_.xSize(), generationPositionCheckerLargeDiameterThreshold, alwaysGenerateLargeParticles);
      std::vector<GenerationPositionChecker::PositionDescription_T> particlesToBeCreated; // locally for this process

      WALBERLA_LOG_INFO_ON_ROOT("Creating particles between z = " << zMin << " and " << zMax);
      // first generate all particle descriptions (position, diameter) of the particles that will be generated
      for (auto ptUnscaled : grid_generator::HCPGrid(creationDomain.getExtended(Vec3(-0.5_r * spacing, -0.5_r * spacing, 0_r)), pointOfReference, spacing))
      {
         Vec3 pt(ptUnscaled[0] * spacingScaling[0], ptUnscaled[1] * spacingScaling[1],
                 ptUnscaled[2] * spacingScaling[2]); // scale back
         auto diameter = diameterGenerator->get();
         if (domainSetup_ == "container")
         {
            auto domainCenter             = simulationDomain_.center();
            auto distanceFromDomainCenter = pt - domainCenter;
            distanceFromDomainCenter[2]   = real_t(0);
            auto distance                 = distanceFromDomainCenter.length();
            real_t containerRadius        = real_t(0.5) * simulationDomain_.xSize();
            if (distance > containerRadius - real_t(0.5) * spacing) continue;
         }

         // for non-spherical particles only insufficient diameter information is available -> scale with safety factor, here max scaling factor
         generationPositionChecker.addParticleIfLarge(pt, diameter * std::sqrt(shapeGenerator->getMaxDiameterScalingFactor())); //std::cbrt(shapeGenerator->getMaxDiameterScalingFactor()));

         if (!particleDomain_->isContainedInLocalSubdomain(pt, real_t(0))) continue;

         particlesToBeCreated.emplace_back(GenerationPositionChecker::PositionDescription_T(pt, diameter));
      }

      // generate, possibly exclude some due to overlaps with others
      for(const auto& particleDescription: particlesToBeCreated)
      {
         auto pt = particleDescription.first;
         auto diameter = particleDescription.second;

         if (!generationPositionChecker.isGenerationAllowed(pt, diameter)) continue;

         // create particle
         auto p = particleStorage_->create();
         p->getPositionRef() = pt;

         shapeGenerator->setShape(diameter, maximumAllowedInteractionRadius, p->getBaseShapeRef(), p->getInteractionRadiusRef());

         p->getBaseShapeRef()->updateMassAndInertia(particleDensity_);

         p->setLinearVelocity( Vec3(0.1_r*math::realRandom(-initialVelocity, initialVelocity, gen_),
                                    0.1_r*math::realRandom(-initialVelocity, initialVelocity, gen_),
                                    -initialVelocity) );

         p->setAngularVelocity( 0.1_r*Vec3(math::realRandom(-initialVelocity, initialVelocity, gen_),
                                           math::realRandom(-initialVelocity, initialVelocity, gen_),
                                           math::realRandom(-initialVelocity, initialVelocity, gen_)) / diameter );

         p->getOwnerRef() = walberla::mpi::MPIManager::instance()->rank();
         p->getTypeRef() = 0;
      }
   }

private:
   std::shared_ptr<data::ParticleStorage> particleStorage_;
   std::shared_ptr<domain::IDomain> particleDomain_;
   AABB simulationDomain_;
   std::string domainSetup_;
   real_t particleDensity_;
   bool scaleGenerationSpacingWithForm_;
   std::mt19937 gen_;
};

void addConfigToDatabase(Config & config,
                         std::map< std::string, walberla::int64_t > & integerProperties,
                         std::map< std::string, double > & realProperties,
                         std::map< std::string, std::string > & stringProperties)
{

   const Config::BlockHandle mainConf = config.getBlock("ParticlePacking");

   Vector3<uint_t> numBlocksPerDirection = mainConf.getParameter< Vector3<uint_t> >("numBlocksPerDirection");
   integerProperties["numBlocksX"] = int64_c(numBlocksPerDirection[0]);
   integerProperties["numBlocksY"] = int64_c(numBlocksPerDirection[1]);
   integerProperties["numBlocksZ"] = int64_c(numBlocksPerDirection[2]);
   integerProperties["useHashGrids"] = (mainConf.getParameter<bool>("useHashGrids")) ? 1 : 0;
   stringProperties["domainSetup"] = mainConf.getParameter<std::string>("domainSetup");
   stringProperties["particleDistribution"] = mainConf.getParameter<std::string>("particleDistribution");
   stringProperties["particleShape"] = mainConf.getParameter<std::string>("particleShape");
   stringProperties["solver"] = mainConf.getParameter<std::string>("solver");
   realProperties["domainWidth"] = mainConf.getParameter<double>("domainWidth");
   realProperties["domainHeight"] = mainConf.getParameter<double>("domainHeight");
   realProperties["particleDensity"] = mainConf.getParameter<double>("particleDensity");
   realProperties["ambientDensity"] = mainConf.getParameter<double>("ambientDensity");
   realProperties["gravitationalAcceleration"] = mainConf.getParameter<double>("gravitationalAcceleration");
   realProperties["limitVelocity"] = mainConf.getParameter<double>("limitVelocity");
   realProperties["totalParticleMass"] = mainConf.getParameter<double>("totalParticleMass");

   const Config::BlockHandle solverConf = config.getBlock("Solver");
   realProperties["dt"] = solverConf.getParameter<double>("dt");
   realProperties["frictionCoefficient"] = solverConf.getParameter<double>("frictionCoefficient");
   realProperties["coefficientOfRestitution"] = solverConf.getParameter<double>("coefficientOfRestitution");
   const Config::BlockHandle solverHCSITSConf = solverConf.getBlock("HCSITS");
   integerProperties["hcsits_numberOfIterations"] = solverHCSITSConf.getParameter<int64_t>("numberOfIterations");
   stringProperties["hcsits_relaxationModel"] = solverHCSITSConf.getParameter<std::string>("relaxationModel");
   realProperties["hcsits_errorReductionParameter"] = solverHCSITSConf.getParameter<double>("errorReductionParameter");
   realProperties["hcsits_relaxationParameter"] = solverHCSITSConf.getParameter<double>("relaxationParameter");
   const Config::BlockHandle solverDEMConf = solverConf.getBlock("DEM");
   realProperties["dem_collisionTimeNonDim"] = solverDEMConf.getParameter<double>("collisionTime") / solverConf.getParameter<double>("dt");
   realProperties["dem_poissonsRatio"] = solverDEMConf.getParameter<double>("poissonsRatio");

   const Config::BlockHandle distributionConf = config.getBlock("Distribution");
   integerProperties["distribution_randomSeed"] = distributionConf.getParameter<int64_t>("randomSeed");
   const Config::BlockHandle uniformConf = distributionConf.getBlock("Uniform");
   realProperties["distribution_uniform_diameter"] = uniformConf.getParameter<real_t>("diameter");
   const Config::BlockHandle logNormalConf = distributionConf.getBlock("LogNormal");
   realProperties["distribution_logNormal_mu"] = logNormalConf.getParameter<real_t>("mu");
   realProperties["distribution_logNormal_variance"] = logNormalConf.getParameter<real_t>("variance");
   const Config::BlockHandle diamMassFracsConf = distributionConf.getBlock("DiameterMassFractions");
   stringProperties["distribution_diamMassFracs_diameters"] = diamMassFracsConf.getParameter<std::string>("diameters");
   stringProperties["distribution_diamMassFracs_massFractions"] = diamMassFracsConf.getParameter<std::string>("massFractions");
   const Config::BlockHandle sievingConf = distributionConf.getBlock("SievingCurve");
   stringProperties["distribution_sievingCurve_sieveSizes"] = sievingConf.getParameter<std::string>("sieveSizes");
   stringProperties["distribution_sievingCurve_massFractions"] = sievingConf.getParameter<std::string>("massFractions");
   integerProperties["distribution_sievingCurve_useDiscreteForm"] = (sievingConf.getParameter<bool>("useDiscreteForm")) ? 1 : 0;

   const Config::BlockHandle shapeConf = config.getBlock("Shape");
   stringProperties["shape_scaleMode"] = shapeConf.getParameter<std::string>("scaleMode");

   const Config::BlockHandle ellipsoidConf = shapeConf.getBlock("Ellipsoid");
   Vec3 ellipsoid_semiAxes = ellipsoidConf.getParameter< Vec3 >("semiAxes");
   realProperties["shape_ellipsoid_semiAxis0"] = double(ellipsoid_semiAxes[0]);
   realProperties["shape_ellipsoid_semiAxis1"] = double(ellipsoid_semiAxes[1]);
   realProperties["shape_ellipsoid_semiAxis2"] = double(ellipsoid_semiAxes[2]);
   const Config::BlockHandle equivalentEllipsoidConf = shapeConf.getBlock("EquivalentEllipsoid");
   stringProperties["shape_equivalentEllipsoid_path"] = equivalentEllipsoidConf.getParameter<std::string>("path");
   const Config::BlockHandle ellipsoidDistributionConf = shapeConf.getBlock("EllipsoidFormDistribution");
   realProperties["shape_ellipsoidFromDistribution_elongationMean"] = ellipsoidDistributionConf.getParameter<real_t>("elongationMean");
   realProperties["shape_ellipsoidFromDistribution_elongationStdDev"] = ellipsoidDistributionConf.getParameter<real_t>("elongationStdDev");
   realProperties["shape_ellipsoidFromDistribution_flatnessMean"] = ellipsoidDistributionConf.getParameter<real_t>("flatnessMean");
   realProperties["shape_ellipsoidFromDistribution_flatnessStdDev"] = ellipsoidDistributionConf.getParameter<real_t>("flatnessStdDev");
   const Config::BlockHandle meshConf = shapeConf.getBlock("Mesh");
   stringProperties["shape_mesh_path"] = meshConf.getParameter<std::string>("path");
   const Config::BlockHandle meshDistributionConf = shapeConf.getBlock("MeshFormDistribution");
   stringProperties["shape_meshFromDistribution_path"] = meshDistributionConf.getParameter<std::string>("path");
   realProperties["shape_meshFromDistribution_elongationMean"] = meshDistributionConf.getParameter<real_t>("elongationMean");
   realProperties["shape_meshFromDistribution_elongationStdDev"] = meshDistributionConf.getParameter<real_t>("elongationStdDev");
   realProperties["shape_meshFromDistribution_flatnessMean"] = meshDistributionConf.getParameter<real_t>("flatnessMean");
   realProperties["shape_meshFromDistribution_flatnessStdDev"] = meshDistributionConf.getParameter<real_t>("flatnessStdDev");
   const Config::BlockHandle meshesUnscaledConf = shapeConf.getBlock("UnscaledMeshesPerFraction");
   stringProperties["shape_meshesUnscaled_folder"] = meshesUnscaledConf.getParameter<std::string>("folder");

   const Config::BlockHandle evaluationConf = config.getBlock("evaluation");
   stringProperties["evaluation_histogramBins"] = evaluationConf.getParameter<std::string>("histogramBins");
   realProperties["evaluation_layerHeight"] = evaluationConf.getParameter<real_t>("layerHeight");

   const Config::BlockHandle generationConf = config.getBlock("Generation");
   realProperties["initialVelocity"] = generationConf.getParameter<double>("initialVelocity");
   realProperties["initialGenerationHeightRatioStart"] = generationConf.getParameter<double>("initialGenerationHeightRatioStart");
   realProperties["initialGenerationHeightRatioEnd"] = generationConf.getParameter<double>("initialGenerationHeightRatioEnd");
   realProperties["generationSpacing"] = generationConf.getParameter<double>("generationSpacing");
   stringProperties["generationPositionCheckerLargeDiameterThreshold"] = generationConf.getParameter<std::string>("generationPositionCheckerLargeDiameterThreshold");
   realProperties["generationHeightRatioStart"] = generationConf.getParameter<double>("generationHeightRatioStart");
   realProperties["generationHeightRatioEnd"] = generationConf.getParameter<double>("generationHeightRatioEnd");
   integerProperties["scaleGenerationSpacingWithForm"] = (generationConf.getParameter<bool>("scaleGenerationSpacingWithForm")) ? 1 : 0;

   integerProperties["shaking"] = (mainConf.getParameter<bool>("shaking")) ? 1 : 0;
   const Config::BlockHandle shakingConf = config.getBlock("Shaking");
   realProperties["shaking_amplitude"] = shakingConf.getParameter<double>("amplitude");
   realProperties["shaking_period"] = shakingConf.getParameter<double>("period");
   realProperties["shaking_duration"] = shakingConf.getParameter<double>("duration");
   integerProperties["shaking_activeFromBeginning"] = (shakingConf.getParameter<bool>("activeFromBeginning")) ? 1 : 0;

   const Config::BlockHandle dampingConf = config.getBlock("Damping");
   stringProperties["damping_method"] = dampingConf.getParameter<std::string>("method");
   realProperties["damping_velocityDampingCoefficient"] = dampingConf.getParameter<real_t>("velocityDampingCoefficient");
   realProperties["damping_forceDampingCoefficient"] = dampingConf.getParameter<real_t>("forceDampingCoefficient");

   const Config::BlockHandle terminationConf = config.getBlock("Termination");
   realProperties["terminalVelocity"] = terminationConf.getParameter<double>("terminalVelocity");
   realProperties["terminalRelativeHeightChange"] = terminationConf.getParameter<double>("terminalRelativeHeightChange");
   realProperties["terminationCheckingSpacing"] = terminationConf.getParameter<double>("terminationCheckingSpacing");
   realProperties["minimalTerminalRunTime"] = terminationConf.getParameter<double>("minimalTerminalRunTime");
   realProperties["maximalTerminalRunTime"] = terminationConf.getParameter<double>("maximalTerminalRunTime");

}

class SelectTensorGlyphForEllipsoids
{
public:
   using return_type = vtk::TensorGlyph;
   return_type operator()(const data::Particle& p) const {
      WALBERLA_CHECK_EQUAL(p.getBaseShape()->getShapeType(), data::Ellipsoid::SHAPE_TYPE);

      auto ellipsoid = static_cast<data::Ellipsoid*>(p.getBaseShape().get());
      return vtk::createTensorGlyph(ellipsoid->getSemiAxes(),p.getRotation());
   }
};



/*
 * Application to generate dense random packings of particles
 *
 * Main features:
 * - two domain setups: cylindrical container, horizontally periodic
 * - two simulation approaches: discrete element method (DEM), hard-contact semi-implicit timestepping solver (HCSITS)
 * - different size distributions
 * - different shapes: spherical, ellipsoidal, polygonal as given by mesh
 * - different generation methods: uniformly, large/small-specialized (for binary cases)
 * - evaluation of vertical porosity profile
 * - VTK visualization
 * - logging of final result and all properties into SQlite database
 * - requires OpenMesh
 * - optional: horizontal forcing to model different deposition conditions
 * - checkpointing
 * - load balancing (use with care in periodic setup as three processes per periodic direction is not enforced
 *   by load balancing but required for correct particle simulation)
 *
 * Simulation process:
 * - Generation phase: continuous generation in upper part of domain and settling due to gravity
 * - Shaking phase (optional, can also be active during generation phase): Shaking in a horizontal direction to compactify packing
 * - Termination phase: Run until converged state is reached
 *
 * See corresponding publication for more infos and cite if using this application:
 * Rettinger, C., Rüde, U., Vollmer, S., Frings, R. - "Effect of sediment form and form distribution on porosity:
 * a simulation study based on the discrete element method". Granular Matter 24, 118 (2022).
 * https://doi.org/10.1007/s10035-022-01275-x
 *
 */
int main(int argc, char **argv) {

   /// Setup
   Environment env(argc, argv);

   /// Config
   auto cfg = env.config();
   WALBERLA_LOG_INFO_ON_ROOT(*cfg);
   if (cfg == nullptr) WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle mainConf = cfg->getBlock("ParticlePacking");

   std::string domainSetup = mainConf.getParameter<std::string>("domainSetup");
   WALBERLA_CHECK(domainSetup == "container" || domainSetup == "periodic");
   real_t domainWidth = mainConf.getParameter<real_t>("domainWidth");
   real_t domainHeight = mainConf.getParameter<real_t>("domainHeight");
   real_t particleDensity = mainConf.getParameter<real_t>("particleDensity");
   real_t ambientDensity = mainConf.getParameter<real_t>("ambientDensity");
   real_t gravitationalAcceleration = mainConf.getParameter<real_t>("gravitationalAcceleration");
   real_t reducedGravitationalAcceleration = (particleDensity - ambientDensity) / particleDensity * gravitationalAcceleration;

   std::string particleDistribution = mainConf.getParameter<std::string>("particleDistribution");
   std::string particleShape = mainConf.getParameter<std::string>("particleShape");

   real_t limitVelocity = mainConf.getParameter<real_t>("limitVelocity");
   real_t totalParticleMass = mainConf.getParameter<real_t>("totalParticleMass");

   uint_t propertyEvaluationSpacing = mainConf.getParameter<uint_t>("propertyEvaluationSpacingInTimeSteps");
   real_t visSpacingInSeconds = mainConf.getParameter<real_t>("visSpacing");
   real_t infoSpacingInSeconds = mainConf.getParameter<real_t>("infoSpacing");
   real_t loggingSpacingInSeconds = mainConf.getParameter<real_t>("loggingSpacing");
   real_t performanceLoggingSpacingInSeconds = mainConf.getParameter<real_t>("performanceLoggingSpacing");
   Vector3<uint_t> numBlocksPerDirection = mainConf.getParameter< Vector3<uint_t> >("numBlocksPerDirection");

   bool useHashGrids = mainConf.getParameter<bool>("useHashGrids");

   std::string solver = mainConf.getParameter<std::string>("solver");

   int particleSortingSpacing = mainConf.getParameter< int >("particleSortingSpacing");


   const Config::BlockHandle solverConf = cfg->getBlock("Solver");
   real_t dt = solverConf.getParameter<real_t>("dt");
   real_t frictionCoefficient = solverConf.getParameter<real_t>("frictionCoefficient");
   real_t coefficientOfRestitution = solverConf.getParameter<real_t>("coefficientOfRestitution");

   uint_t visSpacing = uint_c(visSpacingInSeconds / dt);
   uint_t infoSpacing = uint_c(infoSpacingInSeconds / dt);
   uint_t loggingSpacing = uint_c(loggingSpacingInSeconds / dt);
   uint_t performanceLoggingSpacing = uint_c(performanceLoggingSpacingInSeconds / dt);
   WALBERLA_LOG_INFO_ON_ROOT("Simulation property evaluation spacing = " << propertyEvaluationSpacing);
   WALBERLA_LOG_INFO_ON_ROOT("VTK spacing = " << visSpacing << ", info spacing = " << infoSpacing << ", logging spacing = " << loggingSpacing);
   WALBERLA_CHECK_LESS_EQUAL(propertyEvaluationSpacing, infoSpacing);
   WALBERLA_CHECK_LESS_EQUAL(propertyEvaluationSpacing, loggingSpacing);

   const Config::BlockHandle solverHCSITSConf = solverConf.getBlock("HCSITS");
   real_t hcsits_errorReductionParameter = solverHCSITSConf.getParameter<real_t>("errorReductionParameter");
   real_t hcsits_relaxationParameter = solverHCSITSConf.getParameter<real_t>("relaxationParameter");
   std::string hcsits_relaxationModel = solverHCSITSConf.getParameter<std::string>("relaxationModel");
   uint_t hcsits_numberOfIterations = solverHCSITSConf.getParameter<uint_t>("numberOfIterations");

   const Config::BlockHandle solverDEMConf = solverConf.getBlock("DEM");
   real_t dem_collisionTime = solverDEMConf.getParameter<real_t>("collisionTime");
   //real_t dem_stiffnessNormal = solverDEMConf.getParameter<real_t>("stiffnessNormal");
   real_t dem_poissonsRatio = solverDEMConf.getParameter<real_t>("poissonsRatio");
   real_t dem_kappa = real_t(2) * ( real_t(1) - dem_poissonsRatio ) / ( real_t(2) - dem_poissonsRatio ) ; // from Thornton et al

   const Config::BlockHandle generationConf = cfg->getBlock("Generation");
   real_t initialVelocity = generationConf.getParameter<real_t>("initialVelocity");
   real_t initialGenerationHeightRatioStart = generationConf.getParameter<real_t>("initialGenerationHeightRatioStart");
   real_t initialGenerationHeightRatioEnd = generationConf.getParameter<real_t>("initialGenerationHeightRatioEnd");
   real_t generationSpacing = generationConf.getParameter<real_t>("generationSpacing");
   WALBERLA_CHECK_GREATER(domainWidth, generationSpacing, "Generation Spacing has to be smaller than domain size");
   real_t generationHeightRatioStart = generationConf.getParameter<real_t>("generationHeightRatioStart");
   real_t generationHeightRatioEnd = generationConf.getParameter<real_t>("generationHeightRatioEnd");
   bool scaleGenerationSpacingWithForm = generationConf.getParameter<bool>("scaleGenerationSpacingWithForm");
   std::string generationPositionCheckerLargeDiameterThresholdStr = generationConf.getParameter<std::string>("generationPositionCheckerLargeDiameterThreshold");

   bool shaking = mainConf.getParameter<bool>("shaking");
   const Config::BlockHandle shakingConf = cfg->getBlock("Shaking");
   real_t shaking_amplitude = shakingConf.getParameter<real_t>("amplitude");
   real_t shaking_period = shakingConf.getParameter<real_t>("period");
   real_t shaking_duration = shakingConf.getParameter<real_t>("duration");
   bool shaking_activeFromBeginning = shakingConf.getParameter<bool>("activeFromBeginning");

   const Config::BlockHandle dampingConf = cfg->getBlock("Damping");
   std::string damping_method = dampingConf.getParameter<std::string>("method");
   real_t damping_velocityDampingCoefficient = dampingConf.getParameter<real_t>("velocityDampingCoefficient");
   real_t damping_forceDampingCoefficient = dampingConf.getParameter<real_t>("forceDampingCoefficient");

   const Config::BlockHandle terminationConf = cfg->getBlock("Termination");
   real_t terminalVelocity = terminationConf.getParameter<real_t>("terminalVelocity");
   real_t terminalRelativeHeightChange = terminationConf.getParameter<real_t>("terminalRelativeHeightChange");
   real_t terminationCheckingSpacing = terminationConf.getParameter<real_t>("terminationCheckingSpacing");
   real_t minimalTerminalRunTime = terminationConf.getParameter<real_t>("minimalTerminalRunTime");
   real_t maximalTerminalRunTime = terminationConf.getParameter<real_t>("maximalTerminalRunTime");
   if(!(maximalTerminalRunTime < 0_r)) WALBERLA_CHECK_GREATER_EQUAL(maximalTerminalRunTime, minimalTerminalRunTime, "Maximal terminal run time has to be larger than minimal one!");

   bool useLoadBalancing = mainConf.getParameter<bool>("loadBalancing");
   const Config::BlockHandle loadBalancingConf = cfg->getBlock("LoadBalancing");
   auto loadBalancing_spacingInSeconds = loadBalancingConf.getParameter<real_t>("spacing");
   std::string loadBalancing_method = loadBalancingConf.getParameter<std::string>("method");
   auto loadBalancing_baseWeight = loadBalancingConf.getParameter<real_t>("baseWeight");

   const Config::BlockHandle checkpointRestartConf = cfg->getBlock("CheckpointRestart");
   auto checkPointing_spacingInSeconds = checkpointRestartConf.getParameter<real_t>("spacing");
   std::string checkPointing_existingUID = checkpointRestartConf.getParameter<std::string>("existingUID");
   std::string checkPointing_folder = checkpointRestartConf.getParameter<std::string>("folder");

   uint_t checkPointing_spacing = 0;
   if(checkPointing_spacingInSeconds > 0) checkPointing_spacing = uint_c(checkPointing_spacingInSeconds / dt);

   std::string uniqueFileIdentifier;
   bool restarting = false;
   // state variables that require proper initialization
   uint_t timestep = 0;
   bool isDampingActive = false;
   bool isShakingActive = false;
   real_t oldAvgParticleHeight = real_t(1);
   real_t oldMaxParticleHeight = real_t(1);
   real_t timeLastTerminationCheck = real_t(0);
   real_t timeLastCreation = real_t(0);
   real_t timeBeginShaking = real_t(-1);
   real_t timeEndShaking = real_t(-1);
   real_t timeBeginDamping = real_t(-1);

   if(checkPointing_existingUID == "None")
   {
      // clean start
      uniqueFileIdentifier = std::to_string(std::chrono::system_clock::now().time_since_epoch().count()); // used as hash to identify this run
      walberla::mpi::broadcastObject(uniqueFileIdentifier);
      if(shaking && shaking_activeFromBeginning)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Will use shaking from beginning.");
         isShakingActive = true;
         timeBeginShaking = real_t(0);
      }
   } else
   {
      // restart from existing file
      restarting = true;
      uniqueFileIdentifier = checkPointing_existingUID;

      timestep = checkpointRestartConf.getParameter<uint_t>("timestep");
      isDampingActive = checkpointRestartConf.getParameter<bool>("isDampingActive");
      isShakingActive = checkpointRestartConf.getParameter<bool>("isShakingActive");
      oldAvgParticleHeight = checkpointRestartConf.getParameter<real_t>("oldAvgParticleHeight");
      oldMaxParticleHeight = checkpointRestartConf.getParameter<real_t>("oldMaxParticleHeight");
      timeLastTerminationCheck = checkpointRestartConf.getParameter<real_t>("timeLastTerminationCheck");
      timeLastCreation = checkpointRestartConf.getParameter<real_t>("timeLastCreation");
      timeBeginShaking = checkpointRestartConf.getParameter<real_t>("timeBeginShaking");
      timeEndShaking = checkpointRestartConf.getParameter<real_t>("timeEndShaking");
      timeBeginDamping = checkpointRestartConf.getParameter<real_t>("timeBeginDamping");
   }
   std::string checkPointFileNameBase = checkPointing_folder+"/"+uniqueFileIdentifier + "_checkpoint";
   std::string checkPointFileName_forest = checkPointFileNameBase + "_forest.txt";
   std::string checkPointFileName_mesapd = checkPointFileNameBase + "_mesapd.txt";

   const Config::BlockHandle horizontalForcingConf = cfg->getBlock("HorizontalForcing");
   Vec3 horizontalForcing_acceleration = horizontalForcingConf.getParameter< Vec3 >("accelerationVec");
   real_t horizontalForcing_NonDimFixingHeight = horizontalForcingConf.getParameter<real_t>("nonDimFixingHeight");
   auto reducedHorizontalAcceleration = (particleDensity - ambientDensity) / particleDensity * horizontalForcing_acceleration;

   const Config::BlockHandle evaluationConf = cfg->getBlock("evaluation");
   auto evaluationHistogramBins = parseStringToVector<real_t>(evaluationConf.getParameter<std::string>("histogramBins"));
   //real_t voxelsPerMm = evaluationConf.getParameter<real_t>("voxelsPerMm");
   std::string porosityProfileFolder = evaluationConf.getParameter<std::string>("porosityProfileFolder");
   real_t evaluationLayerHeight = evaluationConf.getParameter<real_t>("layerHeight");
   std::string vtkOutputFolder = evaluationConf.getParameter<std::string>("vtkFolder");
   std::string vtkFinalFolder = evaluationConf.getParameter<std::string>("vtkFinalFolder");
   std::string sqlDBFileName = evaluationConf.getParameter<std::string>("sqlDBFileName");
   bool includeNumberOfContactsInVTK = evaluationConf.getParameter<bool>("includeNumberOfContactsInVTK");

   const Config::BlockHandle shapeConf = cfg->getBlock("Shape");
   ScaleMode shapeScaleMode = str_to_scaleMode(shapeConf.getParameter<std::string>("scaleMode"));
   const Config::BlockHandle distributionConf = cfg->getBlock("Distribution");

   /// BlockForest
   math::AABB simulationDomain(-0.5_r*domainWidth, -0.5_r*domainWidth, 0_r,
                               0.5_r*domainWidth, 0.5_r*domainWidth, domainHeight);
   Vector3<bool> isPeriodic = (domainSetup == "container") ? Vector3<bool>(false) : Vector3<bool>(true, true, false);


   shared_ptr<BlockForest> forest;
   if(restarting)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Reading block forest from file!");

      WALBERLA_MPI_SECTION()
      {
         if (!MPIManager::instance()->rankValid())
            MPIManager::instance()->useWorldComm();
      }
      forest = std::make_shared< BlockForest >( uint_c( MPIManager::instance()->rank() ),
                                               checkPointFileName_forest.c_str(), true, false );
      WALBERLA_LOG_INFO_ON_ROOT("Using domain of size " << forest->getDomain() << " with " <<
                                forest->getXSize() << " x " << forest->getYSize() << " x " << forest->getZSize() << " blocks");
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT("Creating domain of size " << simulationDomain);
      forest = blockforest::createBlockForest(simulationDomain, numBlocksPerDirection, isPeriodic);
      if(checkPointing_spacing != uint_t(0))
      {
         WALBERLA_LOG_INFO_ON_ROOT("Writing block forest to file!");
         forest->saveToFile(checkPointFileName_forest);
      }
   }


   auto domain = std::make_shared<mesa_pd::domain::BlockForestDomain>(forest);

   auto loadBalancingInfoCollection = make_shared<blockforest::InfoCollection>();
   uint_t loadBalancing_spacing = 0;
   if(useLoadBalancing)
   {
      loadBalancing_spacing = uint_c(loadBalancing_spacingInSeconds / dt);

      forest->recalculateBlockLevelsInRefresh( false ); // = only load balancing, no refinement checking
      forest->alwaysRebalanceInRefresh( true ); //load balancing every time refresh is triggered
      forest->reevaluateMinTargetLevelsAfterForcedRefinement( false );
      forest->allowRefreshChangingDepth( false );
      forest->allowMultipleRefreshCycles( false ); // otherwise info collections are invalid

      forest->setRefreshPhantomBlockDataAssignmentFunction( blockforest::WeightAssignmentFunctor( loadBalancingInfoCollection, loadBalancing_baseWeight ) );

      bool curveAllGather = true;
      bool balanceLevelwise = true;
      bool useHilbert = (loadBalancing_method == "Hilbert");

      WALBERLA_LOG_INFO_ON_ROOT("Using load balancing with " << ((useHilbert) ? "Hilbert" : "Morton") << " curve and spacing " << loadBalancing_spacingInSeconds << " s ( = " << loadBalancing_spacing << " time steps)");

      forest->setRefreshPhantomBlockMigrationPreparationFunction( blockforest::DynamicCurveBalance< blockforest::WeightAssignmentFunctor::PhantomBlockWeight >( useHilbert, curveAllGather, balanceLevelwise ) );
   }

   /// MESAPD Data
   auto particleStorage = std::make_shared<data::ParticleStorage>(1);
   data::ParticleAccessorWithBaseShape particleAccessor(particleStorage);

   BlockDataID particleStorageID;
   if(restarting)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Initializing particles from checkpointing file " << checkPointFileName_mesapd );
      particleStorageID = forest->loadBlockData( checkPointFileName_mesapd, mesa_pd::domain::createBlockForestDataHandling(particleStorage), "Particle Storage" );

      mesa_pd::mpi::ClearNextNeighborSync CNNS;
      CNNS(particleAccessor);
      mesa_pd::mpi::ClearGhostOwnerSync CGOS;
      CGOS(particleAccessor);
   } else {
      particleStorageID = forest->addBlockData(mesa_pd::domain::createBlockForestDataHandling(particleStorage), "Particle Storage");
   }

   auto contactStorage = std::make_shared<data::ContactStorage>(1);
   data::ContactAccessor contactAccessor(contactStorage);

   // configure shape creation
   std::unique_ptr<ShapeGenerator> shapeGenerator;
   if(particleShape == "Sphere")
   {
      shapeGenerator = std::make_unique<SphereGenerator>();
   } else if(particleShape == "Ellipsoid")
   {
      auto ellipsoidConfig = shapeConf.getBlock("Ellipsoid");
      std::vector<Vec3> semiAxes = {ellipsoidConfig.getParameter<Vec3>("semiAxes")};

      std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator = std::make_unique<SampleFormGenerator>(semiAxes, shapeScaleMode);
      shapeGenerator = std::make_unique<EllipsoidGenerator>(std::move(normalizedFormGenerator));
   } else if(particleShape == "EquivalentEllipsoid")
   {
      auto ellipsoidConfig = shapeConf.getBlock("EquivalentEllipsoid");
      std::string meshPath = ellipsoidConfig.getParameter<std::string>("path");

      auto meshFileNames = getMeshFilesFromPath(meshPath);
      auto semiAxes = extractSemiAxesFromMeshFiles(meshFileNames);
      std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator = std::make_unique<SampleFormGenerator>(semiAxes, shapeScaleMode);

      shapeGenerator = std::make_unique<EllipsoidGenerator>(std::move(normalizedFormGenerator));
   } else if(particleShape == "EllipsoidFormDistribution")
   {
      auto ellipsoidConfig = shapeConf.getBlock("EllipsoidFormDistribution");
      auto elongationMean = ellipsoidConfig.getParameter<real_t>("elongationMean");
      auto elongationStdDev = ellipsoidConfig.getParameter<real_t>("elongationStdDev");
      auto flatnessMean = ellipsoidConfig.getParameter<real_t>("flatnessMean");
      auto flatnessStdDev = ellipsoidConfig.getParameter<real_t>("flatnessStdDev");

      std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator = std::make_unique<DistributionFormGenerator>(elongationMean, elongationStdDev, flatnessMean, flatnessStdDev, shapeScaleMode);

      shapeGenerator = std::make_unique<EllipsoidGenerator>(std::move(normalizedFormGenerator));
   }
   else if(particleShape == "Mesh")
   {
      auto meshConfig = shapeConf.getBlock("Mesh");
      std::string meshPath = meshConfig.getParameter<std::string>("path");

      auto meshFileNames = getMeshFilesFromPath(meshPath);
      std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator = std::make_unique<ConstFormGenerator>();

      shapeGenerator = std::make_unique<MeshesGenerator>(meshFileNames, shapeScaleMode, std::move(normalizedFormGenerator));
   } else if(particleShape == "MeshFormDistribution")
   {
      auto meshConfig = shapeConf.getBlock("MeshFormDistribution");
      std::string meshPath = meshConfig.getParameter<std::string>("path");
      auto elongationMean = meshConfig.getParameter<real_t>("elongationMean");
      auto elongationStdDev = meshConfig.getParameter<real_t>("elongationStdDev");
      auto flatnessMean = meshConfig.getParameter<real_t>("flatnessMean");
      auto flatnessStdDev = meshConfig.getParameter<real_t>("flatnessStdDev");

      auto meshFileNames = getMeshFilesFromPath(meshPath);
      std::unique_ptr<NormalizedFormGenerator> normalizedFormGenerator = std::make_unique<DistributionFormGenerator>(elongationMean, elongationStdDev, flatnessMean, flatnessStdDev, shapeScaleMode);

      shapeGenerator = std::make_unique<MeshesGenerator>(meshFileNames, shapeScaleMode, std::move(normalizedFormGenerator));
   } else if(particleShape == "UnscaledMeshesPerFraction")
   {
      shapeGenerator = std::make_unique<UnscaledMeshesPerFractionGenerator>(shapeConf,  parseStringToVector<real_t>(distributionConf.getBlock("SievingCurve").getParameter<std::string>("massFractions")));
   } else
   {
      WALBERLA_ABORT("Unknown shape " << particleShape);
   }
   WALBERLA_LOG_INFO_ON_ROOT("Will create particles with ");
   WALBERLA_LOG_INFO_ON_ROOT(" - maximum diameter scaling of " << shapeGenerator->getMaxDiameterScalingFactor());
   WALBERLA_LOG_INFO_ON_ROOT(" - normal volume " << shapeGenerator->getNormalVolume());
   WALBERLA_LOG_INFO_ON_ROOT(" - " << (shapeGenerator->generatesSingleShape() ? "single shape" : "multiple shapes"));

   // configure size creation
   int randomSeedFromConfig = distributionConf.getParameter<int>("randomSeed");
   uint_t randomSeed = (randomSeedFromConfig >= 0) ? uint_c(randomSeedFromConfig) : uint_c(time(nullptr));
   WALBERLA_LOG_INFO_ON_ROOT("Random seed of " << randomSeed);

   std::unique_ptr<DiameterGenerator> diameterGenerator;
   real_t minGenerationParticleDiameter = 0_r;
   real_t maxGenerationParticleDiameter = std::numeric_limits<real_t>::max();
   real_t d50 = 0_r;

   if(particleDistribution == "LogNormal")
   {
      const Config::BlockHandle logNormalConf = distributionConf.getBlock("LogNormal");
      real_t mu = logNormalConf.getParameter<real_t>("mu");
      real_t variance = logNormalConf.getParameter<real_t>("variance");
      diameterGenerator = std::make_unique<LogNormal>(mu, variance, randomSeed);
      // min and max diameter not determinable
      d50 = mu;
      WALBERLA_LOG_INFO_ON_ROOT("Using log-normal distribution with mu = " << mu << ", var = " << variance);
   }
   else if(particleDistribution == "Uniform")
   {
      const Config::BlockHandle uniformConf = distributionConf.getBlock("Uniform");
      real_t diameter = uniformConf.getParameter<real_t>("diameter");
      diameterGenerator = std::make_unique<Uniform>(diameter);
      minGenerationParticleDiameter = diameter;
      maxGenerationParticleDiameter = diameter;
      d50 = diameter;
      WALBERLA_LOG_INFO_ON_ROOT("Using uniform distribution");
   }
   else if(particleDistribution == "DiameterMassFractions")
   {
      const Config::BlockHandle sievingConf = distributionConf.getBlock("DiameterMassFractions");
      auto diameters = parseStringToVector<real_t>(sievingConf.getParameter<std::string>("diameters"));
      auto massFractions = parseStringToVector<real_t>(sievingConf.getParameter<std::string>("massFractions"));
      diameterGenerator = std::make_unique<DiscreteSieving>(diameters, massFractions, randomSeed, shapeGenerator->getNormalVolume(), totalParticleMass, particleDensity);

      maxGenerationParticleDiameter = real_t(0);
      minGenerationParticleDiameter = std::numeric_limits<real_t>::max();
      for(uint_t i = 0; i < diameters.size(); ++i) {
         if(massFractions[i] > real_t(0)) {
            maxGenerationParticleDiameter = std::max(maxGenerationParticleDiameter, diameters[i]);
            minGenerationParticleDiameter = std::min(minGenerationParticleDiameter, diameters[i]);
         }
      }
      d50 = computePercentileFromSieveDistribution(diameters, massFractions, 50_r);
      WALBERLA_LOG_INFO_ON_ROOT("Using diameter - mass fraction distribution");
   }
   else if(particleDistribution == "SievingCurve")
   {
      const Config::BlockHandle sievingConf = distributionConf.getBlock("SievingCurve");
      auto sieveSizes = parseStringToVector<real_t>(sievingConf.getParameter<std::string>("sieveSizes"));
      auto massFractions = parseStringToVector<real_t>(sievingConf.getParameter<std::string>("massFractions"));
      bool useDiscreteForm = sievingConf.getParameter<bool>("useDiscreteForm");

      auto diameters = getMeanDiametersFromSieveSizes(sieveSizes);
      d50 = computePercentileFromSieveDistribution(diameters, massFractions, 50_r);
      real_t d16 = computePercentileFromSieveDistribution(diameters, massFractions, 16_r);
      real_t d84 = computePercentileFromSieveDistribution(diameters, massFractions, 84_r);
      real_t stdDev = std::sqrt(d84/d16);
      WALBERLA_LOG_INFO_ON_ROOT("Curve properties: D50 = " << d50 << ", D16 = " << d16 << ", D84 = " << d84 << ", estimated std. dev. = " << stdDev);

      maxGenerationParticleDiameter = real_t(0);
      minGenerationParticleDiameter = std::numeric_limits<real_t>::max();
      if(useDiscreteForm)
      {
         diameterGenerator = std::make_unique<DiscreteSieving>(diameters, massFractions, randomSeed, shapeGenerator->getNormalVolume(), totalParticleMass, particleDensity);
         for(uint_t i = 0; i < diameters.size(); ++i) {
            if(massFractions[i] > real_t(0)) {
               maxGenerationParticleDiameter = std::max(maxGenerationParticleDiameter, diameters[i]);
               minGenerationParticleDiameter = std::min(minGenerationParticleDiameter, diameters[i]);
            }
         }
         WALBERLA_LOG_INFO_ON_ROOT("Using discrete sieving curve distribution");

      } else {
         diameterGenerator = std::make_unique<ContinuousSieving>(sieveSizes, massFractions, randomSeed, shapeGenerator->getNormalVolume(), totalParticleMass, particleDensity);
         for(uint_t i = 0; i < sieveSizes.size()-1; ++i) {
            if(massFractions[i] > real_t(0)) {
               maxGenerationParticleDiameter = std::max(maxGenerationParticleDiameter, std::max(sieveSizes[i],sieveSizes[i+1]));
               minGenerationParticleDiameter = std::min(minGenerationParticleDiameter, std::min(sieveSizes[i],sieveSizes[i+1]));
            }
         }
         WALBERLA_LOG_INFO_ON_ROOT("Using piece-wise constant / continuous sieving curve distribution");
      }
   }
   else
   {
      WALBERLA_ABORT("Unknown particle distribution specified: " << particleDistribution);
   }

   WALBERLA_LOG_INFO_ON_ROOT("Generate with diameters in range [" << minGenerationParticleDiameter << ", " << maxGenerationParticleDiameter << "] and D50 = " << d50);

   bool useOpenMP = false;

   real_t smallestBlockSize = std::min(simulationDomain.xSize() / real_c(numBlocksPerDirection[0]),
                                       std::min(simulationDomain.ySize() / real_c(numBlocksPerDirection[1]),
                                                simulationDomain.zSize() / real_c(numBlocksPerDirection[2])));

   // plane at top and bottom
   createPlane(particleStorage, Vector3<real_t>(0), Vector3<real_t>(0, 0, 1));
   createPlane(particleStorage, Vector3<real_t>(0_r,0_r,simulationDomain.zMax()), Vector3<real_t>(0, 0, -1));

   real_t domainVolume = simulationDomain.volume();
   if(domainSetup == "container")
   {
      createCylindricalBoundary(particleStorage, Vector3<real_t>(0), Vector3<real_t>(0, 0, 1), 0.5_r*domainWidth);
      domainVolume = math::pi * domainWidth * domainWidth * 0.25_r * simulationDomain.zSize();
   }

   real_t maximumAllowedInteractionRadius = std::numeric_limits<real_t>::infinity();
   if(domainSetup == "periodic")
   {
      // avoid that two large particles are next to each other and would, due to periodic mapping, have 2 different contact points with each other |( p1 () p2 ()| p1  )
      maximumAllowedInteractionRadius = 0.25_r * domainWidth; // max diameter = domainWidth / 2
      WALBERLA_LOG_INFO_ON_ROOT("Periodic case: the maximum interaction radius is restricted to " << maximumAllowedInteractionRadius << " to ensure valid periodic interaction" );
      if(numBlocksPerDirection[0] < 3 || numBlocksPerDirection[1] < 3) WALBERLA_LOG_INFO_ON_ROOT("Warning: At least 3 blocks per periodic direction required for proper simulation!")
   }

   // fill domain with particles initially
   real_t maxGenerationHeight = simulationDomain.zMax() - generationSpacing;
   real_t minGenerationHeight = generationSpacing;

   real_t generationPositionCheckerLargeDiameterThreshold = -1_r;
   if(generationPositionCheckerLargeDiameterThresholdStr == "OFF") generationPositionCheckerLargeDiameterThreshold = -1_r;
   else if(generationPositionCheckerLargeDiameterThresholdStr == "D50") generationPositionCheckerLargeDiameterThreshold = d50;
   else{
      std::istringstream os(generationPositionCheckerLargeDiameterThresholdStr);
      os >> generationPositionCheckerLargeDiameterThreshold;
   }

   WALBERLA_LOG_INFO_ON_ROOT("Generation spacing = " << generationSpacing << ", with large diameter threshold = " << generationPositionCheckerLargeDiameterThreshold);

   ParticleCreator particleCreator(particleStorage, domain, simulationDomain, domainSetup, particleDensity, scaleGenerationSpacingWithForm);
   if(timestep == 0)
   {
      particleCreator.createParticles(std::max(minGenerationHeight, initialGenerationHeightRatioStart * simulationDomain.zMax()),
                                      std::min(maxGenerationHeight, initialGenerationHeightRatioEnd * simulationDomain.zMax()),
                                      generationSpacing, diameterGenerator, shapeGenerator, initialVelocity,
                                      maximumAllowedInteractionRadius, generationPositionCheckerLargeDiameterThreshold );
   }


   math::DistributedSample diameterSample;
   particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                    [&diameterSample](const size_t idx, data::ParticleAccessorWithBaseShape& ac){diameterSample.insert(real_c(2)*ac.getInteractionRadius(idx));}, particleAccessor);

   diameterSample.mpiAllGather();
   WALBERLA_LOG_INFO_ON_ROOT("Statistics of initially created particles' interaction diameters: " << diameterSample.format());

   real_t maxParticleDiameter = maxGenerationParticleDiameter * shapeGenerator->getMaxDiameterScalingFactor();
   if(maxParticleDiameter < diameterSample.max())
   {
      WALBERLA_LOG_INFO_ON_ROOT("Maximum interaction diameter from samples is larger than estimated maximum diameter, will use sampled one instead.")
      maxParticleDiameter = 1.1_r * diameterSample.max(); // 10% safety margin
   }
   if(maxParticleDiameter > 2_r * maximumAllowedInteractionRadius)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Warning: Maximum expected particle interaction diameter ("<< maxParticleDiameter << ") is larger than maximum allowed interaction diameter - check that the generated size & form distributions match the expected ones!");
      maxParticleDiameter = 2_r * maximumAllowedInteractionRadius;
   }

   bool useNextNeighborSync = 2_r * smallestBlockSize > maxParticleDiameter;

   WALBERLA_LOG_INFO_ON_ROOT("Sync info: maximum expected interaction diameter = " << maxParticleDiameter << " and smallest block size = " << smallestBlockSize);

   // sync functionality
   kernel::AssocToBlock associateToBlock(forest);
   std::function<void(void)> syncCall;
   if(useNextNeighborSync)
   {
      WALBERLA_LOG_INFO_ON_ROOT("Using next neighbor sync!");
      syncCall = [&particleStorage,&forest,&domain](){
         mpi::SyncNextNeighborsBlockForest syncParticles;
         syncParticles(*particleStorage, forest, domain);
      };
   } else {
      WALBERLA_LOG_INFO_ON_ROOT("Using ghost owner sync!");
      syncCall = [&particleStorage,&domain](){
         mpi::SyncGhostOwners syncGhostOwnersFunc;
         syncGhostOwnersFunc(*particleStorage, *domain);
      };
   }

   // initial sync
   particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor, associateToBlock, particleAccessor);
   if(useNextNeighborSync){ syncCall(); }
   else { for(uint_t i = 0; i < uint_c(std::ceil(maxParticleDiameter/smallestBlockSize)); ++i) syncCall(); }


   // create linked cells data structure
   real_t linkedCellWidth = 1.01_r * maxParticleDiameter;
   WALBERLA_LOG_INFO_ON_ROOT("Using linked cells with cell width = " << linkedCellWidth);
   auto linkedCells = std::make_unique<data::LinkedCells>(domain->getUnionOfLocalAABBs().getExtended(linkedCellWidth), linkedCellWidth );

   {
      auto info = evaluateParticleInfo(particleAccessor);
      WALBERLA_LOG_INFO_ON_ROOT(info);
   }

   /// VTK Output
   auto vtkDomainOutput = walberla::vtk::createVTKOutput_DomainDecomposition( forest, "domain_decomposition", visSpacing, vtkOutputFolder, "simulation_step" );

   if(visSpacing > 0 && !useLoadBalancing)
   {
      vtkDomainOutput->write();
   }

   // mesapd particle output
   auto particleVtkOutput = make_shared<vtk::ParticleVtkOutput>(particleStorage);
   particleVtkOutput->addOutput<data::SelectParticleUid>("uid");
   particleVtkOutput->addOutput<data::SelectParticleOwner>("owner");
   particleVtkOutput->addOutput<data::SelectParticleInteractionRadius>("interactionRadius");
   if(particleShape.find("Ellipsoid") != std::string::npos)
   {
      particleVtkOutput->addOutput<SelectTensorGlyphForEllipsoids>("tensorGlyph");
   }
   particleVtkOutput->addOutput<data::SelectParticleLinearVelocity>("velocity");
   if(includeNumberOfContactsInVTK) particleVtkOutput->addOutput<data::SelectParticleNumContacts>("numContacts");
   auto vtkParticleSelector = [](const mesa_pd::data::ParticleStorage::iterator &pIt) {
      return (pIt->getBaseShape()->getShapeType() != data::HalfSpace::SHAPE_TYPE &&
              pIt->getBaseShape()->getShapeType() != data::CylindricalBoundary::SHAPE_TYPE &&
              !isSet(pIt->getFlags(), data::particle_flags::GHOST));
   };
   particleVtkOutput->setParticleSelector(vtkParticleSelector);
   auto particleVtkWriter = walberla::vtk::createVTKOutput_PointData(particleVtkOutput, "particles", visSpacing, vtkOutputFolder, "simulation_step");

   mesa_pd::MeshParticleVTKOutput< mesh::PolyMesh > meshParticleVTK(particleStorage, "mesh", visSpacing, vtkOutputFolder);
   meshParticleVTK.addFaceOutput< data::SelectParticleUid >("UID");
   meshParticleVTK.addVertexOutput< data::SelectParticleInteractionRadius >("InteractionRadius");
   meshParticleVTK.addFaceOutput< data::SelectParticleLinearVelocity >("LinearVelocity");
   meshParticleVTK.addVertexOutput< data::SelectParticlePosition >("Position");
   if(includeNumberOfContactsInVTK) meshParticleVTK.addVertexOutput< data::SelectParticleNumContacts >("numContacts");
   auto surfaceVelDataSource = make_shared<mesa_pd::SurfaceVelocityVertexDataSource< mesh::PolyMesh, data::ParticleAccessorWithBaseShape > >("SurfaceVelocity", particleAccessor);
   meshParticleVTK.setParticleSelector(vtkParticleSelector);
   meshParticleVTK.addVertexDataSource(surfaceVelDataSource);


   /// MESAPD kernels

   //collision detection
   data::HashGrids hashGrids;
   kernel::InsertParticleIntoLinkedCells initializeLinkedCells;

   //DEM
   kernel::SemiImplicitEuler dem_integration( dt );
   kernel::LinearSpringDashpot dem_collision(1);
   dem_collision.setFrictionCoefficientStatic(0,0,0_r); // no static friction
   dem_collision.setFrictionCoefficientDynamic(0,0, frictionCoefficient);
   // stiffness and damping depend on effective mass -> have to be calculated for each collision individually

   //HCSITS
   kernel::InitContactsForHCSITS hcsits_initContacts(1);
   hcsits_initContacts.setFriction(0,0,frictionCoefficient);
   hcsits_initContacts.setErp(hcsits_errorReductionParameter);
   kernel::InitParticlesForHCSITS hcsits_initParticles;
   hcsits_initParticles.setGlobalAcceleration(reducedHorizontalAcceleration + Vec3(0,0,-reducedGravitationalAcceleration));
   kernel::HCSITSRelaxationStep hcsits_relaxationStep;
   hcsits_relaxationStep.setRelaxationModel(relaxationModelFromString(hcsits_relaxationModel));
   hcsits_relaxationStep.setCor(coefficientOfRestitution); // Only effective for PGSM
   kernel::IntegrateParticlesHCSITS hcsits_integration;

   // sync
   mpi::ReduceContactHistory reduceAndSwapContactHistory;
   mpi::BroadcastProperty broadcastKernel;
   mpi::ReduceProperty reductionKernel;

   WALBERLA_LOG_INFO_ON_ROOT("Domain of volume " << domainVolume << " m^3.");
   WALBERLA_LOG_INFO_ON_ROOT("Will terminate generation when particle mass is above " << totalParticleMass << " kg.");

   real_t velocityDampingFactor = std::pow(damping_velocityDampingCoefficient, dt);
   if(damping_method == "force")  {WALBERLA_LOG_INFO_ON_ROOT("Once all particles are created, will apply direction-dependent damping of all forces with factor  " << damping_forceDampingCoefficient << ".");}
   else if(damping_method == "velocity")  {WALBERLA_LOG_INFO_ON_ROOT("Once all particles are created, will apply velocity damping of  " << velocityDampingFactor << " per time step.");}
   else {WALBERLA_ABORT("Unknown damping method " << damping_method << ".");}

   real_t maximumTimeBetweenCreation = (generationHeightRatioEnd-generationHeightRatioStart)*simulationDomain.zSize() / initialVelocity; // = time, particles need at max to clear/pass the creation domain
   WALBERLA_LOG_INFO_ON_ROOT("Maximum time between creation steps: " << maximumTimeBetweenCreation);

   if(limitVelocity > 0_r) WALBERLA_LOG_INFO_ON_ROOT("Will apply limiting of translational particle velocity to maximal magnitude of " << limitVelocity);

   SizeEvaluator particleSizeEvaluator(shapeScaleMode);
   std::vector<std::tuple<std::string, std::function<real_t(Vec3)>>> particleShapeEvaluators = {std::make_tuple("flatness",getFlatnessFromSemiAxes),
                                                                                                std::make_tuple("elongation",getElongationFromSemiAxes),
                                                                                                std::make_tuple("equancy",getEquancyFromSemiAxes)};
   uint_t numShapeBins = 17;
   std::vector<std::vector<real_t>> particleShapeBins(particleShapeEvaluators.size());
   for(auto& pSB: particleShapeBins)
   {
      pSB = std::vector<real_t>(numShapeBins, real_t(0));
      real_t binBegin = 0_r;
      real_t binEnd = 1_r;
      real_t inc = (binEnd - binBegin) / real_t(numShapeBins-1);
      real_t val = binBegin;
      for(uint_t i = 0; i < numShapeBins; ++i, val += inc) pSB[i] = val;
   }

   ParticleHistogram particleHistogram(evaluationHistogramBins, particleSizeEvaluator, particleShapeBins, particleShapeEvaluators);

   particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                    particleHistogram, particleAccessor);
   particleHistogram.evaluate();
   WALBERLA_LOG_INFO_ON_ROOT(particleHistogram);

   auto particleInfo = evaluateParticleInfo(particleAccessor);

   if(restarting)
   {
      if(particleInfo.particleVolume * particleDensity < totalParticleMass)
      {
         WALBERLA_ABORT("Cannot use this checkpoint config since further particle generation is necessary but not supported!");
         // note: this restriction is due to the UID generator that has its own, unrecoverable state
         // thus, a new generation of particles would result in the same UIDs being used as the ones from the already
         // present particles stored in the checkpointfile
      }
   }

   PorosityPerHorizontalLayerEvaluator porosityEvaluator(evaluationLayerHeight, simulationDomain, domainSetup);

   std::string loggingFileName = porosityProfileFolder + "/" +  uniqueFileIdentifier + "_logging.txt";
   WALBERLA_LOG_INFO_ON_ROOT("Writing logging file to " << loggingFileName);
   LoggingWriter loggingWriter(loggingFileName, !restarting);

   // horizontal forcing
   real_t horizontalForcing_fixingHeight = horizontalForcing_NonDimFixingHeight * maxParticleDiameter;
   bool fixBottomParticles = (horizontalForcing_fixingHeight > 0_r);
   bool applyHorizontalForcing = (std::abs(horizontalForcing_acceleration[0])>0_r || std::abs(horizontalForcing_acceleration[1])>0_r);
   if(applyHorizontalForcing) WALBERLA_LOG_INFO_ON_ROOT("Applying horizontal forcing with acc = " << horizontalForcing_acceleration);

   WcTimingTree timing;
   WcTimer simulationTimer;

   bool terminateSimulation = false;
   while (!terminateSimulation) {

      real_t currentTime = dt * real_c(timestep);

      timing.start("Sorting");
      if(particleSortingSpacing > 0 && timestep % uint_c(particleSortingSpacing) == 0 && !useHashGrids)
      {
         sorting::LinearizedCompareFunctor linearSorting(linkedCells->domain_, linkedCells->numCellsPerDim_);
         particleStorage->sort(linearSorting);
      }
      timing.stop("Sorting");

      timing.start("VTK");
      if(particleShape.find("Mesh") != std::string::npos) meshParticleVTK(particleAccessor);
      else particleVtkWriter->write();
      if(useLoadBalancing) vtkDomainOutput->write();
      timing.stop("VTK");


      contactStorage->clear();

      if(useHashGrids)
      {
         timing.start("Hash grid");
         hashGrids.clearAll();
         particleStorage->forEachParticle(useOpenMP, kernel::SelectAll(), particleAccessor, hashGrids, particleAccessor);
         timing.stop("Hash grid");

         timing.start("Contact detection");
         hashGrids.forEachParticlePairHalf(useOpenMP, kernel::ExcludeInfiniteInfinite(), particleAccessor,
                                           [&domain, &contactStorage](size_t idx1, size_t idx2, data::ParticleAccessorWithBaseShape &ac){
                                              kernel::DoubleCast double_cast;
                                              mpi::ContactFilter contact_filter;
                                              collision_detection::GeneralContactDetection contactDetection;
                                              //Attention: does not use contact threshold in general case (GJK)

                                              if (double_cast(idx1, idx2, ac, contactDetection, ac)) {
                                                 if (contact_filter(contactDetection.getIdx1(), contactDetection.getIdx2(), ac, contactDetection.getContactPoint(), *domain)) {
                                                    auto c = contactStorage->create();
                                                    c->setId1(contactDetection.getIdx1());
                                                    c->setId2(contactDetection.getIdx2());
                                                    c->setDistance(contactDetection.getPenetrationDepth());
                                                    c->setNormal(contactDetection.getContactNormal());
                                                    c->setPosition(contactDetection.getContactPoint());
                                                 }
                                              }
                                              }, particleAccessor);
         timing.stop("Contact detection");

      } else
      {
         // use linked cells

         timing.start("Linked cells");
         linkedCells->clear();
         particleStorage->forEachParticle(useOpenMP,kernel::SelectAll(),particleAccessor,
                                          initializeLinkedCells, particleAccessor, *linkedCells);
         timing.stop("Linked cells");

         timing.start("Contact detection");
         if(particleShape == "Sphere")
         {
            collision_detection::AnalyticContactDetection contactDetection;
            //acd.getContactThreshold() = contactThreshold;
            kernel::DetectAndStoreContacts detectAndStore(*contactStorage);
            linkedCells->forEachParticlePairHalf(useOpenMP, kernel::ExcludeInfiniteInfinite(), particleAccessor,
                                                detectAndStore, particleAccessor, *domain, contactDetection);

         } else
         {
            linkedCells->forEachParticlePairHalf(useOpenMP, kernel::ExcludeInfiniteInfinite(), particleAccessor,
                                                [&domain, &contactStorage](size_t idx1, size_t idx2, data::ParticleAccessorWithBaseShape &ac){

                                                   collision_detection::GeneralContactDetection contactDetection;
                                                   //Attention: does not use contact threshold in general case (GJK)

                                                   // coarse collision detection via interaction radii
                                                   data::Sphere sp1(ac.getInteractionRadius(idx1));
                                                   data::Sphere sp2(ac.getInteractionRadius(idx2));
                                                   if(contactDetection(idx1, idx2, sp1, sp2, ac)) {
                                                      //NOTE: this works also for infinite particles ( plane, cylindrical boundary) since contact_detection return true
                                                      // and the following contact_filter treats all local-global interactions independently of contact detection result, which would be non-sense for this interaction

                                                      mpi::ContactFilter contact_filter;
                                                      if (contact_filter(contactDetection.getIdx1(), contactDetection.getIdx2(), ac, contactDetection.getContactPoint(), *domain)) {
                                                         //NOTE: usually we first do fine collision detection and then the exact contact location determines the process that handles this contact
                                                         // however, along periodic boundaries, the GJK/EPA for meshes seems to be numerical unstable and yields (sometimes) different contact points for the same interaction pair (but periodically transformed)
                                                         // as a result, the same contact appears twice and potentially handled by two processes simultaneously.
                                                         // thus we change the ordering and do the contact filtering according to the result of the coarse collision detection, i.e. the bounding sphere check
                                                         kernel::DoubleCast double_cast;
                                                         if (double_cast(idx1, idx2, ac, contactDetection, ac)) {
                                                            auto c = contactStorage->create();
                                                            c->setId1(contactDetection.getIdx1());
                                                            c->setId2(contactDetection.getIdx2());
                                                            c->setDistance(contactDetection.getPenetrationDepth());
                                                            c->setNormal(contactDetection.getContactNormal());
                                                            c->setPosition(contactDetection.getContactPoint());
                                                         }
                                                      }
                                                   }}, particleAccessor);
         }

         timing.stop("Contact detection");
      }

      timing.start("Contact eval");
      if(includeNumberOfContactsInVTK)
      {
         particleStorage->forEachParticle(useOpenMP, kernel::SelectAll(), particleAccessor,
                                          [](size_t p_idx, data::ParticleAccessorWithBaseShape& ac){ac.setNumContacts(p_idx,0);}, particleAccessor);
         contactStorage->forEachContact(useOpenMP, kernel::SelectAll(), contactAccessor,
                                        [](size_t c, data::ContactAccessor &ca, data::ParticleAccessorWithBaseShape &pa) {
                                          auto idx1 = ca.getId1(c);
                                          auto idx2 = ca.getId2(c);
                                          pa.getNumContactsRef(idx1)++;
                                          pa.getNumContactsRef(idx2)++;
                                        }, contactAccessor, particleAccessor);

         reductionKernel.operator()<NumContactNotification>(*particleStorage);
      }
      timing.stop("Contact eval");


      timing.start("Shaking");
      if(isShakingActive)
      {
         real_t shaking_common_term = 2_r * math::pi / shaking_period;
         real_t shakingAcceleration = shaking_amplitude * std::sin((currentTime - timeBeginShaking) * shaking_common_term) * shaking_common_term * shaking_common_term;
         if(solver == "DEM")
         {
            particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                             [shakingAcceleration](const size_t idx, auto &ac){addForceAtomic(idx, ac, Vec3(shakingAcceleration,0_r,0_r) * ac.getMass(idx));}, particleAccessor);
         } else
         {
            hcsits_initParticles.setGlobalAcceleration((applyHorizontalForcing ? reducedHorizontalAcceleration : Vec3(0_r)) + Vec3(shakingAcceleration,0_r,-reducedGravitationalAcceleration));
         }
      }
      timing.stop("Shaking");

      if(solver == "HCSITS")
      {
         timing.start("HCSITS");

         timing.start("Init contacts");
         contactStorage->forEachContact(useOpenMP, kernel::SelectAll(), contactAccessor,
                                        hcsits_initContacts, contactAccessor, particleAccessor);
         timing.stop("Init contacts");
         timing.start("Init particles");
         particleStorage->forEachParticle(useOpenMP, kernel::SelectAll(), particleAccessor,
                                          hcsits_initParticles, particleAccessor, dt);
         timing.stop("Init particles");

         timing.start("Velocity update");
         VelocityUpdateNotification::Parameters::relaxationParam = real_t(1.0); // must be set to 1.0 such that dv and dw caused by external forces and torques are not falsely altered
         reductionKernel.operator()<VelocityCorrectionNotification>(*particleStorage);
         broadcastKernel.operator()<VelocityUpdateNotification>(*particleStorage);
         timing.stop("Velocity update");

         VelocityUpdateNotification::Parameters::relaxationParam = hcsits_relaxationParameter;
         for(uint_t i = uint_t(0); i < hcsits_numberOfIterations; i++){
            timing.start("Relaxation step");
            contactStorage->forEachContact(useOpenMP, kernel::SelectAll(), contactAccessor,
                                           hcsits_relaxationStep, contactAccessor, particleAccessor, dt);
            timing.stop("Relaxation step");
            timing.start("Velocity update");
            reductionKernel.operator()<VelocityCorrectionNotification>(*particleStorage);
            broadcastKernel.operator()<VelocityUpdateNotification>(*particleStorage);
            timing.stop("Velocity update");
         }
         timing.start("Integration");
         particleStorage->forEachParticle(useOpenMP, kernel::SelectAll(), particleAccessor,
                                          hcsits_integration, particleAccessor, dt);
         timing.stop("Integration");
         timing.stop("HCSITS");
      }
      else if(solver == "DEM")
      {
         timing.start("DEM");
         timing.start("Collision");
         contactStorage->forEachContact(useOpenMP, kernel::SelectAll(), contactAccessor,
                                        [&dem_collision, coefficientOfRestitution, dem_collisionTime, dem_kappa, dt](size_t c, data::ContactAccessor &ca, data::ParticleAccessorWithBaseShape &pa){
                                           auto idx1 = ca.getId1(c);
                                           auto idx2 = ca.getId2(c);
                                           auto meff = real_t(1) / (pa.getInvMass(idx1) + pa.getInvMass(idx2));

                                           dem_collision.setStiffnessAndDamping(0,0,coefficientOfRestitution, dem_collisionTime, dem_kappa, meff);

                                           dem_collision(idx1, idx2, pa, ca.getPosition(c), ca.getNormal(c), ca.getDistance(c), dt);
                                        },
                                        contactAccessor, particleAccessor);
         timing.stop("Collision");


         timing.start("Apply gravity");
         particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                             [reducedGravitationalAcceleration](const size_t idx, auto &ac){addForceAtomic(idx, ac, Vec3(0_r,0_r,-reducedGravitationalAcceleration) * ac.getMass(idx));}, particleAccessor);
         timing.stop("Apply gravity");

         timing.start("Apply horizontal forcing");
         if(applyHorizontalForcing)
         {
            particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                             [reducedHorizontalAcceleration](const size_t idx, auto &ac){addForceAtomic(idx, ac, reducedHorizontalAcceleration * ac.getMass(idx));}, particleAccessor);
         }
         timing.stop("Apply horizontal forcing");

         timing.start("Reduce");
         reduceAndSwapContactHistory(*particleStorage);
         reductionKernel.operator()<ForceTorqueNotification>(*particleStorage);
         timing.stop("Reduce");

         if(isDampingActive && damping_method == "force")
         {
            // apply damping on (all) forces, see e.g. Zhao, 2017 - "Particle shape effects on fabric of granular random packing", and YADE docu (https://yade-dem.org/doc/formulation.html#numerical-damping)
            particleStorage->forEachParticle(useOpenMP, kernel::SelectAll(), particleAccessor,
                                             [damping_forceDampingCoefficient, dt](size_t idx, data::ParticleAccessorWithBaseShape &ac){
                                                auto force = ac.getForce(idx);
                                                auto velEstimated = ac.getLinearVelocity(idx) + force * ac.getInvMass(idx) * dt; // = velocity integration in SemiImplicitEuler, might need to be changed for other integrator
                                                Vec3 dampingForce(real_c(mesa_pd::sgn(force[0] * velEstimated[0])) * force[0],
                                                                  real_c(mesa_pd::sgn(force[1] * velEstimated[1])) * force[1],
                                                                  real_c(mesa_pd::sgn(force[2] * velEstimated[2])) * force[2]);
                                                dampingForce *= -damping_forceDampingCoefficient;
                                                ac.getForceRef(idx) += dampingForce; },
                                             particleAccessor);
         }

         timing.start("Integration");
         particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                          dem_integration, particleAccessor);
         timing.stop("Integration");

         timing.stop("DEM");

      }

      if(limitVelocity > 0_r)
      {
         timing.start("Velocity limiting");
         particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                          [limitVelocity](const size_t idx, auto &ac){
                                             auto velMagnitude = ac.getLinearVelocity(idx).length();
                                             if(velMagnitude > limitVelocity) ac.getLinearVelocityRef(idx) *= (limitVelocity / velMagnitude );
                                          }, particleAccessor);
         timing.stop("Velocity limiting");
      }


      timing.start("Sync");
      syncCall();
      particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                       associateToBlock, particleAccessor);
      timing.stop("Sync");

      timing.start("Evaluate particles");

      timing.start("Property evaluation");
      if(propertyEvaluationSpacing > 0 && timestep % propertyEvaluationSpacing == 0) particleInfo = evaluateParticleInfo(particleAccessor);

      if(fixBottomParticles){
         porosityEvaluator.clear();
         particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                          porosityEvaluator, particleAccessor);
         porosityEvaluator.evaluate();
         real_t estimatedPackingHeight = porosityEvaluator.estimatePackingHeight();
         if(estimatedPackingHeight > horizontalForcing_fixingHeight)
         {
            WALBERLA_LOG_INFO_ON_ROOT("Fixing particles below " << horizontalForcing_fixingHeight);
            fixBottomParticles = false; // only do it once
            fixParticlesBelowFixingHeight(particleAccessor, horizontalForcing_fixingHeight);
         }
      }
      timing.stop("Property evaluation");


      if(particleInfo.particleVolume * particleDensity < totalParticleMass)
      {

         timing.start("Generation");
         // check if generation
         if(particleInfo.maximumHeight < (generationHeightRatioStart * simulationDomain.zSize() - std::max(maxGenerationParticleDiameter,generationSpacing)) || currentTime - timeLastCreation > maximumTimeBetweenCreation)
         {
            particleCreator.createParticles( std::max(minGenerationHeight, generationHeightRatioStart * simulationDomain.zMax()),
                                             std::min(maxGenerationHeight, generationHeightRatioEnd * simulationDomain.zMax()),
                                             generationSpacing, diameterGenerator, shapeGenerator, initialVelocity,
                                             maximumAllowedInteractionRadius, generationPositionCheckerLargeDiameterThreshold);

            particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                             associateToBlock, particleAccessor);

            if(useNextNeighborSync){ syncCall(); }
            else { for(uint_t i = 0; i < uint_c(std::ceil(maxParticleDiameter/smallestBlockSize)); ++i) syncCall(); }

            timeLastCreation = currentTime;

            // write current particle distribution info
            particleHistogram.clear();
            particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                             particleHistogram, particleAccessor);
            particleHistogram.evaluate();
            WALBERLA_LOG_INFO_ON_ROOT(particleHistogram);

            // reevaluate particle properties
            particleInfo = evaluateParticleInfo(particleAccessor);

         }
         timing.stop("Generation");
      } else if(shaking)
      {
         timing.start("Shaking");
         // apply shaking
         if(timeEndShaking < 0_r){
            if(!isShakingActive)
            {
               isShakingActive = true;
               timeBeginShaking = currentTime;
               timeEndShaking = currentTime + shaking_duration;
               WALBERLA_LOG_INFO_ON_ROOT("Beginning of shaking at time " << currentTime << " s for " << shaking_duration << " s.");
            } else {
               //timeEndShaking = real_c(std::ceil((currentTime + shaking_duration) / shaking_period)) * shaking_period; // make sure to shake only full periods
               timeEndShaking = currentTime + shaking_duration; // since its unclear if full periods are really necessary and actually "improve" results, we skip this here
               WALBERLA_LOG_INFO_ON_ROOT("Continue of shaking at time " << currentTime << " s until time " << timeEndShaking << " s.");
            }
         }

         if(currentTime > timeEndShaking)
         {
            WALBERLA_LOG_INFO_ON_ROOT("Ending of shaking at time " << currentTime << " s.");
            shaking = false;
            isShakingActive = false;
         }
         timing.stop("Shaking");

      } else
      {
         timing.start("Damping");

         if(timeBeginDamping < 0_r){
            timeBeginDamping = currentTime;
            WALBERLA_LOG_INFO_ON_ROOT("Beginning of damping at time " << timeBeginDamping << " s.");
            isDampingActive = true;

            if(applyHorizontalForcing)
            {
               applyHorizontalForcing = false;
               WALBERLA_LOG_INFO_ON_ROOT("Switching off horizontal forcing");
            }
         }

         // apply damping
         if(damping_method == "force")
         {
            // was carried out before DEM integration

         } else if(damping_method == "velocity")
         {
            particleStorage->forEachParticle(useOpenMP, kernel::SelectAll(), particleAccessor,
                                             [velocityDampingFactor](size_t idx, data::ParticleAccessorWithBaseShape &ac){
                                                ac.getLinearVelocityRef(idx) *= velocityDampingFactor;
                                                ac.getAngularVelocityRef(idx) *= velocityDampingFactor;},
                                             particleAccessor);
         }


         // check if termination
         if(currentTime - timeBeginDamping > minimalTerminalRunTime)
         {
            if(currentTime - timeLastTerminationCheck > terminationCheckingSpacing)
            {
               if(particleInfo.maximumVelocity < terminalVelocity)
               {
                  WALBERLA_LOG_INFO_ON_ROOT("Reached terminal max velocity - terminating.");
                  terminateSimulation = true;
               }

               real_t relDiffAvgHeight = std::abs(particleInfo.heightOfMass - oldAvgParticleHeight) / oldAvgParticleHeight;
               real_t relDiffMaxHeight = std::abs(particleInfo.maximumHeight - oldMaxParticleHeight) / oldMaxParticleHeight;
               if(relDiffMaxHeight < 10_r * terminalRelativeHeightChange && relDiffAvgHeight < terminalRelativeHeightChange)
               {
                  // check of max height has to be included to avoid early termination if only little mass is created per generation step
                  WALBERLA_LOG_INFO_ON_ROOT("Reached converged maximum and mass-averaged height - terminating.");
                  terminateSimulation = true;
               }

               oldAvgParticleHeight = particleInfo.heightOfMass;
               oldMaxParticleHeight = particleInfo.maximumHeight;
               timeLastTerminationCheck = currentTime;
            }
         }

         if(maximalTerminalRunTime > 0 && currentTime - timeBeginDamping > maximalTerminalRunTime)
         {
            WALBERLA_LOG_INFO_ON_ROOT("Reached terminal run time - terminating.");
            terminateSimulation = true;
         }

         timing.stop("Damping");
      }
      timing.stop("Evaluate particles");

      if(( infoSpacing > 0 && timestep % infoSpacing == 0) || (loggingSpacing > 0 && timestep % loggingSpacing == 0))
      {
         timing.start("Evaluate infos");
         auto contactInfo = evaluateContactInfo(contactAccessor, particleAccessor);

         porosityEvaluator.clear();
         particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                          porosityEvaluator, particleAccessor);
         porosityEvaluator.evaluate();
         real_t estimatedPorosity = porosityEvaluator.estimateTotalPorosity();
         real_t estimatedPackingHeight = porosityEvaluator.estimatePackingHeight();

         if(loggingSpacing > 0 && timestep % loggingSpacing == 0)
         {
            loggingWriter(currentTime, particleInfo, contactInfo, estimatedPorosity);
         }

         if(infoSpacing > 0 && timestep % infoSpacing == 0) {
            WALBERLA_LOG_INFO_ON_ROOT("t = " << timestep << " = " << currentTime << " s");
            WALBERLA_LOG_INFO_ON_ROOT(particleInfo << " => " << particleInfo.particleVolume * particleDensity << " kg" << ", current porosity = " << estimatedPorosity << ", packing height = " << estimatedPackingHeight);
            real_t ensembleAverageDiameter = diameterFromSphereVolume(particleInfo.particleVolume / real_c(particleInfo.numParticles));
            WALBERLA_LOG_INFO_ON_ROOT(contactInfo << " => " << contactInfo.maximumPenetrationDepth / ensembleAverageDiameter * real_t(100) << "% of avg diameter " << ensembleAverageDiameter);
         }

         timing.stop("Evaluate infos");

      }

      if(performanceLoggingSpacing > 0 && timestep % performanceLoggingSpacing == 0)
      {
         auto reducedTT = timing.getReduced();
         WALBERLA_LOG_INFO_ON_ROOT(reducedTT);
      }

      timing.start("Checkpointing");
      if(checkPointing_spacing > 0 && timestep > 0 && timestep % checkPointing_spacing == 0)
      {
         WALBERLA_LOG_INFO_ON_ROOT("Writing checkpointing file " << checkPointFileName_mesapd << " in time step " << timestep);
         forest->saveBlockData( checkPointFileName_mesapd, particleStorageID );

         std::vector<config::Config::Block*> cfgBlocks;
         cfg->getWritableGlobalBlock().getWritableBlocks("CheckpointRestart",cfgBlocks,1,1);
         auto checkpointRestartConfMod = cfgBlocks[0];
         checkpointRestartConfMod->setOrAddParameter("existingUID", uniqueFileIdentifier);
         checkpointRestartConfMod->setOrAddParameter("timestep", std::to_string(timestep+1));
         checkpointRestartConfMod->setOrAddParameter("isDampingActive", std::to_string(isDampingActive));
         checkpointRestartConfMod->setOrAddParameter("isShakingActive", std::to_string(isShakingActive));
         checkpointRestartConfMod->setOrAddParameter("oldAvgParticleHeight", std::to_string(oldAvgParticleHeight));
         checkpointRestartConfMod->setOrAddParameter("oldMaxParticleHeight", std::to_string(oldMaxParticleHeight));
         checkpointRestartConfMod->setOrAddParameter("timeLastTerminationCheck", std::to_string(timeLastTerminationCheck));
         checkpointRestartConfMod->setOrAddParameter("timeLastCreation", std::to_string(timeLastCreation));
         checkpointRestartConfMod->setOrAddParameter("timeBeginShaking", std::to_string(timeBeginShaking));
         checkpointRestartConfMod->setOrAddParameter("timeEndShaking", std::to_string(timeEndShaking));
         checkpointRestartConfMod->setOrAddParameter("timeBeginDamping", std::to_string(timeBeginDamping));

         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file;
            std::string configFileCopyName = checkPointing_folder + "/" += uniqueFileIdentifier + ".cfg";
            file.open( configFileCopyName.c_str() );
            file << *cfg;
            file.close();
         }

      }
      timing.stop("Checkpointing");

      timing.start("Load balancing");
      if(loadBalancing_spacing > 0 && timestep % loadBalancing_spacing == 0)
      {

         domain::createWithNeighborhood( particleAccessor, *forest, *loadBalancingInfoCollection );

         mesa_pd::mpi::ClearGhostOwnerSync CGOS;
         CGOS(particleAccessor);

         contactStorage->clear();
         forest->refresh();
         domain->refresh();

         linkedCells = std::make_unique<data::LinkedCells>(domain->getUnionOfLocalAABBs().getExtended(linkedCellWidth), linkedCellWidth );

         particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor, associateToBlock, particleAccessor);
         if(useNextNeighborSync){ syncCall(); }
         else { for(uint_t i = 0; i < uint_c(std::ceil(maxParticleDiameter/smallestBlockSize)); ++i) syncCall(); }

         sorting::LinearizedCompareFunctor linearSorting(linkedCells->domain_, linkedCells->numCellsPerDim_);
         particleStorage->sort(linearSorting);

      }
      timing.stop("Load balancing");

      ++timestep;
   }

   simulationTimer.end();

   particleHistogram.clear();
   particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                    particleHistogram, particleAccessor);
   particleHistogram.evaluate();
   WALBERLA_LOG_INFO_ON_ROOT(particleHistogram);

   porosityEvaluator.clear();
   particleStorage->forEachParticle(useOpenMP, kernel::SelectLocal(), particleAccessor,
                                    porosityEvaluator, particleAccessor);
   porosityEvaluator.evaluate();
   real_t estimatedFinalPorosity = porosityEvaluator.estimateTotalPorosity();
   WALBERLA_LOG_INFO_ON_ROOT("Estimated total porosity based on layers = " << estimatedFinalPorosity);

   std::string porosityFileName = porosityProfileFolder + "/" +  uniqueFileIdentifier + "_layers.txt";
   WALBERLA_LOG_INFO_ON_ROOT("Writing porosity profile file to " << porosityFileName);
   porosityEvaluator.printToFile(porosityFileName);

   ContactInfoPerHorizontalLayerEvaluator contactEvaluator(evaluationLayerHeight, simulationDomain);
   contactStorage->forEachContact(useOpenMP, kernel::SelectAll(), particleAccessor,
                                  contactEvaluator, contactAccessor);
   contactEvaluator.evaluate();
   std::string contactInfoFileName = porosityProfileFolder + "/" +  uniqueFileIdentifier + "_contact_layers.txt";
   WALBERLA_LOG_INFO_ON_ROOT("Writing contact info profile file to " << contactInfoFileName);
   contactEvaluator.printToFile(contactInfoFileName);

   auto reducedTT = timing.getReduced();
   WALBERLA_LOG_INFO_ON_ROOT(reducedTT);

   bool logToProcessLocalFiles = false;
   std::string particleInfoFileName = porosityProfileFolder + "/" +  uniqueFileIdentifier + "_particle_info";
   if(logToProcessLocalFiles)
   {
      particleInfoFileName += "_" + std::to_string(walberla::mpi::MPIManager::instance()->rank()) + ".txt";
      WALBERLA_LOG_INFO_ON_ROOT("Writing particle info file to process local files like " << particleInfoFileName);

   } else {
      particleInfoFileName += ".txt";
      WALBERLA_LOG_INFO_ON_ROOT("Writing particle info file to " << particleInfoFileName);
   }
   auto particleInfoString = assembleParticleInformation(*particleStorage, particleSizeEvaluator, 12);
   writeParticleInformationToFile(particleInfoFileName, particleInfoString, logToProcessLocalFiles);

   // write to sqlite data base
   particleInfo = evaluateParticleInfo(particleAccessor);
   auto contactInfo = evaluateContactInfo(contactAccessor, particleAccessor);

   WALBERLA_ROOT_SECTION() {
      std::map<std::string, walberla::int64_t> sql_integerProperties;
      std::map<std::string, double> sql_realProperties;
      std::map<std::string, std::string> sql_stringProperties;
      addConfigToDatabase(*cfg, sql_integerProperties, sql_realProperties, sql_stringProperties);

      // store particle info
      sql_integerProperties["numParticles"] = int64_c(particleInfo.numParticles);
      sql_realProperties["maxParticlePosition"] = double(particleInfo.maximumHeight);
      sql_realProperties["particleVolume"] = double(particleInfo.particleVolume);

      // store contact info
      sql_integerProperties["numContacts"] = int64_c(contactInfo.numContacts);
      sql_realProperties["maxPenetrationDepth"] = double(contactInfo.maximumPenetrationDepth);
      sql_realProperties["avgPenetrationDepth"] = double(contactInfo.averagePenetrationDepth);

      // other info
      //sql_realProperties["simulationTime"] = double(reducedTT["Simulation"].total());
      sql_realProperties["simulationTime"] = double(simulationTimer.total());
      sql_integerProperties["numProcesses"] = int64_c(walberla::mpi::MPIManager::instance()->numProcesses());
      sql_integerProperties["timesteps"] = int64_c(timestep);
      sql_stringProperties["file_identifier"] = uniqueFileIdentifier;
      sql_realProperties["generationSpacing"] = double(generationSpacing);
      std::string histogramData = "";
      for (auto h : particleHistogram.getMassFractionHistogram()) histogramData += std::to_string(h) + " ";
      sql_stringProperties["evaluation_histogramData"] = histogramData;
      std::string numberHistogramData = "";
      for (auto h : particleHistogram.getNumberHistogram()) numberHistogramData += std::to_string(h) + " ";
      sql_stringProperties["evaluation_numberHistogramData"] = numberHistogramData;
      sql_integerProperties["singleShape"] = (shapeGenerator->generatesSingleShape()) ? 1 : 0;
      sql_realProperties["maxAllowedInteractionRadius"] = double(maximumAllowedInteractionRadius);

      for(uint_t i = 0; i < particleHistogram.getNumberOfShapeEvaluators(); ++i)
      {
         std::string shapeBins = "";
         for (auto b : particleHistogram.getShapeBins(i)) shapeBins += std::to_string(b) + " ";
         sql_stringProperties["evaluation_"+std::get<0>(particleHistogram.getShapeEvaluator(i))+"_bins"] = shapeBins;

         std::string shapeHistogram = "";
         for (auto h : particleHistogram.getShapeHistogram(i)) shapeHistogram += std::to_string(h) + " ";
         sql_stringProperties["evaluation_"+std::get<0>(particleHistogram.getShapeEvaluator(i))+"_histogramData"] = shapeHistogram;
      }

      WALBERLA_LOG_INFO_ON_ROOT("Storing run and timing data in sql database file " << sqlDBFileName);
      auto sql_runID = sqlite::storeRunInSqliteDB(sqlDBFileName, sql_integerProperties, sql_stringProperties, sql_realProperties);
      sqlite::storeTimingTreeInSqliteDB(sqlDBFileName, sql_runID, reducedTT, "Timing");
   }

   if(!vtkFinalFolder.empty())
   {

      WALBERLA_LOG_INFO_ON_ROOT("Writing final VTK file to folder " << vtkFinalFolder);
      if(particleShape.find("Mesh") != std::string::npos)
      {
         mesa_pd::MeshParticleVTKOutput< mesh::PolyMesh > finalMeshParticleVTK(particleStorage, uniqueFileIdentifier, uint_t(1), vtkFinalFolder);
         finalMeshParticleVTK.addFaceOutput< data::SelectParticleUid >("UID");
         finalMeshParticleVTK.addVertexOutput< data::SelectParticleInteractionRadius >("InteractionRadius");
         finalMeshParticleVTK.addFaceOutput< data::SelectParticleLinearVelocity >("LinearVelocity");
         finalMeshParticleVTK.addVertexOutput< data::SelectParticlePosition >("Position");
         if(includeNumberOfContactsInVTK) finalMeshParticleVTK.addVertexOutput< data::SelectParticleNumContacts>("numContacts");
         finalMeshParticleVTK.addVertexDataSource(surfaceVelDataSource);
         finalMeshParticleVTK.setParticleSelector(vtkParticleSelector);
         finalMeshParticleVTK(particleAccessor);
      }
      else {
         auto finalParticleVtkWriter = walberla::vtk::createVTKOutput_PointData(particleVtkOutput, uniqueFileIdentifier, uint_t(1), vtkFinalFolder, "final");
         finalParticleVtkWriter->write();
      }

      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         std::string configFileCopyName = vtkFinalFolder + "/" + uniqueFileIdentifier + ".cfg";
         WALBERLA_LOG_INFO_ON_ROOT("Storing config file as " << configFileCopyName);
         file.open( configFileCopyName.c_str() );
         file << *cfg;
         file.close();
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT("Simulation terminated successfully");

   return EXIT_SUCCESS;
}
} // namespace mesa_pd
} // namespace walberla

int main( int argc, char* argv[] ) {
   return walberla::mesa_pd::main( argc, argv );
}