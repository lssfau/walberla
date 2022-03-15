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

#pragma once

#include "mesa_pd/data/ContactAccessor.h"
#include "mesa_pd/data/ParticleStorage.h"

#include "Utility.h"
#include "ShapeGeneration.h"

namespace walberla {
namespace mesa_pd {


struct ParticleInfo
{
   real_t averageVelocity = 0_r;
   real_t maximumVelocity = 0_r;
   uint_t numParticles = 0;
   real_t maximumHeight = 0_r;
   real_t particleVolume = 0_r;
   real_t heightOfMass = 0_r;
   real_t kinEnergy = 0_r; // = m * v^2/2
   real_t meanCoordinationNumber = 0_r;

   void allReduce()
   {
      walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(averageVelocity, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(maximumVelocity, walberla::mpi::MAX);
      walberla::mpi::allReduceInplace(maximumHeight, walberla::mpi::MAX);
      walberla::mpi::allReduceInplace(particleVolume, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(heightOfMass, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(kinEnergy, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(meanCoordinationNumber, walberla::mpi::SUM);

      averageVelocity /= real_c(numParticles);
      heightOfMass /= particleVolume;
      meanCoordinationNumber /= real_c(numParticles);
   }
};

std::ostream &operator<<(std::ostream &os, ParticleInfo const &m) {
   return os << "Particle Info: uAvg = " << m.averageVelocity << ", uMax = " << m.maximumVelocity
             << ", numParticles = " << m.numParticles << ", zMax = " << m.maximumHeight << ", Vp = "
             << m.particleVolume << ", zMass = " << m.heightOfMass << ", kin. energy = " << m.kinEnergy << ", MCN = " << m.meanCoordinationNumber;
}

template< typename Accessor_T>
ParticleInfo evaluateParticleInfo(const Accessor_T & ac)
{

   ParticleInfo info;
   for(uint_t i = 0; i < ac.size(); ++i)
   {
      if (isSet(ac.getFlags(i), data::particle_flags::GHOST)) continue;
      if (isSet(ac.getFlags(i), data::particle_flags::GLOBAL)) continue;

      ++info.numParticles;
      real_t velMagnitude = ac.getLinearVelocity(i).length();
      real_t particleVolume = ac.getBaseShape(i)->getVolume();
      real_t particleMass = ac.getBaseShape(i)->getMass();
      real_t height = ac.getPosition(i)[2];
      info.averageVelocity += velMagnitude;
      info.maximumVelocity = std::max(info.maximumVelocity, velMagnitude);
      info.maximumHeight = std::max(info.maximumHeight, height);
      info.particleVolume += particleVolume;
      info.heightOfMass += particleVolume*height;
      info.kinEnergy += 0.5_r * particleMass * velMagnitude * velMagnitude;
      info.meanCoordinationNumber += real_c(ac.getNumContacts(i));
   }

   info.allReduce();

   return info;
}

struct ContactInfo
{
   real_t averagePenetrationDepth = 0_r;
   real_t maximumPenetrationDepth = 0_r;
   uint_t numContacts = 0;

   void allReduce()
   {
      walberla::mpi::allReduceInplace(numContacts, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(averagePenetrationDepth, walberla::mpi::SUM);
      walberla::mpi::allReduceInplace(maximumPenetrationDepth, walberla::mpi::MAX);

      if(numContacts > 0) averagePenetrationDepth /= real_t(numContacts);
   }
};

std::ostream &operator<<(std::ostream &os, ContactInfo const &m) {
   return os << "Contact Info: numContacts = " << m.numContacts << ", deltaAvg = " << m.averagePenetrationDepth << ", deltaMax = " << m.maximumPenetrationDepth;
}

ContactInfo evaluateContactInfo(const data::ContactAccessor & ca)
{
   ContactInfo info;
   for(uint_t i = 0; i < ca.size(); ++i)
   {
      real_t penetrationDepth = -ca.getDistance(i);
      info.maximumPenetrationDepth = std::max(info.maximumPenetrationDepth, penetrationDepth);
      if(penetrationDepth > 0_r) {
         info.averagePenetrationDepth += penetrationDepth;
         ++info.numContacts;
      }
   }
   info.allReduce();
   return info;
}

class SizeEvaluator
{
public:
   explicit SizeEvaluator(ScaleMode scaleMode) : scaleMode_(scaleMode) { }

   struct SizeInfo
   {
      real_t size;
      Vec3 shapeSemiAxes;
      real_t volume;
   };

   SizeInfo get(data::BaseShape & bs)
   {
      SizeInfo sizeInfo;
      sizeInfo.volume = bs.getVolume();

      auto shapeType = bs.getShapeType();
      if(shapeType == data::ConvexPolyhedron::SHAPE_TYPE) {
         auto convexPolyhedron = static_cast<data::ConvexPolyhedron*>(&bs);
         sizeInfo.shapeSemiAxes = semiAxesFromInertiaTensor(convexPolyhedron->getInertiaBF(), convexPolyhedron->getMass());
      } else if(shapeType == data::Ellipsoid::SHAPE_TYPE) {
         auto ellipsoid = static_cast<data::Ellipsoid*>(&bs);
         sizeInfo.shapeSemiAxes = ellipsoid->getSemiAxes();
      } else if(shapeType == data::Sphere::SHAPE_TYPE) {
         auto sphere = static_cast<data::Sphere*>(&bs);
         sizeInfo.shapeSemiAxes = Vec3(sphere->getRadius());
      } else {
         WALBERLA_ABORT("Shape handling not implemented!");
      }

      switch(scaleMode_)
      {
         case ScaleMode::sphereEquivalent:
         {
            sizeInfo.size = diameterFromSphereVolume(sizeInfo.volume);
            break;
         }
         case ScaleMode::sieveLike:
         {
            sizeInfo.size = sizeFromSemiAxes(sizeInfo.shapeSemiAxes);
            break;
         }
         default: WALBERLA_ABORT("Unknown shape scale mode");
      }
      return sizeInfo;
   }

private:
   ScaleMode scaleMode_;
};

// shape evaluation functions
real_t getFlatnessFromSemiAxes(Vec3 semiAxes)
{
   sortVector(semiAxes);
   return semiAxes[0] / semiAxes[1];
}

real_t getElongationFromSemiAxes(Vec3 semiAxes)
{
   sortVector(semiAxes);
   return semiAxes[1] / semiAxes[2];
}

real_t getEquancyFromSemiAxes(Vec3 semiAxes)
{
   sortVector(semiAxes);
   return semiAxes[0] / semiAxes[2];
}

class ParticleHistogram
{
public:
   ParticleHistogram(std::vector<real_t> sizeBins,
                     SizeEvaluator sizeEvaluator,
                     std::vector<std::vector<real_t>> shapeBins,
                     std::vector<std::tuple<std::string,std::function<real_t(Vec3)>>> shapeEvaluators) :
   sizeBins_(sizeBins), sizeEvaluator_(sizeEvaluator),
   shapeBins_(shapeBins), shapeEvaluators_(shapeEvaluators) {
      bool areBinsAscending = (sizeBins[1] - sizeBins[0]) > real_t(0);
      if(!areBinsAscending) std::reverse(sizeBins_.begin(), sizeBins_.end());

      massFractionHistogram_ = std::vector<real_t>(sizeBins.size() + 1,real_t(0));
      numberHistogram_ = std::vector<uint_t>(sizeBins.size() + 1,0);

      WALBERLA_CHECK_EQUAL(shapeBins.size(), shapeEvaluators.size(), "Different number of shape evaluators and bins!");

      shapeHistograms_ = std::vector<std::vector<uint_t>>(shapeBins.size());
      for(uint_t i = 0; i < shapeBins.size(); ++i)
      {
         shapeHistograms_[i] = std::vector<uint_t>(shapeBins[i].size() + 1,0);
      }
   }

   template<typename Accessor_T>
   void operator()(size_t p_idx, Accessor_T& ac)
   {
      if (isSet(ac.getFlags(p_idx), data::particle_flags::GLOBAL) || isSet(ac.getFlags(p_idx), data::particle_flags::GHOST)) return;

      auto sizeInfo = sizeEvaluator_.get(*ac.getBaseShape(p_idx));

      auto size = sizeInfo.size;
      auto volume = sizeInfo.volume;
      auto semiAxes = sizeInfo.shapeSemiAxes;

      auto sizeHistogramIdx = getHistogramIdx(sizeBins_, size);
      massFractionHistogram_[sizeHistogramIdx] += volume;
      numberHistogram_[sizeHistogramIdx]++;

      for(uint_t i = 0; i < shapeBins_.size(); ++i)
      {
         auto shapeParam = std::get<1>(shapeEvaluators_[i])(semiAxes);
         auto shapeHistogramIdx = getHistogramIdx(shapeBins_[i], shapeParam);
         shapeHistograms_[i][shapeHistogramIdx]++;
      }
   }

   void evaluate()
   {
      walberla::mpi::reduceInplace(numberHistogram_, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(massFractionHistogram_, walberla::mpi::SUM);
      for(uint_t i = 0; i < shapeHistograms_.size(); ++i) walberla::mpi::reduceInplace(shapeHistograms_[i], walberla::mpi::SUM);
      WALBERLA_ROOT_SECTION()
      {
         auto total = std::accumulate(massFractionHistogram_.begin(), massFractionHistogram_.end(), real_t(0));
         for(auto& h: massFractionHistogram_) h /= total;
      }
   }

   void clear()
   {
      for(auto& h: massFractionHistogram_) h = real_t(0);
      for(auto& h: numberHistogram_) h = 0;
      for(uint_t i = 0; i < shapeHistograms_.size(); ++i)
      {
         for(auto& h: shapeHistograms_[i]) h = 0;
      }
   }

   std::vector<real_t> getSizeBins() const {return sizeBins_;}
   std::vector<real_t> getMassFractionHistogram() const {return massFractionHistogram_;}
   std::vector<uint_t> getNumberHistogram() const {return numberHistogram_;}

   uint_t getNumberOfShapeEvaluators() const {return shapeEvaluators_.size();}
   std::vector<real_t> getShapeBins(uint_t idx) const {return shapeBins_[idx];}
   std::vector<uint_t> getShapeHistogram(uint_t idx) const {return shapeHistograms_[idx];}
   std::tuple<std::string, std::function<real_t(Vec3)>> getShapeEvaluator(uint_t idx) const {return shapeEvaluators_[idx];}

private:

   uint_t getHistogramIdx(const std::vector<real_t> & bins, real_t value) const
   {
      auto numBins = bins.size();
      if(value >= bins[numBins-1]) {
         return numBins;
      }
      else {
         for(uint_t b = 0; b < numBins; ++b){
            if(value < bins[b])
            {
               return b;
            }
         }
      }
      WALBERLA_ABORT("No valid histogram bin found for value " << value);
      return 0;
   }


   std::vector<real_t> sizeBins_;
   SizeEvaluator sizeEvaluator_;
   std::vector<std::vector<real_t>> shapeBins_;
   std::vector<std::tuple<std::string,std::function<real_t(Vec3)>>> shapeEvaluators_; // vector of parameter label and function pairs

   std::vector<real_t> massFractionHistogram_;
   std::vector<uint_t> numberHistogram_;
   std::vector<std::vector<uint_t>> shapeHistograms_;
};

std::ostream &operator<<(std::ostream &os, ParticleHistogram const &ph) {
   auto bins = ph.getSizeBins();
   os << "Size bins: | 0 - ";
   for(auto b: bins) os << b*real_t(1e3) << " | " << b*real_t(1e3) << " - ";
   os << "infty | \n";

   auto hist = ph.getMassFractionHistogram();
   os << "Mass Fraction Hist: | ";
   for(auto h: hist) os << h << " | ";

   os << "\n";

   auto numHist = ph.getNumberHistogram();
   os << "Number Hist: | ";
   for(auto h: numHist) os << h << " | ";

   for(uint_t i = 0; i < ph.getNumberOfShapeEvaluators(); ++i)
   {
      os << "\n\nShape parameter " << i << " ("<< std::get<0>(ph.getShapeEvaluator(i)) << ") bins: -infty - ";
      for(auto b: ph.getShapeBins(i)) os << b << " | " << b << " - ";
      os << "infty \n";
      os << "Number Hist: | ";
      for(auto h: ph.getShapeHistogram(i)) os << h << " | ";
   }

   return os;
}

class PorosityPerHorizontalLayerEvaluator
{
public:
   PorosityPerHorizontalLayerEvaluator(real_t layerHeight, const AABB & simulationDomain, std::string domainSetup)
         : layerHeight_(layerHeight)
   {
      if(domainSetup == "container") {
         layerVolume_ = simulationDomain.xSize() * simulationDomain.xSize() * math::pi / real_t(4) * layerHeight;
      }
      else {
         layerVolume_ = simulationDomain.xSize() * simulationDomain.ySize() * layerHeight;
      }
      porosityPerLayer_ = std::vector<real_t>(uint_c(simulationDomain.zSize() / layerHeight), real_t(0));
      radiusPerLayer_ = std::vector<real_t>(uint_c(simulationDomain.zSize() / layerHeight), real_t(0));
   }

   template<typename Accessor_T>
   void operator()(size_t p_idx, Accessor_T& ac)
   {
      if (isSet(ac.getFlags(p_idx), data::particle_flags::GLOBAL)) return;

      real_t volume = ac.getBaseShapeRef(p_idx)->getVolume();
      real_t radius = radiusFromSphereVolume(volume);
      real_t zPos = ac.getPosition(p_idx)[2];
      int heightIdxBegin = int(std::floor((zPos - radius)/layerHeight_));
      int heightIdxEnd = int((zPos + radius)/layerHeight_) + 1;
      for(int i = heightIdxBegin; i < heightIdxEnd; ++i)
      {
         real_t sphereSegmentVolume = calculateSphericalSegmentVolume(real_c(i) * layerHeight_ - zPos,
                                                                      real_c(i+1) * layerHeight_ - zPos,
                                                                      radius);
         uint_t writeIdx = uint_c(std::min(std::max(i, 0), int(porosityPerLayer_.size() - 1) ) );
         // i could be negative if overlap with bottom plane -> add this to idx 0 to not lose volume
         // i could also be too large if evaluated during the simulation so cap its max value to avoid seg faults
         porosityPerLayer_[writeIdx] += sphereSegmentVolume;
         radiusPerLayer_[writeIdx] += sphereSegmentVolume * radius;

      }
   }

   void evaluate()
   {
      walberla::mpi::reduceInplace(porosityPerLayer_, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(radiusPerLayer_, walberla::mpi::SUM);
      WALBERLA_ROOT_SECTION()
      {
         for(uint_t i = 0; i < radiusPerLayer_.size(); ++i) radiusPerLayer_[i] /= porosityPerLayer_[i]; // = volume-averaged radius per layer
         for(auto& p: porosityPerLayer_) p /= layerVolume_;
      }
   }

   void clear()
   {
      for(uint_t i = 0; i < porosityPerLayer_.size(); ++i)
      {
         porosityPerLayer_[i] = 0_r;
         radiusPerLayer_[i] = 0_r;
      }
   }

   void printToFile(std::string fileName)
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName.c_str() );
         file << "#\t z\t porosity\t avgRadius\n";

         for(uint_t i = 0; i < porosityPerLayer_.size(); ++i)
         {
            file << layerHeight_ * (real_t(0.5)+real_c(i)) << " " << porosityPerLayer_[i] << " " << radiusPerLayer_[i] << "\n";
         }
         file.close();
      }
   }

   real_t estimateTotalPorosity()
   {
      /*
      uint_t zMaxIdx = 0;
      // find first empty index
      while(porosityPerLayer_[zMaxIdx] > 0_r)
      {
         ++zMaxIdx;
      }
      zMaxIdx = uint_c(0.95_r * real_c(zMaxIdx)); // remove top 5% according to Schruff, 2018
      return 1_r - std::accumulate(porosityPerLayer_.begin(), std::next(porosityPerLayer_.begin(), zMaxIdx), 0_r) / real_c(zMaxIdx+1);
      */

      const real_t cutOffPorosity = 0.5_r; // some value
      const real_t cutOffPhi = 1_r-cutOffPorosity; // solid volume fraction

      uint_t endEvalIdx = 0;
      uint_t numLayers = porosityPerLayer_.size();
      for(uint_t i = numLayers-1; i > 0; --i)
      {
         if(porosityPerLayer_[i] <= cutOffPhi && porosityPerLayer_[i-1] >= cutOffPhi)
         {
            endEvalIdx = i;
            break;
         }
      }
      if(endEvalIdx > 0) return 1_r - std::accumulate(porosityPerLayer_.begin(), std::next(porosityPerLayer_.begin(), static_cast<long>(endEvalIdx)), 0_r) / real_c(endEvalIdx);
      else return 1_r;

   }

private:

   real_t calculateSphericalSegmentVolume(real_t lowerLayerHeight, real_t upperLayerHeight, real_t radius) const
   {
      real_t r1 = std::max(std::min(lowerLayerHeight,radius),-radius);
      real_t r2 = std::max(std::min(upperLayerHeight,radius),-radius);
      real_t a1 = std::sqrt(std::max(radius*radius - r1*r1, real_t(0)));
      real_t a2 = std::sqrt(std::max(radius*radius - r2*r2, real_t(0)));
      real_t h = r2-r1;
      // wikipedia Kugelschicht
      return math::pi * h / real_t(6) * (real_t(3)*a1*a1 + real_t(3)*a2*a2 + h*h);
   }

   real_t layerHeight_;
   real_t layerVolume_;
   std::vector<real_t> porosityPerLayer_;
   std::vector<real_t> radiusPerLayer_;
};


class ContactInfoPerHorizontalLayerEvaluator
{
public:
   ContactInfoPerHorizontalLayerEvaluator(real_t layerHeight, const AABB & simulationDomain)
         : layerHeight_(layerHeight)
   {
      contactsPerLayer_ = std::vector<uint_t>(uint_c(simulationDomain.zSize() / layerHeight), 0);
      avgPenetrationPerLayer_ = std::vector<real_t>(uint_c(simulationDomain.zSize() / layerHeight), 0_r);
      maxPenetrationPerLayer_ = std::vector<real_t>(uint_c(simulationDomain.zSize() / layerHeight), 0_r);
   }

   template<typename ContactAccessor_T>
   void operator()(size_t c_idx, ContactAccessor_T& cac)
   {
      real_t zPos = cac.getPosition(c_idx)[2];
      real_t penetrationDepth = -cac.getDistance(c_idx);
      uint_t heightIdx= uint_c(std::max(0_r,zPos / layerHeight_));
      contactsPerLayer_[heightIdx]++;
      avgPenetrationPerLayer_[heightIdx] += penetrationDepth;
      maxPenetrationPerLayer_[heightIdx] = std::max(maxPenetrationPerLayer_[heightIdx], penetrationDepth);
   }

   void evaluate()
   {
      walberla::mpi::reduceInplace(contactsPerLayer_, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(avgPenetrationPerLayer_, walberla::mpi::SUM);
      walberla::mpi::reduceInplace(maxPenetrationPerLayer_, walberla::mpi::MAX);
      WALBERLA_ROOT_SECTION()
      {
         for(uint_t i = 0; i < contactsPerLayer_.size(); ++i)
         {
            if(contactsPerLayer_[i] > 0)
            {
               avgPenetrationPerLayer_[i] /= real_c(contactsPerLayer_[i]);
            }
         }
      }
   }

   void clear()
   {
      for(uint_t i = 0; i < contactsPerLayer_.size(); ++i)
      {
         contactsPerLayer_[i] = 0;
         avgPenetrationPerLayer_[i] = 0_r;
         maxPenetrationPerLayer_[i] = 0_r;
      }
   }

   void printToFile(std::string fileName)
   {
      WALBERLA_ROOT_SECTION()
      {
         std::ofstream file;
         file.open( fileName.c_str() );
         file << "#\t z\t numContacts\t avgPenetration\t maxPenetration\n";

         for(uint_t i = 0; i < contactsPerLayer_.size(); ++i)
         {
            file << layerHeight_ * (real_t(0.5)+real_c(i)) << " " << contactsPerLayer_[i]
                 << " " << avgPenetrationPerLayer_[i] << " " << maxPenetrationPerLayer_[i] << "\n";
         }
         file.close();
      }
   }

private:

   real_t layerHeight_;
   std::vector<uint_t> contactsPerLayer_;
   std::vector<real_t> avgPenetrationPerLayer_;
   std::vector<real_t> maxPenetrationPerLayer_;
};


class LoggingWriter
{
public:
   explicit LoggingWriter(std::string fileName) : fileName_(fileName){
      WALBERLA_ROOT_SECTION() {
         std::ofstream file;
         file.open(fileName_.c_str());
         file << "# t numParticles maxVel avgVel maxHeight massHeight kinEnergy numContacts maxPenetration avgPenetration porosity MCN\n";
         file.close();
      }
   }

   void operator()(real_t t, const ParticleInfo & particleInfo, const ContactInfo & contactInfo, real_t porosity)
   {
      WALBERLA_ROOT_SECTION() {
         std::ofstream file;
         file.open(fileName_.c_str(), std::ios_base::app);
         file << t << " " << particleInfo.numParticles << " " << particleInfo.maximumVelocity << " "
              << particleInfo.averageVelocity << " " << particleInfo.maximumHeight << " " << particleInfo.heightOfMass << " "
              << particleInfo.kinEnergy << " " << contactInfo.numContacts << " " << contactInfo.maximumPenetrationDepth << " "
              << contactInfo.averagePenetrationDepth << " " << porosity << " " << particleInfo.meanCoordinationNumber << "\n";
         file.close();
      }
   }

private:
   std::string fileName_;
};


std::string assembleParticleInformation(data::ParticleStorage& ps, SizeEvaluator sizeEvaluator, int precision)
{
   std::ostringstream ossData;

   for (auto pIt : ps)
   {

      using namespace data::particle_flags;
      if( isSet(pIt->getFlags(), GLOBAL) || isSet(pIt->getFlags(), GHOST)) continue;

      auto position = pIt->getPosition();
      auto sizeInfo = sizeEvaluator.get(*pIt->getBaseShape());
      auto rotM = pIt->getRotation().getMatrix();

      ossData << std::setprecision( precision );
      ossData << position[0] << " " << position[1] << " " << position[2] << " "
              << sizeInfo.size << " " << sizeInfo.volume << " "
              << sizeInfo.shapeSemiAxes[0] << " " << sizeInfo.shapeSemiAxes[1] << " " << sizeInfo.shapeSemiAxes[2] << " "
              << pIt->getNumContacts() << " "
              << rotM[0] << " " << rotM[1] << " "<< rotM[2] << " "<< rotM[3] << " "<< rotM[4] << " "<< rotM[5] << " "<< rotM[6] << " "<< rotM[7] << " "<< rotM[8]
              << "\n";
   }

   return ossData.str();
}



} // namespace msa_pd
} // namespace walberla

