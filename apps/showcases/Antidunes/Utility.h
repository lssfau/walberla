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
//! \file   Utility.h
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/shape/Sphere.h"

#include <algorithm>
#include <core/mpi/Broadcast.h>
#include <core/mpi/MPITextFile.h>
#include <core/mpi/Reduce.h>
#include <functional>
#include <iterator>

namespace walberla
{
namespace antidunes
{

struct SphereSelector
{
   template< typename ParticleAccessor_T >
   bool inline operator()(const size_t particleIdx, const ParticleAccessor_T& ac) const
   {
      static_assert(std::is_base_of< mesa_pd::data::IAccessor, ParticleAccessor_T >::value,
                    "Provide a valid accessor as template");
      return ac.getBaseShape(particleIdx)->getShapeType() == mesa_pd::data::Sphere::SHAPE_TYPE;
   }
};

void renameFile(const std::string& oldName, const std::string& newName)
{
   int result = std::rename(oldName.c_str(), newName.c_str());
   if (result != 0)
      WALBERLA_LOG_WARNING_ON_ROOT("Could not rename file " << oldName << " to " << newName << " with error code "
                                                            << result);
}

void write2DVectorToFile(const std::vector< real_t >& vec, uint_t len1, uint_t len2, std::string filename)
{
   std::ofstream file;
   file.open(filename.c_str());
   file.precision(5);

   file << "# " << len1 << " " << len2 << "\n";

   for (uint_t j = uint_t(0); j < len2; ++j)
   {
      for (uint_t i = uint_t(0); i < len1; ++i)
      {
         file << vec[i + j * len1] << "\n";
      }
   }
   file.close();
}

template< typename ParticleAccessor_T >
class BedloadTransportEvaluator
{
 public:
   BedloadTransportEvaluator(const shared_ptr< ParticleAccessor_T >& ac, real_t normalizationFactor,
                             uint_t numParticles)
      : ac_(ac), normalizationFactor_(normalizationFactor), numParticles_(numParticles)
   {}

   void operator()()
   {
      real_t transportRate(real_t(0));
      real_t velocity(real_t(0));

      for (uint_t i = uint_t(0); i < ac_->size(); ++i)
      {
         if (!isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GHOST) &&
             !isSet(ac_->getFlags(i), mesa_pd::data::particle_flags::GLOBAL))
         {
            auto velX = ac_->getLinearVelocity(i)[0];
            transportRate += velX * ac_->getVolume(i);
            velocity += velX;
         }
      }

      // only reduce to root
      WALBERLA_MPI_SECTION()
      {
         mpi::reduceInplace(transportRate, mpi::SUM);
         mpi::reduceInplace(velocity, mpi::SUM);
      }

      avgTransportRate_ = transportRate * normalizationFactor_;
      averageVelocity_  = velocity / real_c(numParticles_);
   }

   // sum_i V_p,i * u_p,i / (L*W)
   real_t getTransportRate() const { return avgTransportRate_; }

   real_t getAverageVelocity() const { return averageVelocity_; }

 private:
   shared_ptr< ParticleAccessor_T > ac_;
   real_t normalizationFactor_;
   uint_t numParticles_;
   real_t averageVelocity_;
   real_t avgTransportRate_;
};

template< typename ParticleAccessor_T >
Vector3< real_t > getTotalHydrodynamicForceOnParticles(const shared_ptr< ParticleAccessor_T >& ac)
{
   Vector3< real_t > totalHydrodynamicForce(0_r);
   for (uint_t i = uint_t(0); i < ac->size(); ++i)
   {
      if (!isSet(ac->getFlags(i), mesa_pd::data::particle_flags::GHOST) &&
          !isSet(ac->getFlags(i), mesa_pd::data::particle_flags::GLOBAL))
      {
         totalHydrodynamicForce += ac->getHydrodynamicForce(i);
      }
   }
   // only reduce to root
   mpi::reduceInplace(totalHydrodynamicForce, mpi::SUM);

   return totalHydrodynamicForce;
}

// evaluates slices of solid volume fraction and fill level
template< typename PdfField_T, typename AntidunesBoundaryHandling_T, typename FlagField_T, typename ScalarField_T >
class AverageDataSliceEvaluator
{
 public:
   AverageDataSliceEvaluator(const shared_ptr< StructuredBlockStorage >& blocks, const ConstBlockDataID& flagFieldID,
                             const ConstBlockDataID& fillFieldID, const ConstBlockDataID& pdfFieldID)
      : blocks_(blocks), flagFieldID_(flagFieldID), fillFieldID_(fillFieldID), pdfFieldID_(pdfFieldID)
   {
      xlen_ = blocks_->getNumberOfXCells();
      ylen_ = blocks_->getNumberOfYCells();
      zlen_ = blocks_->getNumberOfZCells();

      x_z_SolidVolumeFraction_ = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_FillLevel_           = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_VelocityX_           = std::vector< real_t >(xlen_ * zlen_, real_t(0));
      x_z_FluidCellCount_      = std::vector< uint_t >(xlen_ * zlen_, 0);
      maxFluidZPos_            = uint_t(0);
   }

   void operator()()
   {
      // erase data
      std::fill(x_z_SolidVolumeFraction_.begin(), x_z_SolidVolumeFraction_.end(), real_t(0));
      std::fill(x_z_FillLevel_.begin(), x_z_FillLevel_.end(), real_t(0));
      std::fill(x_z_VelocityX_.begin(), x_z_VelocityX_.end(), real_t(0));
      std::fill(x_z_FluidCellCount_.begin(), x_z_FluidCellCount_.end(), 0);

      // fill contributions
      for (auto block = blocks_->begin(); block != blocks_->end(); ++block)
      {
         const PdfField_T* const pdfField     = block->getData< const PdfField_T >(pdfFieldID_);
         const FlagField_T* const flagField   = block->getData< const FlagField_T >(flagFieldID_);
         const ScalarField_T* const fillField = block->getData< const ScalarField_T >(fillFieldID_);

         const auto solidMO     = flagField->getFlag(AntidunesBoundaryHandling_T::movingObstacleFlagID);
         const auto solidNoSlip = flagField->getFlag(AntidunesBoundaryHandling_T::noSlipFlagID);

         CellInterval xyz = flagField->xyzSize();
         Cell globalCell;

         maxFluidZPos_ = uint_t(0);

         // iterate all (inner) cells in the field
         for (auto cell = xyz.begin(); cell != xyz.end(); ++cell)
         {
            blocks_->transformBlockLocalToGlobalCell(globalCell, *block, *cell);
            auto entryIdx = uint_c(globalCell.x()) + uint_c(globalCell.z()) * xlen_;
            if (flagField->isFlagSet(*cell, solidMO) || flagField->isFlagSet(*cell, solidNoSlip))
            {
               x_z_SolidVolumeFraction_[entryIdx] += real_t(1);
               x_z_FillLevel_[entryIdx] += real_t(1);
            }
            else
            {
               auto fillLevel = fillField->get(*cell);
               x_z_FillLevel_[entryIdx] += fillLevel;
               if (fillLevel > 0_r)
               {
                  x_z_VelocityX_[entryIdx] += pdfField->getVelocity(*cell)[0];
                  ++x_z_FluidCellCount_[entryIdx];

                  maxFluidZPos_ = std::max(uint_t(globalCell.z()), maxFluidZPos_);
               }
            }
         }
      }

      // reduce this information to the root process
      mpi::reduceInplace(x_z_SolidVolumeFraction_, mpi::SUM);
      mpi::reduceInplace(x_z_FillLevel_, mpi::SUM);
      mpi::reduceInplace(x_z_VelocityX_, mpi::SUM);
      mpi::reduceInplace(x_z_FluidCellCount_, mpi::SUM);
      mpi::reduceInplace(maxFluidZPos_, mpi::MAX);

      // normalize
      for (uint_t i = 0; i < x_z_VelocityX_.size(); ++i)
      {
         if (x_z_FluidCellCount_[i] > 0) x_z_VelocityX_[i] /= real_c(x_z_FluidCellCount_[i]);
      }
      real_t invNumYCells = 1_r / real_c(ylen_);
      std::for_each(x_z_SolidVolumeFraction_.begin(), x_z_SolidVolumeFraction_.end(),
                    [invNumYCells](real_t& n) { n *= invNumYCells; });
      std::for_each(x_z_FillLevel_.begin(), x_z_FillLevel_.end(), [invNumYCells](real_t& n) { n *= invNumYCells; });

      // note: only root process has the complete information!
   }

   std::vector< real_t >& getSolidVolumeFractionVector() { return x_z_SolidVolumeFraction_; }
   std::vector< real_t >& getFillLevelVector() { return x_z_FillLevel_; }
   std::vector< real_t >& getVelocityXVector() { return x_z_VelocityX_; }
   std::vector< uint_t >& getFluidCellCountVector() { return x_z_FluidCellCount_; }
   uint_t getMaxFluidZPos() { return maxFluidZPos_; }

   uint_t getXLen() const { return xlen_; }
   uint_t getYLen() const { return ylen_; }
   uint_t getZLen() const { return zlen_; }

 private:
   shared_ptr< StructuredBlockStorage > blocks_;
   const ConstBlockDataID flagFieldID_;
   const ConstBlockDataID fillFieldID_;
   const ConstBlockDataID pdfFieldID_;

   uint_t xlen_;
   uint_t ylen_;
   uint_t zlen_;
   uint_t maxFluidZPos_;
   std::vector< real_t > x_z_SolidVolumeFraction_;
   std::vector< real_t > x_z_FillLevel_;
   std::vector< real_t > x_z_VelocityX_;
   std::vector< uint_t > x_z_FluidCellCount_;
};

void writeSphereInformationToFile(const std::string& filename, walberla::mesa_pd::data::ParticleStorage& ps,
                                  Vector3< real_t >& domainSize, int precision = 12)
{
   std::ostringstream ossData;
   ossData << std::setprecision(precision);

   WALBERLA_ROOT_SECTION() { ossData << domainSize[0] << " " << domainSize[1] << " " << domainSize[2] << "\n"; }

   for (auto pIt : ps)
   {
      using namespace walberla::mesa_pd::data;
      if (pIt->getBaseShape()->getShapeType() != Sphere::SHAPE_TYPE) continue;
      using namespace walberla::mesa_pd::data::particle_flags;
      if (isSet(pIt->getFlags(), GHOST)) continue;
      auto sp = static_cast< Sphere* >(pIt->getBaseShape().get());

      auto position = pIt->getPosition();

      ossData << pIt->getUid() << " " << position[0] << " " << position[1] << " " << position[2] << " "
              << sp->getRadius() << '\n';
   }

   walberla::mpi::writeMPITextFile(filename, ossData.str());
}

void getAvgDiameterScalingFactor(const std::string& filename, const Vector3< uint_t >& domainSize,
                                 const uint_t bedCopiesInX, const uint_t bedCopiesInY, real_t& avgDiameter,
                                 real_t& scalingFactor)
{
   using namespace walberla;

   std::string textFile;

   WALBERLA_ROOT_SECTION()
   {
      std::ifstream t(filename.c_str());
      if (!t) { WALBERLA_ABORT("Invalid input file " << filename << "\n"); }
      std::stringstream buffer;
      buffer << t.rdbuf();
      textFile = buffer.str();

      std::istringstream fileIss(textFile);
      std::string line;

      // first line contains generation domain sizes
      std::getline(fileIss, line);
      Vector3< real_t > generationDomainSize_SI(0_r);
      std::istringstream firstLine(line);
      firstLine >> generationDomainSize_SI[0] >> generationDomainSize_SI[1] >> generationDomainSize_SI[2];

      real_t diameter_SI  = 0.0;
      uint_t numParticles = 0;
      while (std::getline(fileIss, line))
      {
         std::istringstream iss(line);

         mesa_pd::data::ParticleStorage::uid_type uID;
         mesa_pd::data::ParticleStorage::position_type pos(0_r);
         walberla::real_t radius = 0;
         iss >> uID >> pos[0] >> pos[1] >> pos[2] >> radius;
         WALBERLA_CHECK_GREATER(radius, 0_r, "Invalid radius of " << radius << " found in input file!")

         diameter_SI += 2_r * radius;

         numParticles++;
      }
      diameter_SI /= real_t(numParticles);

      scalingFactor = real_c(domainSize[0]) / (generationDomainSize_SI[0] * real_c(bedCopiesInX));
      avgDiameter   = diameter_SI * scalingFactor;

      WALBERLA_CHECK_EQUAL(uint_c(scalingFactor * generationDomainSize_SI[1] * real_c(bedCopiesInY)), domainSize[1],
                           "Error: Generated bed with copies and simulation domain do not match!")
   }

   walberla::mpi::broadcastObject(scalingFactor);
   walberla::mpi::broadcastObject(avgDiameter);
}

void initSpheresFromFile(const std::string& filename, walberla::mesa_pd::data::ParticleStorage& ps,
                         const walberla::mesa_pd::domain::IDomain& domain, walberla::real_t density,
                         const Vector3< uint_t >& domainSize,
                         const std::function< bool(walberla::Vector3< real_t >) >& particleCreateFunction,
                         math::AABB simulationDomain, uint_t bedCopiesInX, uint_t bedCopiesInY, uint_t& numParticles,
                         real_t& maxParticleHeight, const real_t& scalingFactor)
{
   using namespace walberla;
   using namespace walberla::mesa_pd;
   using namespace walberla::mesa_pd::data;

   auto rank = walberla::mpi::MPIManager::instance()->rank();

   std::string textFile;

   WALBERLA_ROOT_SECTION()
   {
      std::ifstream t(filename.c_str());
      if (!t) { WALBERLA_ABORT("Invalid input file " << filename << "\n"); }
      std::stringstream buffer;
      buffer << t.rdbuf();
      textFile = buffer.str();
   }

   walberla::mpi::broadcastObject(textFile);

   std::istringstream fileIss(textFile);
   std::string line;

   // first line contains generation domain sizes
   std::getline(fileIss, line);
   Vector3< real_t > generationDomainSize_SI(0_r);
   std::istringstream firstLine(line);
   firstLine >> generationDomainSize_SI[0] >> generationDomainSize_SI[1] >> generationDomainSize_SI[2];
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(generationDomainSize_SI)

   WALBERLA_CHECK_EQUAL(uint_c(scalingFactor * generationDomainSize_SI[0] * real_c(bedCopiesInX)), domainSize[0],
                        "Error: Generated bed with copies and simulation domain do not match in x!")
   WALBERLA_CHECK_EQUAL(uint_c(scalingFactor * generationDomainSize_SI[1] * real_c(bedCopiesInY)), domainSize[1],
                        "Error: Generated bed with copies and simulation domain do not match in y!")

   numParticles      = 0;
   maxParticleHeight = 0_r;

   while (std::getline(fileIss, line))
   {
      std::istringstream iss(line);

      data::ParticleStorage::uid_type uID;
      data::ParticleStorage::position_type pos;
      walberla::real_t radius;
      iss >> uID >> pos[0] >> pos[1] >> pos[2] >> radius;
      radius *= scalingFactor;

      for (uint_t copyInYDir = 0; copyInYDir < bedCopiesInY; ++copyInYDir)
      {
         for (uint_t copyInXDir = 0; copyInXDir < bedCopiesInX; ++copyInXDir)
         {
            auto particlePos = pos;

            particlePos[0] += real_c(copyInXDir) * generationDomainSize_SI[0];
            particlePos[1] += real_c(copyInYDir) * generationDomainSize_SI[1];

            particlePos *= scalingFactor;

            maxParticleHeight = std::max(maxParticleHeight, particlePos[2] + radius);

            if (!particleCreateFunction(particlePos)) continue;

            WALBERLA_CHECK(simulationDomain.contains(particlePos),
                           "Particle read from file is not contained in simulation domain");

            if (!domain.isContainedInProcessSubdomain(uint_c(rank), particlePos)) continue;

            auto pIt = ps.create();
            pIt->setPosition(particlePos);
            pIt->getBaseShapeRef() = std::make_shared< data::Sphere >(radius);
            pIt->getBaseShapeRef()->updateMassAndInertia(density);
            pIt->setInteractionRadius(radius);
            pIt->setOwner(rank);
            pIt->setType(0);

            numParticles++;
         }
      }

      WALBERLA_CHECK_EQUAL(iss.tellg(), -1);
   }
   walberla::mpi::allReduceInplace(maxParticleHeight, walberla::mpi::MAX);
   walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
}

void getAverageVelocity(const mesa_pd::data::ParticleAccessorWithBaseShape& ac, real_t& averageVelocity,
                        real_t& maxVelocity, uint_t& numParticles, real_t& maxHeight)
{
   averageVelocity = real_t(0);
   maxVelocity     = real_t(0);
   numParticles    = uint_t(0);
   maxHeight       = real_t(0);
   for (uint_t i = 0; i < ac.size(); ++i)
   {
      if (isSet(ac.getFlags(i), walberla::mesa_pd::data::particle_flags::GHOST)) continue;
      if (isSet(ac.getFlags(i), walberla::mesa_pd::data::particle_flags::GLOBAL)) continue;

      ++numParticles;
      real_t velMagnitude = ac.getLinearVelocity(i).length();
      averageVelocity += velMagnitude;
      maxVelocity = std::max(maxVelocity, velMagnitude);
      maxHeight   = std::max(maxHeight, ac.getPosition(i)[2]);
   }

   walberla::mpi::allReduceInplace(numParticles, walberla::mpi::SUM);
   walberla::mpi::allReduceInplace(averageVelocity, walberla::mpi::SUM);
   walberla::mpi::allReduceInplace(maxVelocity, walberla::mpi::MAX);
   walberla::mpi::allReduceInplace(maxHeight, walberla::mpi::MAX);

   averageVelocity /= real_t(numParticles);
}

auto createPlane(mesa_pd::data::ParticleStorage& ps, const mesa_pd::Vec3& pos, const mesa_pd::Vec3& normal)
{
   auto p0 = ps.create(true);
   p0->setPosition(pos);
   p0->setBaseShape(std::make_shared< mesa_pd::data::HalfSpace >(normal));
   p0->getBaseShapeRef()->updateMassAndInertia(real_t(1));
   p0->setOwner(walberla::mpi::MPIManager::instance()->rank());
   p0->setType(0);
   p0->setInteractionRadius(std::numeric_limits< real_t >::infinity());
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::GLOBAL);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::INFINITE);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::FIXED);
   mesa_pd::data::particle_flags::set(p0->getFlagsRef(), mesa_pd::data::particle_flags::NON_COMMUNICATING);
   return p0;
}

} // namespace antidunes
} // namespace walberla
