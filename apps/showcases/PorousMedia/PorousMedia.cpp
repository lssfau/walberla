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
//! \file PorousMedia.cpp
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//
//======================================================================================================================

//======================================================================================================================
// This showcase simulates flow through a porous medium that is composed of spherical particles. It must be given a
// parameter file as command-line argument, such as "PorousMedia.prm".
//
// The porous medium is located in a channel with periodic side walls and is read from a file that can be generated with
// "PackedBedCreation.cpp". The flow is driven by a velocity boundary at the top of the channel and a pressure boundary
// (with rho=1) at its bottom. Each spherical particle is mapped to the fluid domain with no-slip boundary conditions
// and the whole fluid domain is initialized with the velocity that is prescribed at the top wall boundary condition.
// The computational grid is statically refined with 3 levels. The lattice Boltzmann kernel features a cumulant LBM
// model and is generated using lbmpy. The simulation is considered stationary, i.e., converged when the L2 norm of the
// relative difference between the current and former time-averaged quantities is less than the specified threshold.
// The showcase also includes various additional features for post-processing, such as time-averaging of density and
// velocity, gradient computation, pressure drop computation, and monitoring of instantaneous quantities at different
// coordinates and slices through the domain.
//
// IMPORTANT REMARK: The setup is intended to be run on large compute clusters with particles being resolved by 80 or
// more cells per diameter, setups with smaller resolution are physically inaccurate and might even become unstable.
//======================================================================================================================

#include "blockforest/BlockForest.h"
#include "blockforest/SetupBlockForest.h"
#include "blockforest/StructuredBlockForest.h"
#include "blockforest/loadbalancing/StaticCurve.h"

#include "boundary/BoundaryHandling.h"

#include "core/Environment.h"
#include "core/SharedFunctor.h"
#include "core/StringUtility.h"
#include "core/cell/CellInterval.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"
#include "core/mpi/Broadcast.h"
#include "core/stringToNum.h"
#include "core/timing/TimingPool.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/AddToStorage.h"
#include "field/FlagField.h"
#include "field/communication/PackInfo.h"
#include "field/vtk/FlagFieldCellFilter.h"
#include "field/vtk/VTKWriter.h"

#include "lbm/PerformanceLogger.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/SimplePressure.h"
#include "lbm/boundary/UBB.h"
#include "lbm/field/AddToStorage.h"
#include "lbm/field/PdfField.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/refinement/all.h"
#include "lbm/vtk/Density.h"
#include "lbm/vtk/Velocity.h"

#include "pe/basic.h"
#include "pe/rigidbody/SetBodyTypeIDs.h"
#include "pe/rigidbody/Sphere.h"
#include "pe/rigidbody/StorageDataHandling.h"
#include "pe/vtk/SphereVtkOutput.h"

#include "pe_coupling/momentum_exchange_method/BodyMapping.h"
#include "pe_coupling/momentum_exchange_method/boundary/SimpleBB.h"

#include "stencil/D3Q6.h"

#include "timeloop/SweepTimeloop.h"

#include "vtk/ChainedFilter.h"
#include "vtk/VTKOutput.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "LbCodeGen_LatticeModel.h"

namespace walberla
{
//////////////////////
// GLOBAL VARIABLES //
//////////////////////

// 4 ghost layers are required when using grid refinement
const uint_t FieldGhostLayers = 4;

///////////
// USING //
///////////

using LatticeModel_T         = lbm::LbCodeGen_LatticeModel;
using Stencil_T              = LatticeModel_T::Stencil;
using CommunicationStencil_T = LatticeModel_T::CommunicationStencil;
using PdfField_T             = lbm::PdfField< LatticeModel_T >;

using flag_t      = uint8_t;
using FlagField_T = FlagField< flag_t >;
using NoSlip_T    = lbm::NoSlip< LatticeModel_T, flag_t >;
using UBB_T       = lbm::UBB< LatticeModel_T, flag_t >;
using Pressure_T  = lbm::SimplePressure< LatticeModel_T, flag_t >;
using MO_SBB_T    = pe_coupling::SimpleBB< LatticeModel_T, FlagField_T >;

using BoundaryHandling_T = BoundaryHandling< FlagField_T, Stencil_T, NoSlip_T, UBB_T, Pressure_T, MO_SBB_T >;

using BodyField_T   = GhostLayerField< pe::BodyID, 1 >;
using BodyTypeTuple = std::tuple< pe::Sphere >; // must contain all PE body types used in this simulation
using VectorField_T = GhostLayerField< real_t, LatticeModel_T::Stencil::D >;

// AvgField_T stores averaged quantities:
// 0: density, 1: velocity_x, 2: velocity_y, 3: velocity_z, 4: velocity_magn
// 5: density_old, 6: velocity_x_old, 7: velocity_y_old, 8: velocity_z_old, 9: velocity_magn_old
using AvgField_T = GhostLayerField< real_t, 10 >;

// number of specific values per cell for defining the field sizes below
constexpr uint_t densValues = 1; // number of values stored for density
constexpr uint_t velValues  = 3; // number of values (components) stored for the velocity
constexpr uint_t gradValues = 9; // number of values (components) stored for the velocity gradient
constexpr uint_t vortValues = 3; // number of values (components) stored for the vorticity
constexpr uint_t instValues =
   densValues + velValues + gradValues + vortValues; // total number of instantaneous values to be stored

// number of snapshots of the whole domain to be cached in InstField_T, i.e., variable controls how often
// instantaneous vtk output is written to disk; this allows to write the output in junks to reduce file-IO
constexpr uint_t instFieldEntries = 1;

// InstField_T stores instantaneous quantities (as float since sufficient for vtk output):
// 0: density,
// 1: velocity_x (U), 2: velocity_y (V), 3: velocity_z (W),
// [gradient] 4: dU/dx, 5: dU/dy, 6: dU/dz, 7: dV/dx, 8: dV/dy, 9: dV/dz, 10: dW/dx, 11: dW/dy, 12: dW/dz
// [vorticity] 13: dW/dy-dV/dz, 14: dU/dz-dW/dx, 15: dV/dx-dU/dy
using InstField_T = GhostLayerField< float, instValues * instFieldEntries >;

// GradientField_T stores gradients (of time-averaged quantities):
// 0: dU/dx, 1: dU/dy, 2: dU/dz, 3: dV/dx, 4: dV/dy, 5: dV/dz, 6: dW/dx, 7: dW/dy, 8: dW/dz
using GradientField_T = GhostLayerField< real_t, gradValues >;

///////////
// FLAGS //
///////////

const FlagUID Fluid_Flag("fluid");
const FlagUID UBB_Flag("velocity bounce back");
const FlagUID NoSlip_Flag("no slip");
const FlagUID SolidNearFluid_Flag("solid near fluid");
const FlagUID Pressure_Flag("pressure boundary");
const FlagUID MO_SBB_Flag("velocity bounce back for particles");

///////////////
// UTILITIES //
///////////////

// convert a number to string with the specified number of digits (precision)
template< typename T >
std::string toStringWithPrecision(const T a_value, const int precision = 12);

// write the content of a std::vector to a line in a file with "timestep" being the line index
template< typename T >
void writeVector(const std::vector< T >& data, const uint_t& timestep, const std::string& filename);

// read the geometry of a packed bed from file to vector "particleCoordinates" (file contains locations of particles
// relative to the particle diameter)
void readPackedBedFromFile(const std::string& filename, std::vector< Vector3< real_t > >& particleCoordinates);

// compute the void fraction in a specified axis-aligned bounding box (AABB)
template< typename BoundaryHandling_T >
real_t getVoidFraction(const shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& boundaryHandlingID,
                       const AABB& computeRegionAABB);

// mark solid cells that have fluid cells in their neighborhood (used for writing first layer of solid cells to vtk
// output => zero-velocity at the solid boundary is also displayed in vtk output)
template< typename BoundaryHandling_T >
void setSolidNearFluidFlag(const std::shared_ptr< StructuredBlockForest >& blocks,
                           const BlockDataID& boundaryHandlingID, FlagUID flag);

// data handling for loading a field of type AvgField_T from file
template< typename AvgField_T >
class AvgFieldHandling : public field::BlockDataHandling< AvgField_T >
{
 public:
   AvgFieldHandling(const weak_ptr< StructuredBlockStorage >& blocks) : blocks_(blocks) {}

 protected:
   AvgField_T* allocate(IBlock* const block) override { return allocateDispatch(block); }

   AvgField_T* reallocate(IBlock* const block) override { return allocateDispatch(block); }

 private:
   weak_ptr< StructuredBlockStorage > blocks_;

   AvgField_T* allocateDispatch(IBlock* const block)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);

      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR(blocks);

      return new AvgField_T(blocks->getNumberOfXCells(*block), blocks->getNumberOfYCells(*block),
                            blocks->getNumberOfZCells(*block), uint_c(0), real_c(0), field::fzyx);
   }
}; // class AvgFieldHandling

// compute the L2 norm of the relative difference between the current and the previous time-averaged quantities
template< typename BoundaryHandling_T, typename AvgField_T >
class AvgDiffEvaluator
{
 public:
   AvgDiffEvaluator(const std::shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& boundaryHandlingID,
                    const ConstBlockDataID& avgFieldID, uint_t interval,
                    const std::shared_ptr< real_t >& avgDensityDiffL2,
                    const std::shared_ptr< Vector3< real_t > >& avgVelocityDiffL2)
      : blocks_(blocks), boundaryHandlingID_(boundaryHandlingID), avgFieldID_(avgFieldID), interval_(interval),
        avgDensityDiffL2_(avgDensityDiffL2), avgVelocityDiffL2_(avgVelocityDiffL2), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % interval_ != uint_c(0)) { return; }

      // compute the L2 norm of the relative difference of density and velocity components in all fluid cells
      computeDifferenceL2();
   }

 private:
   const std::shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID avgFieldID_;
   const uint_t interval_;

   const std::shared_ptr< real_t > avgDensityDiffL2_;
   const std::shared_ptr< Vector3< real_t > > avgVelocityDiffL2_;

   uint_t executionCounter_; // number of times operator() has been called

   // L2 norm as defined in Krueger et al., The Lattice Boltzmann Method, p. 138
   void computeDifferenceL2()
   {
      uint_t numFluidCells = uint_c(0);

      // current time-averaged density and velocity, space-averaged over all fluid cells
      real_t avgDensity                    = real_c(0);
      real_t avgDensitySquared             = real_c(0);
      Vector3< real_t > avgVelocity        = Vector3< real_t >(real_c(0));
      Vector3< real_t > avgVelocitySquared = Vector3< real_t >(real_c(0));

      real_t avgDensityDiffL2             = real_c(0);
      Vector3< real_t > avgVelocityDiffL2 = Vector3< real_t >(real_c(0));

      for (auto blockIterator = blocks_->begin(); blockIterator != blocks_->end(); ++blockIterator)
      {
         BoundaryHandling_T* boundaryHandling =
            blockIterator->template getData< BoundaryHandling_T >(boundaryHandlingID_);
         AvgField_T* velDensAvgField = blockIterator->template getData< AvgField_T >(avgFieldID_);

         // compute space-average over all fluid cells of the time-averaged quantities
         // clang-format off
         WALBERLA_FOR_ALL_CELLS_XYZ(velDensAvgField,
            if (boundaryHandling->isDomain(x, y, z)) {
               ++numFluidCells;

               avgDensity += velDensAvgField->get(x, y, z, uint_c(0));
               avgVelocity[0] += velDensAvgField->get(x, y, z, uint_c(1));
               avgVelocity[1] += velDensAvgField->get(x, y, z, uint_c(2));
               avgVelocity[2] += velDensAvgField->get(x, y, z, uint_c(3));

               avgDensitySquared += velDensAvgField->get(x, y, z, uint_c(0)) * velDensAvgField->get(x, y, z, uint_c(0));
               avgVelocitySquared[0] +=
                  velDensAvgField->get(x, y, z, uint_c(1)) * velDensAvgField->get(x, y, z, uint_c(1));
               avgVelocitySquared[1] +=
                  velDensAvgField->get(x, y, z, uint_c(2)) * velDensAvgField->get(x, y, z, uint_c(2));
               avgVelocitySquared[2] +=
                  velDensAvgField->get(x, y, z, uint_c(3)) * velDensAvgField->get(x, y, z, uint_c(3));

               avgDensityDiffL2 += real_c(
                  std::pow(velDensAvgField->get(x, y, z, uint_c(0)) - velDensAvgField->get(x, y, z, uint_c(5)), 2));
               avgVelocityDiffL2[0] += real_c(
                  std::pow(velDensAvgField->get(x, y, z, uint_c(1)) - velDensAvgField->get(x, y, z, uint_c(6)), 2));
               avgVelocityDiffL2[1] += real_c(
                  std::pow(velDensAvgField->get(x, y, z, uint_c(2)) - velDensAvgField->get(x, y, z, uint_c(7)), 2));
               avgVelocityDiffL2[2] += real_c(
                  std::pow(velDensAvgField->get(x, y, z, uint_c(3)) - velDensAvgField->get(x, y, z, uint_c(8)), 2));
            }
         ) // WALBERLA_FOR_ALL_CELLS
         // clang-format on
      }

      // sum values among all processes
      mpi::allReduceInplace< uint_t >(numFluidCells, mpi::SUM);
      mpi::allReduceInplace< real_t >(avgDensity, mpi::SUM);
      mpi::allReduceInplace< real_t >(avgVelocity, mpi::SUM);

      mpi::allReduceInplace< real_t >(avgDensitySquared, mpi::SUM);
      mpi::allReduceInplace< real_t >(avgVelocitySquared, mpi::SUM);

      mpi::allReduceInplace< real_t >(avgDensityDiffL2, mpi::SUM);
      mpi::allReduceInplace< real_t >(avgVelocityDiffL2, mpi::SUM);

      // compute space-average
      avgDensity /= real_c(numFluidCells);
      avgVelocity[0] /= real_c(numFluidCells);
      avgVelocity[1] /= real_c(numFluidCells);
      avgVelocity[2] /= real_c(numFluidCells);

      // compute L2 norm
      avgDensityDiffL2     = real_c(std::pow(avgDensityDiffL2 / avgDensitySquared, 0.5));
      avgVelocityDiffL2[0] = real_c(std::pow(avgVelocityDiffL2[0] / avgVelocitySquared[0], 0.5));
      avgVelocityDiffL2[1] = real_c(std::pow(avgVelocityDiffL2[1] / avgVelocitySquared[1], 0.5));
      avgVelocityDiffL2[2] = real_c(std::pow(avgVelocityDiffL2[2] / avgVelocitySquared[2], 0.5));

      *avgDensityDiffL2_       = avgDensityDiffL2;
      (*avgVelocityDiffL2_)[0] = avgVelocityDiffL2[0];
      (*avgVelocityDiffL2_)[1] = avgVelocityDiffL2[1];
      (*avgVelocityDiffL2_)[2] = avgVelocityDiffL2[2];

      WALBERLA_LOG_DEVEL_ON_ROOT("numFluidCells= " << numFluidCells);
      WALBERLA_LOG_DEVEL_ON_ROOT("avgDensity= " << std::abs(avgDensity));
      WALBERLA_LOG_DEVEL_ON_ROOT("avgVelocity[0]= " << std::abs(avgVelocity[0]));
      WALBERLA_LOG_DEVEL_ON_ROOT("avgVelocity[1]= " << std::abs(avgVelocity[1]));
      WALBERLA_LOG_DEVEL_ON_ROOT("avgVelocity[2]= " << std::abs(avgVelocity[2]));

      if (math::isnan(avgDensity) || math::isinf(avgDensity) || math::isnan(avgVelocity[0]) ||
          math::isinf(avgVelocity[0]) || math::isnan(avgVelocity[1]) || math::isinf(avgVelocity[1]) ||
          math::isnan(avgVelocity[2]) || math::isinf(avgVelocity[2]))
      {
         WALBERLA_ABORT("Simulation got unstable with nan or inf.")
      }
   }

}; // class AvgDiffEvaluator

// compute the time-average of the density and velocity-components, and store the result in avgFieldID;
// AvgField contains the averages of: density, velocity_x, velocity_y, velocity_z, velocity_magn, density_old,
// velocity_x_old, velocity_y_old, velocity_z_old, velocity_magn_old
template< typename PdfField_T, typename BoundaryHandling_T, typename AvgField_T >
class VelDensAverager
{
 public:
   VelDensAverager(const std::shared_ptr< StructuredBlockForest >& blocks, const ConstBlockDataID& pdfFieldID,
                   const BlockDataID& boundaryHandlingID, const ConstBlockDataID& avgFieldID, uint_t interval,
                   const uint_t compCounterLoaded)
      : blocks_(blocks), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID), avgFieldID_(avgFieldID),
        interval_(interval), compCounterLoaded_(compCounterLoaded), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % interval_ != uint_c(0)) { return; }

      // compute time-average of density and velocity in every fluid cell
      averageDensityVelocity();
   }

 private:
   const std::shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID avgFieldID_;
   const uint_t interval_;
   const uint_t compCounterLoaded_; // number of times this function has been called for a loaded velDensAvgField before

   uint_t executionCounter_; // number of times operator() has been called

   void averageDensityVelocity()
   {
      uint_t compCounter =
         executionCounter_ / interval_ + compCounterLoaded_; // number of times this function has been called

      real_t compCounterInv = real_c(1) / real_c(compCounter);

      for (auto blockIterator = blocks_->begin(); blockIterator != blocks_->end(); ++blockIterator)
      {
         PdfField_T* pdfField = blockIterator->template getData< PdfField_T >(pdfFieldID_);
         BoundaryHandling_T* boundaryHandling =
            blockIterator->template getData< BoundaryHandling_T >(boundaryHandlingID_);
         AvgField_T* velDensAvgField = blockIterator->template getData< AvgField_T >(avgFieldID_);

         // clang-format off
         WALBERLA_FOR_ALL_CELLS_XYZ(velDensAvgField,
            if (boundaryHandling->isDomain(x, y, z)) {
               // store the former average values
               velDensAvgField->get(x, y, z, uint_c(5)) = velDensAvgField->get(x, y, z, uint_c(0));
               velDensAvgField->get(x, y, z, uint_c(6)) = velDensAvgField->get(x, y, z, uint_c(1));
               velDensAvgField->get(x, y, z, uint_c(7)) = velDensAvgField->get(x, y, z, uint_c(2));
               velDensAvgField->get(x, y, z, uint_c(8)) = velDensAvgField->get(x, y, z, uint_c(3));
               velDensAvgField->get(x, y, z, uint_c(9)) = velDensAvgField->get(x, y, z, uint_c(4));

               // get current density and velocity from PDF field
               Vector3< real_t > velocity;
               real_t density = pdfField->getDensityAndVelocity(velocity, x, y, z);

               // calculate the current average using the accumulated values (obtained by multiplication with
               // "compCounter-1", i.e., by reverting the average computation)
               velDensAvgField->get(x, y, z, uint_c(0)) =
                  (velDensAvgField->get(x, y, z, uint_c(0)) * (real_c(compCounter) - real_c(1)) + density) *
                  compCounterInv;
               velDensAvgField->get(x, y, z, uint_c(1)) =
                  (velDensAvgField->get(x, y, z, uint_c(1)) * (real_c(compCounter) - real_c(1)) + velocity[0]) *
                  compCounterInv;
               velDensAvgField->get(x, y, z, uint_c(2)) =
                  (velDensAvgField->get(x, y, z, uint_c(2)) * (real_c(compCounter) - real_c(1)) + velocity[1]) *
                  compCounterInv;
               velDensAvgField->get(x, y, z, uint_c(3)) =
                  (velDensAvgField->get(x, y, z, uint_c(3)) * (real_c(compCounter) - real_c(1)) + velocity[2]) *
                  compCounterInv;
               velDensAvgField->get(x, y, z, uint_c(4)) =
                  (velDensAvgField->get(x, y, z, uint_c(4)) * (real_c(compCounter) - real_c(1)) + velocity.length()) *
                  compCounterInv;
            }
         ) // WALBERLA_FOR_ALL_CELLS
         // clang-format on
      }
   }
}; // class VelDensAverager

// compute the gradient of the velocity and store it in gradientFieldID
template< typename BoundaryHandling_T, typename AvgField_T, typename GradientField_T >
class GradientComputer
{
 public:
   GradientComputer(const std::shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& boundaryHandlingID,
                    const ConstBlockDataID& avgFieldID, const ConstBlockDataID& gradientFieldID)
      : blocks_(blocks), boundaryHandlingID_(boundaryHandlingID), avgFieldID_(avgFieldID),
        gradientFieldID_(gradientFieldID)
   {}

   void operator()()
   {
      // sync avgField field, as ghost layer values are needed in computeVelocityGradient()
      blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > avgFieldSync(blocks_);
      avgFieldSync.addPackInfo(make_shared< field::refinement::PackInfo< AvgField_T, stencil::D3Q27 > >(avgFieldID_));
      avgFieldSync();

      // compute the gradient of the velocity
      computeVelocityGradient();
   }

 private:
   const std::shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID avgFieldID_;
   const BlockDataID gradientFieldID_;

   void computeVelocityGradient()
   {
      const real_t dx        = real_t(1); // dx on finest level
      const uint_t numLevels = blocks_->getNumberOfLevels();

      for (auto blockIterator = blocks_->begin(); blockIterator != blocks_->end(); ++blockIterator)
      {
         GradientField_T* gradientField = blockIterator->template getData< GradientField_T >(gradientFieldID_);
         AvgField_T* velDensAvgField    = blockIterator->template getData< AvgField_T >(avgFieldID_);
         BoundaryHandling_T* boundaryHandling =
            blockIterator->template getData< BoundaryHandling_T >(boundaryHandlingID_);

         // compute dx for the level that the block resides on
         const uint_t blockLevel         = blocks_->getLevel(*blockIterator);
         const uint_t levelScalingFactor = (uint_t(1) << (numLevels - uint_t(1) - blockLevel));
         const real_t dxOnLevel          = dx * real_c(levelScalingFactor);
         const real_t invDelta           = real_c(1) / (dxOnLevel * real_c(2));

         // clang-format off
         WALBERLA_FOR_ALL_CELLS_XYZ(gradientField,
            if (boundaryHandling->isDomain(x, y, z)) {
               // dU/dx
               gradientField->get(x, y, z, uint_c(0)) = (velDensAvgField->get(x + cell_idx_c(1), y, z, uint_c(1)) -
                                                         velDensAvgField->get(x - cell_idx_c(1), y, z, uint_c(1))) *
                                                        invDelta;

               // dU/dy
               gradientField->get(x, y, z, uint_c(1)) = (velDensAvgField->get(x, y + cell_idx_c(1), z, uint_c(1)) -
                                                         velDensAvgField->get(x, y - cell_idx_c(1), z, uint_c(1))) *
                                                        invDelta;

               // dU/dz
               gradientField->get(x, y, z, uint_c(2)) = (velDensAvgField->get(x, y, z + cell_idx_c(1), uint_c(1)) -
                                                         velDensAvgField->get(x, y, z - cell_idx_c(1), uint_c(1))) *
                                                        invDelta;

               // dV/dx
               gradientField->get(x, y, z, uint_c(3)) = (velDensAvgField->get(x + cell_idx_c(1), y, z, uint_c(2)) -
                                                         velDensAvgField->get(x - cell_idx_c(1), y, z, uint_c(2))) *
                                                        invDelta;

               // dV/dy
               gradientField->get(x, y, z, uint_c(4)) = (velDensAvgField->get(x, y + cell_idx_c(1), z, uint_c(2)) -
                                                         velDensAvgField->get(x, y - cell_idx_c(1), z, uint_c(2))) *
                                                        invDelta;

               // dV/dz
               gradientField->get(x, y, z, uint_c(5)) = (velDensAvgField->get(x, y, z + cell_idx_c(1), uint_c(2)) -
                                                         velDensAvgField->get(x, y, z - cell_idx_c(1), uint_c(2))) *
                                                        invDelta;

               // dW/dx
               gradientField->get(x, y, z, uint_c(6)) = (velDensAvgField->get(x + cell_idx_c(1), y, z, uint_c(3)) -
                                                         velDensAvgField->get(x - cell_idx_c(1), y, z, uint_c(3))) *
                                                        invDelta;

               // dW/dy
               gradientField->get(x, y, z, uint_c(7)) = (velDensAvgField->get(x, y + cell_idx_c(1), z, uint_c(3)) -
                                                         velDensAvgField->get(x, y - cell_idx_c(1), z, uint_c(3))) *
                                                        invDelta;

               // dW/dz
               gradientField->get(x, y, z, uint_c(8)) = (velDensAvgField->get(x, y, z + cell_idx_c(1), uint_c(3)) -
                                                         velDensAvgField->get(x, y, z - cell_idx_c(1), uint_c(3))) *
                                                        invDelta;
         }) // WALBERLA_FOR_ALL_CELLS
         // clang-format on
      }
   }
}; // class GradientComputer

// store "instFieldEntries" instantaneous snapshots of density and velocity-components into instFieldID
template< typename PdfField_T, typename BoundaryHandling_T, typename InstField_T >
class InstStorer
{
 public:
   InstStorer(const std::shared_ptr< StructuredBlockForest >& blocks, const ConstBlockDataID& pdfFieldID,
              const BlockDataID& boundaryHandlingID, const ConstBlockDataID& instFieldID, uint_t interval)
      : blocks_(blocks), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID), instFieldID_(instFieldID),
        interval_(interval), entryCounter_(uint_c(0)), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % interval_ != uint_c(0)) { return; }

      // update ghost layer in PDF field as its values are read in computeGradientVorticity()
      blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync(blocks_);
      pdfGhostLayerSync.addPackInfo(
         make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >(pdfFieldID_));
      pdfGhostLayerSync();

      // store density and velocity, and compute gradient and vorticity
      computeGradientVorticity();

      ++entryCounter_;
      if (entryCounter_ >= instFieldEntries) { entryCounter_ = uint_c(0); }
   }

 private:
   const std::shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const BlockDataID instFieldID_;
   const uint_t interval_;

   uint_t entryCounter_;     // number of current entries in instFieldID_
   uint_t executionCounter_; // number of times operator() has been called

   void computeGradientVorticity()
   {
      // position in the field at which data is to be written
      uint_t entryPosition = entryCounter_ * instValues;

      const real_t dx        = real_t(1); // dx on finest level
      const uint_t numLevels = blocks_->getNumberOfLevels();

      for (auto blockIterator = blocks_->begin(); blockIterator != blocks_->end(); ++blockIterator)
      {
         PdfField_T* pdfField = blockIterator->template getData< PdfField_T >(pdfFieldID_);
         BoundaryHandling_T* boundaryHandling =
            blockIterator->template getData< BoundaryHandling_T >(boundaryHandlingID_);
         InstField_T* instField = blockIterator->template getData< InstField_T >(instFieldID_);

         // compute dx for the level that the block resides on
         const uint_t blockLevel         = blocks_->getLevel(*blockIterator);
         const uint_t levelScalingFactor = (uint_t(1) << (numLevels - uint_t(1) - blockLevel));
         const real_t dxOnLevel          = dx * real_c(levelScalingFactor);
         const real_t invDelta           = real_c(1) / (dxOnLevel * real_c(2));

         // clang-format off
         WALBERLA_FOR_ALL_CELLS_XYZ(pdfField,

            if (boundaryHandling->isDomain(x, y, z)) {
               // get density and velocity from PDF field
               Vector3< real_t > velocity;
               const float_t density = float_c(pdfField->getDensityAndVelocity(velocity, x, y, z));

               instField->get(x, y, z, uint_c(0) + entryPosition) = density;
               instField->get(x, y, z, uint_c(1) + entryPosition) = float_c(velocity[0]);
               instField->get(x, y, z, uint_c(2) + entryPosition) = float_c(velocity[1]);
               instField->get(x, y, z, uint_c(3) + entryPosition) = float_c(velocity[2]);

               // get velocities from PDF field (not from instField, to avoid requiring ghost layers in
               // instField which would tremendously increase the memory consumption)
               const Vector3< real_t > velPlusX  = pdfField->getVelocity(x + cell_idx_c(1), y, z);
               const Vector3< real_t > velMinusX = pdfField->getVelocity(x - cell_idx_c(1), y, z);
               const Vector3< real_t > velPlusY  = pdfField->getVelocity(x, y + cell_idx_c(1), z);
               const Vector3< real_t > velMinusY = pdfField->getVelocity(x, y - cell_idx_c(1), z);
               const Vector3< real_t > velPlusZ  = pdfField->getVelocity(x, y, z + cell_idx_c(1));
               const Vector3< real_t > velMinusZ = pdfField->getVelocity(x, y, z - cell_idx_c(1));

               // compute gradients (not stored directly in field to avoid float-conversion before vorticity
               // computation)
               real_t dudx = real_c((velPlusX[0] - velMinusX[0]) * invDelta);
               real_t dudy = real_c((velPlusY[0] - velMinusY[0]) * invDelta);
               real_t dudz = real_c((velPlusZ[0] - velMinusZ[0]) * invDelta);
               real_t dvdx = real_c((velPlusX[1] - velMinusX[1]) * invDelta);
               real_t dvdy = real_c((velPlusY[1] - velMinusY[1]) * invDelta);
               real_t dvdz = real_c((velPlusZ[1] - velMinusZ[1]) * invDelta);
               real_t dwdx = real_c((velPlusX[2] - velMinusX[2]) * invDelta);
               real_t dwdy = real_c((velPlusY[2] - velMinusY[2]) * invDelta);
               real_t dwdz = real_c((velPlusZ[2] - velMinusZ[2]) * invDelta);

               // store gradients
               instField->get(x, y, z, uint_c(4) + entryPosition)  = float_c(dudx);
               instField->get(x, y, z, uint_c(5) + entryPosition)  = float_c(dudy);
               instField->get(x, y, z, uint_c(6) + entryPosition)  = float_c(dudz);
               instField->get(x, y, z, uint_c(7) + entryPosition)  = float_c(dvdx);
               instField->get(x, y, z, uint_c(8) + entryPosition)  = float_c(dvdy);
               instField->get(x, y, z, uint_c(9) + entryPosition)  = float_c(dvdz);
               instField->get(x, y, z, uint_c(10) + entryPosition) = float_c(dwdx);
               instField->get(x, y, z, uint_c(11) + entryPosition) = float_c(dwdy);
               instField->get(x, y, z, uint_c(12) + entryPosition) = float_c(dwdz);

               // compute and store vorticity
               instField->get(x, y, z, uint_c(13) + entryPosition) = float_c(dwdy - dvdz);
               instField->get(x, y, z, uint_c(14) + entryPosition) = float_c(dudz - dwdx);
               instField->get(x, y, z, uint_c(15) + entryPosition) = float_c(dvdx - dudy);
            }
         ) // WALBERLA_FOR_ALL_CELLS
         // clang-format on
      }
   }
}; // class InstStorer

// get the current velocity and density at the given coordinate "evalCoord"
template< typename PdfField_T >
class VelDensPointEvaluator
{
 public:
   VelDensPointEvaluator(const std::shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& pdfFieldID,
                         const Vector3< real_t >& evalCoord, uint_t evalInterval,
                         const std::shared_ptr< std::vector< real_t > >& densityVelocity)
      : blocks_(blocks), pdfFieldID_(pdfFieldID), evalCoord_(evalCoord), evalInterval_(evalInterval),
        densityVelocity_(densityVelocity), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % evalInterval_ != uint_c(0)) { return; }

      // update density and velocity at evalCoord_
      updateDensityVelocity();
   }

 private:
   const std::shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID pdfFieldID_;
   const Vector3< real_t > evalCoord_;
   const uint_t evalInterval_;
   const std::shared_ptr< std::vector< real_t > > densityVelocity_;

   uint_t executionCounter_; // number of times operator() has been called

   void updateDensityVelocity()
   {
      real_t density = real_c(0);
      Vector3< real_t > velocity(real_c(0));

      for (auto blockIterator = blocks_->begin(); blockIterator != blocks_->end(); ++blockIterator)
      {
         if (blockIterator->getAABB().contains(evalCoord_))
         {
            Cell localCell       = blocks_->getBlockLocalCell(*blockIterator, evalCoord_);
            PdfField_T* pdfField = blockIterator->template getData< PdfField_T >(pdfFieldID_);

            velocity = pdfField->getVelocity(localCell);
            density  = pdfField->getDensity(localCell);

            *densityVelocity_ = { density, velocity[0], velocity[1], velocity[2], length(velocity) };
         }
      }
   }
}; // class VelDensPointEvaluator

// get the pressure drop in z-direction bounded by the specified AABB
template< typename PdfField_T, typename BoundaryHandling_T >
class PressureDropEvaluator
{
 public:
   PressureDropEvaluator(const shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& pdfFieldID,
                         const BlockDataID& boundaryHandlingID, const AABB& evalAABB, const uint_t evalInterval,
                         const std::shared_ptr< real_t >& pressureDrop)
      : blocks_(blocks), pdfFieldID_(pdfFieldID), boundaryHandlingID_(boundaryHandlingID), evalAABB_(evalAABB),
        evalInterval_(evalInterval), pressureDrop_(pressureDrop), executionCounter_(uint_c(0))
   {}

   void operator()()
   {
      ++executionCounter_;

      // only evaluate in given intervals
      if (executionCounter_ % evalInterval_ != uint_c(0)) { return; }

      // update pressure drop
      updatePressureDrop();
   }

 private:
   const shared_ptr< StructuredBlockForest > blocks_;
   const BlockDataID pdfFieldID_;
   const BlockDataID boundaryHandlingID_;
   const AABB evalAABB_;
   const uint_t evalInterval_;
   const std::shared_ptr< real_t > pressureDrop_;

   uint_t executionCounter_; // number of times operator() has been called

   void updatePressureDrop()
   {
      real_t lowerDensity       = real_c(0);
      real_t upperDensity       = real_c(0);
      uint_t numCellsLowerSlice = uint_c(0);
      uint_t numCellsUpperSlice = uint_c(0);

      for (auto blockIterator = blocks_->begin(); blockIterator != blocks_->end(); ++blockIterator)
      {
         const uint_t level = blocks_->getLevel(*blockIterator);

         // global lower cell BB slice for density measurement (shifted 1 cell away from porous medium in
         // z-direction)
         CellInterval lowerCellBB = blocks_->getCellBBFromAABB(evalAABB_, level);
         lowerCellBB.zMax()       = lowerCellBB.zMin() - cell_idx_c(1);
         lowerCellBB.zMin()       = lowerCellBB.zMax();
         numCellsLowerSlice       = lowerCellBB.numCells();

         // global upper cell BB slice for density measurement (shifted 1 cell away from porous medium in
         // z-direction)
         CellInterval upperCellBB = blocks_->getCellBBFromAABB(evalAABB_, level);
         upperCellBB.zMin()       = upperCellBB.zMax() + cell_idx_c(1);
         upperCellBB.zMax()       = upperCellBB.zMin();
         numCellsUpperSlice       = upperCellBB.numCells();

         // block's cell BB that intersects the global slices (still in global coordinates)
         CellInterval lowerBlockCellBB = blocks_->getBlockCellBB(*blockIterator);
         lowerBlockCellBB.intersect(lowerCellBB);
         CellInterval upperBlockCellBB = blocks_->getBlockCellBB(*blockIterator);
         upperBlockCellBB.intersect(upperCellBB);

         // transform the global coordinates of relevant cells to block local coordinates
         CellInterval lowerBlockLocalCellBB;
         CellInterval upperBlockLocalCellBB;
         blocks_->transformGlobalToBlockLocalCellInterval(lowerBlockLocalCellBB, *blockIterator, lowerBlockCellBB);
         blocks_->transformGlobalToBlockLocalCellInterval(upperBlockLocalCellBB, *blockIterator, upperBlockCellBB);

         PdfField_T* pdfField = blockIterator->template getData< PdfField_T >(pdfFieldID_);
         BoundaryHandling_T* boundaryHandling =
            blockIterator->template getData< BoundaryHandling_T >(boundaryHandlingID_);

         // sum density of relevant cells
         // clang-format off
         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(lowerBlockLocalCellBB,
                                                if (boundaryHandling->isDomain(x, y, z))
                                                {
                                                   lowerDensity += pdfField->getDensity(x,y,z);
                                                }
         ) // WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ

         WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(upperBlockLocalCellBB,
                                                if (boundaryHandling->isDomain(x, y, z))
                                                {
                                                   upperDensity += pdfField->getDensity(x,y,z);
                                                }
         ) // WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ
         // clang-format on
      }

      mpi::allReduceInplace< real_t >(lowerDensity, mpi::SUM);
      mpi::allReduceInplace< real_t >(upperDensity, mpi::SUM);

      // average density
      lowerDensity /= real_c(numCellsLowerSlice);
      upperDensity /= real_c(numCellsUpperSlice);

      *pressureDrop_ = (upperDensity - lowerDensity) / real_c(3);
   }
}; // class PressureDropEvaluator

// write vector-fields to vtk
template< typename OutputType /*= float */, typename FieldType, uint_t VEC_SIZE_ARG >
class FieldVectorVTKWriter : public vtk::BlockCellDataWriter< OutputType, VEC_SIZE_ARG >
{
 public:
   FieldVectorVTKWriter(const ConstBlockDataID& fieldID, const std::string& id, const cell_idx_t f)
      : vtk::BlockCellDataWriter< OutputType, VEC_SIZE_ARG >(id), bdid_(fieldID), field_(nullptr), f_(f)
   {}

 protected:
   void configure() override
   {
      WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
      field_ = this->block_->template getData< FieldType >(bdid_);
   }

   OutputType evaluate(const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const cell_idx_t f) override
   {
      WALBERLA_ASSERT_NOT_NULLPTR(field_);

      return numeric_cast< OutputType >(field_->get(x, y, z, uint_c(f_) + uint_c(f)));
   }

   const ConstBlockDataID bdid_;
   const FieldType* field_;
   const cell_idx_t f_;
}; // class FieldVectorVTKWriter

////////////////
// PARAMETERS //
////////////////

struct Setup
{
   std::string packedBedGeometryFile; // filename of the file that defines the porous geometry

   Vector3< uint_t > numCoarseBlocks; // number of coarse blocks in x,y,z, direction
   Vector3< uint_t > cellsPerBlock;   // cells inside each block in x,y,z direction (uniform distribution)
   uint_t levels;                     // number of refinement levels

   Vector3< uint_t > domainSize; // domain size in x-, y-, z-direction
   real_t bedShiftZ;             // shift of the packed bed in z-direction (to extend outflow region)
   uint_t particleDiameter;      // cells per diameter of the particles (WITHOUT radial shrinking)
   real_t radiusScale;           // scaling factor for the particles' radius (for radial shrinking)

   real_t omega;                      // relaxation rate
   real_t nu;                         // kinematic viscosity
   real_t Re;                         // particle Reynolds number
   Vector3< real_t > topwallVelocity; // top wall velocity

   uint_t maxTimesteps; // maximum number of time steps to simulate
   uint_t minTimesteps; // minimum number of time steps to simulate (to avoid premature convergence)

   std::string checkpointFilename; // name of the file which stores the blockforest, PdfField and AvgField
   bool storePdfField;             // store the PDF field
   bool loadPdfField;              // load the PDF field
   bool storeBlockForest;          // store the blockforest
   bool loadBlockForest;           // load the blockforest
   bool storeAvgField;             // store the field that contains the time-averaged quantities
   bool loadAvgField;              // load the field that contains the time-averaged quantities
   uint_t avgCounterLoaded;        // number of times the average computation was performed in a loaded averaged field
   uint_t fieldStoreInterval;      // time step interval for storing PdfField and AvgField

   bool vtkFluidField;           // fluid field vtk-output
   uint_t vtkFluidFieldInterval; // time step interval of fluid field vtk-output
   real_t vtkFieldResolution;    // resolution of the fluid and flag field vtk-output
   bool vtkFlagField;            // flag field vtk-output
   bool vtkDomainDecomposition;  // domain decomposition vtk-output
   bool vtkBodies;               // particles vtk-output
   uint_t vtkAvgFieldInterval;   // time step interval of time-averaged field vtk-output

   uint_t pressureDropInterval; // time step interval for evaluation of the pressure drop

   uint_t evalInterval;           // time step interval for evaluation at the specified coordinates and planes
   bool evalVtk;                  // xz-, and xy-planes at P1, P2, P3 vtk-output
   uint_t evalVtkInterval;        // time step interval for writing vtk output at the specified coordinate
   Vector3< real_t > evalCoordP1; // global coordinate of the evaluation point P1
   Vector3< real_t > evalCoordP2; // global coordinate of the evaluation point P2
   Vector3< real_t > evalCoordP3; // global coordinate of the evaluation point P3

   real_t convThreshold; // threshold for considering the simulation converged

   uint_t perfLogInterval; // time step interval in which the performance is logged, i.e., written to output

   void printSetup()
   {
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(packedBedGeometryFile);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(numCoarseBlocks);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(cellsPerBlock);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(levels);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(particleDiameter);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(radiusScale);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(domainSize);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(bedShiftZ);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(omega);
      WALBERLA_LOG_DEVEL_ON_ROOT("tau = " << real_c(1) / omega);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(Re);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(nu);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(topwallVelocity);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(maxTimesteps);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(minTimesteps);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(storePdfField);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(loadPdfField);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(storeBlockForest);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(loadBlockForest);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(storeAvgField);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(loadAvgField);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(fieldStoreInterval);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(avgCounterLoaded);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkFluidField);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkFluidFieldInterval);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkFieldResolution);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkFlagField);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkDomainDecomposition);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(vtkAvgFieldInterval);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pressureDropInterval);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evalInterval);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evalVtk);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evalVtkInterval);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evalCoordP1);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evalCoordP2);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(evalCoordP3);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(convThreshold);
      WALBERLA_LOG_DEVEL_VAR_ON_ROOT(perfLogInterval);
   }
};

///////////////////////
// BOUNDARY HANDLING //
///////////////////////

class MyBoundaryHandling
{
 public:
   MyBoundaryHandling(const weak_ptr< StructuredBlockForest >& forest, const BlockDataID& flagFieldID,
                      const BlockDataID& pdfFieldID, const BlockDataID& bodyFieldID,
                      const Vector3< real_t >& topWallVelocity)
      : forest_(forest), flagFieldID_(flagFieldID), pdfFieldID_(pdfFieldID), bodyFieldID_(bodyFieldID),
        topWallVelocity_(topWallVelocity)
   {}

   BoundaryHandling_T* operator()(IBlock* const block, const StructuredBlockStorage* const storage) const;

 private:
   weak_ptr< StructuredBlockForest > forest_;

   const BlockDataID flagFieldID_;
   const BlockDataID pdfFieldID_;
   const BlockDataID bodyFieldID_;

   const Vector3< real_t > topWallVelocity_;
};

BoundaryHandling_T* MyBoundaryHandling::operator()(IBlock* const block,
                                                   const StructuredBlockStorage* const storage) const
{
   WALBERLA_ASSERT_NOT_NULLPTR(block);
   WALBERLA_ASSERT_NOT_NULLPTR(storage);

   FlagField_T* flagField = block->getData< FlagField_T >(flagFieldID_);
   PdfField_T* pdfField   = block->getData< PdfField_T >(pdfFieldID_);
   BodyField_T* bodyField = block->getData< BodyField_T >(bodyFieldID_);

   const flag_t fluid = flagField->registerFlag(Fluid_Flag);

   BoundaryHandling_T* handling =
      new BoundaryHandling_T("boundary handling", flagField, fluid, NoSlip_T("no slip", NoSlip_Flag, pdfField),
                             UBB_T("velocity bounce back", UBB_Flag, pdfField),
                             Pressure_T("pressure boundary", Pressure_Flag, pdfField, real_c(1)),
                             MO_SBB_T("MO_SBB", MO_SBB_Flag, pdfField, flagField, bodyField, fluid, *storage, *block));

   auto forest = forest_.lock();
   WALBERLA_CHECK_NOT_NULLPTR(forest);

   const uint_t level = forest->getLevel(*block);

   CellInterval domainBB = storage->getDomainCellBB(level);

   // extend the domain in x,y,z direction with the ghost layers at any domain boundary
   domainBB.xMin() -= cell_idx_c(FieldGhostLayers);
   domainBB.xMax() += cell_idx_c(FieldGhostLayers);

   domainBB.yMin() -= cell_idx_c(FieldGhostLayers);
   domainBB.yMax() += cell_idx_c(FieldGhostLayers);

   domainBB.zMin() -= cell_idx_c(FieldGhostLayers);
   domainBB.zMax() += cell_idx_c(FieldGhostLayers);

   // set boundary conditions (unspecified boundaries are periodic implicitly)
   // velocity boundary TOP
   CellInterval top(domainBB.xMin(), domainBB.yMin(), domainBB.zMax() - cell_idx_c(FieldGhostLayers), domainBB.xMax(),
                    domainBB.yMax(), domainBB.zMax());
   storage->transformGlobalToBlockLocalCellInterval(top, *block);
   handling->forceBoundary(UBB_Flag, top, typename UBB_T::Velocity(topWallVelocity_));

   // pressure boundary BOTTOM
   CellInterval bottom(domainBB.xMin(), domainBB.yMin(), domainBB.zMin(), domainBB.xMax(), domainBB.yMax(),
                       domainBB.zMin() + cell_idx_c(FieldGhostLayers));
   storage->transformGlobalToBlockLocalCellInterval(bottom, *block);
   handling->forceBoundary(Pressure_Flag, bottom); // set density to 1

   // mark all remaining cells as fluid cells
   handling->fillWithDomain(domainBB);

   return handling;
}

//////////////////////
// REFINEMENT SETUP //
//////////////////////

// define the refinement region, i.e., assign the appropriate refinement level to each block
static void refinementSelection(SetupBlockForest& forest, const uint_t levels, AABB refinementAABB1,
                                AABB refinementAABB2)
{
   for (auto block = forest.begin(); block != forest.end(); ++block)
   {
      uint_t blockLevel = block->getLevel();
      AABB blockAABB    = block->getAABB();

      // the following needs to be changed when NOT using 3 levels of refinement
      if (blockAABB.intersects(refinementAABB1) && blockLevel < (levels - uint_t(1))) { block->setMarker(true); }
      if (blockAABB.intersects(refinementAABB2) && blockLevel < (levels - uint_t(2))) { block->setMarker(true); }
   }
}

// assign the expected workload and memory load to each block (based on its refinement level)
static void workloadAndMemoryAssignment(SetupBlockForest& forest)
{
   for (auto block = forest.begin(); block != forest.end(); ++block)
   {
      block->setWorkload(numeric_cast< workload_t >(real_c(uint_t(1) << block->getLevel())));
      block->setMemory(numeric_cast< memory_t >(1));
   }
}

// create the block structure with static refinement
static shared_ptr< StructuredBlockForest > createBlockStructure(const AABB& domainAABB, const AABB& refinementAABB1,
                                                                const AABB& refinementAABB2,
                                                                const Vector3< uint_t >& cellsPerBlock,
                                                                uint_t numberOfLevels, bool keepGlobalBlockInformation)
{
   SetupBlockForest sforest;

   // compute the number of fine blocks
   const Vector3< uint_t > numberOfFineBlocksPerDirection(uint_c(domainAABB.size(0)) / cellsPerBlock[0],
                                                          uint_c(domainAABB.size(1)) / cellsPerBlock[1],
                                                          uint_c(domainAABB.size(2)) / cellsPerBlock[2]);

   for (uint_t i = uint_c(0); i < uint_c(3); ++i)
   {
      WALBERLA_CHECK_EQUAL(numberOfFineBlocksPerDirection[i] * cellsPerBlock[i], uint_c(domainAABB.size(i)),
                           "Domain can not be decomposed in direction " << i << " into fine blocks of size "
                                                                        << cellsPerBlock[i]);
   }

   // check validity of number of fine blocks by comparing to number of specified coarse blocks
   const uint_t levelScalingFactor = (uint_c(1) << (numberOfLevels - uint_c(1)));
   const Vector3< uint_t > numberOfCoarseBlocksPerDirection(numberOfFineBlocksPerDirection / levelScalingFactor);

   for (uint_t i = uint_c(0); i < uint_c(3); ++i)
   {
      WALBERLA_CHECK_EQUAL(numberOfCoarseBlocksPerDirection[i] * levelScalingFactor, numberOfFineBlocksPerDirection[i],
                           "Domain can not be refined in direction "
                              << i << " according to the specified number of levels!");
   }

   // mark the region for static refinement
   sforest.addRefinementSelectionFunction(
      std::bind(refinementSelection, std::placeholders::_1, numberOfLevels, refinementAABB1, refinementAABB2));

   // assign workload and memory to each blocks according to the refinement level
   sforest.addWorkloadMemorySUIDAssignmentFunction(std::bind(workloadAndMemoryAssignment, std::placeholders::_1));

   // initialize the blockforest setup with periodicity in x- and y-direction
   sforest.init(domainAABB, numberOfCoarseBlocksPerDirection[0], numberOfCoarseBlocksPerDirection[1],
                numberOfCoarseBlocksPerDirection[2], true, true, false);

   // calculate process distribution
   const memory_t memoryLimit = math::Limits< memory_t >::inf();

   // perform load balancing
   sforest.balanceLoad(blockforest::StaticLevelwiseCurveBalance(true), uint_c(MPIManager::instance()->numProcesses()),
                       real_t(0), memoryLimit, true);

   WALBERLA_LOG_INFO_ON_ROOT(sforest);

   MPIManager::instance()->useWorldComm();

   // create StructuredBlockForest (encapsulates a newly created BlockForest)
   shared_ptr< StructuredBlockForest > sbf = make_shared< StructuredBlockForest >(
      make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()), sforest, keepGlobalBlockInformation),
      cellsPerBlock[0], cellsPerBlock[1], cellsPerBlock[2]);
   sbf->createCellBoundingBoxes();

   return sbf;
}

//////////
// MAIN //
//////////

int main(int argc, char** argv)
{
   Environment walberlaEnv(argc, argv);

   if (argc < 2) { WALBERLA_ABORT("Please specify a parameter file as input argument."); }

   Setup setup;

   ///////////////////////////
   // SIMULATION PROPERTIES //
   ///////////////////////////

   auto domainParameters = walberlaEnv.config()->getOneBlock("DomainParameters");

   setup.packedBedGeometryFile = domainParameters.getParameter< std::string >("packedBedGeometryFile");

   std::ifstream geomFile(setup.packedBedGeometryFile);
   if (!geomFile) { WALBERLA_ABORT("Could not find the geometry file " << setup.packedBedGeometryFile); }
   else
   {
      geomFile.close();
   }

   // cells per diameter (WITHOUT radial shrinking of particles)
   setup.particleDiameter = domainParameters.getParameter< uint_t >("cellsPerDiameter");
   setup.radiusScale      = domainParameters.getParameter< real_t >("radiusScale");

   Vector3< uint_t > relDomainSize = domainParameters.getParameter< Vector3< uint_t > >("relDomainSize");
   // shift the packed bed in z-direction to compensate domain extension (particle coordinates originally created for
   // relDomainSize[2] == 16)
   setup.bedShiftZ  = real_c(relDomainSize[2]) - real_c(16);
   setup.domainSize = relDomainSize * setup.particleDiameter;
   const AABB simulationDomain(real_t(0), real_t(0), real_t(0), real_c(setup.domainSize[0]),
                               real_c(setup.domainSize[1]), real_c(setup.domainSize[2]));

   setup.levels                    = uint_c(3);
   const uint_t finestLevel        = setup.levels - uint_t(1);
   const uint_t levelScalingFactor = (uint_t(1) << finestLevel);

   setup.numCoarseBlocks = domainParameters.getParameter< Vector3< uint_t > >("numCoarseBlocks");

   setup.cellsPerBlock[0] =
      uint_c(ceil(real_c(setup.domainSize[0]) / real_c(setup.numCoarseBlocks[0] * levelScalingFactor)));
   setup.cellsPerBlock[1] =
      uint_c(ceil(real_c(setup.domainSize[1]) / real_c(setup.numCoarseBlocks[1] * levelScalingFactor)));
   setup.cellsPerBlock[2] =
      uint_c(ceil(real_c(setup.domainSize[2]) / real_c(setup.numCoarseBlocks[2] * levelScalingFactor)));

   const Vector3< uint_t > remainingCells(
      setup.cellsPerBlock[0] * setup.numCoarseBlocks[0] * levelScalingFactor - setup.domainSize[0],
      setup.cellsPerBlock[1] * setup.numCoarseBlocks[1] * levelScalingFactor - setup.domainSize[1],
      setup.cellsPerBlock[2] * setup.numCoarseBlocks[2] * levelScalingFactor - setup.domainSize[2]);

   if (setup.numCoarseBlocks[0] < uint_c(3) || setup.numCoarseBlocks[1] < uint_c(3))
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
         "The number of blocks in x- and y-direction must not be smaller than 3. This is a strict "
         "requirement of the PE with periodic boundaries in these directions.");
   }

   if (remainingCells[0] != uint_c(0) || remainingCells[1] != uint_c(0) || remainingCells[2] != uint_c(0))
   {
      WALBERLA_ABORT("The domain size can not be divided by the number of blocks in x- or y-direction.");
   }

   auto simulationParameters = walberlaEnv.config()->getOneBlock("SimulationParameters");

   // read the porous geometry from file
   std::vector< Vector3< real_t > > particleCoordinates;
   readPackedBedFromFile(setup.packedBedGeometryFile, particleCoordinates);

   setup.Re = simulationParameters.getParameter< real_t >(
      "Re"); // Reynolds number (defined using the inlet-velocity and particleDiameter)

   setup.omega = simulationParameters.getParameter< real_t >("omega"); // on finest level
   const real_t omegaLevel0 =
      lbm::collision_model::levelDependentRelaxationParameter(uint_c(0), setup.omega, finestLevel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(omegaLevel0);

   const real_t omegaLevel1 =
      lbm::collision_model::levelDependentRelaxationParameter(uint_c(1), setup.omega, finestLevel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(omegaLevel1);

   const real_t omegaLevel2 =
      lbm::collision_model::levelDependentRelaxationParameter(uint_c(2), setup.omega, finestLevel);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(omegaLevel2);

   setup.nu = lbm::collision_model::viscosityFromOmega(setup.omega);

   const real_t inletVelocity = setup.nu * setup.Re / real_c(setup.particleDiameter);
   setup.topwallVelocity      = Vector3< real_t >(real_c(0), real_c(0), -inletVelocity);

   setup.maxTimesteps = uint_c< real_t >(
      simulationParameters.getParameter< real_t >("maxTimesteps")); // real_t for input in scientific notation

   setup.minTimesteps = uint_c< real_t >(
      simulationParameters.getParameter< real_t >("minTimesteps")); // real_t for input in scientific notation

   auto outputParameters    = walberlaEnv.config()->getOneBlock("OutputParameters");
   setup.checkpointFilename = std::to_string(setup.particleDiameter) + "_" + std::to_string(setup.levels) + "_" +
                              std::to_string(setup.numCoarseBlocks[0]) + "x" +
                              std::to_string(setup.numCoarseBlocks[1]) + "x" + std::to_string(setup.numCoarseBlocks[2]);
   const std::string pdfFieldFile    = setup.checkpointFilename + ".pdfField";
   const std::string blockforestFile = setup.checkpointFilename + ".blockforest";
   const std::string avgFieldFile    = setup.checkpointFilename + ".avgField";

   setup.storePdfField    = outputParameters.getParameter< bool >("storePdfField");
   setup.loadPdfField     = outputParameters.getParameter< bool >("loadPdfField");
   setup.storeBlockForest = outputParameters.getParameter< bool >("storeBlockForest");
   setup.loadBlockForest  = outputParameters.getParameter< bool >("loadBlockForest");

   setup.storeAvgField      = outputParameters.getParameter< bool >("storeAvgField");
   setup.loadAvgField       = outputParameters.getParameter< bool >("loadAvgField");
   setup.avgCounterLoaded   = uint_c< real_t >(outputParameters.getParameter< real_t >("avgCounterLoaded"));
   setup.fieldStoreInterval = uint_c< real_t >(outputParameters.getParameter< real_t >("fieldStoreInterval"));

   setup.vtkFluidField          = outputParameters.getParameter< bool >("vtkFluidField");
   setup.vtkFluidFieldInterval  = uint_c< real_t >(outputParameters.getParameter< real_t >("vtkFluidFieldInterval"));
   setup.vtkFieldResolution     = outputParameters.getParameter< real_t >("vtkFieldResolution");
   setup.vtkFlagField           = outputParameters.getParameter< bool >("vtkFlagField");
   setup.vtkDomainDecomposition = outputParameters.getParameter< bool >("vtkDomainDecomposition");
   setup.vtkBodies              = outputParameters.getParameter< bool >("vtkBodies");
   setup.vtkAvgFieldInterval    = uint_c< real_t >(outputParameters.getParameter< real_t >("vtkAvgFieldInterval"));

   setup.pressureDropInterval = uint_c< real_t >(outputParameters.getParameter< real_t >("pressureDropInterval"));

   setup.convThreshold = outputParameters.getParameter< real_t >("convThreshold");

   setup.evalInterval       = uint_c< real_t >(outputParameters.getParameter< real_t >("evalInterval"));
   setup.evalVtk            = outputParameters.getParameter< bool >("evalVtk");
   setup.evalVtkInterval    = uint_c(instFieldEntries) * setup.evalInterval;
   Vector3< real_t > evalP1 = outputParameters.getParameter< Vector3< real_t > >("evalP1");
   Vector3< real_t > evalP2 = outputParameters.getParameter< Vector3< real_t > >("evalP2");
   Vector3< real_t > evalP3 = outputParameters.getParameter< Vector3< real_t > >("evalP3");

   evalP1[2] += setup.bedShiftZ;
   evalP2[2] += setup.bedShiftZ;
   evalP3[2] += setup.bedShiftZ;

   for (uint_t i = uint_c(0); i != uint_c(3); ++i)
   {
      setup.evalCoordP1[i] = evalP1[i] * real_c(setup.particleDiameter);
      setup.evalCoordP2[i] = evalP2[i] * real_c(setup.particleDiameter);
      setup.evalCoordP3[i] = evalP3[i] * real_c(setup.particleDiameter);
   }

   // definition of slices through the domain at which velocity and density will be averaged
   const AABB xyPlaneP1AABB(real_c(0), real_c(0), setup.evalCoordP1[2], real_c(setup.domainSize[0]),
                            real_c(setup.domainSize[1]), setup.evalCoordP1[2] + real_c(0.1));
   const AABB xyPlaneP2AABB(real_c(0), real_c(0), setup.evalCoordP2[2], real_c(setup.domainSize[0]),
                            real_c(setup.domainSize[1]), setup.evalCoordP2[2] + real_c(0.1));
   const AABB xyPlaneP3AABB(real_c(0), real_c(0), setup.evalCoordP3[2], real_c(setup.domainSize[0]),
                            real_c(setup.domainSize[1]), setup.evalCoordP3[2] + real_c(0.1));
   const AABB xzPlaneAABB(real_c(0), real_c(2) * real_c(setup.particleDiameter),
                          setup.bedShiftZ * real_c(setup.particleDiameter), real_c(setup.domainSize[0]),
                          real_c(2) * real_c(setup.particleDiameter) + real_c(0.1), real_c(setup.domainSize[2]));

   setup.perfLogInterval = uint_c< real_t >(outputParameters.getParameter< real_t >("perfLogInterval"));

   setup.printSetup();
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(instFieldEntries);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(blockforestFile);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(avgFieldFile);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(pdfFieldFile);

   //////////////////////
   // SIMULATION SETUP //
   //////////////////////

   // refinement region in and around the packed bed
   const AABB porousMediaRefinementAABB(real_c(0), real_c(0), real_c(5.5) * real_c(setup.particleDiameter),
                                        real_c(setup.domainSize[0]), real_c(setup.domainSize[1]),
                                        real_c(setup.domainSize[2]) - real_c(4.5) * real_c(setup.particleDiameter));

   // refinement region near the outflow
   const AABB outflowRefinementAABB(real_c(0), real_c(0), real_c(3) * real_c(setup.particleDiameter),
                                    real_c(setup.domainSize[0]), real_c(setup.domainSize[1]),
                                    real_c(5.5) * real_c(setup.particleDiameter));

   // create the uniform block grid
   shared_ptr< StructuredBlockForest > blocks;

   if (setup.loadBlockForest)
   {
      // load blockforest from file
      MPIManager::instance()->useWorldComm();

      blocks = make_shared< StructuredBlockForest >(
         make_shared< BlockForest >(uint_c(MPIManager::instance()->rank()), blockforestFile.c_str(), true, false),
         setup.cellsPerBlock[0], setup.cellsPerBlock[1], setup.cellsPerBlock[2]);
      blocks->createCellBoundingBoxes();
   }
   else
   {
      blocks = createBlockStructure(simulationDomain, porousMediaRefinementAABB, outflowRefinementAABB,
                                    setup.cellsPerBlock, setup.levels, false);

      // write blockforest to file
      if (setup.storeBlockForest) { blocks->getBlockForest().saveToFile(blockforestFile); }
   }

   // define the lattice model
   LatticeModel_T latticeModel = LatticeModel_T(omegaLevel0);

   // add the PDF field
   BlockDataID pdfFieldID;

   // layout of PDF field MUST be fzyx when using lbmpy codegen and WALBERLA_OPTIMIZE_FOR_LOCALHOST
   if (setup.loadPdfField)
   {
      // load PDF field from file
      shared_ptr< lbm::internal::PdfFieldHandling< LatticeModel_T > > dataHandling =
         make_shared< lbm::internal::PdfFieldHandling< LatticeModel_T > >(
            blocks, latticeModel, false, Vector3< real_t >(real_c(0)), real_c(1), FieldGhostLayers, field::fzyx);
      pdfFieldID = (blocks->getBlockStorage()).loadBlockData(pdfFieldFile, dataHandling, "pdf field");
   }
   else
   {
      pdfFieldID = lbm::addPdfFieldToStorage(blocks, "pdf field", latticeModel, setup.topwallVelocity, real_c(1),
                                             FieldGhostLayers, field::fzyx);
   }

   // add avgField that stores the time-averaged density and velocity of all fluid cells:
   // density, velocity_x, velocity_y, velocity_z, velocity_magn, density_old, velocity_x_old, velocity_y_old,
   // velocity_z_old, velocity_magn_old
   BlockDataID avgFieldID;

   if (setup.loadAvgField)
   {
      // load avgField from file
      std::shared_ptr< AvgFieldHandling< AvgField_T > > dataHandling =
         std::make_shared< AvgFieldHandling< AvgField_T > >(blocks);
      avgFieldID = (blocks->getBlockStorage()).loadBlockData(avgFieldFile, dataHandling, "avg field");
   }
   else
   {
      avgFieldID =
         field::addToStorage< AvgField_T >(blocks, "avg field", real_c(0), field::fzyx, FieldGhostLayers, false);
   }

   // add field for caching instantaneous data that will be written in junks to decrease file-IO (explicitly without
   // ghost layers)
   BlockDataID instFieldID;
   instFieldID = field::addToStorage< InstField_T >(blocks, "inst field", float_c(0), field::fzyx, uint_c(0), false);

   // add field that contains the gradient of the averaged velocity (explicitly without ghost layers)
   BlockDataID gradientFieldID =
      field::addToStorage< GradientField_T >(blocks, "gradient field", real_c(0), field::fzyx, uint_c(0), false);

   // add flag field
   BlockDataID flagFieldID = field::addFlagFieldToStorage< FlagField_T >(blocks, "flag field", FieldGhostLayers, false);

   // add field for particle creation (bodies)
   BlockDataID bodyFieldID =
      field::addToStorage< BodyField_T >(blocks, "body field", nullptr, field::fzyx, FieldGhostLayers, false);

   // add boundary handling
   BlockDataID boundaryHandlingID = blocks->addStructuredBlockData< BoundaryHandling_T >(
      MyBoundaryHandling(blocks, flagFieldID, pdfFieldID, bodyFieldID, setup.topwallVelocity), "boundary handling");

   /////////////////
   // PE COUPLING //
   /////////////////

   // generate IDs of specified PE body types
   pe::SetBodyTypeIDs< BodyTypeTuple >::execute();

   // add global body storage for PE bodies
   shared_ptr< pe::BodyStorage > globalBodyStorage = make_shared< pe::BodyStorage >();

   // add block-local body storage
   auto bodyStorageID = blocks->addBlockData(pe::createStorageDataHandling< BodyTypeTuple >(), "pe body storage");

   // get AABB of the packed bed (tangential to the outermost spherical particles in z-direction)
   AABB packedBedAABB;
   real_t sphereMinZ = real_c(setup.domainSize[2]);
   real_t sphereMaxZ = real_c(0);

   // create the porous medium with spherical particles
   for (auto it = particleCoordinates.begin(); it != particleCoordinates.end(); ++it)
   {
      Vector3< real_t > coord;

      coord = ((*it) * real_c(setup.particleDiameter));

      // shift the packed bed bedShiftZ*dp upwards in z-direction to extend outflow
      coord[2] += setup.bedShiftZ * real_c(setup.particleDiameter);

      // create spherical particle and shrink them radially
      pe::createSphere(*globalBodyStorage, blocks->getBlockStorage(), bodyStorageID, 0, coord,
                       real_c(setup.particleDiameter) * real_c(0.5) * setup.radiusScale);

      if (coord[2] < sphereMinZ) { sphereMinZ = coord[2]; }
      if (coord[2] > sphereMaxZ) { sphereMaxZ = coord[2]; }
   }

   // shift coordinates such that they are tangential to the spheres
   sphereMinZ -= real_c(setup.particleDiameter) * real_c(0.5) * setup.radiusScale;
   sphereMaxZ += real_c(setup.particleDiameter) * real_c(0.5) * setup.radiusScale;

   // update AABB of the packed bed
   packedBedAABB =
      AABB(real_c(0), real_c(0), sphereMinZ, real_c(setup.domainSize[0]), real_c(setup.domainSize[1]), sphereMaxZ);

   // synchronize each particle among all relevant blocks; needs to be called multiple times if a body spans more than
   // one block because information is only transferred to direct neighboring blocks in each call
   const real_t overlap = real_c(1.5);
   for (uint_t i = uint_c(0); i <= uint_c(real_c(0.5) * real_c(setup.particleDiameter)); ++i)
   {
      pe::syncShadowOwners< BodyTypeTuple >(blocks->getBlockForest(), bodyStorageID, nullptr, overlap);
   }

   // assign boundary conditions and map the particles into the domain
   pe_coupling::mapMovingBodies< BoundaryHandling_T >(*blocks, boundaryHandlingID, bodyStorageID, *globalBodyStorage,
                                                      bodyFieldID, MO_SBB_Flag);

   // mark solid cells that have fluid cell flags in their neighborhood (for writing them to vtk)
   setSolidNearFluidFlag< BoundaryHandling_T >(blocks, boundaryHandlingID, SolidNearFluid_Flag);

   real_t voidFraction = getVoidFraction< BoundaryHandling_T >(blocks, boundaryHandlingID, packedBedAABB);
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(voidFraction);

   //////////////////////////
   // VTK OUTPUT - GENERAL //
   //////////////////////////

   decltype(vtk::createVTKOutput_BlockData(blocks)) fluidFieldVTK = nullptr;

   // vtk fluid field output
   if (setup.vtkFluidField)
   {
      // ghost layers MUST NOT be written when sampling resolution is adjusted (see below)
      fluidFieldVTK = vtk::createVTKOutput_BlockData(blocks, "fluid_field", setup.vtkFluidFieldInterval, uint_c(0));

      blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync(blocks);
      pdfGhostLayerSync.addPackInfo(make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >(pdfFieldID));
      fluidFieldVTK->addBeforeFunction(pdfGhostLayerSync);
      fluidFieldVTK->setSamplingResolution(setup.vtkFieldResolution);

      field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldID);
      fluidFilter.addFlag(Fluid_Flag);
      fluidFilter.addFlag(SolidNearFluid_Flag); // write first layer of solid cells to vtk
      fluidFieldVTK->addCellInclusionFilter(fluidFilter);

      auto velocityWriter =
         make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "VelocityFromPDF");
      auto densityWriter = make_shared< lbm::DensityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "DensityFromPDF");
      fluidFieldVTK->addCellDataWriter(velocityWriter);
      fluidFieldVTK->addCellDataWriter(densityWriter);

      fluidFieldVTK->write(true, 0);
   }

   // vtk flag field output
   if (setup.vtkFlagField)
   {
      auto flagFieldVTK = vtk::createVTKOutput_BlockData(blocks, "flag_field", uint_c(1), uint_c(0));

      flagFieldVTK->setSamplingResolution(setup.vtkFieldResolution);
      flagFieldVTK->addCellDataWriter(make_shared< field::VTKWriter< FlagField_T > >(flagFieldID, "FlagField"));

      flagFieldVTK->write();
   }

   // vtk domain decomposition output
   if (setup.vtkDomainDecomposition)
   {
      vtk::writeDomainDecomposition(blocks, "domain_decomposition", "vtk_out", "write_call", true, true, 0);
   }

   // vtk PE bodies output
   if (setup.vtkBodies)
   {
      auto bodiesVTK = vtk::createVTKOutput_PointData(
         make_shared< pe::SphereVtkOutput >(bodyStorageID, blocks->getBlockStorage()), "bodies");
      bodiesVTK->write();
   }

   ////////////////////////////////////////////////////////////
   // VTK OUTPUT - TIME-AVERAGE AND INSTANTANEOUS QUANTITIES //
   ////////////////////////////////////////////////////////////

   decltype(vtk::createVTKOutput_BlockData(blocks)) avgVelDensVTK   = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) avgXZPlaneVTK   = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) avgXYPlaneP1VTK = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) avgXYPlaneP2VTK = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) avgXYPlaneP3VTK = nullptr;

   decltype(vtk::createVTKOutput_BlockData(blocks)) instXZPlaneVTK   = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) instXYPlaneP1VTK = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) instXYPlaneP2VTK = nullptr;
   decltype(vtk::createVTKOutput_BlockData(blocks)) instXYPlaneP3VTK = nullptr;

   if (setup.evalVtk)
   {
      blockforest::communication::NonUniformBufferedScheme< stencil::D3Q27 > pdfGhostLayerSync(blocks);
      pdfGhostLayerSync.addPackInfo(make_shared< lbm::refinement::PdfFieldSyncPackInfo< LatticeModel_T > >(pdfFieldID));
      auto velocityWriter =
         make_shared< lbm::VelocityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "VelocityFromPDF");
      auto densityWriter = make_shared< lbm::DensityVTKWriter< LatticeModel_T, float > >(pdfFieldID, "DensityFromPDF");

      // filters for xz-plane
      field::FlagFieldCellFilter< FlagField_T > fluidFilter(flagFieldID);
      fluidFilter.addFlag(Fluid_Flag);
      fluidFilter.addFlag(SolidNearFluid_Flag);
      vtk::ChainedFilter combinedSliceFilterXZPlane;
      vtk::AABBCellFilter aabbSliceFilterXZPlane(xzPlaneAABB);
      combinedSliceFilterXZPlane.addFilter(fluidFilter);
      combinedSliceFilterXZPlane.addFilter(aabbSliceFilterXZPlane);

      // filters for xy-plane at P1
      vtk::ChainedFilter combinedSliceFilterXYPlaneP1;
      vtk::AABBCellFilter aabbSliceFilterXYPlaneP1(xyPlaneP1AABB);
      combinedSliceFilterXYPlaneP1.addFilter(fluidFilter);
      combinedSliceFilterXYPlaneP1.addFilter(aabbSliceFilterXYPlaneP1);

      // filters for xy-plane at P2
      vtk::ChainedFilter combinedSliceFilterXYPlaneP2;
      vtk::AABBCellFilter aabbSliceFilterXYPlaneP2(xyPlaneP2AABB);
      combinedSliceFilterXYPlaneP2.addFilter(fluidFilter);
      combinedSliceFilterXYPlaneP2.addFilter(aabbSliceFilterXYPlaneP2);

      // filters for xy-plane at P3
      vtk::ChainedFilter combinedSliceFilterXYPlaneP3;
      vtk::AABBCellFilter aabbSliceFilterXYPlaneP3(xyPlaneP3AABB);
      combinedSliceFilterXYPlaneP3.addFilter(fluidFilter);
      combinedSliceFilterXYPlaneP3.addFilter(aabbSliceFilterXYPlaneP3);

      // vtk output for instantaneous quantities
      instXZPlaneVTK = vtk::createVTKOutput_BlockData(blocks, "inst-xz-plane", uint_c(1), uint_c(0));
      instXZPlaneVTK->addCellInclusionFilter(combinedSliceFilterXZPlane);

      instXYPlaneP1VTK = vtk::createVTKOutput_BlockData(blocks, "inst-xy-plane-p1", uint_c(1), uint_c(0));
      instXYPlaneP1VTK->addCellInclusionFilter(combinedSliceFilterXYPlaneP1);

      instXYPlaneP2VTK = vtk::createVTKOutput_BlockData(blocks, "inst-xy-plane-p2", uint_c(1), uint_c(0));
      instXYPlaneP2VTK->addCellInclusionFilter(combinedSliceFilterXYPlaneP2);

      instXYPlaneP3VTK = vtk::createVTKOutput_BlockData(blocks, "inst-xy-plane-p3", uint_c(1), uint_c(0));
      instXYPlaneP3VTK->addCellInclusionFilter(combinedSliceFilterXYPlaneP3);

      // iterate over whole field for writing all cached values
      for (uint_t i = uint_c(0); i != uint_c(instValues * instFieldEntries); i += uint_c(instValues))
      {
         std::string writeTimestep = std::to_string(setup.evalInterval * uint_c(real_c(i / instValues)));

         instXZPlaneVTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, densValues > >(
            instFieldID, "density_" + writeTimestep, i + uint_c(0)));
         instXZPlaneVTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, velValues > >(
            instFieldID, "velocity_" + writeTimestep, i + densValues));
         instXZPlaneVTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, gradValues > >(
            instFieldID, "gradient_" + writeTimestep, i + densValues + velValues));
         instXZPlaneVTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, vortValues > >(
            instFieldID, "vorticity_" + writeTimestep, i + densValues + velValues + gradValues));

         instXYPlaneP1VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, densValues > >(
            instFieldID, "density_" + writeTimestep, i + uint_c(0)));
         instXYPlaneP1VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, velValues > >(
            instFieldID, "velocity_" + writeTimestep, i + densValues));
         instXYPlaneP1VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, gradValues > >(
            instFieldID, "gradient_" + writeTimestep, i + densValues + velValues));
         instXYPlaneP1VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, vortValues > >(
            instFieldID, "vorticity_" + writeTimestep, i + densValues + velValues + gradValues));

         instXYPlaneP2VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, densValues > >(
            instFieldID, "density_" + writeTimestep, i + uint_c(0)));
         instXYPlaneP2VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, velValues > >(
            instFieldID, "velocity_" + writeTimestep, i + densValues));
         instXYPlaneP2VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, gradValues > >(
            instFieldID, "gradient_" + writeTimestep, i + densValues + velValues));
         instXYPlaneP2VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, vortValues > >(
            instFieldID, "vorticity_" + writeTimestep, i + densValues + velValues + gradValues));

         instXYPlaneP3VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, densValues > >(
            instFieldID, "density_" + writeTimestep, i + uint_c(0)));
         instXYPlaneP3VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, velValues > >(
            instFieldID, "velocity_" + writeTimestep, i + densValues));
         instXYPlaneP3VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, gradValues > >(
            instFieldID, "gradient_" + writeTimestep, i + densValues + velValues));
         instXYPlaneP3VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, InstField_T, vortValues > >(
            instFieldID, "vorticity_" + writeTimestep, i + densValues + velValues + gradValues));
      }

      // vtk output for time-averaged quantities in the whole domain
      avgVelDensVTK = vtk::createVTKOutput_BlockData(blocks, "avgVelDens", uint_c(1), uint_c(0));
      avgVelDensVTK->addCellInclusionFilter(fluidFilter);
      avgVelDensVTK->addBeforeFunction(GradientComputer< BoundaryHandling_T, AvgField_T, GradientField_T >(
         blocks, boundaryHandlingID, avgFieldID, gradientFieldID));
      avgVelDensVTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, densValues > >(avgFieldID, "density", uint_c(0)));
      avgVelDensVTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, velValues > >(avgFieldID, "velocity", densValues));
      avgVelDensVTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, GradientField_T, gradValues > >(
         gradientFieldID, "gradient", uint_c(0)));

      // vtk output for time-averaged xz-plane
      avgXZPlaneVTK = vtk::createVTKOutput_BlockData(blocks, "avg-xz-plane", uint_c(1), uint_c(0));
      avgXZPlaneVTK->addBeforeFunction(GradientComputer< BoundaryHandling_T, AvgField_T, GradientField_T >(
         blocks, boundaryHandlingID, avgFieldID, gradientFieldID));
      avgXZPlaneVTK->addCellInclusionFilter(combinedSliceFilterXZPlane);
      avgXZPlaneVTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, densValues > >(avgFieldID, "density", uint_c(0)));
      avgXZPlaneVTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, velValues > >(avgFieldID, "velocity", densValues));
      avgXZPlaneVTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, GradientField_T, gradValues > >(
         gradientFieldID, "gradient", uint_c(0)));

      // vtk output for time-averaged xy-plane at P1
      avgXYPlaneP1VTK = vtk::createVTKOutput_BlockData(blocks, "avg-xy-plane-p1", uint_c(1), uint_c(0));
      avgXYPlaneP1VTK->addBeforeFunction(GradientComputer< BoundaryHandling_T, AvgField_T, GradientField_T >(
         blocks, boundaryHandlingID, avgFieldID, gradientFieldID));
      avgXYPlaneP1VTK->addCellInclusionFilter(combinedSliceFilterXYPlaneP1);
      avgXYPlaneP1VTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, densValues > >(avgFieldID, "density", uint_c(0)));
      avgXYPlaneP1VTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, velValues > >(avgFieldID, "velocity", densValues));
      avgXYPlaneP1VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, GradientField_T, gradValues > >(
         gradientFieldID, "gradient", uint_c(0)));

      // vtk output for time-averaged xy-plane at P2
      avgXYPlaneP2VTK = vtk::createVTKOutput_BlockData(blocks, "avg-xy-plane-p2", uint_c(1), uint_c(0));
      avgXYPlaneP2VTK->addBeforeFunction(GradientComputer< BoundaryHandling_T, AvgField_T, GradientField_T >(
         blocks, boundaryHandlingID, avgFieldID, gradientFieldID));
      avgXYPlaneP2VTK->addCellInclusionFilter(combinedSliceFilterXYPlaneP2);
      avgXYPlaneP2VTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, densValues > >(avgFieldID, "density", uint_c(0)));
      avgXYPlaneP2VTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, velValues > >(avgFieldID, "velocity", densValues));
      avgXYPlaneP2VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, GradientField_T, gradValues > >(
         gradientFieldID, "gradient", uint_c(0)));

      // vtk output for time-averaged xy-plane at P3
      avgXYPlaneP3VTK = vtk::createVTKOutput_BlockData(blocks, "avg-xy-plane-p3", uint_c(1), uint_c(0));
      avgXYPlaneP3VTK->addBeforeFunction(GradientComputer< BoundaryHandling_T, AvgField_T, GradientField_T >(
         blocks, boundaryHandlingID, avgFieldID, gradientFieldID));
      avgXYPlaneP3VTK->addCellInclusionFilter(combinedSliceFilterXYPlaneP3);
      avgXYPlaneP3VTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, densValues > >(avgFieldID, "density", uint_c(0)));
      avgXYPlaneP3VTK->addCellDataWriter(
         make_shared< FieldVectorVTKWriter< float, AvgField_T, velValues > >(avgFieldID, "velocity", densValues));
      avgXYPlaneP3VTK->addCellDataWriter(make_shared< FieldVectorVTKWriter< float, GradientField_T, gradValues > >(
         gradientFieldID, "gradient", uint_c(0)));
   }

   ///////////////
   // TIME LOOP //
   ///////////////

   // initialize time loop
   SweepTimeloop timeloop(blocks->getBlockStorage(), setup.maxTimesteps);

   // add sweep for LBM with refinement to timeloop
   auto lbmSweep           = make_shared< LatticeModel_T::Sweep >(pdfFieldID);
   auto refinementTimestep = lbm::refinement::makeTimeStep< LatticeModel_T, BoundaryHandling_T, LatticeModel_T::Sweep >(
      blocks, lbmSweep, pdfFieldID, boundaryHandlingID);
   timeloop.addFuncBeforeTimeStep(
      makeSharedFunctor< lbm::refinement::TimeStep< LatticeModel_T, LatticeModel_T::Sweep, BoundaryHandling_T > >(
         refinementTimestep),
      "LBM sweep (with refinement)");

   // add performance logger to time loop
   lbm::PerformanceLogger< FlagField_T > perfLogger(blocks, flagFieldID, Fluid_Flag, setup.perfLogInterval);
   timeloop.addFuncAfterTimeStep(perfLogger, "Performance Logging");

   // add vtk output to time loop
   if (fluidFieldVTK != nullptr) { timeloop.addFuncAfterTimeStep(vtk::writeFiles(fluidFieldVTK), "VTK - fluid field"); }

   // add sweep for getting velocity and density at predefined coordinates
   auto densityVelocityP1 = std::make_shared< std::vector< real_t > >();
   auto densityVelocityP2 = std::make_shared< std::vector< real_t > >();
   auto densityVelocityP3 = std::make_shared< std::vector< real_t > >();
   timeloop.addFuncAfterTimeStep(
      VelDensPointEvaluator< PdfField_T >(blocks, pdfFieldID, setup.evalCoordP1, setup.evalInterval, densityVelocityP1),
      "Evaluation at P1");
   timeloop.addFuncAfterTimeStep(
      VelDensPointEvaluator< PdfField_T >(blocks, pdfFieldID, setup.evalCoordP2, setup.evalInterval, densityVelocityP2),
      "Evaluation at P2");
   timeloop.addFuncAfterTimeStep(
      VelDensPointEvaluator< PdfField_T >(blocks, pdfFieldID, setup.evalCoordP3, setup.evalInterval, densityVelocityP3),
      "Evaluation at P3");

   const std::string fileResultP1 = std::string("result-p1-") + std::to_string(uint_c(setup.evalCoordP1[0])) + "_" +
                                    std::to_string(uint_c(setup.evalCoordP1[1])) + "_" +
                                    std::to_string(uint_c(setup.evalCoordP1[2])) + ".txt";
   const std::string fileResultP2 = std::string("result-p2-") + std::to_string(uint_c(setup.evalCoordP2[0])) + "_" +
                                    std::to_string(uint_c(setup.evalCoordP2[1])) + "_" +
                                    std::to_string(uint_c(setup.evalCoordP2[2])) + ".txt";
   const std::string fileResultP3 = std::string("result-p3-") + std::to_string(uint_c(setup.evalCoordP3[0])) + "_" +
                                    std::to_string(uint_c(setup.evalCoordP3[1])) + "_" +
                                    std::to_string(uint_c(setup.evalCoordP3[2])) + ".txt";

   // add sweep for pressure drop evaluation to time loop
   auto pressureDrop = std::make_shared< real_t >(real_c(0));
   timeloop.addFuncAfterTimeStep(
      PressureDropEvaluator< PdfField_T, BoundaryHandling_T >(blocks, pdfFieldID, boundaryHandlingID, packedBedAABB,
                                                              setup.pressureDropInterval, pressureDrop),
      "Pressure drop evaluation");

   // add sweep for average computation to time loop
   timeloop.addFuncAfterTimeStep(
      VelDensAverager< PdfField_T, BoundaryHandling_T, AvgField_T >(blocks, pdfFieldID, boundaryHandlingID, avgFieldID,
                                                                    setup.evalInterval, setup.avgCounterLoaded),
      "Time-averaging");

   // add sweep for storing of instantaneous quantities to time loop
   timeloop.addFuncAfterTimeStep(InstStorer< PdfField_T, BoundaryHandling_T, InstField_T >(
                                    blocks, pdfFieldID, boundaryHandlingID, instFieldID, setup.evalInterval),
                                 "Instantaneous quantity evaluation");

   // add L2 computation of difference between two successive field average values
   auto avgDensityDiffL2  = std::make_shared< real_t >(real_c(1e15)); // initialize with clearly wrong values
   auto avgVelocityDiffL2 = std::make_shared< Vector3< real_t > >(real_c(1e15));
   timeloop.addFuncAfterTimeStep(
      AvgDiffEvaluator< BoundaryHandling_T, AvgField_T >(blocks, boundaryHandlingID, avgFieldID, setup.evalInterval,
                                                         avgDensityDiffL2, avgVelocityDiffL2),
      "Average difference norm computation");

   WcTimingPool timeloopTiming;

   // run time loop
   for (uint_t timestep = uint_c(1); timestep <= setup.maxTimesteps; ++timestep)
   {
      // perform time step
      timeloop.singleStep(timeloopTiming, true);

      // print current L2 norm of relative difference between current and former time-averaged density, and velocity
      if (timestep % setup.evalInterval == uint_c(0))
      {
         WALBERLA_LOG_DEVEL_ON_ROOT(
            "L2 norm of relative difference between current and former average density =" << *avgDensityDiffL2);
         WALBERLA_LOG_DEVEL_ON_ROOT(
            "L2 norm of relative difference between current and former average velocity =" << *avgVelocityDiffL2);
      }

      // check convergence
      if (timestep % setup.evalInterval == uint_c(0) && timestep > setup.minTimesteps &&
          *avgDensityDiffL2 < setup.convThreshold && (*avgVelocityDiffL2)[0] < setup.convThreshold &&
          (*avgVelocityDiffL2)[1] < setup.convThreshold && (*avgVelocityDiffL2)[2] < setup.convThreshold)
      {
         // write the final time-averaged quantities to vtk output
         if (avgVelDensVTK != nullptr) { avgVelDensVTK->write(true, 0); }
         if (avgXZPlaneVTK != nullptr) { avgXZPlaneVTK->write(true, 0); }
         if (avgXYPlaneP1VTK != nullptr) { avgXYPlaneP1VTK->write(true, 0); }
         if (avgXYPlaneP2VTK != nullptr) { avgXYPlaneP2VTK->write(true, 0); }
         if (avgXYPlaneP3VTK != nullptr) { avgXYPlaneP3VTK->write(true, 0); }

         // store the final PDF field
         if (setup.storePdfField)
         {
            blocks->saveBlockData("final_" + std::to_string(timestep) + "_" + pdfFieldFile, pdfFieldID);
         }

         // store the final time-averaged velocity-density field
         if (setup.storeAvgField)
         {
            blocks->saveBlockData("final_" + std::to_string(timestep) + "_" + avgFieldFile, avgFieldID);
         }

         WALBERLA_LOG_DEVEL_ON_ROOT(
            "Time-averaging has reached a stationary value. Simulation is considered converged.")
         break;
      }

      // write current density and velocity at P1, P2, and P3 to files
      if (timestep % setup.evalInterval == uint_c(0))
      {
         if (!(*densityVelocityP1).empty()) { writeVector(*densityVelocityP1, timestep, fileResultP1); }
         if (!(*densityVelocityP2).empty()) { writeVector(*densityVelocityP2, timestep, fileResultP2); }
         if (!(*densityVelocityP3).empty()) { writeVector(*densityVelocityP3, timestep, fileResultP3); }
      }

      // print current pressure drop
      if (timestep % setup.pressureDropInterval == uint_c(0))
      {
         WALBERLA_LOG_DEVEL_ON_ROOT("Pressure drop =" << *pressureDrop);
      }

      // print time loop statistics
      if (timestep % setup.perfLogInterval == uint_c(0)) { timeloopTiming.logResultOnRoot(); }

      // store PDF field
      if (setup.storePdfField && timestep % setup.fieldStoreInterval == uint_c(0))
      {
         blocks->saveBlockData(std::to_string(timestep) + "_" + pdfFieldFile, pdfFieldID);
      }

      // store time-averaged field
      if (setup.storeAvgField && timestep % setup.fieldStoreInterval == uint_c(0))
      {
         blocks->saveBlockData(std::to_string(timestep) + "_" + avgFieldFile, avgFieldID);
      }

      // write the time-averaged quantities to vtk output
      if (timestep % setup.vtkAvgFieldInterval == uint_c(0))
      {
         if (avgVelDensVTK != nullptr) { avgVelDensVTK->write(true, 0); }
         if (avgXZPlaneVTK != nullptr) { avgXZPlaneVTK->write(true, 0); }
         if (avgXYPlaneP1VTK != nullptr) { avgXYPlaneP1VTK->write(true, 0); }
         if (avgXYPlaneP2VTK != nullptr) { avgXYPlaneP2VTK->write(true, 0); }
         if (avgXYPlaneP3VTK != nullptr) { avgXYPlaneP3VTK->write(true, 0); }
      }

      // write instantaneous quantities to vtk output
      if (timestep % setup.evalVtkInterval == uint_c(0))
      {
         if (instXZPlaneVTK != nullptr) { instXZPlaneVTK->write(true, 0); }
         if (instXYPlaneP1VTK != nullptr) { instXYPlaneP1VTK->write(true, 0); }
         if (instXYPlaneP2VTK != nullptr) { instXYPlaneP2VTK->write(true, 0); }
         if (instXYPlaneP3VTK != nullptr) { instXYPlaneP3VTK->write(true, 0); }
      }
   }

   timeloopTiming.logResultOnRoot();

   return EXIT_SUCCESS;
}

template< typename BoundaryHandling_T >
void setSolidNearFluidFlag(const std::shared_ptr< StructuredBlockForest >& blocks,
                           const BlockDataID& boundaryHandlingID, FlagUID flag)
{
   for (auto blockIterator = blocks->begin(); blockIterator != blocks->end(); ++blockIterator)
   {
      BoundaryHandling_T* boundaryHandling = blockIterator->getData< BoundaryHandling_T >(boundaryHandlingID);
      FlagField_T* flagField               = boundaryHandling->getFlagField();

      auto solidNearFluidFlag = flagField->getOrRegisterFlag(flag);
      auto solidSphereFlag    = flagField->getFlag(MO_SBB_Flag);
      auto fluidFlag          = flagField->getFlag(Fluid_Flag);

      // clang-format off
      WALBERLA_FOR_ALL_CELLS_XYZ(flagField,
                                 if (flagField->isFlagSet(x, y, z, solidSphereFlag)) {
                                    for (auto dir = stencil::D3Q6::beginNoCenter(); dir != stencil::D3Q6::end(); ++dir)
                                    {
                                       // check if neighboring cell is fluid cell
                                       if (flagField->isFlagSet(x + cell_idx_c(dir.cx()), y + cell_idx_c(dir.cy()),
                                                                z + cell_idx_c(dir.cz()), fluidFlag))
                                       {
                                          // set flag in this cell
                                          flagField->addFlag(x, y, z, solidNearFluidFlag);
                                          break;
                                       }
                                    }
                                 }
      ) // WALBERLA_FOR_ALL_CELLS_XYZ
      //clang-format on
   }
}

template< typename BoundaryHandling_T >
real_t getVoidFraction(const shared_ptr< StructuredBlockForest >& blocks, const BlockDataID& boundaryHandlingID,
                       const AABB& computeRegionAABB)
{
   WALBERLA_LOG_DEVEL_VAR_ON_ROOT(computeRegionAABB);
   uint_t numTotalCells = uint_c(0);
   uint_t numFluidCells = uint_c(0);

   for (auto blockIterator = blocks->begin(); blockIterator != blocks->end(); ++blockIterator)
   {
      const uint_t level           = blocks->getLevel(*blockIterator);
      CellInterval computeRegionBB = blocks->getCellBBFromAABB(computeRegionAABB, level);

      // block's cell BB that intersects the computeRegionBB (in global coordinates)
      CellInterval blockCellBB = blocks->getBlockCellBB(*blockIterator);
      blockCellBB.intersect(computeRegionBB);

      // transform the global coordinates of relevant cells to block local coordinates
      CellInterval blockLocalCellBB;
      blocks->transformGlobalToBlockLocalCellInterval(blockLocalCellBB, *blockIterator, blockCellBB);

      BoundaryHandling_T* boundaryHandling = blockIterator->getData< BoundaryHandling_T >(boundaryHandlingID);

      // count total number of cells and number of fluid cells
      // clang-format off
      WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ(blockLocalCellBB,
                                             ++numTotalCells;
                                             if (boundaryHandling->isDomain(x, y, z)) { ++numFluidCells; }
      ) // WALBERLA_FOR_ALL_CELLS_IN_INTERVAL_XYZ
      // clang-format on
   }

   // sum values over all processes
   mpi::allReduceInplace< uint_t >(numTotalCells, mpi::SUM);
   mpi::allReduceInplace< uint_t >(numFluidCells, mpi::SUM);

   // compute void fraction
   return real_c(numFluidCells) / real_c(numTotalCells);
}

void readPackedBedFromFile(const std::string& filename, std::vector< Vector3< real_t > >& particleCoordinates)
{
   // read location of particles from file
   std::ifstream file;
   file.open(filename);
   std::string line;

   std::getline(file, line); // skip file header
   while (std::getline(file, line))
   {
      std::vector< std::string > v = string_split(line, "\t");
      particleCoordinates.emplace_back(stringToNum< real_t >((v[0])), stringToNum< real_t >((v[1])),
                                       stringToNum< real_t >((v[2])));
   }
   file.close();
}

template< typename T >
void writeVector(const std::vector< T >& data, const uint_t& timestep, const std::string& filename)
{
   std::fstream file;
   file.open(filename, std::fstream::app);

   file << timestep;
   for (const auto i : data)
   {
      file << "\t" << toStringWithPrecision(i, 12);
   }
   file << "\n";
   file.close();
}

template< typename T >
std::string toStringWithPrecision(const T a_value, const int n /*= 12*/)
{
   std::ostringstream out;
   out.precision(n);

   out << std::fixed << a_value;
   return out.str();
}

} // namespace walberla

int main(int argc, char** argv) { return walberla::main(argc, argv); }
