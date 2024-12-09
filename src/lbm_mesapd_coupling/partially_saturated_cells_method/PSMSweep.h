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
//! \file PSMSweep.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Modification of pe_coupling/partially_saturated_cells_method/PSMSweep.h
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

#include "lbm/lattice_model/all.h"
#include "lbm/sweeps/StreamPull.h"
#include "lbm/sweeps/SweepBase.h"

#include "lbm_mesapd_coupling/DataTypes.h"
#include "lbm_mesapd_coupling/utility/ParticleFunctions.h"

#include "mesa_pd/common/ParticleFunctions.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

// functions to calculate the cell coverage weighting B
template< int Weighting_T >
real_t calculateWeighting(const real_t& /*epsilon*/, const real_t& /*tau*/)
{
   WALBERLA_STATIC_ASSERT(Weighting_T == 1 || Weighting_T == 2);
   exit(EXIT_FAILURE);
}
template<>
inline real_t calculateWeighting< 1 >(const real_t& epsilon, const real_t& /*tau*/)
{
   return epsilon;
}
template<>
inline real_t calculateWeighting< 2 >(const real_t& epsilon, const real_t& tau)
{
   return epsilon * (tau - real_c(0.5)) / ((real_c(1) - epsilon) + (tau - real_c(0.5)));
}

/*!\brief LBM sweep for the partially saturated cells method
 *
 * Literature:
 * Original idea:  Noble et al. - A lattice-Boltzmann method for partially saturated computational cells, 1998
 * Overview paper: Owen et al. - An efficient framework for fluid-structure interaction using the lattice Boltzmann
 * method and immersed moving boundaries, 2010
 *
 * Based on the overview paper, at least three different versions for the calculation of the solid collision term
 * (omega_S) and two different versions for the cell coverage weighting function (B) are available. They can be selected
 * via the template parameters SolidCollision_T ( = {1,2,3} ) and Weighting_T ( = {1,2} ). SolidCollision_T = 1: Eq. 28
 * Weighting_T = 1: B = solid volume fraction SolidCollision_T = 2: Eq. 29                     Weighting_T = 2: Eq. 30
 * SolidCollision_T = 3: Eq. 33
 *
 * For the calculation of the force acting on the particle, an additional minus sign has to be added compared to Eqs. 31
 * and 32 in the paper.
 *
 */
template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
class PSMSweep : public lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >
{
 public:
   using PdfField_T =
      typename lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::PdfField_T;
   using Stencil_T = typename LatticeModel_T::Stencil;

   PSMSweep(const BlockDataID& pdfFieldID, const BlockDataID& particleAndVolumeFractionFieldID,
            const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
            const Filter_T& _filter                         = walberla::field::DefaultEvaluationFilter(),
            const DensityVelocityIn_T& _densityVelocityIn   = lbm::DefaultDensityEquilibriumVelocityCalculation(),
            const DensityVelocityOut_T& _densityVelocityOut = lbm::DefaultDensityVelocityCallback())
      : lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >(
           pdfFieldID, _filter, _densityVelocityIn, _densityVelocityOut),
        particleAndVolumeFractionFieldID_(particleAndVolumeFractionFieldID), blockStorage_(blockStorage), accessor_(ac)
   {}

   PSMSweep(const BlockDataID& src, const BlockDataID& dst, const BlockDataID& particleAndVolumeFractionFieldID,
            const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
            const Filter_T& _filter                         = walberla::field::DefaultEvaluationFilter(),
            const DensityVelocityIn_T& _densityVelocityIn   = lbm::DefaultDensityEquilibriumVelocityCalculation(),
            const DensityVelocityOut_T& _densityVelocityOut = lbm::DefaultDensityVelocityCallback())
      : lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >(
           src, dst, _filter, _densityVelocityIn, _densityVelocityOut),
        particleAndVolumeFractionFieldID_(particleAndVolumeFractionFieldID), blockStorage_(blockStorage), accessor_(ac)
   {}

   void operator()(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0))
   {
      streamCollide(block, numberOfGhostLayersToInclude);
   }

   void streamCollide(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0));

   void stream(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0));
   void collide(IBlock* const block, const uint_t numberOfGhostLayersToInclude = uint_t(0));

   inline ParticleAndVolumeFractionField_T* getParticleAndVolumeFractionField(IBlock* const block) const
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);
      ParticleAndVolumeFractionField_T* particleAndVolumeFractionField =
         block->getData< ParticleAndVolumeFractionField_T >(particleAndVolumeFractionFieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(particleAndVolumeFractionField);
      return particleAndVolumeFractionField;
   }

   inline void getFields(IBlock* const block, PdfField_T*& src, PdfField_T*& dst,
                         ParticleAndVolumeFractionField_T*& particleAndVolumeFractionField)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(block);

      src                            = this->getSrcField(block);
      dst                            = this->getDstField(block, src);
      particleAndVolumeFractionField = getParticleAndVolumeFractionField(block);

      WALBERLA_ASSERT_NOT_NULLPTR(src);
      WALBERLA_ASSERT_NOT_NULLPTR(dst);
      WALBERLA_ASSERT_NOT_NULLPTR(particleAndVolumeFractionField);

      WALBERLA_ASSERT_EQUAL(src->xyzSize(), dst->xyzSize());
      WALBERLA_ASSERT_EQUAL(src->xyzSize(), particleAndVolumeFractionField->xyzSize());
   }

 private:
   const BlockDataID particleAndVolumeFractionFieldID_;
   shared_ptr< StructuredBlockStorage > blockStorage_;
   const shared_ptr< ParticleAccessor_T > accessor_;
};

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
void PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T, Weighting_T,
               ParticleAccessor_T >::streamCollide(IBlock* const block, const uint_t numberOfGhostLayersToInclude)
{
   PdfField_T* src(nullptr);
   PdfField_T* dst(nullptr);
   ParticleAndVolumeFractionField_T* particleAndVolumeFractionField(nullptr);

   getFields(block, src, dst, particleAndVolumeFractionField);

   const real_t dxCurrentLevel      = blockStorage_->dx(blockStorage_->getLevel(*block));
   const real_t lengthScalingFactor = dxCurrentLevel;
   const real_t forceScalingFactor  = lengthScalingFactor * lengthScalingFactor;

   WALBERLA_ASSERT_GREATER(src->nrOfGhostLayers(), numberOfGhostLayersToInclude);
   WALBERLA_ASSERT_GREATER_EQUAL(dst->nrOfGhostLayers(), numberOfGhostLayersToInclude);

   const auto& lm = src->latticeModel();
   dst->resetLatticeModel(
      lm); /* required so that member functions for getting density and equilibrium velocity can be called for dst! */

   this->filter(*block);
   this->densityVelocityIn(*block);
   this->densityVelocityOut(*block);

   const auto& collisionModel = lm.collisionModel();

   const real_t omega = collisionModel.omega();
   const real_t tau   = real_c(1) / omega;

   const cell_idx_t xSize = cell_idx_c(src->xSize());
   const cell_idx_t ySize = cell_idx_c(src->ySize());
   const cell_idx_t zSize = cell_idx_c(src->zSize());
   const cell_idx_t gl    = cell_idx_c(numberOfGhostLayersToInclude);
   for (cell_idx_t z = -gl; z < (zSize + gl); ++z)
   {
      for (cell_idx_t y = -gl; y < (ySize + gl); ++y)
      {
         for (cell_idx_t x = -gl; x < (xSize + gl); ++x)
         {
            if (this->filter(x, y, z))
            {
               using namespace stencil;

               real_t pdfs[Stencil_T::Size];

               // stream pull & temporal storage of PDFs
               for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
               {
                  dst->get(x, y, z, d.toIdx()) = src->get(x - d.cx(), y - d.cy(), z - d.cz(), d.toIdx());
                  pdfs[d.toIdx()]              = dst->get(x, y, z, d.toIdx());
               }

               Vector3< real_t > velocity;
               real_t rho = this->densityVelocityIn(velocity, dst, x, y, z);

               this->densityVelocityOut(x, y, z, lm, velocity, rho);

               // equilibrium distributions
               auto pdfs_equ = lbm::EquilibriumDistribution< LatticeModel_T >::get(velocity, rho);

               // possible external forcing on fluid
               const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >(
                  x, y, z, velocity, rho, omega, collisionModel.omega_bulk(), collisionModel.omega_odd());

               // check if particle is present
               if (particleAndVolumeFractionField->get(x, y, z).size() != size_t(0))
               {
                  // total coverage ratio in the cell
                  real_t Bn = real_t(0);

                  // averaged solid collision operator for all intersecting particles s
                  // = \sum_s B_s * \Omega_s_i
                  std::vector< real_t > omega_n(Stencil_T::Size, real_t(0));

                  // get center of cell
                  Vector3< real_t > cellCenter = blockStorage_->getBlockLocalCellCenter(*block, Cell(x, y, z));

                  for (auto particleFracIt = particleAndVolumeFractionField->get(x, y, z).begin();
                       particleFracIt != particleAndVolumeFractionField->get(x, y, z).end(); ++particleFracIt)
                  {
                     real_t omega_s(real_c(0));
                     Vector3< real_t > forceOnParticle(real_c(0));

                     real_t Bs = calculateWeighting< Weighting_T >((*particleFracIt).second, tau);
                     Bn += Bs;

                     // particle velocity at cell center
                     const size_t idx = accessor_->uidToIdx(particleFracIt->first);
                     WALBERLA_ASSERT_UNEQUAL(idx, accessor_->getInvalidIdx(), "Index of particle is invalid!");
                     const auto particleVelocity = mesa_pd::getVelocityAtWFPoint(idx, *accessor_, cellCenter);

                     // equilibrium distributions with solid velocity
                     auto pdfs_equ_solid = lbm::EquilibriumDistribution< LatticeModel_T >::get(particleVelocity, rho);

                     for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
                     {
                        // Different solid collision operators available
                        if (SolidCollision_T == 1)
                        {
                           omega_s =
                              pdfs[d.toInvIdx()] - pdfs_equ[d.toInvIdx()] + pdfs_equ_solid[d.toIdx()] - pdfs[d.toIdx()];
                        }
                        else if (SolidCollision_T == 2)
                        {
                           omega_s = pdfs_equ_solid[d.toIdx()] - pdfs[d.toIdx()] +
                                     (real_c(1) - omega) * (pdfs[d.toIdx()] - pdfs_equ[d.toIdx()]);
                        }
                        else if (SolidCollision_T == 3)
                        {
                           omega_s = pdfs[d.toInvIdx()] - pdfs_equ_solid[d.toInvIdx()] + pdfs_equ_solid[d.toIdx()] -
                                     pdfs[d.toIdx()];
                        }
                        else
                        {
                           WALBERLA_STATIC_ASSERT(SolidCollision_T >= 1 && SolidCollision_T <= 3);
                           exit(EXIT_FAILURE);
                        }
                        real_t BsOmegaS = Bs * omega_s;

                        omega_n[d.toIdx()] += BsOmegaS;

                        forceOnParticle[0] -= BsOmegaS * real_c(d.cx());
                        forceOnParticle[1] -= BsOmegaS * real_c(d.cy());
                        forceOnParticle[2] -= BsOmegaS * real_c(d.cz());
                     }

                     // scale force when using refinement with (dx)^3 / dt
                     forceOnParticle *= forceScalingFactor;

                     // only if cell inside inner domain
                     if (dst->isInInnerPart(Cell(x, y, z)))
                     {
                        // apply force (and automatically torque) on particle
                        lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic(idx, *accessor_, forceOnParticle,
                                                                               cellCenter);
                     }
                  }

                  // collide step
                  for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
                  {
                     // external forcing
                     const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >(
                        x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[d.toIdx()], real_c(d.cx()),
                        real_c(d.cy()), real_c(d.cz()), omega, collisionModel.omega_bulk(), collisionModel.omega_odd());

                     dst->get(x, y, z, d.toIdx()) =
                        pdfs[d.toIdx()] - omega * (real_c(1) - Bn) * (pdfs[d.toIdx()] - pdfs_equ[d.toIdx()]) // SRT
                        + omega_n[d.toIdx()] + (real_c(1) - Bn) * forceTerm;
                  }
               }
               else
               {
                  // SRT collide step
                  for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
                  {
                     // external forcing
                     const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >(
                        x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[d.toIdx()], real_c(d.cx()),
                        real_c(d.cy()), real_c(d.cz()), omega, collisionModel.omega_bulk(), collisionModel.omega_odd());

                     dst->get(x, y, z, d.toIdx()) =
                        pdfs[d.toIdx()] - omega * (pdfs[d.toIdx()] - pdfs_equ[d.toIdx()]) + forceTerm;
                  }
               }
            }
         }
      }
   }
   src->swapDataPointers(dst);
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
void PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T, Weighting_T,
               ParticleAccessor_T >::stream(IBlock* const block, const uint_t numberOfGhostLayersToInclude)
{
   PdfField_T* src(nullptr);
   PdfField_T* dst(nullptr);
   lbm::SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::getFields(block, src, dst);
   lbm::StreamPull< LatticeModel_T >::execute(src, dst, block, this->filter_, numberOfGhostLayersToInclude);
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
void PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T, Weighting_T,
               ParticleAccessor_T >::collide(IBlock* const block, const uint_t numberOfGhostLayersToInclude)
{
   PdfField_T* src                                                  = this->getSrcField(block);
   ParticleAndVolumeFractionField_T* particleAndVolumeFractionField = getParticleAndVolumeFractionField(block);

   WALBERLA_ASSERT_GREATER(src->nrOfGhostLayers(), numberOfGhostLayersToInclude);

   const real_t dxCurrentLevel      = blockStorage_->dx(blockStorage_->getLevel(*block));
   const real_t lengthScalingFactor = dxCurrentLevel;
   const real_t forceScalingFactor  = lengthScalingFactor * lengthScalingFactor;

   this->filter(*block);
   this->densityVelocityIn(*block);
   this->densityVelocityOut(*block);

   const auto& lm             = src->latticeModel();
   const auto& collisionModel = lm.collisionModel();

   const real_t omega = collisionModel.omega();
   const real_t tau   = real_c(1) / omega;

   const cell_idx_t xSize = cell_idx_c(src->xSize());
   const cell_idx_t ySize = cell_idx_c(src->ySize());
   const cell_idx_t zSize = cell_idx_c(src->zSize());
   const cell_idx_t gl    = cell_idx_c(numberOfGhostLayersToInclude);
   for (cell_idx_t z = -gl; z < (zSize + gl); ++z)
   {
      for (cell_idx_t y = -gl; y < (ySize + gl); ++y)
      {
         for (cell_idx_t x = -gl; x < (xSize + gl); ++x)
         {
            if (this->filter(x, y, z))
            {
               using namespace stencil;

               real_t pdfs[Stencil_T::Size];

               // temporal storage of PDFs
               for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
               {
                  pdfs[d.toIdx()] = src->get(x, y, z, d.toIdx());
               }

               Vector3< real_t > velocity;
               real_t rho = this->densityVelocityIn(velocity, src, x, y, z);

               this->densityVelocityOut(x, y, z, lm, velocity, rho);

               // equilibrium distributions
               auto pdfs_equ = lbm::EquilibriumDistribution< LatticeModel_T >::get(velocity, rho);

               // possible external forcing on fluid
               const auto commonForceTerms = lm.forceModel().template directionIndependentTerms< LatticeModel_T >(
                  x, y, z, velocity, rho, omega, collisionModel.omega_bulk(), collisionModel.omega_odd());

               // check if a particle is present
               if (particleAndVolumeFractionField->get(x, y, z).size() != size_t(0))
               {
                  // total coverage ratio in the cell
                  real_t Bn = real_t(0);

                  // averaged solid collision operator for all intersecting particles s
                  // = \sum_s B_s * \Omega_s_i
                  std::vector< real_t > omega_n(Stencil_T::Size, real_t(0));

                  // get center of cell
                  Vector3< real_t > cellCenter = blockStorage_->getBlockLocalCellCenter(*block, Cell(x, y, z));

                  for (auto particleFracIt = particleAndVolumeFractionField->get(x, y, z).begin();
                       particleFracIt != particleAndVolumeFractionField->get(x, y, z).end(); ++particleFracIt)
                  {
                     real_t omega_s(real_c(0));
                     Vector3< real_t > forceOnParticle(real_c(0));

                     real_t Bs = calculateWeighting< Weighting_T >((*particleFracIt).second, tau);
                     Bn += Bs;

                     // particle velocity at cell center
                     const size_t idx = accessor_->uidToIdx(particleFracIt->first);
                     WALBERLA_ASSERT_UNEQUAL(idx, accessor_->getInvalidIdx(), "Index of particle is invalid!");
                     const auto particleVelocity = mesa_pd::getVelocityAtWFPoint(idx, *accessor_, cellCenter);

                     // equilibrium distributions with solid velocity
                     auto pdfs_equ_solid = lbm::EquilibriumDistribution< LatticeModel_T >::get(particleVelocity, rho);

                     for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
                     {
                        // Different solid collision operators available
                        if (SolidCollision_T == 1)
                        {
                           omega_s =
                              pdfs[d.toInvIdx()] - pdfs_equ[d.toInvIdx()] + pdfs_equ_solid[d.toIdx()] - pdfs[d.toIdx()];
                        }
                        else if (SolidCollision_T == 2)
                        {
                           omega_s = pdfs_equ_solid[d.toIdx()] - pdfs[d.toIdx()] +
                                     (real_c(1) - omega) * (pdfs[d.toIdx()] - pdfs_equ[d.toIdx()]);
                        }
                        else if (SolidCollision_T == 3)
                        {
                           omega_s = pdfs[d.toInvIdx()] - pdfs_equ_solid[d.toInvIdx()] + pdfs_equ_solid[d.toIdx()] -
                                     pdfs[d.toIdx()];
                        }
                        else
                        {
                           WALBERLA_STATIC_ASSERT(SolidCollision_T >= 1 && SolidCollision_T <= 3);
                           exit(EXIT_FAILURE);
                        }
                        real_t BsOmegaS = Bs * omega_s;

                        omega_n[d.toIdx()] += BsOmegaS;

                        forceOnParticle[0] -= BsOmegaS * real_c(d.cx());
                        forceOnParticle[1] -= BsOmegaS * real_c(d.cy());
                        forceOnParticle[2] -= BsOmegaS * real_c(d.cz());
                     }

                     // scale force when using refinement with (dx)^3 / dt
                     forceOnParticle *= forceScalingFactor;

                     // only if cell inside inner domain
                     if (src->isInInnerPart(Cell(x, y, z)))
                     {
                        // apply force (and automatically torque) on particle
                        lbm_mesapd_coupling::addHydrodynamicForceAtWFPosAtomic(idx, *accessor_, forceOnParticle,
                                                                               cellCenter);
                     }
                  }

                  // collide step
                  for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
                  {
                     // external forcing
                     const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >(
                        x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[d.toIdx()], real_c(d.cx()),
                        real_c(d.cy()), real_c(d.cz()), omega, collisionModel.omega_bulk(), collisionModel.omega_odd());

                     src->get(x, y, z, d.toIdx()) =
                        pdfs[d.toIdx()] - omega * (real_c(1) - Bn) * (pdfs[d.toIdx()] - pdfs_equ[d.toIdx()]) // SRT
                        + omega_n[d.toIdx()] + (real_c(1) - Bn) * forceTerm;
                  }
               }
               else
               {
                  // SRT collide step
                  for (auto d = Stencil_T::begin(); d != Stencil_T::end(); ++d)
                  {
                     // external forcing
                     const real_t forceTerm = lm.forceModel().template forceTerm< LatticeModel_T >(
                        x, y, z, velocity, rho, commonForceTerms, LatticeModel_T::w[d.toIdx()], real_c(d.cx()),
                        real_c(d.cy()), real_c(d.cz()), omega, collisionModel.omega_bulk(), collisionModel.omega_odd());

                     src->get(x, y, z, d.toIdx()) =
                        pdfs[d.toIdx()] - omega * (pdfs[d.toIdx()] - pdfs_equ[d.toIdx()]) + forceTerm;
                  }
               }
            }
         }
      }
   }
}

/////////////////////////////
// makePSMSweep FUNCTIONS //
////////////////////////////

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T,
                      Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& pdfFieldID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const Filter_T& filter, const DensityVelocityIn_T& densityVelocityIn,
                const DensityVelocityOut_T& densityVelocityOut)
{
   using PSMS_T = PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T,
                            Weighting_T, ParticleAccessor_T >;
   return shared_ptr< PSMS_T >(new PSMS_T(pdfFieldID, particleAndVolumeFractionFieldID, blockStorage, ac, filter,
                                          densityVelocityIn, densityVelocityOut));
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T,
                      Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& srcID, const BlockDataID& dstID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const Filter_T& filter, const DensityVelocityIn_T& densityVelocityIn,
                const DensityVelocityOut_T& densityVelocityOut)
{
   using PSMS_T = PSMSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, SolidCollision_T,
                            Weighting_T, ParticleAccessor_T >;
   return shared_ptr< PSMS_T >(new PSMS_T(srcID, dstID, particleAndVolumeFractionFieldID, blockStorage, ac, filter,
                                          densityVelocityIn, densityVelocityOut));
}

// only block data IDs of PDF data

template< typename LatticeModel_T, int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter,
                      lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                      SolidCollision_T, Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& pdfFieldID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac)
{
   return makePSMSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter,
                        lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                        SolidCollision_T, Weighting_T >(
      pdfFieldID, particleAndVolumeFractionFieldID, blockStorage, ac, walberla::field::DefaultEvaluationFilter(),
      lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback());
}

template< typename LatticeModel_T, int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter,
                      lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                      SolidCollision_T, Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& srcID, const BlockDataID& dstID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac)
{
   return makePSMSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter,
                        lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                        SolidCollision_T, Weighting_T >(
      srcID, dstID, particleAndVolumeFractionFieldID, blockStorage, ac, walberla::field::DefaultEvaluationFilter(),
      lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback());
}

// block data IDs of PDF data + flag field as filter

template< typename LatticeModel_T, typename FlagField_T, int SolidCollision_T, int Weighting_T,
          typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                      lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                      SolidCollision_T, Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& pdfFieldID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const ConstBlockDataID& flagFieldID, const Set< FlagUID >& cellsToEvaluate)
{
   return makePSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                        lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                        SolidCollision_T, Weighting_T >(
      pdfFieldID, particleAndVolumeFractionFieldID, blockStorage, ac,
      walberla::field::FlagFieldEvaluationFilter< FlagField_T >(flagFieldID, cellsToEvaluate),
      lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback());
}

template< typename LatticeModel_T, typename FlagField_T, int SolidCollision_T, int Weighting_T,
          typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                      lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                      SolidCollision_T, Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& srcID, const BlockDataID& dstID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const ConstBlockDataID& flagFieldID, const Set< FlagUID >& cellsToEvaluate)
{
   return makePSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                        lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::DefaultDensityVelocityCallback,
                        SolidCollision_T, Weighting_T >(
      srcID, dstID, particleAndVolumeFractionFieldID, blockStorage, ac,
      walberla::field::FlagFieldEvaluationFilter< FlagField_T >(flagFieldID, cellsToEvaluate),
      lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::DefaultDensityVelocityCallback());
}

// block data IDs of PDF data + flag field as filter + block data ID of velocity field (out)

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, int SolidCollision_T,
          int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                      lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::VelocityCallback< VelocityField_T >,
                      SolidCollision_T, Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& pdfFieldID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const ConstBlockDataID& flagFieldID, const Set< FlagUID >& cellsToEvaluate,
                const BlockDataID& velocityFieldID)
{
   return makePSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                        lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::VelocityCallback< VelocityField_T >,
                        SolidCollision_T, Weighting_T >(
      pdfFieldID, particleAndVolumeFractionFieldID, blockStorage, ac,
      walberla::field::FlagFieldEvaluationFilter< FlagField_T >(flagFieldID, cellsToEvaluate),
      lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::VelocityCallback< VelocityField_T >(velocityFieldID));
}

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, int SolidCollision_T,
          int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                      lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::VelocityCallback< VelocityField_T >,
                      SolidCollision_T, Weighting_T, ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& srcID, const BlockDataID& dstID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const ConstBlockDataID& flagFieldID, const Set< FlagUID >& cellsToEvaluate,
                const BlockDataID& velocityFieldID)
{
   return makePSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                        lbm::DefaultDensityEquilibriumVelocityCalculation, lbm::VelocityCallback< VelocityField_T >,
                        SolidCollision_T, Weighting_T >(
      srcID, dstID, particleAndVolumeFractionFieldID, blockStorage, ac,
      walberla::field::FlagFieldEvaluationFilter< FlagField_T >(flagFieldID, cellsToEvaluate),
      lbm::DefaultDensityEquilibriumVelocityCalculation(), lbm::VelocityCallback< VelocityField_T >(velocityFieldID));
}

// block data IDs of PDF data + flag field as filter + block data IDs of velocity and density field (out)

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, typename DensityField_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                      lbm::DefaultDensityEquilibriumVelocityCalculation,
                      lbm::DensityVelocityCallback< VelocityField_T, DensityField_T >, SolidCollision_T, Weighting_T,
                      ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& pdfFieldID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const ConstBlockDataID& flagFieldID, const Set< FlagUID >& cellsToEvaluate,
                const BlockDataID& velocityFieldID, const BlockDataID& densityFieldID)
{
   return makePSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                        lbm::DefaultDensityEquilibriumVelocityCalculation,
                        lbm::DensityVelocityCallback< VelocityField_T, DensityField_T >, SolidCollision_T,
                        Weighting_T >(
      pdfFieldID, particleAndVolumeFractionFieldID, blockStorage, ac,
      walberla::field::FlagFieldEvaluationFilter< FlagField_T >(flagFieldID, cellsToEvaluate),
      lbm::DefaultDensityEquilibriumVelocityCalculation(),
      lbm::DensityVelocityCallback< VelocityField_T, DensityField_T >(velocityFieldID, densityFieldID));
}

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, typename DensityField_T,
          int SolidCollision_T, int Weighting_T, typename ParticleAccessor_T >
shared_ptr< PSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                      lbm::DefaultDensityEquilibriumVelocityCalculation,
                      lbm::DensityVelocityCallback< VelocityField_T, DensityField_T >, SolidCollision_T, Weighting_T,
                      ParticleAccessor_T > >
   makePSMSweep(const BlockDataID& srcID, const BlockDataID& dstID, const BlockDataID& particleAndVolumeFractionFieldID,
                const shared_ptr< StructuredBlockStorage >& blockStorage, const shared_ptr< ParticleAccessor_T >& ac,
                const ConstBlockDataID& flagFieldID, const Set< FlagUID >& cellsToEvaluate,
                const BlockDataID& velocityFieldID, const BlockDataID& densityFieldID)
{
   return makePSMSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter< FlagField_T >,
                        lbm::DefaultDensityEquilibriumVelocityCalculation,
                        lbm::DensityVelocityCallback< VelocityField_T, DensityField_T >, SolidCollision_T,
                        Weighting_T >(
      srcID, dstID, particleAndVolumeFractionFieldID, blockStorage, ac,
      walberla::field::FlagFieldEvaluationFilter< FlagField_T >(flagFieldID, cellsToEvaluate),
      lbm::DefaultDensityEquilibriumVelocityCalculation(),
      lbm::DensityVelocityCallback< VelocityField_T, DensityField_T >(velocityFieldID, densityFieldID));
}

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
