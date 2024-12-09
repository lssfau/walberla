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
//! \file PSMUtility.h
//! \ingroup lbm_mesapd_coupling
//! \author Samuel Kemmler <samuel.kemmler@fau.de>
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//! \brief Modification of pe_coupling/partially_saturated_cells_method/PSMUtility.h
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/Field.h"
#include "field/GhostLayerField.h"

#include "lbm/field/PdfField.h"

namespace walberla
{
namespace lbm_mesapd_coupling
{
namespace psm
{

/*!\brief Evaluates the macroscopic velocity for a given cell when using the PSM method.
 *
 * returns the cell local macroscopic velocity for the PSM method
 * cell local velocity = \f$ (1-\sum_s B_s) * u + \sum_s B_s * v_s \f$
 * u = fluid velocity
 * v_s = velocity of object s at cell center
 *
 * lbm::PdfField< LatticeModel_T > is typically typedef'ed as PdfField_T;
 * GhostLayerField< std::vector< std::pair< walberla::id_t, real_t > >, 1 > is typically typedef'ed as
 * ParticleAndVolumeFractionField_T;
 *
 * Weighting_T is like in the PSMSweep.
 */
template< typename LatticeModel_T, int Weighting_T, typename ParticleAccessor_T >
Vector3< real_t > getPSMMacroscopicVelocity(const IBlock& block, lbm::PdfField< LatticeModel_T >* pdfField,
                                            ParticleAndVolumeFractionField_T* particleAndVolumeFractionField,
                                            StructuredBlockStorage& blockStorage, const Cell& cell,
                                            const ParticleAccessor_T& ac)
{
   static_assert(LatticeModel_T::compressible == false, "Only works with incompressible models!");
   WALBERLA_ASSERT_NOT_NULLPTR(pdfField);
   WALBERLA_ASSERT_NOT_NULLPTR(particleAndVolumeFractionField);

   const real_t ralaxationTime =
      real_c(1) / pdfField->latticeModel().collisionModel().omega(cell.x(), cell.y(), cell.z());

   Vector3< real_t > velocity(real_t(0));
   real_t totalSolidWeightingInCell = real_t(0);

   for (auto particleFracIt = particleAndVolumeFractionField->get(cell).begin();
        particleFracIt != particleAndVolumeFractionField->get(cell).end(); ++particleFracIt)
   {
      const Vector3< real_t > coordsCellCenter = blockStorage.getBlockLocalCellCenter(block, cell);

      const real_t eps = (*particleFracIt).second;

      const real_t Bs = calculateWeighting< Weighting_T >(eps, ralaxationTime);

      const size_t idx = ac.uidToIdx(particleFracIt->first);
      WALBERLA_ASSERT_UNEQUAL(idx, ac.getInvalidIdx(), "Index of particle is invalid!");
      velocity += Bs * mesa_pd::getVelocityAtWFPoint(idx, ac, coordsCellCenter);

      totalSolidWeightingInCell += Bs;
   }

   velocity += pdfField->getVelocity(cell) * (real_c(1) - totalSolidWeightingInCell);

   return velocity;
}

/*!\brief Initializes the PDF field inside the particles according to the velocities of the particles.
 *
 * As the Partially Saturated Cells method relies on executing the LBM sweep also inside the particles, it is good
 * practice (and for some PSM variants also required) to initialize the PDF field ( i.e. the velocity ) in agreement
 * with possible initial velocities of the particles. This is also the case in the presence of external forces acting on
 * the fluid, as these will often shift the macroscopic velocities during the initialization of the PDF field.
 *
 * Note, that the ParticleAndVolumeFractionMapping for PSM has to be called first to have a valid field.
 *
 * Only the velocity of cells intersecting with particles is set, pure fluid cells remain unchanged.
 */
template< typename LatticeModel_T, int Weighting_T, typename ParticleAccessor_T >
void initializeDomainForPSM(StructuredBlockStorage& blockStorage, const BlockDataID& pdfFieldID,
                            const BlockDataID& particleAndVolumeFractionFieldID, const ParticleAccessor_T& ac)
{
   using PdfField_T = lbm::PdfField< LatticeModel_T >;

   // iterate all blocks with an iterator 'block'
   for (auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt)
   {
      // get the field data out of the block
      PdfField_T* pdfField = blockIt->getData< PdfField_T >(pdfFieldID);
      ParticleAndVolumeFractionField_T* particleAndVolumeFractionField =
         blockIt->getData< ParticleAndVolumeFractionField_T >(particleAndVolumeFractionFieldID);

      auto xyzFieldSize = pdfField->xyzSize();
      for (auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt)
      {
         Cell cell(*fieldIt);

         const real_t ralaxationTime =
            real_c(1) / pdfField->latticeModel().collisionModel().omega(cell.x(), cell.y(), cell.z());

         Vector3< real_t > weightedAverageParticleVelocityInCell(real_t(0));
         real_t totalSolidWeightingInCell = real_t(0);

         for (auto particleFracIt = particleAndVolumeFractionField->get(cell).begin();
              particleFracIt != particleAndVolumeFractionField->get(cell).end(); ++particleFracIt)
         {
            const Vector3< real_t > coordsCellCenter = blockStorage.getBlockLocalCellCenter(*blockIt, cell);

            const real_t eps = (*particleFracIt).second;

            const real_t Bs = calculateWeighting< Weighting_T >(eps, ralaxationTime);

            const size_t idx = ac.uidToIdx(particleFracIt->first);
            WALBERLA_ASSERT_UNEQUAL(idx, ac.getInvalidIdx(), "Index of particle is invalid!");
            weightedAverageParticleVelocityInCell += Bs * mesa_pd::getVelocityAtWFPoint(idx, ac, coordsCellCenter);

            totalSolidWeightingInCell += Bs;
         }

         if (totalSolidWeightingInCell > real_t(0))
         {
            Vector3< real_t > fluidVelocityInCell(real_t(0));
            const real_t rho = pdfField->getDensityAndVelocity(fluidVelocityInCell, cell);

            // set the PDFs to equilibrium with the density rho and the average velocity of all intersecting particles
            pdfField->setToEquilibrium(cell,
                                       fluidVelocityInCell * (real_c(1) - totalSolidWeightingInCell) +
                                          weightedAverageParticleVelocityInCell,
                                       rho);
         }
      }
   }
}

} // namespace psm
} // namespace lbm_mesapd_coupling
} // namespace walberla
