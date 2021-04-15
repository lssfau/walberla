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
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/Field.h"
#include "field/GhostLayerField.h"
#include "lbm/field/PdfField.h"
#include "pe/Types.h"

namespace walberla {
namespace pe_coupling {

/*!\brief Evaluates the macroscopic velocity for a given cell when using the PSM method.
 *
 * returns the cell local macroscopic velocity for the PSM method
 * cell local velocity = \f$ (1-\sum_s B_s) * u + \sum_s B_s * v_s \f$
 * u = fluid velocity
 * v_s = velocity of object s at cell center
 *
 * lbm::PdfField< LatticeModel_T > is typically typedef'ed as PdfField_T;
 * GhostLayerField< std::vector< std::pair< pe::BodyID, real_t > >, 1 > is typically typedef'ed as BodyAndVolumeFractionField_T;
 *
 * Weighting_T is like in the PSMSweep.
 */
template < typename LatticeModel_T, int Weighting_T >
Vector3<real_t> getPSMMacroscopicVelocity( const IBlock & block,
                                           lbm::PdfField< LatticeModel_T > * pdfField,
                                           GhostLayerField< std::vector< std::pair< pe::BodyID, real_t > >, 1 > * bodyAndVolumeFractionField,
                                           StructuredBlockStorage & blockStorage,
                                           const Cell & cell )
{
   static_assert( LatticeModel_T::compressible == false, "Only works with incompressible models!" );
   WALBERLA_ASSERT_NOT_NULLPTR( pdfField );
   WALBERLA_ASSERT_NOT_NULLPTR( bodyAndVolumeFractionField );

   const real_t ralaxationTime = real_c(1) / pdfField->latticeModel().collisionModel().omega( cell.x(), cell.y(), cell.z() );

   Vector3<real_t> velocity( real_t(0) );
   real_t totalSolidWeightingInCell = real_t(0);

   for( auto bodyFracIt = bodyAndVolumeFractionField->get( cell ).begin(); bodyFracIt != bodyAndVolumeFractionField->get( cell ).end(); ++bodyFracIt )
   {
      const Vector3< real_t > coordsCellCenter = blockStorage.getBlockLocalCellCenter( block, cell );

      const real_t eps = (*bodyFracIt).second;

      const real_t Bs = calculateWeighting< Weighting_T >( eps, ralaxationTime );

      velocity += Bs * (*bodyFracIt).first->velFromWF( coordsCellCenter );

      totalSolidWeightingInCell += Bs;
   }

   velocity += pdfField->getVelocity( cell ) * ( real_c(1) - totalSolidWeightingInCell );

   return velocity;
}


/*!\brief Initializes the PDF field inside the bodies according to the velocities of the bodies.
 *
 * As the Partially Saturated Cells method relies on executing the LBM sweep also inside the bodies, it is good practice (and for some PSM variants also required)
 * to initialize the PDF field ( i.e. the velocity ) in agreement with possible initial velocities of the bodies.
 * This is also the case in the presence of external forces acting on the fluid, as these will often shift the macroscopic velocities during the
 * initialization of the PDF field.
 *
 * Note, that the BodyAndVolumeFractionMapping for PSM has to be called first to have a valid field.
 *
 * Only the velocity of cells intersecting with bodies is set, pure fluid cells remain unchanged.
 */
template < typename LatticeModel_T, int Weighting_T >
void initializeDomainForPSM( StructuredBlockStorage & blockStorage,
                             const BlockDataID & pdfFieldID, const BlockDataID & bodyAndVolumeFractionFieldID )
{
   using PdfField_T = lbm::PdfField<LatticeModel_T>;
   using BodyAndVolumeFraction_T = std::pair<pe::BodyID, real_t>;
   using BodyAndVolumeFractionField_T = GhostLayerField<std::vector<BodyAndVolumeFraction_T>, 1>;

   // iterate all blocks with an iterator 'block'
   for( auto blockIt = blockStorage.begin(); blockIt != blockStorage.end(); ++blockIt )
   {
      // get the field data out of the block
      PdfField_T* pdfField = blockIt->getData< PdfField_T > ( pdfFieldID );
      BodyAndVolumeFractionField_T* bodyAndVolumeFractionField = blockIt->getData< BodyAndVolumeFractionField_T > ( bodyAndVolumeFractionFieldID );

      auto xyzFieldSize = pdfField->xyzSize();
      for( auto fieldIt = xyzFieldSize.begin(); fieldIt != xyzFieldSize.end(); ++fieldIt )
      {
         Cell cell( *fieldIt );

         const real_t ralaxationTime = real_c(1) / pdfField->latticeModel().collisionModel().omega( cell.x(), cell.y(), cell.z() );

         Vector3<real_t> weightedAverageBodyVelocityInCell( real_t(0) );
         real_t totalSolidWeightingInCell = real_t(0);

         for( auto bodyFracIt = bodyAndVolumeFractionField->get(cell).begin(); bodyFracIt != bodyAndVolumeFractionField->get(cell).end(); ++bodyFracIt )
         {
            const Vector3< real_t > coordsCellCenter = blockStorage.getBlockLocalCellCenter( *blockIt, cell );

            const real_t eps = (*bodyFracIt).second;

            const real_t Bs = calculateWeighting< Weighting_T >( eps, ralaxationTime );

            weightedAverageBodyVelocityInCell += Bs * (*bodyFracIt).first->velFromWF( coordsCellCenter );

            totalSolidWeightingInCell += Bs;
         }

         if( totalSolidWeightingInCell > real_t(0) )
         {
            Vector3< real_t > fluidVelocityInCell( real_t(0) );
            const real_t rho = pdfField->getDensityAndVelocity( fluidVelocityInCell, cell );

            // set the PDFs to equilibrium with the density rho and the average velocity of all intersecting bodies
            pdfField->setToEquilibrium( cell, fluidVelocityInCell * ( real_c(1) - totalSolidWeightingInCell ) + weightedAverageBodyVelocityInCell, rho );
         }
      }
   }
}

} // namespace pe_coupling
} // namespace walberla

