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
//! \file GNSVelocityFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Evaluator for the fluid-phase velocity, given a (GNS) PDf field.
 *
 * Evaluator for the fluid-phase velocity when using GNS-LBM (in contrast to the volume-averaged velocity that is
 * obtained by calculating the first order moment of the PDFs).
 *
 * See: Z. Guo, T. S. Zhao - "Lattice Boltzmann model for incompressible flows through porous media",
 * Phys. Rev. E 66 (2002)036304. doi:10.1103/PhysRevE.66.036304.
 *
 * Since in LBM, the macroscopic velocity depends on the local fluid forcing, one has to pay attention to the currently
 * set forces in the corresponding force field before evaluating the velocity.
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class GNSVelocityFieldEvaluator
{
public:
   using PdfField_T = lbm::PdfField<LatticeModel_T>;
   using VelocityField_T = GhostLayerField<Vector3<real_t>, 1>;
   using ScalarField_T = GhostLayerField<real_t, 1>;

   GNSVelocityFieldEvaluator( const BlockDataID & velocityFieldID, const ConstBlockDataID & pdfFieldID, const ConstBlockDataID & solidVolumeFractionFieldID,
                              const ConstBlockDataID & boundaryHandlingID )
      : velocityFieldID_( velocityFieldID ), pdfFieldID_( pdfFieldID ), solidVolumeFractionFieldID_( solidVolumeFractionFieldID ),
        boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()(IBlock * const block)
   {
      const PdfField_T* pdfField                  = block->getData< PdfField_T >( pdfFieldID_ );
      VelocityField_T* velocityField              = block->getData< VelocityField_T >( velocityFieldID_ );
      const ScalarField_T* svfField               = block->getData< ScalarField_T >( solidVolumeFractionFieldID_ );
      const BoundaryHandling_T* boundaryHandling  = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            real_t fluidVolumeFraction = real_t(1) - svfField->get(x,y,z);
            Vector3<real_t> volumeAverageVelocity = pdfField->getVelocity(x,y,z);
            velocityField->get(x,y,z) = volumeAverageVelocity / fluidVolumeFraction;
         }
      );
   }

private:
   const BlockDataID velocityFieldID_;
   const ConstBlockDataID pdfFieldID_;
   const ConstBlockDataID solidVolumeFractionFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
