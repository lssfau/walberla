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
//! \file VelocityFieldEvaluator.h
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

/*!\brief Evaluator of the velocity field, given a PDF field.
 *
 * Calculates the velocity in each domain cell, given a PDF field.
 * Since in LBM, the macroscopic velocity depends on the local fluid forcing, one has to pay attention to the currently
 * set forces in the corresponding force field before evaluating the velocity.
 *
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class VelocityFieldEvaluator
{
public:
   using PdfField_T = lbm::PdfField<LatticeModel_T>;
   using VelocityField_T = GhostLayerField<Vector3<real_t>, 1>;

   VelocityFieldEvaluator( const BlockDataID & velocityFieldID,
                           const ConstBlockDataID & pdfFieldID, const ConstBlockDataID & boundaryHandlingID )
      : velocityFieldID_( velocityFieldID ), pdfFieldID_( pdfFieldID ), boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()(IBlock * const block)
   {
      const PdfField_T* pdfField                   = block->getData< PdfField_T >( pdfFieldID_ );
      VelocityField_T* velocityField         = block->getData< VelocityField_T >( velocityFieldID_ );
      const BoundaryHandling_T * boundaryHandling  = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            WALBERLA_ASSERT( !math::isnan(pdfField->getVelocity(x,y,z)[0]) &&
                             !math::isnan(pdfField->getVelocity(x,y,z)[1]) &&
                             !math::isnan(pdfField->getVelocity(x,y,z)[2]), "NaN found when evaluating velocity in cell " << Cell(x,y,z) );
            velocityField->get(x,y,z) = pdfField->getVelocity(x,y,z);
         }
      );
   }

private:
   const BlockDataID velocityFieldID_;
   const ConstBlockDataID pdfFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
