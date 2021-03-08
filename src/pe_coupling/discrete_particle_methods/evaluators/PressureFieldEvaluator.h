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
//! \file PressureFieldEvaluator.h
//! \ingroup pe_coupling
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/GhostLayerField.h"

namespace walberla {
namespace pe_coupling {
namespace discrete_particle_methods {

using math::Vector3;

/*!\brief Evaluator of the pressure field, given a PDF field.
 *
 * Uses p = rho / 3 = sum_q f_q / 3, to compute the pressure field from a given PDF field.
 */
template <typename LatticeModel_T, typename BoundaryHandling_T>
class PressureFieldEvaluator
{
public:
   using PdfField_T = lbm::PdfField<LatticeModel_T>;
   using ScalarField_T = GhostLayerField<real_t, 1>;

   PressureFieldEvaluator( const BlockDataID & pressureFieldID, const ConstBlockDataID & pdfFieldID, const ConstBlockDataID & boundaryHandlingID )
      : pressureFieldID_( pressureFieldID ), pdfFieldID_( pdfFieldID ), boundaryHandlingID_( boundaryHandlingID )
   { }

   void operator()(IBlock * const block)
   {
      const PdfField_T* pdfField                   = block->getData< PdfField_T >( pdfFieldID_ );
      ScalarField_T* pressureField                 = block->getData< ScalarField_T >( pressureFieldID_ );
      const BoundaryHandling_T * boundaryHandling  = block->getData< BoundaryHandling_T >( boundaryHandlingID_ );

      const real_t c_s_sqr = real_t(1)/real_t(3);
      WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,
         if( boundaryHandling->isDomain(x,y,z) )
         {
            real_t density = pdfField->getDensity(x,y,z);
            WALBERLA_ASSERT( !math::isnan(density), "NaN found when evaluating density in cell " << Cell(x,y,z) );
            pressureField->get(x,y,z) = c_s_sqr * density;
         }
      );
   }

private:
   const BlockDataID pressureFieldID_;
   const ConstBlockDataID pdfFieldID_;
   const ConstBlockDataID boundaryHandlingID_;
};


} // namespace discrete_particle_methods
} // namespace pe_coupling
} // namespace walberla
