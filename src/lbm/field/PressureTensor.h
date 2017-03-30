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
//! \file PressureTensor.h
//! \ingroup lbm
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Matrix3.h"

// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {

namespace internal {

template< typename LatticeModel_T, typename FieldPtrOrIterator >
void getPressureTensor( Matrix3< real_t > & pressureTensor, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   /* This function calculates the pressure tensor, also known as momentum flux tensor. Note that this is
      not identical to the traceless viscous stress tensor. They are connected as follows:
         PressureTensor == Pressure * IdentityMatrix + Density * tensorProduct(Velocity,Velocity) - StressTensor
      The (non-traceless) stress tensor in the Navier-Stokes equation is
         Pressure * IdentityMatrix - StressTensor == PressureTensor - Density * tensorProduct(Velocity,Velocity)
      For details, see sections 2.2.4 and 2.3.1 in Ulf Schiller's dissertation.
      Also note that the pressure is only defined up to a constant offset. Physically however, only pressure
      differences are relevant anyway. */

   auto d = LatticeModel_T::Stencil::begin();

   const auto & firstPdf = it[ d.toIdx() ];
   auto c = Vector3<real_t>(d.cx(), d.cy(), d.cz());

   pressureTensor = tensorProduct(c,c) * firstPdf;

   ++d;

   while( d != LatticeModel_T::Stencil::end() )
   {
      const auto & pdfValue = it[ d.toIdx() ];
      c = Vector3<real_t>(d.cx(), d.cy(), d.cz());

      pressureTensor += tensorProduct(c,c) * pdfValue;

      ++d;
   }

   /* Nonequilibrium components of the pressure tensor need to be scaled with (1-1/2omega), see eq. 16 in
      appendix of Hou et al., Simulation of Cavity Flow by the Lattice-Boltzmann Method, 1995.
      The factor below is correct for SRT. An exact expression for TRT is not known, but the approximation
      is still good. For MRT, this approximation can usually also be used, but an exact expression in mode
      space is available. */
   real_t fac = latticeModel.collisionModel().viscosity(it.x(),it.y(),it.z()) / ( latticeModel.collisionModel().viscosity(it.x(),it.y(),it.z()) + real_t(1.0/6.0) );
   pressureTensor[1] *= fac; // xy
   pressureTensor[2] *= fac; // xz
   pressureTensor[3] *= fac; // yx
   pressureTensor[5] *= fac; // yz
   pressureTensor[6] *= fac; // zx
   pressureTensor[7] *= fac; // zy
}

template< typename LatticeModel_T, typename PdfField_T >
void getPressureTensor( Matrix3< real_t > & pressureTensor, const LatticeModel_T & latticeModel, const PdfField_T & pdf,
                         const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
{
   const auto & xyz0 = pdf(x,y,z,0);

   auto d = LatticeModel_T::Stencil::begin();

   const auto & firstPdf = pdf.getF( &xyz0, d.toIdx() );
   auto c = Vector3<real_t>(d.cx(), d.cy(), d.cz());

   pressureTensor = tensorProduct(c,c) * firstPdf;

   ++d;

   while( d != LatticeModel_T::Stencil::end() )
   {
      const auto & pdfValue = pdf.getF( &xyz0, d.toIdx() );
      c = Vector3<real_t>(d.cx(), d.cy(), d.cz());

      pressureTensor += tensorProduct(c,c) * pdfValue;

      ++d;
   }

   real_t fac = latticeModel.collisionModel().viscosity(x,y,z) / ( latticeModel.collisionModel().viscosity(x,y,z) + real_t(1.0/6.0) );
   pressureTensor[1] *= fac; // xy
   pressureTensor[2] *= fac; // xz
   pressureTensor[3] *= fac; // yx
   pressureTensor[5] *= fac; // yz
   pressureTensor[6] *= fac; // zx
   pressureTensor[7] *= fac; // zy
}

} // namespace internal


template< typename LatticeModel_T >
struct PressureTensor
{
   template< typename FieldPtrOrIterator >
   static void get( Matrix3< real_t > & pressureTensor, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
   {
      internal::getPressureTensor<LatticeModel_T, FieldPtrOrIterator>( pressureTensor, latticeModel, it );
   }

   template< typename PdfField_T >
   static void get( Matrix3< real_t > & pressureTensor, const LatticeModel_T & latticeModel, const PdfField_T & pdf,
                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      internal::getPressureTensor<LatticeModel_T, PdfField_T>( pressureTensor, latticeModel, pdf, x, y, z );
   }
};

} // namespace lbm
} // namespace walberla
