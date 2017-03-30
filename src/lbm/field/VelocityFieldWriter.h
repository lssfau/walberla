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
//! \file VelocityFieldWriter.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "field/iterators/IteratorMacros.h"



namespace walberla {
namespace lbm {



namespace internal {

template< typename PdfField_T, typename VelocityField_T, typename T >
struct VelocityFieldWriterHelper
{
   static void store( const PdfField_T * pdfField, VelocityField_T * velocityField,
                      const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      const Vector3<real_t> velocity = pdfField->getVelocity(x,y,z);
      velocityField->get(x,y,z,0) = velocity[0];
      velocityField->get(x,y,z,1) = velocity[1];
      velocityField->get(x,y,z,2) = velocity[2];

   }
};

template< typename PdfField_T, typename VelocityField_T >
struct VelocityFieldWriterHelper< PdfField_T, VelocityField_T, Vector3<real_t> >
{
   static void store( const PdfField_T * pdfField, VelocityField_T * velocityField,
                      const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      velocityField->get(x,y,z) = pdfField->getVelocity(x,y,z);
   }
};

}



template< typename PdfField_T, typename VelocityField_T >
class VelocityFieldWriter
{
public:

   VelocityFieldWriter( const ConstBlockDataID & pdfFieldId, const BlockDataID & velocityFieldId ) :
      pdfFieldId_( pdfFieldId ), velocityFieldId_( velocityFieldId )
   {}

   void operator()( IBlock * const block )
   {
      const PdfField_T * pdfField = block->template getData< PdfField_T >( pdfFieldId_ );
      VelocityField_T * velocityField = block->template getData< VelocityField_T >( velocityFieldId_ );

      WALBERLA_ASSERT_EQUAL( pdfField->xyzSize(), velocityField->xyzSize() );

      WALBERLA_FOR_ALL_CELLS_XYZ( pdfField,

         store( pdfField, velocityField, x, y, z );
      )
   }

private:

   void store( const PdfField_T * pdfField, VelocityField_T * velocityField, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
      internal::VelocityFieldWriterHelper< PdfField_T, VelocityField_T, typename VelocityField_T::value_type >::store( pdfField, velocityField, x, y, z );
   }

   ConstBlockDataID pdfFieldId_;
   BlockDataID velocityFieldId_;
};



} // namespace lbm
} // namespace walberla
