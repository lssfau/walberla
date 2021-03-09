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
//! \file Destroyer.h
//! \ingroup pe_coupling
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================


#pragma once

#include <limits>


namespace walberla {
namespace pe_coupling {

template< typename LatticeModel_T >
class NaNDestroyer
{
public:

   using PdfField_T = lbm::PdfField<LatticeModel_T>;

   void operator()( const cell_idx_t & x, const cell_idx_t & y, const cell_idx_t & z, IBlock * const /*block*/, PdfField_T * const pdfField ) {
      for (auto d = uint_t(0); d < LatticeModel_T::Stencil::Size; ++d)
         pdfField->get(x, y, z, d) = std::numeric_limits<typename PdfField_T::value_type>::quiet_NaN();
   }
};

} // namespace pe_coupling
} // namespace walberla
