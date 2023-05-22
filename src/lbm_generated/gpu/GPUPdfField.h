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
//! \file GPUPdfField.h
//! \ingroup lbm_generated
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "gpu/GPUField.h"

using namespace walberla::gpu;

namespace walberla::lbm_generated {

template< typename LatticeStorageSpecification_T >
class GPUPdfField : public GPUField< real_t >
{
 public:

   //** Type Definitions  **********************************************************************************************
   /*! \name Type Definitions */
   //@{
   using LatticeStorageSpecification = LatticeStorageSpecification_T;
   using Stencil = typename LatticeStorageSpecification_T::Stencil;

   using value_type = typename GPUField<real_t>::value_type;
   //@}
   //*******************************************************************************************************************

   GPUPdfField( uint_t _xSize, uint_t _ySize, uint_t _zSize,
               const LatticeStorageSpecification_T & storageSpecification,
               uint_t _nrOfGhostLayers, const Layout & _layout = zyxf, bool usePitchedMem = true );


   ~GPUPdfField() = default;

 protected:
   LatticeStorageSpecification_T storageSpecification_;
};



template< typename LatticeStorageSpecification_T >
GPUPdfField< LatticeStorageSpecification_T >::GPUPdfField( uint_t _xSize, uint_t _ySize, uint_t _zSize,
                                                          const LatticeStorageSpecification_T & storageSpecification,
                                                          uint_t ghostLayers, const Layout & layout, bool usePitchedMem) :
                    GPUField< real_t>( _xSize, _ySize, _zSize, LatticeStorageSpecification_T::Stencil::Size, ghostLayers, layout, usePitchedMem ), storageSpecification_( storageSpecification )
{
}

} // namespace lbm