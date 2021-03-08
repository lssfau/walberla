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
//! \file ComponentExtractionAdaptor.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "GhostLayerFieldAdaptor.h"


namespace walberla {
namespace field {


   template< typename Field_T, uint_t fComponent, uint_t vectorComponent >
   class ComponentExtractionFunction
   {
   public:
      using basefield_t = Field_T;
      using basefield_iterator = typename Field_T::const_base_iterator;

      using OutputTrait = VectorTrait<typename basefield_t::value_type>;

      using value_type = typename OutputTrait::OutputType;

      static const uint_t F_SIZE = 1;


      static_assert( fComponent < basefield_t::F_SIZE, "Invalid fComponent"  );
      static_assert( vectorComponent < OutputTrait::F_SIZE, "Invalid vectorComponent" );

      value_type operator() ( const basefield_t & baseField,
                              cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t  ) const
      {
         return OutputTrait::get( baseField.get(x,y,z,fComponent), vectorComponent ) ;
      }

      value_type operator() ( const basefield_iterator & it ) const {
         return OutputTrait::get(  it.getF( fComponent ), vectorComponent );
      }
   };


   template<typename field_t, uint_t component, uint_t vectorComponent=0>
   class ComponentExtractionAdaptor : public GhostLayerFieldAdaptor< ComponentExtractionFunction<field_t,component,vectorComponent>, 0 >
   {
   public:
      using Func = ComponentExtractionFunction<field_t, component, vectorComponent>;
      using baseclass = GhostLayerFieldAdaptor<Func, 0>;

      ComponentExtractionAdaptor( const field_t & field, const Func & func = Func() )
         : baseclass( field, func)
      {}
   };



} // namespace field
} // namespace walberla


