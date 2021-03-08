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
//! \file Adaptors.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "MacroscopicValueCalculation.h"
#include "PdfField.h"
#include "core/math/Vector3.h"
#include "field/adaptors/GhostLayerFieldAdaptor.h"
#include "field/iterators/FieldNeighborPointer.h"


namespace walberla {
namespace lbm {



template< typename LatticeModel_T >
class DensityAdaptionFunction
{
public:
   using basefield_t = PdfField<LatticeModel_T>;
   using basefield_iterator = typename basefield_t::const_base_iterator;
   using value_type = real_t;

   static const uint_t F_SIZE = 1u;

   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t /*f*/ = 0 ) const
   {
      return baseField.getDensity(x,y,z);
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      const basefield_t * baseFieldPtr = dynamic_cast<const basefield_t *>( it.getField() );
      return lbm::getDensity( baseFieldPtr->latticeModel(), it );
   }
};



template< typename LatticeModel_T >
class VelocityVectorAdaptionFunction
{
public:
   using basefield_t = PdfField<LatticeModel_T>;
   using basefield_iterator = typename basefield_t::const_base_iterator;
   using value_type = Vector3<real_t>;

   static const uint_t F_SIZE = 1u;

   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t /*f*/ = 0 ) const
   {
      return baseField.getVelocity(x,y,z);
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      const basefield_t * baseFieldPtr = dynamic_cast<const basefield_t *>( it.getField() );
      return lbm::getVelocity( baseFieldPtr->latticeModel(), it );
   }
};


template< typename LatticeModel_T, bool stream = false >
class MomentumDensityVectorAdaptionFunction
{
public:
   using basefield_t = PdfField<LatticeModel_T>;
   using basefield_iterator = typename basefield_t::const_base_iterator;
   using value_type = Vector3<real_t>;

   using Stencil = typename LatticeModel_T::Stencil;

   static const uint_t F_SIZE = 1u;

   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t /*f*/ = 0 ) const
   {
      if( stream )
         return lbm::getMomentumDensity<LatticeModel_T>( baseField.latticeModel(), field::FieldNeighborPointer< basefield_t, const basefield_t, const real_t >( baseField, x, y, z ) );
      else
         return baseField.getMomentumDensity(x,y,z);
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      return operator()( static_cast<const basefield_t&>(*it.getField()), it.x(), it.y(), it.z(), it.f() );
   }
};


template< typename LatticeModel_T >
class VelocityAdaptionFunction
{
public:
   using basefield_t = PdfField<LatticeModel_T>;
   using basefield_iterator = typename basefield_t::const_base_iterator;
   using value_type = real_t;

   static const uint_t F_SIZE = 3u;

   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f = 0 ) const
   {
      WALBERLA_ASSERT_LESS( f, 3 );
      Vector3<real_t> res = baseField.getVelocity(x,y,z);
      return res[ uint_c(f) ];
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      WALBERLA_ASSERT_LESS( it.f(), 3 );
      const basefield_t * baseFieldPtr = dynamic_cast<const basefield_t *>( it.getField() );
      Vector3<real_t> res = lbm::getVelocity( baseFieldPtr->latticeModel(), it );
      return res[ it.f() ];
   }
};


template< typename LatticeModel_T, bool stream = false >
class MomentumDensityAdaptionFunction
{
public:
   using basefield_t = PdfField<LatticeModel_T>;
   using basefield_iterator = typename basefield_t::const_base_iterator;
   using value_type = real_t;

   using Stencil = typename LatticeModel_T::Stencil;

   static const uint_t F_SIZE = 3u;

   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f = 0 ) const
   {
      WALBERLA_ASSERT_LESS( f, 3 );
      Vector3<real_t> res;
      if( stream )
         res = lbm::getMomentumDensity<LatticeModel_T>( baseField.latticeModel(), field::FieldNeighborPointer< basefield_t, basefield_t, real_t >( baseField, x, y, z ) );
      else
         res = baseField.getMomentumDensity(x,y,z);
      return res[ uint_c(f) ];
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      return operator()( static_cast<const basefield_t&>(*it.getField()), it.x(), it.y(), it.z(), it.f() );
   }
};


template< typename LatticeModel_T >
class ShearRateAdaptionFunction
{
public:
   using basefield_t = PdfField<LatticeModel_T>;
   using basefield_iterator = typename basefield_t::const_base_iterator;
   using value_type = real_t;

   static const uint_t F_SIZE = 1u;

   value_type operator() ( const basefield_t & baseField,
                           cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t /*f*/ = 0 ) const
   {
      return baseField.getShearRate(x,y,z);
   }

   value_type operator() ( const basefield_iterator & it ) const
   {
      Vector3<real_t> velocity;
      auto baseFieldPtr = it->getField();
      real_t rho = getDensityAndVelocity( velocity, baseFieldPtr->latticeModel(), it );
      return ShearRate< LatticeModel_T >::get( baseFieldPtr->latticeModel(), it, velocity, rho );
   }
};


// Struct is a workaround for missing templated typedefs
template< typename LatticeModel_T >
struct Adaptor
{
   using Density = field::GhostLayerFieldAdaptor<DensityAdaptionFunction<LatticeModel_T>, 0>;
   using VelocityVector = field::GhostLayerFieldAdaptor<VelocityVectorAdaptionFunction<LatticeModel_T>, 0>;
   using Velocity = field::GhostLayerFieldAdaptor<VelocityAdaptionFunction<LatticeModel_T>, 0>;
   using ShearRate = field::GhostLayerFieldAdaptor<ShearRateAdaptionFunction<LatticeModel_T>, 0>;

   using StreamMomentumDensityVector = field::GhostLayerFieldAdaptor<MomentumDensityVectorAdaptionFunction<LatticeModel_T, true>, 1>;
   using StreamMomentumDensity = field::GhostLayerFieldAdaptor<MomentumDensityAdaptionFunction<LatticeModel_T, true>, 1>;
   using CollideMomentumDensityVector = field::GhostLayerFieldAdaptor<MomentumDensityVectorAdaptionFunction<LatticeModel_T, false>, 0>;
   using CollideMomentumDensity = field::GhostLayerFieldAdaptor<MomentumDensityAdaptionFunction<LatticeModel_T, false>, 0>;
};



} // namespace lbm
} // namespace walberla
