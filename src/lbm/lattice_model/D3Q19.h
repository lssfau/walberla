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
//! \file D3Q19.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LatticeModelBase.h"
#include "stencil/D3Q19.h"

#include <type_traits>


namespace walberla {
namespace lbm {



template< typename CollisionModel_T, bool Compressible = false, typename ForceModel_T = force_model::None, int EquilibriumAccuracyOrder = 2 >
class D3Q19 : public LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >
{
public:

   using CollisionModel = typename LatticeModelBase<CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder>::CollisionModel;
   using ForceModel = typename LatticeModelBase<CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder>::ForceModel;

   using Stencil = stencil::D3Q19;
   using CommunicationStencil = stencil::D3Q19;

   static const char * NAME;

   static const real_t w_0;
   static const real_t w_1;
   static const real_t w_2;
   static const real_t w[19];    // Stencil::Size !
   static const real_t wInv[19]; // Stencil::Size !

   D3Q19( const CollisionModel_T & cm, const ForceModel_T & fm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, fm ) {}

   // available only if the force model == force_model::None
   D3Q19( const CollisionModel_T & cm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, force_model::None() )
   {
      static_assert( (std::is_same< ForceModel_T, force_model::None >::value), "This constructor is only available if the force model is equal to force_model::None!" );
   }

   ~D3Q19() override = default;

protected:

   void config( IBlock & /*block*/, StructuredBlockStorage & /*sbs*/ ) override {}
};

template< typename CM, bool C, typename FM, int EAO > const char*  D3Q19<CM,C,FM,EAO>::NAME = "D3Q19";


template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w_0 = real_t(1.0) / real_t( 3.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w_1 = real_t(1.0) / real_t(18.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w_2 = real_t(1.0) / real_t(36.0);

// must match with the static array 'dir' in stencil::D3Q19
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::w[19] = { real_t(1.0) / real_t( 3.0),   // C
                                                                                                 real_t(1.0) / real_t(18.0),   // N
                                                                                                 real_t(1.0) / real_t(18.0),   // S
                                                                                                 real_t(1.0) / real_t(18.0),   // W
                                                                                                 real_t(1.0) / real_t(18.0),   // E
                                                                                                 real_t(1.0) / real_t(18.0),   // T
                                                                                                 real_t(1.0) / real_t(18.0),   // B
                                                                                                 real_t(1.0) / real_t(36.0),   // NW
                                                                                                 real_t(1.0) / real_t(36.0),   // NE
                                                                                                 real_t(1.0) / real_t(36.0),   // SW
                                                                                                 real_t(1.0) / real_t(36.0),   // SE
                                                                                                 real_t(1.0) / real_t(36.0),   // TN
                                                                                                 real_t(1.0) / real_t(36.0),   // TS
                                                                                                 real_t(1.0) / real_t(36.0),   // TW
                                                                                                 real_t(1.0) / real_t(36.0),   // TE
                                                                                                 real_t(1.0) / real_t(36.0),   // BN
                                                                                                 real_t(1.0) / real_t(36.0),   // BS
                                                                                                 real_t(1.0) / real_t(36.0),   // BW
                                                                                                 real_t(1.0) / real_t(36.0) }; // BE

// must match with the static array 'dir' in stencil::D3Q19
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q19<CM,C,FM,EAO>::wInv[19] = { real_t( 3.0),   // C
                                                                                                    real_t(18.0),   // N
                                                                                                    real_t(18.0),   // S
                                                                                                    real_t(18.0),   // W
                                                                                                    real_t(18.0),   // E
                                                                                                    real_t(18.0),   // T
                                                                                                    real_t(18.0),   // B
                                                                                                    real_t(36.0),   // NW
                                                                                                    real_t(36.0),   // NE
                                                                                                    real_t(36.0),   // SW
                                                                                                    real_t(36.0),   // SE
                                                                                                    real_t(36.0),   // TN
                                                                                                    real_t(36.0),   // TS
                                                                                                    real_t(36.0),   // TW
                                                                                                    real_t(36.0),   // TE
                                                                                                    real_t(36.0),   // BN
                                                                                                    real_t(36.0),   // BS
                                                                                                    real_t(36.0),   // BW
                                                                                                    real_t(36.0) }; // BE



} // namespace lbm
} // namespace walberla
