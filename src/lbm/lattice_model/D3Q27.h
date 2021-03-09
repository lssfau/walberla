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
//! \file D3Q27.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LatticeModelBase.h"
#include "stencil/D3Q27.h"

#include <type_traits>


namespace walberla {
namespace lbm {



template< typename CollisionModel_T, bool Compressible = false, typename ForceModel_T = force_model::None, int EquilibriumAccuracyOrder = 2 >
class D3Q27 : public LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >
{
public:

   static_assert( ( ! std::is_same< CollisionModel_T, collision_model::D3Q19MRT >::value), "D3Q19MRT only works with D3Q19!" );

   using CollisionModel = typename LatticeModelBase<CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder>::CollisionModel;
   using ForceModel = typename LatticeModelBase<CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder>::ForceModel;

   using Stencil = stencil::D3Q27;
   using CommunicationStencil = stencil::D3Q27;

   static const char * NAME;

   static const real_t w_0;
   static const real_t w_1;
   static const real_t w_2;
   static const real_t w_3;
   static const real_t w[27];    // Stencil::Size !
   static const real_t wInv[27]; // Stencil::Size !

   D3Q27( const CollisionModel_T & cm, const ForceModel_T & fm  ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, fm ) {}

   // available only if the force model == force_model::None
   D3Q27( const CollisionModel_T & cm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, force_model::None() )
   {
      static_assert( (std::is_same< ForceModel_T, force_model::None >::value), "This constructor is only available if the force model is equal to force_model::None!" );
   }

   ~D3Q27() override = default;

protected:

   void config( IBlock & /*block*/, StructuredBlockStorage & /*sbs*/ ) override {}
};

template< typename CM, bool C, typename FM, int EAO > const char*  D3Q27<CM,C,FM,EAO>::NAME = "D3Q27";

template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_0 = real_t(8.0) / real_t( 27.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_1 = real_t(2.0) / real_t( 27.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_2 = real_t(1.0) / real_t( 54.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w_3 = real_t(1.0) / real_t(216.0);

// must match with the static array 'dir' in stencil::D3Q27
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::w[27] = { real_t(8.0) / real_t( 27.0),   // C
                                                                                                 real_t(2.0) / real_t( 27.0),   // N
                                                                                                 real_t(2.0) / real_t( 27.0),   // S
                                                                                                 real_t(2.0) / real_t( 27.0),   // W
                                                                                                 real_t(2.0) / real_t( 27.0),   // E
                                                                                                 real_t(2.0) / real_t( 27.0),   // T
                                                                                                 real_t(2.0) / real_t( 27.0),   // B
                                                                                                 real_t(1.0) / real_t( 54.0),   // NW
                                                                                                 real_t(1.0) / real_t( 54.0),   // NE
                                                                                                 real_t(1.0) / real_t( 54.0),   // SW
                                                                                                 real_t(1.0) / real_t( 54.0),   // SE
                                                                                                 real_t(1.0) / real_t( 54.0),   // TN
                                                                                                 real_t(1.0) / real_t( 54.0),   // TS
                                                                                                 real_t(1.0) / real_t( 54.0),   // TW
                                                                                                 real_t(1.0) / real_t( 54.0),   // TE
                                                                                                 real_t(1.0) / real_t( 54.0),   // BN
                                                                                                 real_t(1.0) / real_t( 54.0),   // BS
                                                                                                 real_t(1.0) / real_t( 54.0),   // BW
                                                                                                 real_t(1.0) / real_t( 54.0),   // BE
                                                                                                 real_t(1.0) / real_t(216.0),   // TNE
                                                                                                 real_t(1.0) / real_t(216.0),   // TNW
                                                                                                 real_t(1.0) / real_t(216.0),   // TSE
                                                                                                 real_t(1.0) / real_t(216.0),   // TSW
                                                                                                 real_t(1.0) / real_t(216.0),   // BNE
                                                                                                 real_t(1.0) / real_t(216.0),   // BNW
                                                                                                 real_t(1.0) / real_t(216.0),   // BSE
                                                                                                 real_t(1.0) / real_t(216.0) }; // BSW

// must match with the static array 'dir' in stencil::D3Q27
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q27<CM,C,FM,EAO>::wInv[27] = { real_t(27.0) / real_t( 8.0),   // C
                                                                                                    real_t(27.0) / real_t( 2.0),   // N
                                                                                                    real_t(27.0) / real_t( 2.0),   // S
                                                                                                    real_t(27.0) / real_t( 2.0),   // W
                                                                                                    real_t(27.0) / real_t( 2.0),   // E
                                                                                                    real_t(27.0) / real_t( 2.0),   // T
                                                                                                    real_t(27.0) / real_t( 2.0),   // B
                                                                                                    real_t( 54.0),                 // NW
                                                                                                    real_t( 54.0),                 // NE
                                                                                                    real_t( 54.0),                 // SW
                                                                                                    real_t( 54.0),                 // SE
                                                                                                    real_t( 54.0),                 // TN
                                                                                                    real_t( 54.0),                 // TS
                                                                                                    real_t( 54.0),                 // TW
                                                                                                    real_t( 54.0),                 // TE
                                                                                                    real_t( 54.0),                 // BN
                                                                                                    real_t( 54.0),                 // BS
                                                                                                    real_t( 54.0),                 // BW
                                                                                                    real_t( 54.0),                 // BE
                                                                                                    real_t(216.0),                 // TNE
                                                                                                    real_t(216.0),                 // TNW
                                                                                                    real_t(216.0),                 // TSE
                                                                                                    real_t(216.0),                 // TSW
                                                                                                    real_t(216.0),                 // BNE
                                                                                                    real_t(216.0),                 // BNW
                                                                                                    real_t(216.0),                 // BSE
                                                                                                    real_t(216.0) };               // BSW

} // namespace lbm
} // namespace walberla
