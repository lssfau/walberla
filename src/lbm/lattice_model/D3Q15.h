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
//! \file D3Q15.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "LatticeModelBase.h"
#include "stencil/D3Q15.h"
#include "stencil/D3Q27.h"

#include <type_traits>


namespace walberla {
namespace lbm {



template< typename CollisionModel_T, bool Compressible = false, typename ForceModel_T = force_model::None, int EquilibriumAccuracyOrder = 2 >
class D3Q15 : public LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >
{
public:

   typedef typename LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >::CollisionModel  CollisionModel;
   typedef typename LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >::ForceModel      ForceModel;

   typedef stencil::D3Q15 Stencil;
   typedef stencil::D3Q27 CommunicationStencil;

   static const char * NAME;

   static const real_t w_0;
   static const real_t w_1;
   static const real_t w_2;
   static const real_t w[15];    // Stencil::Size !
   static const real_t wInv[15]; // Stencil::Size !

   D3Q15( const CollisionModel_T & cm, const ForceModel_T & fm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, fm ) {}

   // available only if the force model == force_model::None
   D3Q15( const CollisionModel_T & cm ) :
      LatticeModelBase< CollisionModel_T, Compressible, ForceModel_T, EquilibriumAccuracyOrder >( cm, force_model::None() )
   {
      static_assert( (std::is_same< ForceModel_T, force_model::None >::value), "This constructor is only available if the force model is equal to force_model::None!" );
   }

   virtual ~D3Q15() {}

protected:

   virtual void config( IBlock & /*block*/, StructuredBlockStorage & /*sbs*/ ) {}
};

template< typename CM, bool C, typename FM, int EAO > const char*  D3Q15<CM,C,FM,EAO>::NAME = "D3Q15";

template< typename CM, bool C, typename FM, int EAO > const real_t D3Q15<CM,C,FM,EAO>::w_0 = real_t(2.0) / real_t( 9.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q15<CM,C,FM,EAO>::w_1 = real_t(1.0) / real_t( 9.0);
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q15<CM,C,FM,EAO>::w_2 = real_t(1.0) / real_t(72.0);

// must match with the static array 'dir' in stencil::D3Q15
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q15<CM,C,FM,EAO>::w[15] = { real_t(2.0) / real_t( 9.0),   // C
                                                                                                 real_t(1.0) / real_t( 9.0),   // N
                                                                                                 real_t(1.0) / real_t( 9.0),   // S
                                                                                                 real_t(1.0) / real_t( 9.0),   // W
                                                                                                 real_t(1.0) / real_t( 9.0),   // E
                                                                                                 real_t(1.0) / real_t( 9.0),   // T
                                                                                                 real_t(1.0) / real_t( 9.0),   // B
                                                                                                 real_t(1.0) / real_t(72.0),   // TNE
                                                                                                 real_t(1.0) / real_t(72.0),   // TNW
                                                                                                 real_t(1.0) / real_t(72.0),   // TSE
                                                                                                 real_t(1.0) / real_t(72.0),   // TSW
                                                                                                 real_t(1.0) / real_t(72.0),   // BNE
                                                                                                 real_t(1.0) / real_t(72.0),   // BNW
                                                                                                 real_t(1.0) / real_t(72.0),   // BSE
                                                                                                 real_t(1.0) / real_t(72.0) }; // BSW

// must match with the static array 'dir' in stencil::D3Q15
template< typename CM, bool C, typename FM, int EAO > const real_t D3Q15<CM,C,FM,EAO>::wInv[15] = { real_t( 4.5),   // C
                                                                                                    real_t( 9.0),   // N
                                                                                                    real_t( 9.0),   // S
                                                                                                    real_t( 9.0),   // W
                                                                                                    real_t( 9.0),   // E
                                                                                                    real_t( 9.0),   // T
                                                                                                    real_t( 9.0),   // B
                                                                                                    real_t(72.0),   // TNE
                                                                                                    real_t(72.0),   // TNW
                                                                                                    real_t(72.0),   // TSE
                                                                                                    real_t(72.0),   // TSW
                                                                                                    real_t(72.0),   // BNE
                                                                                                    real_t(72.0),   // BNW
                                                                                                    real_t(72.0),   // BSE
                                                                                                    real_t(72.0) }; // BSW



} // namespace lbm
} // namespace walberla
