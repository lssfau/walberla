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
//! \file EquilibriumDistribution.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include "stencil/Directions.h"

#include <type_traits>


namespace walberla {
namespace lbm {


//////////////////////////////////
// get equilibrium distribution //
//////////////////////////////////



namespace internal {

inline real_t multiplyVelocityDirection( const real_t cx, const real_t cy, const real_t cz, const Vector3< real_t > & velocity )
{
   return cx * velocity[0] + cy * velocity[1] + cz * velocity[2];
}

inline real_t multiplyVelocityDirection( const stencil::Direction & direction, const Vector3< real_t > & velocity )
{
   using namespace stencil;
   return multiplyVelocityDirection( real_c(cx[direction]), real_c(cy[direction]), real_c(cz[direction]), velocity );
}

} // namespace internal



template< typename LatticeModel_T, class Enable = void >
class EquilibriumDistribution;



template< typename LatticeModel_T >
class EquilibriumDistribution< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
	                                                                     LatticeModel_T::equilibriumAccuracyOrder == 2
	                                                                     >::type >
{
public:

   static_assert( LatticeModel_T::compressible == false,         "Only works with incompressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef typename LatticeModel_T::Stencil Stencil;

   static real_t get( const real_t cx, const real_t cy, const real_t cz, const real_t w,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t vel = internal::multiplyVelocityDirection( cx, cy, cz, velocity );
      return w * ( (rho - real_t(1.0)) - real_t(1.5) * velocity.sqrLength() + real_t(4.5) * vel * vel + real_t(3.0) * vel );
   }

   static real_t get( const stencil::Direction direction,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;
      return get( real_c(cx[direction]), real_c(cy[direction]), real_c(cz[direction]), LatticeModel_T::w[ Stencil::idx[direction] ], velocity, rho );
   }

   static real_t getSymmetricPart( const stencil::Direction direction,
                                   const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t vel = internal::multiplyVelocityDirection( direction, velocity );
      return LatticeModel_T::w[ Stencil::idx[direction] ] * ( (rho - real_t(1.0)) - real_t(1.5) * velocity.sqrLength() + real_t(4.5) * vel * vel );
   }

   static real_t getAsymmetricPart( const stencil::Direction direction,
                                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t /*rho*/ = real_t(1.0) )
   {
      return LatticeModel_T::w[ Stencil::idx[direction] ] * real_t(3.0) * internal::multiplyVelocityDirection( direction, velocity );
   }

   static std::vector< real_t > get( const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > equilibrium( Stencil::Size );

      const real_t dirIndependent = (rho - real_t(1.0)) - real_t(1.5) * velocity.sqrLength();
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const real_t vel = internal::multiplyVelocityDirection( *d, velocity );
         equilibrium[d.toIdx()] = LatticeModel_T::w[ d.toIdx() ] * ( dirIndependent + real_t(4.5) * vel * vel + real_t(3.0) * vel );
      }

      return equilibrium;
   }
};



template< typename LatticeModel_T >
class EquilibriumDistribution< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                        LatticeModel_T::equilibriumAccuracyOrder == 1
                                                                        >::type >
{
public:

   static_assert( LatticeModel_T::compressible == false,         "Only works with incompressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 1, "Only works for lattice models that require the equilibrium distribution to be order 1 accurate!" );

   typedef typename LatticeModel_T::Stencil Stencil;

   static real_t get( const real_t cx, const real_t cy, const real_t cz, const real_t w,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      return w * ( (rho - real_t(1.0)) + real_t(3.0) * internal::multiplyVelocityDirection( cx, cy, cz, velocity ) );
   }

   static real_t get( const stencil::Direction direction,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;
      return get( real_c(cx[direction]), real_c(cy[direction]), real_c(cz[direction]), LatticeModel_T::w[ Stencil::idx[direction] ], velocity, rho );
   }

   static real_t getSymmetricPart( const stencil::Direction direction,
                                   const Vector3< real_t > & /*velocity*/ = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      return LatticeModel_T::w[ Stencil::idx[direction] ] * ( rho - real_t(1.0) );
   }

   static real_t getAsymmetricPart( const stencil::Direction direction,
                                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t /*rho*/ = real_t(1.0) )
   {
      return LatticeModel_T::w[ Stencil::idx[direction] ] * real_t(3.0) * internal::multiplyVelocityDirection( direction, velocity );
   }

   static std::vector< real_t > get( const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > equilibrium( Stencil::Size );

      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         equilibrium[d.toIdx()] = LatticeModel_T::w[ d.toIdx() ] * ( (rho - real_t(1.0)) + real_t(3.0) * internal::multiplyVelocityDirection( *d, velocity ) );

      return equilibrium;
   }
};



template< typename LatticeModel_T >
class EquilibriumDistribution< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                        LatticeModel_T::equilibriumAccuracyOrder == 2 >::type >
{
public:

   static_assert( LatticeModel_T::compressible,                  "Only works with compressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 2, "Only works for lattice models that require the equilibrium distribution to be order 2 accurate!" );

   typedef typename LatticeModel_T::Stencil Stencil;

   static real_t get( const real_t cx, const real_t cy, const real_t cz, const real_t w,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t vel = internal::multiplyVelocityDirection( cx, cy, cz, velocity );
      return w * rho * ( real_t(1.0) - real_t(1.5) * velocity.sqrLength() + real_t(4.5) * vel * vel + real_t(3.0) * vel );
   }

   static real_t get( const stencil::Direction direction,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;
      return get( real_c(cx[direction]), real_c(cy[direction]), real_c(cz[direction]), LatticeModel_T::w[ Stencil::idx[direction] ], velocity, rho );
   }

   static real_t getSymmetricPart( const stencil::Direction direction,
                                   const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      const real_t vel = internal::multiplyVelocityDirection( direction, velocity );
      return LatticeModel_T::w[ Stencil::idx[direction] ] * rho * ( real_t(1.0) - real_t(1.5) * velocity.sqrLength() + real_t(4.5) * vel * vel );
   }

   static real_t getAsymmetricPart( const stencil::Direction direction,
                                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      return LatticeModel_T::w[ Stencil::idx[direction] ] * rho * real_t(3.0) * internal::multiplyVelocityDirection( direction, velocity );
   }

   static std::vector< real_t > get( const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > equilibrium( Stencil::Size );

      const real_t dirIndependent = real_t(1.0) - real_t(1.5) * velocity.sqrLength();
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         const real_t vel = internal::multiplyVelocityDirection( *d, velocity );
         equilibrium[d.toIdx()] = LatticeModel_T::w[ d.toIdx() ] * rho * ( dirIndependent + real_t(4.5) * vel * vel + real_t(3.0) * vel );
      }

      return equilibrium;
   }
};



template< typename LatticeModel_T >
class EquilibriumDistribution< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                        LatticeModel_T::equilibriumAccuracyOrder == 1 >::type >
{
public:

   static_assert( LatticeModel_T::compressible,                  "Only works with compressible models!" );
   static_assert( LatticeModel_T::equilibriumAccuracyOrder == 1, "Only works for lattice models that require the equilibrium distribution to be order 1 accurate!" );

   typedef typename LatticeModel_T::Stencil Stencil;

   static real_t get( const real_t cx, const real_t cy, const real_t cz, const real_t w,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      return w * rho * ( real_t(1.0) + real_t(3.0) * internal::multiplyVelocityDirection( cx, cy, cz, velocity ) );
   }

   static real_t get( const stencil::Direction direction,
                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      using namespace stencil;
      return get( real_c(cx[direction]), real_c(cy[direction]), real_c(cz[direction]), LatticeModel_T::w[ Stencil::idx[direction] ], velocity, rho );
   }

   static real_t getSymmetricPart( const stencil::Direction direction,
                                   const Vector3< real_t > & /*velocity*/ = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      return LatticeModel_T::w[ Stencil::idx[direction] ] * rho;
   }

   static real_t getAsymmetricPart( const stencil::Direction direction,
                                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      return LatticeModel_T::w[ Stencil::idx[direction] ] * rho * real_t(3.0) * internal::multiplyVelocityDirection( direction, velocity );
   }

   static std::vector< real_t > get( const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      std::vector< real_t > equilibrium( Stencil::Size );

      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
         equilibrium[d.toIdx()] = LatticeModel_T::w[ d.toIdx() ] * rho * ( real_t(1.0) + real_t(3.0) * internal::multiplyVelocityDirection( *d, velocity ) );

      return equilibrium;
   }
};



} // namespace lbm
} // namespace walberla
