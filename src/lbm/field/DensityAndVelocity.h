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
//! \file DensityAndVelocity.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Equilibrium.h"
#include "core/DataTypes.h"
#include "core/math/Vector3.h"

#include <type_traits>


// Back-end for calculating macroscopic values
// You should never use these functions directly, always refer to the member functions
// of PdfField or the free functions that can be found in MacroscopicValueCalculation.h

namespace walberla {
namespace lbm {

namespace internal {

template< typename LatticeModel_T, class Enable = void >
struct AdaptVelocityToForce
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::internal::AdaptVelocityToForce' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};

template< typename LatticeModel_T >
struct AdaptVelocityToForce< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                      LatticeModel_T::ForceModel::constant &&
                                                                      LatticeModel_T::ForceModel::shiftMacVel
                                                                      >::type >
{
   static Vector3<real_t> get( const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t rho )
   {
      return velocity - latticeModel.forceModel().forceDensity() * real_t(0.5) / rho;
   }

   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator &, const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t rho )
   {
      return get( latticeModel, velocity, rho );
   }

   static Vector3<real_t> get( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T & latticeModel,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      return get( latticeModel, velocity, rho );
   }
};

template< typename LatticeModel_T >
struct AdaptVelocityToForce< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                      LatticeModel_T::ForceModel::constant &&
                                                                      LatticeModel_T::ForceModel::shiftMacVel
                                                                      >::type >
{
   static Vector3<real_t> get( const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t )
   {
      return velocity - latticeModel.forceModel().forceDensity() * real_t(0.5);
   }

   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator &, const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t rho )
   {
      return get( latticeModel, velocity, rho );
   }

   static Vector3<real_t> get( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T & latticeModel,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      return get( latticeModel, velocity, rho );
   }
};

template< typename LatticeModel_T >
struct AdaptVelocityToForce< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                      ! LatticeModel_T::ForceModel::constant &&
                                                                      LatticeModel_T::ForceModel::shiftMacVel
                                                                      >::type >
{
   /*
   static Vector3<real_t> get( const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t rho )
   {
      return velocity - latticeModel.forceModel().forceDensity() * real_t(0.5) / rho;
   }
   */

   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t rho )
   {
      return velocity - latticeModel.forceModel().forceDensity(it.x(),it.y(),it.z()) * real_t(0.5) / rho;
   }

   static Vector3<real_t> get( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & latticeModel,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      return velocity - latticeModel.forceModel().forceDensity(x,y,z) * real_t(0.5) / rho;
   }
};

template< typename LatticeModel_T >
struct AdaptVelocityToForce< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                      ! LatticeModel_T::ForceModel::constant &&
                                                                      LatticeModel_T::ForceModel::shiftMacVel
                                                                      >::type >
{
   /*
   static Vector3<real_t> get( const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t )
   {
      return velocity - latticeModel.forceModel().forceDensity() * real_t(0.5);
   }
   */

   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t )
   {
      return velocity - latticeModel.forceModel().forceDensity(it.x(),it.y(),it.z()) * real_t(0.5);
   }

   static Vector3<real_t> get( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & latticeModel,
                               const Vector3< real_t > & velocity, const real_t )
   {
      return velocity - latticeModel.forceModel().forceDensity(x,y,z) * real_t(0.5);
   }
};

template< typename LatticeModel_T >
struct AdaptVelocityToForce< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::shiftMacVel >::type >
{
   static Vector3<real_t> get( const LatticeModel_T &, const Vector3< real_t > & velocity, const real_t )
   {
      return velocity;
   }

   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator &, const LatticeModel_T &, const Vector3< real_t > & velocity, const real_t )
   {
      return velocity;
   }

   static Vector3<real_t> get( const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T &,
                               const Vector3< real_t > & velocity, const real_t )
   {
      return velocity;
   }
};

} // namespace internal



//////////////////////////////////////
// set density and velocity (x,y,z) //
//////////////////////////////////////

template< typename LatticeModel_T >
struct DensityAndVelocity
{
   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      Vector3< real_t > velAdaptedToForce = internal::AdaptVelocityToForce<LatticeModel_T>::get( it, latticeModel, velocity, rho );
      Equilibrium< LatticeModel_T >::set( it, velAdaptedToForce, rho );
   }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & latticeModel,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      Vector3< real_t > velAdaptedToForce = internal::AdaptVelocityToForce<LatticeModel_T>::get( x, y, z, latticeModel, velocity, rho );
      Equilibrium< LatticeModel_T >::set( pdf, x, y, z, velAdaptedToForce, rho );
   }
};

//////////////////////////////////////
// set density and velocity (range) //
//////////////////////////////////////

template< typename LatticeModel_T, typename FieldIteratorXYZ, class Enable = void >
struct DensityAndVelocityRange
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::DensityAndVelocityRange' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};



template< typename LatticeModel_T, typename FieldIteratorXYZ >
struct DensityAndVelocityRange< LatticeModel_T, FieldIteratorXYZ, typename std::enable_if< LatticeModel_T::ForceModel::constant >::type >
{
   static_assert( LatticeModel_T::ForceModel::constant, "Only works with constant forces!" );

   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end, const LatticeModel_T & latticeModel,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      Vector3< real_t > velAdaptedToForce = internal::AdaptVelocityToForce<LatticeModel_T>::get( latticeModel, velocity, rho );
      EquilibriumRange< LatticeModel_T, FieldIteratorXYZ >::set( begin, end, velAdaptedToForce, rho );
   }
};



template< typename LatticeModel_T, typename FieldIteratorXYZ >
struct DensityAndVelocityRange< LatticeModel_T, FieldIteratorXYZ, typename std::enable_if< ! LatticeModel_T::ForceModel::constant >::type >
{
   static_assert( LatticeModel_T::ForceModel::constant == false, "Does not work with constant forces!" );

   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end, const LatticeModel_T & latticeModel,
                    const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) )
   {
      for( auto cell = begin; cell != end; ++cell )
         DensityAndVelocity< LatticeModel_T >::set( cell, latticeModel, velocity, rho );
   }
};



} // namespace lbm
} // namespace walberla
