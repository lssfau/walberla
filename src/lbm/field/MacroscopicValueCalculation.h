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
//! \file MacroscopicValueCalculation.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "Density.h"
#include "DensityAndMomentumDensity.h"
#include "DensityAndVelocity.h"
#include "Equilibrium.h"
#include "ShearRate.h"
#include "PressureTensor.h"
#include "QCriterion.h"
#include "Vorticity.h"

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/math/Matrix3.h"


namespace walberla {
namespace lbm {



////////////////////////////
// SET DENSITY & VELOCITY //
////////////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void setDensityAndVelocity( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel,
                                   const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );

/////////////////
// EQUILIBRIUM //
/////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void setToEquilibrium( FieldPtrOrIterator & it,
                              const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );

/////////////
// DENSITY //
/////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

/////////////////////////
// DENSITY IN SI UNITS //
/////////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensitySI( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it, const real_t rho_SI );

//////////////////////
// MOMENTUM DENSITY //
//////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getMomentumDensity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getEquilibriumMomentumDensity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );


//////////////
// VELOCITY //
//////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getVelocity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getEquilibriumVelocity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getEquilibriumVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

//////////////////////////
// VELOCITY IN SI UNITS //
//////////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getVelocitySI( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it, const real_t dx_SI, const real_t dt_SI );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel,const FieldPtrOrIterator & it,
                           const real_t dx_SI, const real_t dt_SI );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getVelocitySI( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it, const real_t dxDividedByDt_SI );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                           const real_t dxDividedByDt_SI );

////////////////////////////////
// DENSITY & MOMENTUM DENSITY //
////////////////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );


////////////////////////
// DENSITY & VELOCITY //
////////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndEquilibriumVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

////////////////////////////////////
// DENSITY & VELOCITY IN SI UNITS //
////////////////////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                                       const real_t rho_SI, const real_t dx_SI, const real_t dt_SI );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                                       const real_t rho_SI, const real_t dxDividedByDt_SI );

////////////////
// SHEAR RATE //
////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getShearRate( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

/////////////////////
// PRESSURE TENSOR //
/////////////////////

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Matrix3< real_t > getPressureTensor( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );

template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getPressureTensor( Matrix3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it );



/////////////////////
// IMPLEMENTATIONS //
/////////////////////



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void setDensityAndVelocity( FieldPtrOrIterator & it, const LatticeModel_T & latticeModel, const Vector3< real_t > & velocity, const real_t rho )
{
   DensityAndVelocity< LatticeModel_T >::set( it, latticeModel, velocity, rho );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void setToEquilibrium( FieldPtrOrIterator & it, const Vector3< real_t > & velocity, const real_t rho )
{
   Equilibrium< LatticeModel_T >::set( it, velocity, rho );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   return Density< LatticeModel_T >::get( latticeModel, it );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensitySI( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it, const real_t rho_SI )
{
   return getDensity( latticeModel, it ) * rho_SI;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getMomentumDensity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   Vector3< real_t > momentumDensity;
   getMomentumDensity< LatticeModel_T >( momentumDensity, latticeModel, it );
   return momentumDensity;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   MomentumDensity< LatticeModel_T >::get( momentumDensity, latticeModel, it );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getEquilibriumMomentumDensity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   Vector3< real_t > momentumDensity;
   getEquilibriumMomentumDensity< LatticeModel_T >( momentumDensity, latticeModel, it );
   return momentumDensity;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   MomentumDensity< LatticeModel_T >::getEquilibrium( momentumDensity, latticeModel, it );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getVelocity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   Vector3< real_t > velocity;
   getVelocity< LatticeModel_T >( velocity, latticeModel, it );
   return velocity;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   if( LatticeModel_T::compressible )
   {
      const real_t rho = getDensityAndMomentumDensity< LatticeModel_T >( velocity, latticeModel, it );
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   else
   {
      MomentumDensity< LatticeModel_T >::get( velocity, latticeModel, it );
   }
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getEquilibriumVelocity( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   Vector3< real_t > velocity;
   getEquilibriumVelocity< LatticeModel_T >( velocity, latticeModel, it );
   return velocity;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getEquilibriumVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   if( LatticeModel_T::compressible )
   {
      const real_t rho = getDensityAndEquilibriumMomentumDensity< LatticeModel_T >( velocity, latticeModel, it );
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   else
   {
      MomentumDensity< LatticeModel_T >::getEquilibrium( velocity, latticeModel, it );
   }
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getVelocitySI( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it, const real_t dx_SI, const real_t dt_SI )
{
   Vector3< real_t > velocity;
   getVelocity< LatticeModel_T >( velocity, latticeModel, it );
   return velocity * ( dx_SI / dt_SI );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                           const real_t dx_SI, const real_t dt_SI )
{
   getVelocity< LatticeModel_T >( velocity, latticeModel, it );
   velocity *= ( dx_SI / dt_SI );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Vector3< real_t > getVelocitySI( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it, const real_t dxDividedByDt_SI )
{
   Vector3< real_t > velocity;
   getVelocity< LatticeModel_T >( velocity, latticeModel, it );
   return velocity * dxDividedByDt_SI;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                           const real_t dxDividedByDt_SI )
{
   getVelocity< LatticeModel_T >( velocity, latticeModel, it );
   velocity *= dxDividedByDt_SI;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   return DensityAndMomentumDensity< LatticeModel_T >::get( momentumDensity, latticeModel, it );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   return DensityAndMomentumDensity< LatticeModel_T >::getEquilibrium( momentumDensity, latticeModel, it );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   const real_t rho = getDensityAndMomentumDensity< LatticeModel_T >( velocity, latticeModel, it );
   if( LatticeModel_T::compressible )
   {
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   return rho;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndEquilibriumVelocity( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   const real_t rho = getDensityAndEquilibriumMomentumDensity< LatticeModel_T >( velocity, latticeModel, it );
   if( LatticeModel_T::compressible )
   {
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   return rho;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                                       const real_t rho_SI, const real_t dx_SI, const real_t dt_SI )
{
   const real_t rho = getDensityAndVelocity< LatticeModel_T >( velocity, latticeModel, it );
   velocity *= ( dx_SI / dt_SI );
   return rho * rho_SI;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it,
                                       const real_t rho_SI, const real_t dxDividedByDt_SI )
{
   const real_t rho = getDensityAndVelocity< LatticeModel_T >( velocity, latticeModel, it );
   velocity *= dxDividedByDt_SI;
   return rho * rho_SI;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline real_t getShearRate( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   Vector3< real_t > velocity;
   const real_t rho = getDensityAndVelocity( velocity, latticeModel, it );

   return ShearRate< LatticeModel_T >::get( latticeModel, it, velocity, rho );
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline Matrix3< real_t > getPressureTensor( const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   Matrix3< real_t > pressureTensor;
   getPressureTensor< LatticeModel_T >( pressureTensor, latticeModel, it );
   return pressureTensor;
}



template< typename LatticeModel_T, typename FieldPtrOrIterator >
inline void getPressureTensor( Matrix3< real_t > & pressureTensor, const LatticeModel_T & latticeModel, const FieldPtrOrIterator & it )
{
   PressureTensor<LatticeModel_T>::get( pressureTensor, it );
}


template< typename VelocityField_T, typename Filter_T >
inline real_t getQCriterion(const VelocityField_T &velocityField, const Filter_T &filter,
                            const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                            real_t dx = real_t(1), real_t dy = real_t(1), real_t dz = real_t(1)) {
   return QCriterion::get(velocityField, filter, x, y, z, dx, dy, dz);
}

template< typename VelocityField_T, typename Filter_T >
inline Vector3<real_t> getVorticity(const VelocityField_T &velocityField, const Filter_T &filter,
                                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                    real_t dx = real_t(1), real_t dy = real_t(1), real_t dz = real_t(1)) {
   return Vorticity::get(velocityField, filter, x, y, z, dx, dy, dz);
}



} // namespace lbm
} // namespace walberla
