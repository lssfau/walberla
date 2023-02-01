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
//! \file ForceModel.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "core/math/Matrix3.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/adaptors/VectorFieldAccessor.h"
#include "lbm/lattice_model/CollisionModel.h"
#include "lbm/lattice_model/EquilibriumDistribution.h"

#include <type_traits>



namespace walberla {
namespace lbm {
namespace force_model {



//**********************************************************************************************************************
/*!
*   \section docForceModel Force Model for the LBM
*
*   Every force model must implement the following concept:
*
*   1.  Every force model must have a "tag". See force_model::None_tag and class force_model::None. If you write your
*       own force model and none of the existing tags fit, just define your own empty struct and use it as tag.
*   2.  Every force model must define a type "DirectionIndependentTerms_T" which must be returned by function
*       "directionIndependentTerms" (see 9) and which is passed to function "forceTerm" (see 10) during the collision
*       step. [If your force model does not have or does not need to calculate any common, direction-independent terms,
*       you can use the empty struct force_model::NoDirectionIndependentTerms - see class force_model::None as an example]
*   3.  static const bool shiftMacVel: specifies whether or not during the evaluation of the macroscopic velocity there
*                                      is a shift by F/2
*   4.  static const bool shiftEquVel: specifies whether or not during the evaluation of the equilibrium velocity there
*                                      is a shift by F/2
*   5.  static const bool constant: true, if the body force density is constant for all cells throughout the entire simulation
*                                   (false otherwise)
*   6.  void configure( IBlock & block, StructuredBlockStorage & sbs ): this function is called when a force density model
*                                                                       instance is assigned to a block - can be used to
*                                                                       scale stored values to a different grid level
*   7.  const Vector3<real_t> forceDensity() const: must return the body force density (this function must only exist if constant == true)
*   8.  const Vector3<real_t> forceDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const:
*        -> must return the force density for block local cell (x,y,z) - this function must also exist for constant force densities!
*   9.  "template< typename LatticeModel_T >
*       DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
*                                                              const Vector3<real_t> & velocity, const real_t rho,
*                                                              const real_t omega ) const":
*        -> can be used to pre-calculate common, direction-independent terms that are needed for every lattice direction
*           of cell (x,y,z)
*           Parameters: (x,y,z) = the cell that is currently processed, velocity/rho = the velocity and density of the
*                       equilibrium distribution, omega = relaxation parameter that corresponds to the lattice viscosity
*   10. "template< typename LatticeModel_T >
*       real_t forceTerm( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
*                         const Vector3<real_t> & velocity, const real_t rho,
*                         const DirectionIndependentTerms_T & commonTerms,
*                         const real_t w, const real_t cx, const real_t cy, const real_t cz, const real_t omega ) const":
*        -> must evaluate the force term of the LBM collision step - called from the LB compute kernel for every lattice
*           direction of cell (x,y,z)
*           Parameters: (x,y,z) = the cell that is currently processed, velocity/rho = the velocity and density of the
*                       equilibrium distribution, commonTerms = direction-independent terms calculated by function
*                       "directionIndependentTerms" (see 9), w = lattice model weighting factor, (cx,cy,cz) = lattice
*                       direction, omega = relaxation parameter that corresponds to the lattice viscosity
*   11. "bool setConstantBodyForceIfPossible( const Vector3<real_t> & forceDensity, const uint_t level = uint_t(0) )":
*        -> Must return true if a constant force density can be set. In that case, the force density passed in must be set.
*           'level' is the grid level 'forceDensity' corresponds to.
*
*   ATTENTION: Both forceDensity functions (7+8) and the functions calculating common terms and the force density term (9+10) may be
*              called in parallel by multiple threads and therefore must be thread-safe!
*
*   For an example see classes SimpleConstant or LuoField.
*/
//**********************************************************************************************************************

// For documentation about lattice models and how to assemble a lattice model that consists of
// a collision and a force model see documentation in LatticeModelBase.h

struct None_tag {};
struct Simple_tag {}; // Refers to a simple force model which introduces the following additional force density term in the
                      // collision process: 3 * (lattice_weight)_i * (e)_i * (force)_i [often: (force)_i = rho * acceleration]
                      // Should only be used with constant force densities!
                      // Shifts the macroscopic velocity by F/2, but does not change the equilibrium velocity.
struct EDM_tag {};    // Refers to the exact difference method by Kupershtokh
                      // Shifts the macroscopic velocity by F/2, but does not change the equilibrium velocity.
struct Luo_tag {};    // Refers to the force model by Luo which introduces a more complex force density term in the collision process:
                      // F_i = w_i * [ (c_i - u) / (c_s^2) + (c_i * (c_i * u)) / (c_s^4) ] * F
                      // Shifts the macroscopic velocity by F/2, but does not change the equilibrium velocity.
struct Guo_tag {};    // Refers to the force model by Guo which introduces a more complex force density term in the collision process:
                      // F_i = w_i * [ 1 - 1 / (2 * tau) ] * [ (c_i - u) / (c_s^2) + (c_i * (c_i * u)) / (c_s^4) ] * F
                      // Adapts the calculation of the macroscopic velocity as well as the equilibrium velocity (both shifted by F/2)!

struct Correction_tag{}; // refers to correction terms proposed in Jonas Latt's thesis


/// Used as return type if there are no common, direction-independent terms for the force density term
struct NoDirectionIndependentTerms {};


/// Returns the body force density for level 'targetLevel' given the body force density 'forceDensity_level' on level 'level'
inline Vector3< real_t > levelDependentBodyForce( const uint_t targetLevel, const Vector3< real_t > & forceDensity_level, const uint_t level )
{
   const real_t powTwoTarget = real_c( uint_t(1) << targetLevel );
   const real_t powTwoLevel  = real_c( uint_t(1) << level );
   return forceDensity_level * ( powTwoLevel / powTwoTarget );
}

/// Returns the acceleration for level 'targetLevel' given the acceleration 'acceleration_level' on level 'level'
inline Vector3< real_t > levelDependentAcceleration( const uint_t targetLevel, const Vector3< real_t > & acceleration_level, const uint_t level )
{
   const real_t powTwoTarget = real_c( uint_t(1) << targetLevel );
   const real_t powTwoLevel  = real_c( uint_t(1) << level );
   return acceleration_level * ( powTwoLevel / powTwoTarget );
}



class None
{
public:

   using tag = None_tag;
   using DirectionIndependentTerms_T = NoDirectionIndependentTerms;

   static const bool shiftMacVel = false;
   static const bool shiftEquVel = false;

   static const bool constant = true;

   void pack( mpi::SendBuffer & ) const {}
   void unpack( mpi::RecvBuffer & ) {}

   void configure( IBlock & /*block*/, StructuredBlockStorage & /*sbs*/ ) const {}

   const Vector3<real_t> forceDensity() const { return Vector3<real_t>(); }
   const Vector3<real_t> forceDensity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return Vector3<real_t>(); }

   template< typename LatticeModel_T >
   DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                                                          const Vector3<real_t> & /*velocity*/, const real_t /*rho*/, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   { return DirectionIndependentTerms_T(); }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const Vector3<real_t> & /*velocity*/, const real_t /*rho*/,
                     const DirectionIndependentTerms_T & /*commonTerms*/, const real_t /*w*/,
                     const real_t /*cx*/, const real_t /*cy*/, const real_t /*cz*/, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const { return real_t(0); }

   bool setConstantBodyForceIfPossible( const Vector3<real_t> &, const uint_t = uint_t(0) ) { return false; }
};



/// For the incompressible LBM, in lattice units, the body force density is equal to the acceleration, since
/// bodyForceDensity_lattice = density_lattice * acceleration_lattice and density_lattice is implicitly set to 1
class SimpleConstant
{
public:

   using tag = Simple_tag;
   using DirectionIndependentTerms_T = NoDirectionIndependentTerms;

   static const bool shiftMacVel = true;
   static const bool shiftEquVel = false;

   static const bool constant = true;

   SimpleConstant( const Vector3< real_t > & bodyForceDensity, const uint_t level = uint_t(0) ) :
      bodyForceDensity_( bodyForceDensity ), level_( level ) {}
   SimpleConstant( const real_t vx, const real_t vy, const real_t vz, const uint_t level = uint_t(0) ) :
      bodyForceDensity_(vx,vy,vz), level_( level ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << bodyForceDensity_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> bodyForceDensity_ >> level_; }

   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t level = level_;
      level_ = sbs.getLevel( block );
      bodyForceDensity_ = levelDependentBodyForce( level_, bodyForceDensity_, level );
   }

   /// "forceDensity_level" is the level that corresponds to "acceleration"
   void reset( const Vector3< real_t > & bodyForceDensity, const uint_t forceDensity_level = uint_t(0) )
   {
      bodyForceDensity_ = levelDependentBodyForce( level_, bodyForceDensity, forceDensity_level );
   }

   const Vector3< real_t > & forceDensity() const { return bodyForceDensity_; }
   const Vector3< real_t > & forceDensity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return bodyForceDensity_; }

   template< typename LatticeModel_T >
   DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                                                          const Vector3<real_t> & /*velocity*/, const real_t /*rho*/, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   { return DirectionIndependentTerms_T(); }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const Vector3<real_t> & /*velocity*/, const real_t /*rho*/,
                     const DirectionIndependentTerms_T & /*commonTerms*/, const real_t w,
                     const real_t cx, const real_t cy, const real_t cz, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      return real_t(3.0) * w * ( cx * bodyForceDensity_[0] + cy * bodyForceDensity_[1] + cz * bodyForceDensity_[2] );
   }

   /// "forceDensity_level" is the level that corresponds to "acceleration"
   bool setConstantBodyForceIfPossible( const Vector3< real_t > & bodyForceDensity, const uint_t forceDensity_level = uint_t(0) )
   {
      reset( bodyForceDensity, forceDensity_level );
      return true;
   }

private:

   Vector3< real_t > bodyForceDensity_;

   uint_t level_;
};



/// \cite kupershtokh2003calculations, \cite kupershtokh2004incorporating
template< typename ForceField_T > // ForceField_T: Field with fSize = 1 and data type = Vector3< real_t >
class EDMField
{
private:

   template< typename LatticeModel_T, class Enable = void >
   struct DirectionIndependentTerm
   {
      static real_t get( const real_t ) { return real_t(0); }
   };
   template< typename LatticeModel_T >
   struct DirectionIndependentTerm< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible >::type >
   {
      static real_t get( const real_t rho ) { return real_t(1) / rho; }
   };

   template< typename LatticeModel_T, class Enable = void >
   struct ShiftedVelocity
   {
      static Vector3<real_t> get( const Vector3<real_t> & velocity, const Vector3<real_t> & forceDensity, const real_t )
      {
         return velocity + forceDensity;
      }
   };
   template< typename LatticeModel_T >
   struct ShiftedVelocity< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible >::type >
   {
      static Vector3<real_t> get( const Vector3<real_t> & velocity, const Vector3<real_t> & forceDensity, const real_t rhoInv )
      {
         return velocity + forceDensity * rhoInv;
      }
   };

public:

   using tag = EDM_tag;
   using DirectionIndependentTerms_T = real_t;
   using ForceField = ForceField_T;

   static const bool shiftMacVel = true;
   static const bool shiftEquVel = false;

   static const bool constant = false;

   EDMField( const BlockDataID & forceDensityFieldId ) :
      forceDensityFieldId_( forceDensityFieldId ), forceDensityField_( NULL ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << forceDensityFieldId_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> forceDensityFieldId_; }

   void configure( IBlock & block, StructuredBlockStorage & /*sbs*/ )
   {
      forceDensityField_ = block.getData< ForceField_T >( forceDensityFieldId_ );
   }

   typename field::VectorFieldAccessor<ForceField_T>::vector_or_constRefVector forceDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      return field::VectorFieldAccessor<ForceField_T>::get( forceDensityField_, x,y,z);
   }

   template< typename LatticeModel_T >
   real_t directionIndependentTerms( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                                     const Vector3<real_t> & /*velocity*/, const real_t rho, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      return DirectionIndependentTerm< LatticeModel_T >::get( rho );
   }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const Vector3<real_t> & velocity, const real_t rho,
                     const real_t commonTerm, const real_t w,
                     const real_t cx, const real_t cy, const real_t cz, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      const auto shiftedVelocity = ShiftedVelocity< LatticeModel_T >::get( velocity, forceDensity(x,y,z), commonTerm );

      return EquilibriumDistribution< LatticeModel_T >::get( cx, cy, cz, w, shiftedVelocity, rho )
               - EquilibriumDistribution< LatticeModel_T >::get( cx, cy, cz, w, velocity, rho );
   }

   bool setConstantBodyForceIfPossible( const Vector3< real_t > &, const uint_t = uint_t(0) ) { return false; }

private:

   BlockDataID forceDensityFieldId_;
   ForceField_T * forceDensityField_;
};



/// For the incompressible LBM, in lattice units, the body force is equal to the acceleration, since
/// bodyForceDensity_lattice = density_lattice * acceleration_lattice and density_lattice is implicitly set to 1.
/// \cite luo1993lattice, \cite shan1994simulation
class LuoConstant
{
public:

   using tag = Luo_tag;
   using DirectionIndependentTerms_T = NoDirectionIndependentTerms;

   static const bool shiftMacVel = true;
   static const bool shiftEquVel = false;

   static const bool constant = true;

   LuoConstant( const Vector3< real_t > & bodyForceDensity, const uint_t level = uint_t(0) ) :
      bodyForceDensity_( bodyForceDensity ), level_( level ) {}
   LuoConstant( const real_t vx, const real_t vy, const real_t vz, const uint_t level = uint_t(0) ) :
      bodyForceDensity_(vx,vy,vz), level_( level ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << bodyForceDensity_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> bodyForceDensity_ >> level_; }

   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t level = level_;
      level_ = sbs.getLevel( block );
      bodyForceDensity_ = levelDependentBodyForce( level_, bodyForceDensity_, level );
   }

   /// "forceDensity_level" is the level that corresponds to "bodyForceDensity"
   void reset( const Vector3< real_t > & bodyForceDensity, const uint_t forceDensity_level = uint_t(0) )
   {
      bodyForceDensity_ = levelDependentBodyForce( level_, bodyForceDensity, forceDensity_level );
   }

   const Vector3< real_t > & forceDensity() const { return bodyForceDensity_; }
   const Vector3< real_t > & forceDensity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return bodyForceDensity_; }

   template< typename LatticeModel_T >
   DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                                                          const Vector3<real_t> & /*velocity*/, const real_t /*rho*/, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   { return DirectionIndependentTerms_T(); }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const Vector3<real_t> & velocity, const real_t /*rho*/,
                     const DirectionIndependentTerms_T & /*commonTerms*/, const real_t w,
                     const real_t cx, const real_t cy, const real_t cz, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      const Vector3<real_t> c( cx, cy, cz );
      return real_t(3) * w * ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * bodyForceDensity_ );
   }

   /// "forceDensity_level" is the level that corresponds to "bodyForceDensity"
   bool setConstantBodyForceIfPossible( const Vector3< real_t > & bodyForceDensity, const uint_t forceDensity_level = uint_t(0) )
   {
      reset( bodyForceDensity, forceDensity_level );
      return true;
   }

private:

   Vector3< real_t > bodyForceDensity_;

   uint_t level_;
};



template< typename ForceField_T > // ForceField_T: Field with fSize = 1 and data type = Vector3< real_t >
class LuoField
{
public:

   using tag = Luo_tag;
   using DirectionIndependentTerms_T = NoDirectionIndependentTerms;
   using ForceField = ForceField_T;

   static const bool shiftMacVel = true;
   static const bool shiftEquVel = false;

   static const bool constant = false;

   LuoField( const BlockDataID & forceDensityFieldId ) :
      forceDensityFieldId_( forceDensityFieldId ), forceDensityField_( NULL ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << forceDensityFieldId_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> forceDensityFieldId_; }

   void configure( IBlock & block, StructuredBlockStorage & /*sbs*/ )
   {
      forceDensityField_ = block.getData< ForceField_T >( forceDensityFieldId_ );
   }

   typename field::VectorFieldAccessor<ForceField_T>::vector_or_constRefVector forceDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      return field::VectorFieldAccessor<ForceField_T>::get( forceDensityField_, x,y,z);
   }

   template< typename LatticeModel_T >
   DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                                                          const Vector3<real_t> & /*velocity*/, const real_t /*rho*/, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   { return DirectionIndependentTerms_T(); }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const Vector3<real_t> & velocity, const real_t /*rho*/,
                     const DirectionIndependentTerms_T & /*commonTerms*/, const real_t w,
                     const real_t cx, const real_t cy, const real_t cz, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      const Vector3<real_t> c( cx, cy, cz );
      return real_t(3) * w * ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * forceDensity(x,y,z) );
   }

   bool setConstantBodyForceIfPossible( const Vector3< real_t > &, const uint_t = uint_t(0) ) { return false; }

private:

   BlockDataID forceDensityFieldId_;
   ForceField_T * forceDensityField_;
};



/// For the incompressible LBM, in lattice units, the body force is equal to the acceleration, since
/// bodyForceDensity_lattice = density_lattice * acceleration_lattice and density_lattice is implicitly set to 1.
/// \cite guo2002discrete, \cite schiller2008thermal
class GuoConstant
{
public:

   using tag = Guo_tag;
   using DirectionIndependentTerms_T = Matrix3<real_t>;

   static const bool shiftMacVel = true;
   static const bool shiftEquVel = true;

   static const bool constant = true;

   GuoConstant( const Vector3< real_t > & bodyForceDensity, const uint_t level = uint_t(0) ) :
      bodyForceDensity_( bodyForceDensity ), level_( level ) {}
   GuoConstant( const real_t vx, const real_t vy, const real_t vz, const uint_t level = uint_t(0) ) :
      bodyForceDensity_(vx,vy,vz), level_( level ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << bodyForceDensity_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> bodyForceDensity_ >> level_; }

   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t level = level_;
      level_ = sbs.getLevel( block );
      bodyForceDensity_ = levelDependentBodyForce( level_, bodyForceDensity_, level );
   }

   /// "forceDensity_level" is the level that corresponds to "bodyForceDensity"
   void reset( const Vector3< real_t > & bodyForceDensity, const uint_t forceDensity_level = uint_t(0) )
   {
      bodyForceDensity_ = levelDependentBodyForce( level_, bodyForceDensity, forceDensity_level );
   }

   const Vector3< real_t > & forceDensity() const { return bodyForceDensity_; }
   const Vector3< real_t > & forceDensity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return bodyForceDensity_; }

   template< typename LatticeModel_T >
   DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                                                          const Vector3<real_t> & velocity, const real_t /*rho*/, const real_t omega, const real_t omega_bulk, const real_t /*omega_odd*/ ) const
   {
      if (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::MRT_tag >::value)
      {
         const real_t one_over_d  = real_t(1) / real_t(LatticeModel_T::Stencil::D);

         const auto common = Matrix3<real_t>::makeDiagonalMatrix( velocity * bodyForceDensity_ );
         return (tensorProduct( velocity, bodyForceDensity_ ) +
                 tensorProduct( bodyForceDensity_, velocity ) -
                 common * (real_t(2)*one_over_d) ) * real_t(0.5) * ( real_t(2) - omega )
                + common * ( one_over_d * ( real_t(2) - omega_bulk ) );
      }
      else
      {
         WALBERLA_ASSERT_FLOAT_EQUAL( omega, omega_bulk );
         return DirectionIndependentTerms_T();
      }
   }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const Vector3<real_t> & velocity, const real_t /*rho*/,
                     const DirectionIndependentTerms_T & commonTerms, const real_t w,
                     const real_t cx, const real_t cy, const real_t cz, const real_t omega, const real_t /*omega_bulk*/, const real_t omega_odd ) const
   {
      const Vector3<real_t> c( cx, cy, cz );
      if (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::MRT_tag >::value)
      {
         const real_t one_third  = real_t(1) / real_t(3);

         const real_t common = (commonTerms * ( tensorProduct(c,c) - Matrix3<real_t>::makeDiagonalMatrix(one_third) )).trace();
         return real_t(3.0) * w * ( bodyForceDensity_ * c + real_t(1.5) * common);
      }
      else if (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value)
      {
         return real_t(3.0) * w * ( ( real_t(1) - real_t(0.5) * omega ) *
                                    ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * bodyForceDensity_ ) +
                                    ( omega - omega_odd ) * real_t(0.5) * (c * bodyForceDensity_) );
      }
      else
      {
         return real_t(3.0) * w * ( real_t(1) - real_t(0.5) * omega ) *
                                  ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * bodyForceDensity_ );
      }
   }

   /// "forceDensity_level" is the level that corresponds to "bodyForceDensity"
   bool setConstantBodyForceIfPossible( const Vector3< real_t > & bodyForceDensity, const uint_t forceDensity_level = uint_t(0) )
   {
      reset( bodyForceDensity, forceDensity_level );
      return true;
   }

private:

   Vector3< real_t > bodyForceDensity_;

   uint_t level_;
};



template< typename ForceField_T > // ForceField_T: Field with fSize = 1 and data type = Vector3< real_t >
class GuoField
{
public:

   using tag = Guo_tag;
   using DirectionIndependentTerms_T = Matrix3<real_t>;
   using ForceField = ForceField_T;

   static const bool shiftMacVel = true;
   static const bool shiftEquVel = true;

   static const bool constant = false;

   GuoField( const BlockDataID & forceDensityFieldId ) :
      forceDensityFieldId_( forceDensityFieldId ), forceDensityField_( nullptr ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << forceDensityFieldId_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> forceDensityFieldId_; }

   void configure( IBlock & block, StructuredBlockStorage & /*sbs*/ )
   {
      forceDensityField_ = block.getData< ForceField_T >( forceDensityFieldId_ );
   }

   typename field::VectorFieldAccessor<ForceField_T>::vector_or_constRefVector forceDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      return field::VectorFieldAccessor<ForceField_T>::get( forceDensityField_, x,y,z);
   }

   template< typename LatticeModel_T >
   DirectionIndependentTerms_T directionIndependentTerms( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                          const Vector3<real_t> & velocity, const real_t /*rho*/, const real_t omega, const real_t omega_bulk, const real_t /*omega_odd*/ ) const
   {
      if (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::MRT_tag >::value)
      {
         const real_t one_over_d  = real_t(1) / real_t(LatticeModel_T::Stencil::D);

         const auto common = Matrix3<real_t>::makeDiagonalMatrix( velocity * forceDensity(x,y,z) );
         return (tensorProduct( velocity, forceDensity(x,y,z) ) +
                 tensorProduct( forceDensity(x,y,z), velocity ) -
                 common * (real_t(2)*one_over_d) ) * real_t(0.5) * ( real_t(2) - omega )
                + common * ( one_over_d * ( real_t(2) - omega_bulk ) );
      }
      else
      {
         WALBERLA_ASSERT_FLOAT_EQUAL( omega, omega_bulk );
         return DirectionIndependentTerms_T();
      }
   }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const Vector3<real_t> & velocity, const real_t /*rho*/,
                     const DirectionIndependentTerms_T & commonTerms, const real_t w,
                     const real_t cx, const real_t cy, const real_t cz, const real_t omega, const real_t /*omega_bulk*/, const real_t omega_odd ) const
   {
      const Vector3<real_t> c( cx, cy, cz );
      if (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::MRT_tag >::value)
      {
         const real_t one_third  = real_t(1) / real_t(3);

         const real_t common = (commonTerms * ( tensorProduct(c,c) - Matrix3<real_t>::makeDiagonalMatrix(one_third) )).trace();
         return real_t(3.0) * w * ( forceDensity(x,y,z) * c + real_t(1.5) * common);
      }
      else if (std::is_same< typename LatticeModel_T::CollisionModel::tag, collision_model::TRT_tag >::value)
      {
         return real_t(3.0) * w * ( ( real_t(1) - real_t(0.5) * omega ) *
                                    ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * forceDensity(x,y,z) ) +
                                    ( omega - omega_odd ) * real_t(0.5) * (c * forceDensity(x,y,z)) );
      }
      else
      {
         return real_t(3.0) * w * ( real_t(1) - real_t(0.5) * omega ) *
                                  ( ( c - velocity + ( real_t(3) * ( c * velocity ) * c ) ) * forceDensity(x,y,z) );
      }
   }

   bool setConstantBodyForceIfPossible( const Vector3< real_t > &, const uint_t = uint_t(0) ) { return false; }

private:

   BlockDataID forceDensityFieldId_;
   ForceField_T * forceDensityField_;
};



/// \cite latt2007hydrodynamic
template< typename MomentumDensityField_T >
class Correction
{
public:

   using tag = Correction_tag;
   using DirectionIndependentTerms_T = Vector3<real_t>;

   static const bool shiftMacVel = false;
   static const bool shiftEquVel = false;

   static const bool constant = true;

   Correction( const BlockDataID & previousRhoVelocityId ) :
      forceDensity_( real_t(0) ), previousRhoVelocityId_( previousRhoVelocityId ), previousRhoVelocity_(nullptr) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << forceDensity_ << previousRhoVelocityId_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> forceDensity_ >> previousRhoVelocityId_; }

   void configure( IBlock & block, StructuredBlockStorage & /*sbs*/ )
   {
      previousRhoVelocity_ = block.getData< MomentumDensityField_T >( previousRhoVelocityId_ );
   }

   const Vector3< real_t > & forceDensity() const { return forceDensity_; }
   const Vector3< real_t > & forceDensity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return forceDensity_; }

   template< typename LatticeModel_T >
   Vector3<real_t> directionIndependentTerms( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                              const Vector3<real_t> & velocity, const real_t rho, const real_t omega, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      const auto rhoVelocity = rho * velocity;
      Vector3<real_t> & previousRhoVelocity = previousRhoVelocity_->get(x,y,z);
      Vector3<real_t> result = ( real_c(3) - real_c(1.5) * omega ) * ( rhoVelocity - previousRhoVelocity );
      previousRhoVelocity = rhoVelocity;
      return result;
   }

   template< typename LatticeModel_T >
   real_t forceTerm( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const Vector3<real_t> & /*velocity*/, const real_t /*rho*/,
                     const Vector3<real_t> & commonTerm, const real_t w, const real_t cx, const real_t cy, const real_t cz, const real_t /*omega*/, const real_t /*omega_bulk*/, const real_t /*omega_odd*/ ) const
   {
      return w * ( cx * commonTerm[0] + cy * commonTerm[1] + cz * commonTerm[2] );
   }

   bool setConstantBodyForceIfPossible( const Vector3<real_t> &, const uint_t = uint_t(0) ) { return false; }

private:

   Vector3< real_t > forceDensity_;

   BlockDataID previousRhoVelocityId_;
   MomentumDensityField_T * previousRhoVelocity_;
};



} // namespace force_model
} // namespace lbm
} // namespace walberla
