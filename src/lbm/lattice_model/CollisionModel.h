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
//! \file CollisionModel.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/mpi/BufferSizeTrait.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"
#include "domain_decomposition/StructuredBlockStorage.h"


namespace walberla {
namespace lbm {
namespace collision_model {



// For documentation about lattice models and how to assemble a lattice model that consists of
// a collision and a force model see documentation in LatticeModelBase.h

struct SRT_tag {};
struct TRT_tag {};
struct MRT_tag {};
struct Cumulant_tag{};



/// Returns the relaxation parameter for level 'targetLevel' given the relaxation parameter 'parameter_level' on level 'level'
inline real_t levelDependentRelaxationParameter( const uint_t targetLevel, const real_t parameter_level, const uint_t level )
{
   real_t powFactor(1);
   if( level < targetLevel )
      powFactor = real_c( uint_t(1) << ( targetLevel - level ) );
   else if( level == targetLevel )
      return parameter_level;
   else
      powFactor = real_t(1) / real_c( uint_t(1) << ( level - targetLevel ) );

   const real_t parameter_level_half = real_c(0.5) * parameter_level;
   return parameter_level / ( parameter_level_half + powFactor * ( real_t(1) - parameter_level_half ) );
}

inline real_t viscosityFromOmega( const real_t omega )
{
   static const real_t one_third = real_t(1) / real_t(3);
   return ( real_t(1) / omega - real_t(0.5) ) * one_third;
}

inline real_t omegaFromViscosity( const real_t viscosity )
{
   return real_t(1) / ( real_t(0.5) + real_t(3) * viscosity );
}



//**********************************************************************************************************************
/*! Implementation of LB collision model SRT \cite bhatnagar1954model
*/
//**********************************************************************************************************************

class SRT
{
public:

   typedef SRT_tag tag;
   static const bool constant = true;
   
   SRT( const real_t _omega, const uint_t _level = uint_t(0) ) :
      omega_( _omega ), viscosity_( viscosityFromOmega( _omega ) ), level_( _level ) {}

   void pack( mpi::SendBuffer & buffer ) const { buffer << omega_ << viscosity_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> omega_ >> viscosity_ >> level_; }

   /// Omega is adapted to the "right", level-dependent omega once "configure" is called
   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t _level = level_;
      level_ = sbs.getLevel( block );
      omega_ = levelDependentRelaxationParameter( level_, omega_, _level );
      viscosity_ = viscosityFromOmega( omega_ );
   }

   /// Only call this function if you know what you're doing (changes the viscosity!)
   /// "_omega_level" is the level that corresponds to "_omega"
   void reset( const real_t _omega, const uint_t _omega_level = uint_t(0) )
   {
      omega_ = levelDependentRelaxationParameter( level_, _omega, _omega_level );
      viscosity_ = viscosityFromOmega( omega_ );
   }

   real_t omega() const { return omega_; }
   inline real_t omega_bulk() const { return omega(); }
   real_t viscosity() const { return viscosity_; }

   real_t omega( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                 const Vector3<real_t> & /*velocity*/ = Vector3<real_t>(), const real_t /*rho*/ = real_t(1) ) const { return omega_; }
   real_t viscosity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return viscosity_; }

   real_t viscosity( const uint_t _level ) const
   {
      const real_t _omega = levelDependentRelaxationParameter( _level, omega_, level_ );
      return viscosityFromOmega( _omega );
   }

   uint_t level() const { return level_; }

private:

   real_t omega_;
   real_t viscosity_;

   uint_t level_;
};



template< typename OmegaField_T > // OmegaField_T: Field with fSize = 1 and data type = real_t
class SRTField
{
public:

   typedef SRT_tag tag;
   static const bool constant = false;

   SRTField() :
      omegaFieldId_(), omegaField_( NULL ), level_( uint_t(0) ) {}

   SRTField( const BlockDataID & omegaFieldId, const uint_t _level = uint_t(0) ) :
      omegaFieldId_( omegaFieldId ), omegaField_( nullptr ), level_( _level ) {}

   void setFieldId( const BlockDataID & omegaFieldId, const uint_t _level = uint_t(0) )
   {
      omegaFieldId_ = omegaFieldId;
      level_ = _level;
   }

   void pack( mpi::SendBuffer & buffer ) const { buffer << omegaFieldId_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> omegaFieldId_ >> level_; }

   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t _level = level_;
      level_ = sbs.getLevel( block );
      omegaField_ = block.getData< OmegaField_T >( omegaFieldId_ );
      for( auto it = omegaField_->begin(); it != omegaField_->end(); ++it )
         reset( it.x(), it.y(), it.z(), *it, _level );
   }

   /// Only call this function if you know what you're doing (changes the viscosity!)
   /// "_omega_level" is the level that corresponds to "_omega"
   void reset( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const real_t _omega, const uint_t _omega_level = uint_t(0) )
   {
      omegaField_->get(x,y,z) = levelDependentRelaxationParameter( level_, _omega, _omega_level );
   }

   real_t omega( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                 const Vector3<real_t> & /*velocity*/ = Vector3<real_t>(), const real_t /*rho*/ = real_t(1) ) const
   {
      return omegaField_->get(x,y,z);
   }

   real_t viscosity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const { return viscosityFromOmega( omega(x,y,z) ); }

   uint_t level() const { return level_; }

private:

   BlockDataID omegaFieldId_;
   OmegaField_T * omegaField_;

   uint_t level_;
};



//**********************************************************************************************************************
/*! Implementation of LB collision model TRT \cite ginzburg2008two
*   
*   Instead of providing two relaxation parameters, you can also provide one relaxation parameter omega and a "magic"
*   number 'x'. The two relaxation parameters are then calculated as follows:
*      \f[
*      \lambda_e = \omega
*      \f]
*      \f[
*      x = ( \frac{1}{\lambda_e} - 0.5 )( \frac{1}{\lambda_d} - 0.5 )
*      \f]
*   Good choices for 'x' are 3 / 16 (exact solution for Poiseuille flow with wall distance = 0.5 \cite ginzburg2008study) and
*   around 1 / 4 (best stability).
*/
//**********************************************************************************************************************

class TRT
{
public:

   typedef TRT_tag tag;
   
   static const real_t threeSixteenth;

   TRT( const real_t _lambda_e, const real_t _lambda_d, const uint_t _level = uint_t(0) ) :
      lambda_e_( _lambda_e ), lambda_d_( _lambda_d ),
      magicNumber_( magicNumber( _lambda_e, _lambda_d ) ),
      viscosity_( viscosityFromOmega( _lambda_e ) ), level_( _level )
   {}

   static TRT constructWithMagicNumber( const real_t _omega, const real_t _magicNumber = threeSixteenth,
                                        const uint_t _level = uint_t(0) )
   {
      TRT trt;
      trt.initWithMagicNumber( _omega, _magicNumber, _level );
      return trt;
   }

   void pack( mpi::SendBuffer & buffer ) const { buffer << lambda_e_ << lambda_d_ << magicNumber_ << viscosity_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { buffer >> lambda_e_ >> lambda_d_ >> magicNumber_ >> viscosity_ >> level_; }

   /// Adapts the two relaxation parameters to the "right", level-dependent parameters once "configure" is called
   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t _level = level_;
      level_ = sbs.getLevel( block );
      lambda_e_ = levelDependentRelaxationParameter( level_, lambda_e_, _level );
      lambda_d_ = lambda_d( lambda_e_, magicNumber_ );
      viscosity_ = viscosityFromOmega( lambda_e_ );
   }

   /// Only call this function if you know what you're doing (changes the viscosity!)
   /// "_lambda_level" is the level that corresponds to "_lambda_e" and "_lambda_d"
   void reset( const real_t _lambda_e, const real_t _lambda_d, const uint_t _lambda_level = uint_t(0) )
   {
      lambda_e_ = levelDependentRelaxationParameter( level_, _lambda_e, _lambda_level );
      magicNumber_ = magicNumber( _lambda_e, _lambda_d );
      lambda_d_ = lambda_d( lambda_e_, magicNumber_ );
      viscosity_ = viscosityFromOmega( lambda_e_ );
   }

   /// Only call this function if you know what you're doing (changes the viscosity!)
   /// "omega_level" is the level that corresponds to "omega"
   void resetWithMagicNumber( const real_t _omega, const real_t _magicNumber = threeSixteenth, const uint_t omega_level = uint_t(0) )
   {
      lambda_e_ = levelDependentRelaxationParameter( level_, lambda_e( _omega ), omega_level );
      lambda_d_ = lambda_d( lambda_e_, _magicNumber );
      magicNumber_ = _magicNumber;
      viscosity_ = viscosityFromOmega( lambda_e_ );
   }

   real_t lambda_e() const { return lambda_e_; }
   real_t lambda_d() const { return lambda_d_; }

   real_t viscosity() const { return viscosity_; }
   real_t viscosity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return viscosity_; }

   real_t omega() const { return lambda_e_; }
   real_t omega( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                 const Vector3<real_t> & /*velocity*/ = Vector3<real_t>(), const real_t /*rho*/ = real_t(1) ) const { return omega(); }

   inline real_t omega_bulk() const { return omega(); }

   static real_t lambda_e( const real_t _omega ) { return _omega; }
   static real_t lambda_d( const real_t _omega, const real_t _magicNumber = threeSixteenth )
   {
      return ( real_t(4) - real_t(2) * _omega ) / ( real_t(4) * _magicNumber * _omega + real_t(2) - _omega );
   }
   static real_t magicNumber( const real_t _lambda_e, const real_t _lambda_d )
   {
      return ( ( real_t(2) - _lambda_e ) * ( real_t(2) - _lambda_d ) ) / ( real_t(4) * _lambda_e * _lambda_d );
   }

   real_t viscosity( const uint_t _level ) const
   {
      const real_t _lambda_e = levelDependentRelaxationParameter( _level, lambda_e_, level_ );
      return viscosityFromOmega( _lambda_e );
   }

   uint_t level() const { return level_; }

private:

   TRT() : lambda_e_( real_t(0) ), lambda_d_( real_t(0) ), magicNumber_( real_t(0) ), viscosity_( real_t(0) ), level_( uint_t(0) ) {}
   
   void initWithMagicNumber( const real_t _omega, const real_t _magicNumber, const uint_t _level )
   {
      lambda_e_    = lambda_e( _omega );
      lambda_d_    = lambda_d( _omega, _magicNumber );
      magicNumber_ = _magicNumber;
      viscosity_   = viscosityFromOmega( _omega );
      level_       = _level;
   }
   
   

   real_t lambda_e_;
   real_t lambda_d_;

   real_t magicNumber_; // "magic" number

   real_t viscosity_;

   uint_t level_;
};



//**********************************************************************************************************************
/*! Implementation of LB collision model MRT \cite dhumieres2002multiple
*   
*   Instead of providing all relaxation parameters, you can also provide just two. There are two variants, one that matches TRT \cite ginzburg2008two and one by Pan et al. \cite pan2006evaluation.
*/
//**********************************************************************************************************************

/// Only works with D3Q19!
class D3Q19MRT
{
public:

   typedef MRT_tag tag;
   
   static const real_t threeSixteenth;
   
   /// { 0, _s1, _s2,   0, _s4, 0, _s4, 0, _s4, _s9,  _s10, _s9,  _s10, _s9,  _s9,  _s9,  _s16, _s16, _s16 }
   /// { 0, s_e, s_eps, 0, s_q, 0, s_q, 0, s_q, s_nu, s_pi, s_nu, s_pi, s_nu, s_nu, s_nu, s_m,  s_m,  s_m  }
   D3Q19MRT( const real_t _s1, const real_t _s2, const real_t _s4, const real_t _s9, const real_t _s10, const real_t _s16,
             const uint_t _level = uint_t(0) ) :
      viscosity_( viscosityFromOmega( _s9 ) ), level_( _level )
   {
      s_[0]  = real_t(0);
      s_[1]  = _s1;
      s_[2]  = _s2;
      s_[3]  = real_t(0);
      s_[4]  = _s4;
      s_[5]  = real_t(0);
      s_[6]  = _s4;
      s_[7]  = real_t(0);
      s_[8]  = _s4;
      s_[9]  = _s9;
      s_[10] = _s10;
      s_[11] = _s9;
      s_[12] = _s10;
      s_[13] = _s9;
      s_[14] = _s9;
      s_[15] = _s9;
      s_[16] = _s16;
      s_[17] = _s16;
      s_[18] = _s16;
   }

   static D3Q19MRT constructPan( const real_t lambda_e, const real_t lambda_d, const uint_t _level = uint_t(0) )
   {
      D3Q19MRT mrt;
      mrt.initPan( lambda_e, lambda_d, _level );
      return mrt;
   }
   
   /// Model by Pan et al., An evaluation of lattice Boltzmann schemes for porous medium flow simulation. http://dx.doi.org/10.1016/j.compfluid.2005.03.008
   static D3Q19MRT constructPanWithMagicNumber( const real_t omega, const real_t magicNumber = threeSixteenth, const uint_t _level = uint_t(0) )
   {
      real_t lambda_e = TRT::lambda_e( omega );
      real_t lambda_d = TRT::lambda_d( omega, magicNumber );
      return constructPan( lambda_e, lambda_d, _level );
   }

   /// Supposed to be identical to TRT !
   static D3Q19MRT constructTRT( const real_t lambda_e, const real_t lambda_d, const uint_t _level = uint_t(0) )
   {
      D3Q19MRT mrt;
      mrt.initTRT( lambda_e, lambda_d, _level );
      return mrt;
   }
   static D3Q19MRT constructTRTWithMagicNumber( const real_t omega, const real_t magicNumber = threeSixteenth, const uint_t _level = uint_t(0) )
   {
      real_t lambda_e = TRT::lambda_e( omega );
      real_t lambda_d = TRT::lambda_d( omega, magicNumber );
      return constructTRT( lambda_e, lambda_d, _level );
   }

   void pack( mpi::SendBuffer & buffer ) const { for( uint_t i = uint_t(0); i < uint_t(19); ++i ) buffer << s_[i]; buffer << viscosity_ << level_; }
   void unpack( mpi::RecvBuffer & buffer ) { for( uint_t i = uint_t(0); i < uint_t(19); ++i ) buffer >> s_[i]; buffer >> viscosity_ >> level_; }

   /// Adapts the relaxation parameters to the "right", level-dependent parameters once "configure" is called
   void configure( IBlock & block, StructuredBlockStorage & sbs )
   {
      const uint_t _level = level_;
      level_ = sbs.getLevel( block );

      s_[1]  = levelDependentRelaxationParameter( level_, s_[1],  _level );
      s_[2]  = levelDependentRelaxationParameter( level_, s_[2],  _level );
      s_[4]  = levelDependentRelaxationParameter( level_, s_[4],  _level );
      s_[6]  = levelDependentRelaxationParameter( level_, s_[6],  _level );
      s_[8]  = levelDependentRelaxationParameter( level_, s_[8],  _level );
      s_[9]  = levelDependentRelaxationParameter( level_, s_[9],  _level );
      s_[10] = levelDependentRelaxationParameter( level_, s_[10], _level );
      s_[11] = levelDependentRelaxationParameter( level_, s_[11], _level );
      s_[12] = levelDependentRelaxationParameter( level_, s_[12], _level );
      s_[13] = levelDependentRelaxationParameter( level_, s_[13], _level );
      s_[14] = levelDependentRelaxationParameter( level_, s_[14], _level );
      s_[15] = levelDependentRelaxationParameter( level_, s_[15], _level );
      s_[16] = levelDependentRelaxationParameter( level_, s_[16], _level );
      s_[17] = levelDependentRelaxationParameter( level_, s_[17], _level );
      s_[18] = levelDependentRelaxationParameter( level_, s_[18], _level );
      
      viscosity_ = viscosityFromOmega( s_[9] );

      /*
      const uint_t _level = level_;
      level_ = sbs.getLevel( block );

      const real_t s_nu_scaled = levelDependentRelaxationParameter( level_, s_[9],  _level );
      const real_t linearFactor = s_nu_scaled / s_[9];

      s_[1]  = linearFactor * s_[1];
      s_[2]  = linearFactor * s_[2];
      s_[4]  = linearFactor * s_[4];
      s_[6]  = linearFactor * s_[6];
      s_[8]  = linearFactor * s_[8];
      s_[9]  = s_nu_scaled;
      s_[10] = linearFactor * s_[10];
      s_[11] = s_nu_scaled;
      s_[12] = linearFactor * s_[12];
      s_[13] = s_nu_scaled;
      s_[14] = s_nu_scaled;
      s_[15] = s_nu_scaled;
      s_[16] = linearFactor * s_[16];
      s_[17] = linearFactor * s_[17];
      s_[18] = linearFactor * s_[18];
      
      viscosity_ = viscosityFromOmega( s_[9] );
      */
   }

   real_t s0()  const { return real_t(0); }
   real_t s1()  const { return s_[1];     }
   real_t s2()  const { return s_[2];     }
   real_t s3()  const { return real_t(0); }
   real_t s4()  const { return s_[4];     }
   real_t s5()  const { return real_t(0); }
   real_t s6()  const { return s_[6];     }
   real_t s7()  const { return real_t(0); }
   real_t s8()  const { return s_[8];     }
   real_t s9()  const { return s_[9];     }
   real_t s10() const { return s_[10];    }
   real_t s11() const { return s_[11];    }
   real_t s12() const { return s_[12];    }
   real_t s13() const { return s_[13];    }
   real_t s14() const { return s_[14];    }
   real_t s15() const { return s_[15];    }
   real_t s16() const { return s_[16];    }
   real_t s17() const { return s_[17];    }
   real_t s18() const { return s_[18];    }

   real_t s( const uint_t index ) const
   {
      WALBERLA_ASSERT_LESS( index, 19 );
      return s_[ index ];
   }

   real_t s_e()   const { return s_[1];  }
   real_t s_eps() const { return s_[2];  }
   real_t s_q()   const { return s_[4];  }
   real_t s_nu()  const { return s_[9];  }
   real_t s_pi()  const { return s_[10]; }
   real_t s_m()   const { return s_[16]; }

   real_t omega() const { return s_[9]; }
   real_t omega( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                 const Vector3<real_t> & /*velocity*/ = Vector3<real_t>(), const real_t /*rho*/ = real_t(1) ) const { return omega(); }

   real_t omega_bulk() const { return s_[1]; }

   real_t viscosity() const { return viscosity_; }
   real_t viscosity( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/ ) const { return viscosity_; }

   real_t viscosity( const uint_t _level ) const
   {
      const real_t _s9  = levelDependentRelaxationParameter( _level, s_[9], level_ );
      return viscosityFromOmega( _s9 );
   }

   uint_t level() const { return level_; }

private:

   D3Q19MRT() : viscosity_( real_t(0) ), level_( uint_t(0) )
   {
      s_[0]  = real_t(0); s_[1]  = real_t(0); s_[2]  = real_t(0); s_[3]  = real_t(0); s_[4]  = real_t(0);
      s_[5]  = real_t(0); s_[6]  = real_t(0); s_[7]  = real_t(0); s_[8]  = real_t(0); s_[9]  = real_t(0);
      s_[10] = real_t(0); s_[11] = real_t(0); s_[12] = real_t(0); s_[13] = real_t(0); s_[14] = real_t(0);
      s_[15] = real_t(0); s_[16] = real_t(0); s_[17] = real_t(0); s_[18] = real_t(0);
   }

   void initTRT( const real_t lambda_e, const real_t lambda_d, const uint_t _level = uint_t(0) )
   {
      s_[0]  = real_t(0);
      s_[1]  = lambda_e;
      s_[2]  = lambda_e;
      s_[3]  = real_t(0);
      s_[4]  = lambda_d;
      s_[5]  = real_t(0);
      s_[6]  = lambda_d;
      s_[7]  = real_t(0);
      s_[8]  = lambda_d;
      s_[9]  = lambda_e;
      s_[10] = lambda_e;
      s_[11] = lambda_e;
      s_[12] = lambda_e;
      s_[13] = lambda_e;
      s_[14] = lambda_e;
      s_[15] = lambda_e;
      s_[16] = lambda_d;
      s_[17] = lambda_d;
      s_[18] = lambda_d;
      
      viscosity_ = viscosityFromOmega( lambda_e );
      level_     = _level;
   }
   
   void initPan( const real_t lambda_e, const real_t lambda_d, const uint_t _level = uint_t(0) )
   {
      s_[0]  = real_t(0);
      s_[1]  = lambda_d;
      s_[2]  = lambda_d;
      s_[3]  = real_t(0);
      s_[4]  = lambda_d;
      s_[5]  = real_t(0);
      s_[6]  = lambda_d;
      s_[7]  = real_t(0);
      s_[8]  = lambda_d;
      s_[9]  = lambda_e;
      s_[10] = lambda_d;
      s_[11] = lambda_e;
      s_[12] = lambda_d;
      s_[13] = lambda_e;
      s_[14] = lambda_e;
      s_[15] = lambda_e;
      s_[16] = lambda_d;
      s_[17] = lambda_d;
      s_[18] = lambda_d;
      
      viscosity_ = viscosityFromOmega( lambda_e );
      level_     = _level;
   }
   
   
   
   real_t s_[19];

   real_t viscosity_;

   uint_t level_;
};


//**********************************************************************************************************************
/*! Implementation of Cumulant LB collision operator \cite geier2015cumulant
*   
*   Only applicable for 3D simulation with D3Q27 stencil.
*
*   Here, omega[0] is the most important omega and is treated as omega1 from literature
*/
//**********************************************************************************************************************

class D3Q27Cumulant
{
public:
  
   typedef Cumulant_tag tag;

   /// Initializes all omegas to one except omega1
   D3Q27Cumulant( const real_t _omega1, const uint_t _level = uint_t(0) ) :
      viscosity_( viscosityFromOmega(_omega1) ), level_( _level )
   {
      omega_[0] = _omega1;
      for( uint_t i = uint_t(1); i < uint_t(10); ++i )
         omega_[i] = real_t(1);
   }

   /// Initializes all omegas separately
   D3Q27Cumulant( const real_t _omega1, const real_t _omega2, const real_t _omega3, const real_t _omega4, const real_t _omega5,
                  const real_t _omega6, const real_t _omega7, const real_t _omega8, const real_t _omega9, const real_t _omega10,
                  const uint_t _level = uint_t(0) ) :
      viscosity_( viscosityFromOmega(_omega1) ), level_( _level )
   {
      omega_[0]= _omega1;
      omega_[1]= _omega2;
      omega_[2]= _omega3;
      omega_[3]= _omega4;
      omega_[4]= _omega5;
      omega_[5]= _omega6;
      omega_[6]= _omega7;
      omega_[7]= _omega8;
      omega_[8]= _omega9;
      omega_[9]= _omega10;
   }

   void pack( mpi::SendBuffer &buffer ) const {
      for( uint_t i = uint_t(0); i < uint_t(10); ++i )
         buffer << omega_[i];
      buffer << viscosity_ << level_ ;
   }

   void unpack( mpi::RecvBuffer &buffer ) {
      for( uint_t i = uint_t(0); i < uint_t(10); ++i )
         buffer >> omega_[i] ;
      buffer >> viscosity_ >> level_ ;
   }

   /// Adapts only the main omega, i.e. omega0, to the "right", level-dependent parameters once "configure" is called
   void configure( IBlock &block, StructuredBlockStorage &sbs )
   {
      const uint_t _level = level_;
      level_ = sbs.getLevel( block );
      omega_[0] = levelDependentRelaxationParameter( level_, omega_[0], _level );
      viscosity_= viscosityFromOmega( omega_[0] );
   }

   real_t omega1()  const { return omega_[0]; }
   real_t omega2()  const { return omega_[1]; }
   real_t omega3()  const { return omega_[2]; }
   real_t omega4()  const { return omega_[3]; }
   real_t omega5()  const { return omega_[4]; }
   real_t omega6()  const { return omega_[5]; }
   real_t omega7()  const { return omega_[6]; }
   real_t omega8()  const { return omega_[7]; }
   real_t omega9()  const { return omega_[8]; }
   real_t omega10() const { return omega_[9]; }

   real_t omega() const { return omega_[0]; }
   real_t omega( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/,
                 const Vector3<real_t> & /*velocity*/ = Vector3<real_t>(), const real_t /*rho*/ = real_t(1) ) const { return omega(); }

   real_t omega( uint_t idx ) const { return omega_[idx]; }

   inline real_t viscosity() const { return viscosity_; }

   inline real_t viscosity( const uint_t _level ) const {
      const real_t omegaZero = levelDependentRelaxationParameter( _level, omega_[0], level_ );
      return viscosityFromOmega( omegaZero );
   }

   inline real_t viscosity( const cell_idx_t, const cell_idx_t, const cell_idx_t ) const { return viscosity_; }

   uint_t level() const { return level_; }

private:

   real_t omega_[10];
   real_t viscosity_;
   uint_t level_;
};
  
} // namespace collision_model
} // namespace lbm
} // namespace walberla
