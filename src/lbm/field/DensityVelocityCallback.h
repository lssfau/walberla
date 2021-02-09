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
//! \file DensityVelocityCallback.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/debug/Debug.h"
#include "core/math/Vector3.h"

#include "domain_decomposition/IBlock.h"

#include "field/AddToStorage.h"
#include "field/GhostLayerField.h"

#include <type_traits>



namespace walberla {
namespace lbm {



//template< typename LatticeModel_T > 
//struct VelocityField
//{
//   typedef walberla::field::GhostLayerField< Vector3<real_t>, uint_t(1) > type;
//};
//
//
//
//template< typename LatticeModel_T >
//BlockDataID addVelocityFieldToStorage( const shared_ptr< StructuredBlockStorage > & blocks, const std::string & identifier = std::string(),
//                                       const Vector3<real_t> & initValue = Vector3<real_t>(), uint_t nrOfGhostLayers = 1 )
//{
//   return field::addToStorage< typename VelocityField< LatticeModel_T >::type >( blocks, identifier, initValue, field::fzyx, nrOfGhostLayers );
//}


namespace internal {

template< typename LatticeModel_T, class Enable = void >
struct VelocityCallbackCorrection
{
   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
                                                     "For every lattice model, a fitting specialization of class 'lbm::internal::VelocityCallbackCorrection' is supposed to exist!\n"
                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                            LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t, const cell_idx_t, const cell_idx_t,
                      const LatticeModel_T & latticeModel, const real_t rho )
   {
      velocity += latticeModel.forceModel().force() * real_t(0.5) / rho;
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                            LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t, const cell_idx_t, const cell_idx_t,
                      const LatticeModel_T & latticeModel, const real_t )
   {
      velocity += latticeModel.forceModel().force() * real_t(0.5);
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                            ! LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                      const LatticeModel_T & latticeModel, const real_t rho )
   {
      velocity += latticeModel.forceModel().force(x,y,z) * real_t(0.5) / rho;
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                            ! LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                      const LatticeModel_T & latticeModel, const real_t )
   {
      velocity += latticeModel.forceModel().force(x,y,z) * real_t(0.5);
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                            LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t, const cell_idx_t, const cell_idx_t,
                      const LatticeModel_T & latticeModel, const real_t rho )
   {
      velocity -= latticeModel.forceModel().force() * real_t(0.5) / rho;
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                            LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t, const cell_idx_t, const cell_idx_t,
                      const LatticeModel_T & latticeModel, const real_t )
   {
      velocity -= latticeModel.forceModel().force() * real_t(0.5);
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
                                                                            ! LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                      const LatticeModel_T & latticeModel, const real_t rho )
   {
      velocity -= latticeModel.forceModel().force(x,y,z) * real_t(0.5) / rho;
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
                                                                            ! LatticeModel_T::ForceModel::constant &&
                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                      const LatticeModel_T & latticeModel, const real_t )
   {
      velocity -= latticeModel.forceModel().force(x,y,z) * real_t(0.5);
   }
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::ForceModel::shiftMacVel &&
                                                                            LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > &, const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T &, const real_t ) {}
};

template< typename LatticeModel_T >
struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::shiftMacVel &&
                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
                                                                            >::type >
{
   static void apply( Vector3< real_t > &, const cell_idx_t, const cell_idx_t, const cell_idx_t, const LatticeModel_T &, const real_t ) {}
};

} // namespace internal



//**********************************************************************************************************************
/*!
*   \section docDensityVelocityCalculation Density/Velocity Calculation
*
*   These functors are called after the streaming step in order to calculate the density and velocity which are
*   required for calculating the equilibrium distribution needed in the collision step.
*   The density and velocity returned by these functors are passed to a density/velocity callback functor
*   (see \ref docDensityVelocityCallback) prior to the collision step.
*   These callback functors always work on block local cell coordinates.
*
*   The concept for a density/velocity calculation functor looks like as follows (class with two member functions):
*
*   1. void operator()( IBlock & block )
*      -> called every time a new block is processed
*
*   2. real_t operator()( Vector3<real_t> & velocity, const PdfField_T * const field,
*                         const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
*      -> called every time a new cell is processed and the density and equilibrium velocity must be calculated for the
*         following collision step.
*         Must return the density via the return value of the function and the velocity via the reference parameter
*         'velocity'. The argument 'field' is the PDF field that is currently processed.
*         Might be called in parallel, i.e, must be threat-safe!
*/
//**********************************************************************************************************************



class DefaultDensityEquilibriumVelocityCalculation
{
public:

   void operator()( IBlock & ) const {}

   template< typename PdfField_T >
   real_t operator()( Vector3<real_t> & velocity, const PdfField_T * const field,
                      const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      return field->getDensityAndEquilibriumVelocity( velocity, x, y, z );
   }
};



template< typename VelocityField_T >
class AdvectionDiffusionDensityEquilibriumVelocityCalculation
{
public:

   AdvectionDiffusionDensityEquilibriumVelocityCalculation( const ConstBlockDataID & velocityFieldId ) :
      velocityFieldId_( velocityFieldId ), velocityField_( nullptr ) {}

   void operator()( IBlock & block )
   {
      velocityField_ = block.template getData< VelocityField_T >( velocityFieldId_ );
   }

   /// If you want to write an optimized specialization of an advection diffusion sweep, you probably
   /// need direct access to the external velocity.
   typename VelocityField_T::value_type getVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      return velocityField_->get( x, y, z );
   }

   template< typename PdfField_T >
   real_t operator()( Vector3<real_t> & velocity, const PdfField_T * const field,
                      const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
   {
      velocity = getVelocity( x, y, z );
      return field->getDensity( x, y, z );
   }

private:

   ConstBlockDataID velocityFieldId_;
   VelocityField_T const * velocityField_;

};



//**********************************************************************************************************************
/*!
*   \section docDensityVelocityCallback Density/Velocity Callback
*
*   These callback functors are called after the calculation of the density and equilibrium velocity. The density and
*   velocity passed to these callback functors are the same density and velocity that are used to calculate the
*   equilibrium distribution for the collision step.
*   These callback functors always work on block local cell coordinates.
*
*   The concept for a density/velocity callback functor looks like as follows (class with two member functions):
*
*   1. void operator()( IBlock & block )
*      -> called every time a new block is processed
*
*   2. void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & latticeModel,
*                       const Vector3<real_t> & velocity, const real_t rho )
*      -> called every time a new cell is processed and the density and equilibrium velocity have been calculated prior
*         to being used for performing the collision step of the LBM.
*         Might be called in parallel, i.e, must be threat-safe!
*/
//**********************************************************************************************************************



class DefaultDensityVelocityCallback
{
public:

   void operator()( IBlock & ) const {}

   template< typename LatticeModel_T >
   void operator()( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const LatticeModel_T & /*lm*/,
                    const Vector3<real_t> & /*velocity*/, const real_t /*rho*/ ) const
   {}
};



template< typename VelocityField_T, typename Enable=void >
class VelocityCallback
{
public:

   VelocityCallback( const BlockDataID & fieldId ) : fieldId_( fieldId ), field_( NULL ) {}

   void operator()( IBlock & block )
   {
      field_ = block.template getData< VelocityField_T >( fieldId_ );
   }

   template< typename LatticeModel_T >
   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & lm,
                    Vector3<real_t> velocity, const real_t rho )
   {
      static_assert( VelocityField_T::F_SIZE == 3, "Only valid for Fields with 3 components (F_SIZE==3)" );
      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
      internal::VelocityCallbackCorrection< LatticeModel_T >::apply( velocity, x, y, z, lm, rho );
      field_->get(x,y,z,0) = velocity[0];
      field_->get(x,y,z,1) = velocity[1];
      field_->get(x,y,z,2) = velocity[2];
   }

private:
   BlockDataID fieldId_;
   VelocityField_T * field_;
};


template< typename VelocityField_T  >
class VelocityCallback<VelocityField_T, typename std::enable_if< std::is_same< typename VelocityField_T::value_type, Vector3<real_t> >::value >::type >
{
public:
   VelocityCallback( const BlockDataID & fieldId ) : fieldId_( fieldId ), field_( NULL ) {}

   void operator()( IBlock & block )
   {
      field_ = block.template getData< VelocityField_T >( fieldId_ );
   }

   template< typename LatticeModel_T >
   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & lm,
                    Vector3<real_t> velocity, const real_t rho )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
      internal::VelocityCallbackCorrection< LatticeModel_T >::apply( velocity, x, y, z, lm, rho );
      field_->get(x,y,z) = velocity;
   }

private:
   BlockDataID fieldId_;
   VelocityField_T * field_;
};



template< typename DensityField_T >
class DensityCallback
{
public:

   DensityCallback( const BlockDataID & fieldId ) : fieldId_( fieldId ), field_( NULL ) {}

   void operator()( IBlock & block )
   {
      field_ = block.template getData< DensityField_T >( fieldId_ );
   }

   template< typename LatticeModel_T >
   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & /*lm*/,
                    Vector3<real_t> /*velocity*/, const real_t rho )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
      field_->get(x,y,z) = rho;
   }

private:

   BlockDataID fieldId_;
   DensityField_T * field_;

};



template< typename VelocityField_T, typename DensityField_T >
class DensityVelocityCallback
{
public:

   DensityVelocityCallback( const BlockDataID & velocityFieldId, const BlockDataID & densityFieldId ) :
      vId_( velocityFieldId ), dId_( densityFieldId ), vfield_( NULL ), dfield_( NULL ) {}

   void operator()( IBlock & block )
   {
      vfield_ = block.template getData< VelocityField_T >( vId_ );
      dfield_ = block.template getData< DensityField_T >( dId_ );
   }

   template< typename LatticeModel_T >
   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & lm,
                    Vector3<real_t> velocity, const real_t rho )
   {
      WALBERLA_ASSERT_NOT_NULLPTR( vfield_ );
      WALBERLA_ASSERT_NOT_NULLPTR( dfield_ );
      internal::VelocityCallbackCorrection< LatticeModel_T >::apply( velocity, x, y, z, lm, rho );
      vfield_->get(x,y,z) = velocity;
      dfield_->get(x,y,z) = rho;
   }

private:

   BlockDataID vId_;
   BlockDataID dId_;

   VelocityField_T * vfield_;
   DensityField_T * dfield_;

};



} // namespace lbm
} // namespace walberla



//namespace walberla {
//namespace lbm {
//
//
//
//
//template< typename LatticeModel_T, class Enable = void > 
//struct VelocityField
//{
//   typedef walberla::field::GhostLayerField< Vector3<real_t>, uint_t(1) > type;
//};
//
//template< typename LatticeModel_T > 
//struct VelocityField< LatticeModel_T, typename std::enable_if< LatticeModel_T::ForceModel::shiftMacVel >::type >
//{
//   typedef walberla::field::GhostLayerField< Vector3<real_t>, uint_t(2) > type;
//};
//
//
//
//template< typename LatticeModel_T >
//BlockDataID addVelocityFieldToStorage( const shared_ptr< StructuredBlockStorage > & blocks, const std::string & identifier = std::string(),
//                                       const Vector3<real_t> & initValue = Vector3<real_t>(), uint_t nrOfGhostLayers = 1 )
//{
//   return field::addToStorage< typename VelocityField< LatticeModel_T >::type >( blocks, identifier, initValue, field::fzyx, nrOfGhostLayers );
//}
//
//
//
//namespace internal {
//
//template< typename LatticeModel_T, class Enable = void >
//struct VelocityCallbackCorrection
//{
//   static_assert( never_true<LatticeModel_T>::value, "This static error message is never supposed to be triggered!\n"
//                                                     "For every lattice model, a fitting specialization of class 'lbm::internal::VelocityCallbackCorrection' is supposed to exist!\n"
//                                                     "If you see this message during compilation, please report to the developers of waLBerla." );
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
//                                                                            LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t rho )
//   {
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity + storage;
//      storage = latticeModel.forceModel().force() * real_t(0.5) / rho;
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
//                                                                            LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t )
//   {
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity + storage;
//      storage = latticeModel.forceModel().force() * real_t(0.5);
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
//                                                                            ! LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t rho )
//   {
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity + storage;
//      storage = latticeModel.forceModel().force(x,y,z) * real_t(0.5) / rho;      
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
//                                                                            ! LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t )
//   {
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity + storage;
//      storage = latticeModel.forceModel().force(x,y,z) * real_t(0.5);
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
//                                                                            LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t rho )
//   {
//      field->get(x,y,z) = velocity - latticeModel.forceModel().force() * real_t(0.5) / rho;
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
//                                                                            LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t )
//   {
//      field->get(x,y,z) = velocity - latticeModel.forceModel().force() * real_t(0.5);   
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
//                                                                            ! LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t rho )
//   {
//      field->get(x,y,z) = velocity - latticeModel.forceModel().force(x,y,z) * real_t(0.5) / rho;   
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
//                                                                            ! LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t )
//   {
//      field->get(x,y,z) = velocity - latticeModel.forceModel().force(x,y,z) * real_t(0.5);   
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
//                                                                            LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t rho )
//   {
//      const auto current = latticeModel.forceModel().force() * real_t(0.5) / rho;
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity - current + storage;
//      storage = current;
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
//                                                                            LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t )
//   {
//      const auto current = latticeModel.forceModel().force() * real_t(0.5);
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity - current + storage;
//      storage = current;  
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< LatticeModel_T::compressible &&
//                                                                            ! LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t rho )
//   {
//      const auto current = latticeModel.forceModel().force(x,y,z) * real_t(0.5) / rho;
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity - current + storage;
//      storage = current;
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::compressible &&
//                                                                            ! LatticeModel_T::ForceModel::constant &&
//                                                                            LatticeModel_T::ForceModel::shiftEquVel &&
//                                                                            LatticeModel_T::ForceModel::shiftMacVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
//                      const LatticeModel_T & latticeModel, const real_t )
//   {
//      const auto current = latticeModel.forceModel().force(x,y,z) * real_t(0.5);
//      auto & storage = field->get(x,y,z,1);
//      field->get(x,y,z) = velocity - current + storage;
//      storage = current;
//   }
//};
//
//template< typename LatticeModel_T >
//struct VelocityCallbackCorrection< LatticeModel_T, typename std::enable_if< ! LatticeModel_T::ForceModel::shiftMacVel &&
//                                                                            ! LatticeModel_T::ForceModel::shiftEquVel
//                                                                            >::type >
//{
//   template< typename VelocityField_T >
//   static void apply( VelocityField_T * const & field, const Vector3< real_t > &, const cell_idx_t, const cell_idx_t, const cell_idx_t,
//                      const LatticeModel_T &, const real_t ) {}
//};
//
//} // namespace internal
//
//
//
//class DefaultDensityVelocityCallback
//{
//public:
//
//   void operator()( IBlock & ) const {}
//
//   template< typename LatticeModel_T >
//   void operator()( const cell_idx_t /*x*/, const cell_idx_t /*y*/, const cell_idx_t /*z*/, const LatticeModel_T & /*lm*/,
//                    const Vector3<real_t> & /*velocity*/, const real_t /*rho*/ ) const
//   {}
//};
//
//
//
//template< typename VelocityField_T >
//class VelocityCallback
//{
//public:
//
//   VelocityCallback( const BlockDataID & fieldId ) : fieldId_( fieldId ), field_( NULL ) {}
//
//   void operator()( IBlock & block )
//   {
//      field_ = block.template getData< VelocityField_T >( fieldId_ );
//   }
//
//   template< typename LatticeModel_T >
//   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & lm,
//                    const Vector3<real_t> & velocity, const real_t rho )
//   {
//      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
//      internal::VelocityCallbackCorrection< LatticeModel_T >::apply( field_, velocity, x, y, z, lm, rho );
//   }
//
//private:
//
//   BlockDataID fieldId_;
//   VelocityField_T * field_;
//
//};
//
//
//
//template< typename DensityField_T >
//class DensityCallback
//{
//public:
//
//   DensityCallback( const BlockDataID & fieldId ) : fieldId_( fieldId ), field_( NULL ) {}
//
//   void operator()( IBlock & block )
//   {
//      field_ = block.template getData< DensityField_T >( fieldId_ );
//   }
//
//   template< typename LatticeModel_T >
//   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & /*lm*/,
//                    Vector3<real_t> /*velocity*/, const real_t rho )
//   {
//      WALBERLA_ASSERT_NOT_NULLPTR( field_ );
//      field_->get(x,y,z) = rho;
//   }
//
//private:
//
//   BlockDataID fieldId_;
//   DensityField_T * field_;
//
//};
//
//
//
//template< typename VelocityField_T, typename DensityField_T >
//class DensityVelocityCallback
//{
//public:
//
//   DensityVelocityCallback( const BlockDataID & velocityFieldId, const BlockDataID & densityFieldId ) :
//      vId_( velocityFieldId ), dId_( densityFieldId ), vfield_( NULL ), dfield_( NULL ) {}
//
//   void operator()( IBlock & block )
//   {
//      vfield_ = block.template getData< VelocityField_T >( vId_ );
//      dfield_ = block.template getData< DensityField_T >( dId_ );
//   }
//
//   template< typename LatticeModel_T >
//   void operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LatticeModel_T & lm,
//                    const Vector3<real_t> & velocity, const real_t rho )
//   {
//      WALBERLA_ASSERT_NOT_NULLPTR( vfield_ );
//      internal::VelocityCallbackCorrection< LatticeModel_T >::apply( vfield_, velocity, x, y, z, lm, rho );
//      WALBERLA_ASSERT_NOT_NULLPTR( dfield_ );
//      dfield_->get(x,y,z) = rho;
//   }
//
//private:
//
//   BlockDataID vId_;
//   BlockDataID dId_;
//
//   VelocityField_T * vfield_;
//   DensityField_T * dfield_;
//
//};
//
//
//
//} // namespace lbm
//} // namespace walberla
