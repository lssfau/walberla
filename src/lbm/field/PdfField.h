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
//! \file PdfField.h
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

#include "core/math/Vector3.h"
#include "core/math/Matrix3.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"


namespace walberla {
namespace lbm {



//**********************************************************************************************************************
/*!
*   \brief Specialized field class for PDF fields (fields containing multiple particle distribution functions per cell)
*
*   In addition to a generic GhostLayerField, each PdfField contains a lattice model which - among other things -
*   determines how macroscopic values like density and velocity have to be calculated. For this purpose, a set of member
*   functions exist that can be used to set single cells (or the whole field) to equilibrium or to evaluate the density
*   and/or velocity for specific positions. If you happen to have an iterator to a PdfField, you should use the free
*   functions located in "MacroscopicValueCalculation.h" - you do not have to convert your iterator into x-, y-, and
*   z-coordinates!
*   Also, particle distribution functions (i.e., the values stored in the field) can be accessed using stencil
*   directions, e.g. "pdfField( x, y, z, stencil::NE )".
*
*   Note that the behavior is different for compressible and incompressible lattice models.
*   In the compressible case, the behavior and formulas are as expected from common LBM literature and the density
*   is taken as the 0-th order moment of the PDFs and this value is used.
*   In order to make LBM quasi-incompressible, it was suggested, e.g. in
*     Q. Zou, S. Hou, S. Chen, G.D. Doolen, J. Stat. Phys. 81(1–2), 35 (1995)
*     X. He, L.S. Luo, J. Stat. Phys. 88(3–4), 927 (1997)
*   that the density is implicitly assumed to have a constant value of 1 in many cases (e.g. in getVelocity()),
*   and only the deviation from this value enters some of the formulas, like the equilibrium distribution functions.
*   Additionally, the PDFs are normalized around 0 in the incompressible case to increase the numerical accuracy,
*   i.e. only the deviation of the PDF values from their respective lattice weight is stored.
*   As a result, manually summing of the PDF values will yield the density deviation in this case.
*   But the getDensity() function reverts this normalization (by adding 1) and will yield the physical density.
*   This normalization, however, usually doesn't affect the implementation of functions like LBM sweeps.
*
*/
//**********************************************************************************************************************

template< typename LatticeModel_T >
class PdfField : public GhostLayerField< real_t, LatticeModel_T::Stencil::Size >
{
public:

   //** Type Definitions  **********************************************************************************************
   /*! \name Type Definitions */
   //@{
   using LatticeModel = LatticeModel_T;
   using Stencil = typename LatticeModel_T::Stencil;

   using value_type = typename GhostLayerField<real_t, Stencil::Size>::value_type;

   using iterator = typename GhostLayerField<real_t, Stencil::Size>::iterator;
   using const_iterator = typename GhostLayerField<real_t, Stencil::Size>::const_iterator;

   using reverse_iterator = typename GhostLayerField<real_t, Stencil::Size>::reverse_iterator;
   using const_reverse_iterator = typename GhostLayerField<real_t, Stencil::Size>::const_reverse_iterator;

   using base_iterator = typename GhostLayerField<real_t, Stencil::Size>::base_iterator;
   using const_base_iterator = typename GhostLayerField<real_t, Stencil::Size>::const_base_iterator;

   using Ptr = typename GhostLayerField<real_t, Stencil::Size>::Ptr;
   using ConstPtr = typename GhostLayerField<real_t, Stencil::Size>::ConstPtr;
   //@}
   //*******************************************************************************************************************

   PdfField( const uint_t _xSize, const uint_t _ySize, const uint_t _zSize,
             const LatticeModel_T & _latticeModel,
             const bool initialize = true, const Vector3< real_t > & initialVelocity = Vector3< real_t >( real_t(0.0) ),
             const real_t initialDensity = real_t(1.0),
             const uint_t ghostLayers = uint_t(1), const field::Layout & _layout = field::fzyx,
             const shared_ptr< field::FieldAllocator<real_t> > & alloc = shared_ptr< field::FieldAllocator<real_t> >() );

   ~PdfField() override = default;



   //inline bool operator==( const PdfField & rhs ) const; // TODO! -> ticket
   //inline bool operator!=( const PdfField & rhs ) const { return !operator==( rhs ); }

   inline PdfField * clone()              const;
   inline PdfField * cloneUninitialized() const;
   inline PdfField * cloneShallowCopy()   const;

   const LatticeModel_T & latticeModel() const { return latticeModel_; }
         LatticeModel_T & latticeModel()       { return latticeModel_; }

   void resetLatticeModel( const LatticeModel_T & lm ) { latticeModel_ = lm; }

   /////////////////////////////////////////////////
   // Access functions (with stencil::Direction!) //
   /////////////////////////////////////////////////

   using GhostLayerField< real_t, Stencil::Size >::get;

         real_t & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )       { return get( x, y, z, Stencil::idx[d] ); }
   const real_t & get( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const { return get( x, y, z, Stencil::idx[d] ); }
         real_t & get( const Cell & c, stencil::Direction d )       { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }
   const real_t & get( const Cell & c, stencil::Direction d ) const { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }

   using GhostLayerField< real_t, Stencil::Size >::operator();

         real_t & operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )       { return get( x, y, z, Stencil::idx[d] ); }
   const real_t & operator()( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const { return get( x, y, z, Stencil::idx[d] ); }
         real_t & operator()( const Cell & c, stencil::Direction d )       { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }
   const real_t & operator()( const Cell & c, stencil::Direction d ) const { return get( c.x(), c.y(), c.z(), Stencil::idx[d] ); }

   //////////////////////////////
   // set density and velocity //
   //////////////////////////////

   inline void setDensityAndVelocity( const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );

   inline void setDensityAndVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );
   inline void setDensityAndVelocity( const Cell & cell,
                                      const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );

   /////////////////////
   // set equilibrium //
   /////////////////////

   inline void setToEquilibrium( const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );

   inline void setToEquilibrium( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                 const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );
   inline void setToEquilibrium( const Cell & cell,
                                 const Vector3< real_t > & velocity = Vector3< real_t >( real_t(0.0) ), const real_t rho = real_t(1.0) );

   /////////////////
   // get density //
   /////////////////

   inline real_t getDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline real_t getDensity( const Cell & cell ) const;

   /////////////////////////////
   // get density in SI units //
   /////////////////////////////

   inline real_t getDensitySI( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const real_t rho_SI ) const;
   inline real_t getDensitySI( const Cell & cell,                                          const real_t rho_SI ) const;

   //////////////////////////
   // get momentum density //
   //////////////////////////

   inline Vector3< real_t > getMomentumDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline Vector3< real_t > getMomentumDensity( const Cell & cell ) const;

   inline void getMomentumDensity( Vector3< real_t > & momentumDensity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline void getMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const;

   inline Vector3< real_t > getEquilibriumMomentumDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline Vector3< real_t > getEquilibriumMomentumDensity( const Cell & cell ) const;

   inline void getEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline void getEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const;

   //////////////////
   // get velocity //
   //////////////////

   inline Vector3< real_t > getVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline Vector3< real_t > getVelocity( const Cell & cell ) const;

   inline void getVelocity( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline void getVelocity( Vector3< real_t > & velocity, const Cell & cell ) const;

   inline Vector3< real_t > getEquilibriumVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline Vector3< real_t > getEquilibriumVelocity( const Cell & cell ) const;

   inline void getEquilibriumVelocity( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline void getEquilibriumVelocity( Vector3< real_t > & velocity, const Cell & cell ) const;

   //////////////////////////////
   // get velocity in SI units //
   //////////////////////////////

   inline Vector3< real_t > getVelocitySI( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                           const real_t dx_SI, const real_t dt_SI ) const;
   inline Vector3< real_t > getVelocitySI( const Cell & cell, const real_t dx_SI, const real_t dt_SI ) const;

   inline void getVelocitySI( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                            const real_t dx_SI, const real_t dt_SI ) const;
   inline void getVelocitySI( Vector3< real_t > & velocity, const Cell & cell, const real_t dx_SI, const real_t dt_SI ) const;

   inline Vector3< real_t > getVelocitySI( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const real_t dxDividedByDt_SI ) const;
   inline Vector3< real_t > getVelocitySI( const Cell & cell,                                          const real_t dxDividedByDt_SI ) const;

   inline void getVelocitySI( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const real_t dxDividedByDt_SI ) const;
   inline void getVelocitySI( Vector3< real_t > & velocity, const Cell & cell,                                          const real_t dxDividedByDt_SI ) const;

   //////////////////////////////////////
   // get density and momentum density //
   //////////////////////////////////////

   inline real_t getDensityAndMomentumDensity( Vector3< real_t > & momentumDensity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline real_t getDensityAndMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const;

   inline real_t getDensityAndEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline real_t getDensityAndEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const;

   //////////////////////////////
   // get density and velocity //
   //////////////////////////////

   inline real_t getDensityAndVelocity( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline real_t getDensityAndVelocity( Vector3< real_t > & velocity, const Cell & cell ) const;

   inline real_t getDensityAndEquilibriumVelocity( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline real_t getDensityAndEquilibriumVelocity( Vector3< real_t > & velocity, const Cell & cell ) const;

   //////////////////////////////////////////
   // get density and velocity in SI units //
   //////////////////////////////////////////

   inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                          const real_t rho_SI, const real_t dx_SI, const real_t dt_SI ) const;
   inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const Cell & cell,
                                          const real_t rho_SI, const real_t dx_SI, const real_t dt_SI ) const;

   inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                          const real_t rho_SI, const real_t dxDividedByDt_SI ) const;
   inline real_t getDensityAndVelocitySI( Vector3< real_t > & velocity, const Cell & cell,
                                          const real_t rho_SI, const real_t dxDividedByDt_SI ) const;

   ////////////////////
   // get shear rate //
   ////////////////////

   inline real_t getShearRate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline real_t getShearRate( const Cell & cell ) const;

   /////////////////////////
   // get pressure tensor //
   /////////////////////////

   inline Matrix3< real_t > getPressureTensor( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline Matrix3< real_t > getPressureTensor( const Cell & cell ) const;

   inline void getPressureTensor( Matrix3< real_t > & pressureTensor, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const;
   inline void getPressureTensor( Matrix3< real_t > & pressureTensor, const Cell & cell ) const;


protected:

   //** Shallow Copy ***************************************************************************************************
   /*! \name Shallow Copy */
   //@{
   inline PdfField( const PdfField< LatticeModel_T > & other );
   Field< real_t, Stencil::Size > * cloneShallowCopyInternal() const override { return new PdfField< LatticeModel_T >( *this ); }
   //@}
   //*******************************************************************************************************************

   LatticeModel_T latticeModel_;
};



template< typename LatticeModel_T >
PdfField< LatticeModel_T >::PdfField( const uint_t _xSize, const uint_t _ySize, const uint_t _zSize,
                                      const LatticeModel_T & _latticeModel,
                                      const bool initialize, const Vector3< real_t > & initialVelocity, const real_t initialDensity,
                                      const uint_t ghostLayers, const field::Layout & _layout,
                                      const shared_ptr< field::FieldAllocator<real_t> > & alloc ) :

   GhostLayerField< real_t, Stencil::Size >( _xSize, _ySize, _zSize, ghostLayers, _layout, alloc ),
   latticeModel_( _latticeModel )
{
#ifdef _OPENMP
   // take care of proper thread<->memory assignment (first-touch allocation policy !)
   this->setWithGhostLayer( real_t(0) );
#endif

   if( initialize )
      setDensityAndVelocity( initialVelocity, initialDensity );
}



template< typename LatticeModel_T >
inline PdfField< LatticeModel_T > * PdfField< LatticeModel_T >::clone() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< real_t, Stencil::Size >::clone() );
}

template< typename LatticeModel_T >
inline PdfField< LatticeModel_T > * PdfField< LatticeModel_T >::cloneUninitialized() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< real_t, Stencil::Size >::cloneUninitialized() );
}

template< typename LatticeModel_T >
inline PdfField< LatticeModel_T > * PdfField< LatticeModel_T >::cloneShallowCopy() const
{
   return dynamic_cast< PdfField * >( GhostLayerField< real_t, Stencil::Size >::cloneShallowCopy() );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::setDensityAndVelocity( const Vector3< real_t > & velocity, const real_t rho )
{
   auto beginIterator = this->beginWithGhostLayerXYZ();
   DensityAndVelocityRange< LatticeModel_T, iterator >::set( beginIterator, this->end(), latticeModel_, velocity, rho );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::setDensityAndVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                               const Vector3< real_t > & velocity, const real_t rho )
{
   DensityAndVelocity< LatticeModel_T >::set( *this, x, y, z, latticeModel_, velocity, rho );
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::setDensityAndVelocity( const Cell & cell, const Vector3< real_t > & velocity, const real_t rho )
{
   setDensityAndVelocity( cell.x(), cell.y(), cell.z(), velocity, rho );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::setToEquilibrium( const Vector3< real_t > & velocity, const real_t rho )
{
   auto beginIterator = this->beginWithGhostLayerXYZ();
   EquilibriumRange< LatticeModel_T, iterator >::set( beginIterator, this->end(), velocity, rho );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::setToEquilibrium( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                          const Vector3< real_t > & velocity, const real_t rho )
{
   Equilibrium< LatticeModel_T >::set( *this, x, y, z, velocity, rho );
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::setToEquilibrium( const Cell & cell, const Vector3< real_t > & velocity, const real_t rho )
{
   setToEquilibrium( cell.x(), cell.y(), cell.z(), velocity, rho );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return Density< LatticeModel_T >::get( latticeModel_, *this, x, y, z );
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensity( const Cell & cell ) const
{
   return getDensity( cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensitySI( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const real_t rho_SI ) const
{
   return getDensity(x,y,z) * rho_SI;
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensitySI( const Cell & cell, const real_t rho_SI ) const
{
   return getDensitySI( cell.x(), cell.y(), cell.z(), rho_SI );
}



template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getMomentumDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   Vector3< real_t > momentumDensity;
   getMomentumDensity( momentumDensity, x, y, z );
   return momentumDensity;
}

template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getMomentumDensity( const Cell & cell ) const
{
   return getMomentumDensity( cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getMomentumDensity( Vector3< real_t > & momentumDensity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   MomentumDensity< LatticeModel_T >::get( momentumDensity, latticeModel_, *this, x, y, z );
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const
{
   getMomentumDensity( momentumDensity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getEquilibriumMomentumDensity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   Vector3< real_t > momentumDensity;
   getEquilibriumMomentumDensity( momentumDensity, x, y, z );
   return momentumDensity;
}

template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getEquilibriumMomentumDensity( const Cell & cell ) const
{
   return getEquilibriumMomentumDensity( cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   MomentumDensity< LatticeModel_T >::getEquilibrium( momentumDensity, latticeModel_, *this, x, y, z );
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const
{
   getEquilibriumMomentumDensity( momentumDensity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   Vector3< real_t > velocity;
   getVelocity( velocity, x, y, z );
   return velocity;
}

template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getVelocity( const Cell & cell ) const
{
   return getVelocity( cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getVelocity( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   if( LatticeModel_T::compressible )
   {
      const real_t rho = DensityAndMomentumDensity< LatticeModel_T >::get( velocity, latticeModel_, *this, x, y, z );
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   else
   {
      MomentumDensity< LatticeModel_T >::get( velocity, latticeModel_, *this, x, y, z );
   }
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getVelocity( Vector3< real_t > & velocity, const Cell & cell ) const
{
   getVelocity( velocity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getEquilibriumVelocity( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   Vector3< real_t > velocity;
   getEquilibriumVelocity( velocity, x, y, z );
   return velocity;
}

template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getEquilibriumVelocity( const Cell & cell ) const
{
   return getEquilibriumVelocity( cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getEquilibriumVelocity( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   if( LatticeModel_T::compressible )
   {
      const real_t rho = DensityAndMomentumDensity< LatticeModel_T >::getEquilibrium( velocity, latticeModel_, *this, x, y, z );
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   else
   {
      MomentumDensity< LatticeModel_T >::getEquilibrium( velocity, latticeModel_, *this, x, y, z );
   }
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getEquilibriumVelocity( Vector3< real_t > & velocity, const Cell & cell ) const
{
   getEquilibriumVelocity( velocity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getVelocitySI( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                    const real_t dx_SI, const real_t dt_SI ) const
{
   Vector3< real_t > velocity;
   getVelocitySI( velocity, x, y, z, dx_SI, dt_SI );
   return velocity;
}

template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getVelocitySI( const Cell & cell, const real_t dx_SI, const real_t dt_SI ) const
{
   return getVelocitySI( cell.x(), cell.y(), cell.z(), dx_SI, dt_SI );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getVelocitySI( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                       const real_t dx_SI, const real_t dt_SI ) const
{
   getVelocitySI( velocity, x, y, z, dx_SI / dt_SI );
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getVelocitySI( Vector3< real_t > & velocity, const Cell & cell, const real_t dx_SI, const real_t dt_SI ) const
{
   getVelocitySI( velocity, cell.x(), cell.y(), cell.z(), dx_SI, dt_SI );
}



template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getVelocitySI( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                    const real_t dxDividedByDt_SI ) const
{
   Vector3< real_t > velocity;
   getVelocitySI( velocity, x, y, z, dxDividedByDt_SI );
   return velocity;
}

template< typename LatticeModel_T >
inline Vector3< real_t > PdfField< LatticeModel_T >::getVelocitySI( const Cell & cell, const real_t dxDividedByDt_SI ) const
{
   return getVelocitySI( cell.x(), cell.y(), cell.z(), dxDividedByDt_SI );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getVelocitySI( Vector3< real_t > & velocity, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                       const real_t dxDividedByDt_SI ) const
{
   getVelocity( velocity, x, y, z );
   velocity *= dxDividedByDt_SI;
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getVelocitySI( Vector3< real_t > & velocity, const Cell & cell, const real_t dxDividedByDt_SI ) const
{
   getVelocitySI( velocity, cell.x(), cell.y(), cell.z(), dxDividedByDt_SI );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndMomentumDensity( Vector3< real_t > & momentumDensity,
                                                                        const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return DensityAndMomentumDensity< LatticeModel_T >::get( momentumDensity, latticeModel_, *this, x, y, z );
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const
{
   return getDensityAndMomentumDensity( momentumDensity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity,
                                                                                   const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   return DensityAndMomentumDensity< LatticeModel_T >::getEquilibrium( momentumDensity, latticeModel_, *this, x, y, z );
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndEquilibriumMomentumDensity( Vector3< real_t > & momentumDensity, const Cell & cell ) const
{
   return getDensityAndEquilibriumMomentumDensity( momentumDensity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndVelocity( Vector3< real_t > & velocity,
                                                                 const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   const real_t rho = DensityAndMomentumDensity< LatticeModel_T >::get( velocity, latticeModel_, *this, x, y, z );
   if( LatticeModel_T::compressible )
   {
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   return rho;
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndVelocity( Vector3< real_t > & velocity, const Cell & cell ) const
{
   return getDensityAndVelocity( velocity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndEquilibriumVelocity( Vector3< real_t > & velocity,
                                                                            const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   const real_t rho = DensityAndMomentumDensity< LatticeModel_T >::getEquilibrium( velocity, latticeModel_, *this, x, y, z );
   if( LatticeModel_T::compressible )
   {
      const real_t invRho = real_t(1.0) / rho;
      velocity *= invRho;
   }
   return rho;
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndEquilibriumVelocity( Vector3< real_t > & velocity, const Cell & cell ) const
{
   return getDensityAndEquilibriumVelocity( velocity, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndVelocitySI( Vector3< real_t > & velocity,
                                                                   const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                   const real_t rho_SI, const real_t dx_SI, const real_t dt_SI ) const
{
   return getDensityAndVelocitySI( velocity, x, y, z, rho_SI, dx_SI / dt_SI );
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndVelocitySI( Vector3< real_t > & velocity, const Cell & cell,
                                                                   const real_t rho_SI, const real_t dx_SI, const real_t dt_SI ) const
{
   return getDensityAndVelocitySI( velocity, cell.x(), cell.y(), cell.z(), rho_SI, dx_SI, dt_SI );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndVelocitySI( Vector3< real_t > & velocity,
                                                                   const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                                                                   const real_t rho_SI, const real_t dxDividedByDt_SI ) const
{
   const real_t rho = getDensityAndVelocity( velocity, x, y, z );
   velocity *= dxDividedByDt_SI;
   return rho * rho_SI;
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getDensityAndVelocitySI( Vector3< real_t > & velocity, const Cell & cell,
                                                                   const real_t rho_SI, const real_t dxDividedByDt_SI ) const
{
   return getDensityAndVelocitySI( velocity, cell.x(), cell.y(), cell.z(), rho_SI, dxDividedByDt_SI );
}



template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getShearRate( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   Vector3< real_t > velocity;
   const real_t rho = getDensityAndVelocity( velocity, x, y, z );

   return ShearRate< LatticeModel_T >::get( latticeModel_, *this, x, y, z, velocity, rho );
}

template< typename LatticeModel_T >
inline real_t PdfField< LatticeModel_T >::getShearRate( const Cell & cell ) const
{
   return getShearRate( cell.x(), cell.y(), cell.z() );
}

template< typename LatticeModel_T >
inline Matrix3< real_t > PdfField< LatticeModel_T >::getPressureTensor( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   Matrix3< real_t > pressureTensor;
   getPressureTensor( pressureTensor, x, y, z );
   return pressureTensor;
}

template< typename LatticeModel_T >
inline Matrix3< real_t > PdfField< LatticeModel_T >::getPressureTensor( const Cell & cell ) const
{
   return getPressureTensor( cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getPressureTensor( Matrix3< real_t > & pressureTensor, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const
{
   PressureTensor<LatticeModel_T>::get( pressureTensor, latticeModel_, *this, x, y, z );
}

template< typename LatticeModel_T >
inline void PdfField< LatticeModel_T >::getPressureTensor( Matrix3< real_t > & pressureTensor, const Cell & cell ) const
{
   getPressureTensor( pressureTensor, cell.x(), cell.y(), cell.z() );
}



template< typename LatticeModel_T >
inline PdfField< LatticeModel_T >::PdfField( const PdfField< LatticeModel_T > & other )
   : GhostLayerField< real_t, Stencil::Size >::GhostLayerField( other ),
     latticeModel_( other.latticeModel_ )
{
}



} // namespace lbm
} // namespace walberla
