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
//! \file AdaptorIterator.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "field/iterators/FieldIterator.h"


namespace walberla {
namespace field {


//**********************************************************************************************************************
/*! Iterator class for Adaptors
*
*  Implements same functions as walberla::field::FieldIterator by wrapping a FieldIterator object.
*  When dereferenced the value is first transformed by a functor object provided by the Adaptor.
*
*  Template Parameters
*     AdaptedIterator: The iterator type of the base field
*     FieldAdaptor:    The Field adaptor class
*/
//**********************************************************************************************************************
template < typename AdaptedIterator, typename FieldAdaptor>
class AdaptorIterator
{
public:
   using T = typename FieldAdaptor::value_type;
   using Functor = typename FieldAdaptor::functor_t;
   using OwnType = AdaptorIterator<AdaptedIterator, FieldAdaptor>;

   AdaptorIterator( const AdaptedIterator & baseIterator, const FieldAdaptor * adaptedField )
      : base_( baseIterator ), adaptedField_( adaptedField), functor_( adaptedField->getFunctor() )
   {}

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
   inline bool operator==( const OwnType & it ) const { return base_ == it.base_; }
   inline bool operator!=( const OwnType & it ) const { return base_ != it.base_; }
   //@}
   //*******************************************************************************************************************

   //**Access Functions*************************************************************************************************
   /*!\name Access Functions */
   //@{

   inline T operator*()    const { return functor_( base_ ); }
   inline T operator->()   const { return functor_( base_ ); }

   inline T getF( cell_idx_t cf ) const { return functor_( adaptedField_->getAdaptedField(), x(), y(), z(), f() + cf ); }
   inline T getF( uint_t cf )     const { return functor_( adaptedField_->getAdaptedField(), x(), y(), z(), f() + cell_idx_c(cf) ); }

   inline T neighbor( stencil::Direction d, cell_idx_t cf = 0 ) const {
      using namespace stencil;
      return functor_( x() + cx[d], y() + cy[d], z() + cz[d], f() + cf );
   }

   inline T neighbor( stencil::Direction d, uint_t cf )  const {
      return neighbor( d, cell_idx_c( cf) );
   }

   inline T neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, cell_idx_t cf = 0 ) const {
      return functor_( adaptedField_->getAdaptedField(), x() + cx, y() + cy, z() + cz, f() + cf );
   }

   inline T neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, uint_t cf ) const {
      return neighbor( cx,cy,cz, cell_idx_c( cf ) );
   }

   //@}
   //*******************************************************************************************************************


   //** Coordinates of current position ********************************************************************************
   /*! \name Coordinates of current position */
   //@{
   inline cell_idx_t x() const { return base_.x(); }
   inline cell_idx_t y() const { return base_.y(); }
   inline cell_idx_t z() const { return base_.z(); }
   inline cell_idx_t f() const { return base_.f(); }

   inline Cell cell()    const { return base_.cell(); }
   //@}
   //*******************************************************************************************************************


   //**Utility Functions ***********************************************************************************************
   /*!\name Utility Functions */
   //@{
   void print(std::ostream & str) const { str << "Adapted Iterator of "; base_.print( str );  }

   const FieldAdaptor * getField() const   { return adaptedField_; }
   //@}
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*!\name Operators */
   //@{
   inline AdaptorIterator& operator++() { ++base_; return *this; }
   inline AdaptorIterator& operator--() { --base_; return *this; }
   //@}
   //*******************************************************************************************************************

   //**Fast Iteration **************************************************************************************************
   /*!\name Fast Iteration */
   //@{
   inline void incrOuter()        { base_.incrOuter(); }
   inline void incrInner()        { base_.incrInner(); }
   inline bool testInner() const  { return base_.testInner(); }
   //@}
   //*******************************************************************************************************************


protected:
   AdaptedIterator      base_;
   const FieldAdaptor * adaptedField_;
   const Functor  &     functor_;
};




} // namespace field
} // namespace walberla


