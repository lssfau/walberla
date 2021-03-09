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
//! \file FieldIterator.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "IteratorMacros.h"

#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/debug/Debug.h"
#include "field/Layout.h"

#include "stencil/Directions.h"

#include <iostream>
#include <iterator>
#include <type_traits>

namespace walberla {
namespace field {


   template<typename T, uint_t fSize_> class Field; // forward for friend declaration



   //*******************************************************************************************************************
   /*!
    * \brief Iterator over a 4D array, stored linearized in memory
    *
    * \ingroup field
    *
    *  The iterator is created by specifying a field slice.
    *  It supports neighbor access and knows its own x,y,z,f coordinates.
    *
    *  This is NOT a lightweight object as most iterators are, since it has to store
    *  begin and end for each of the four dimensions plus some data to speed up the calculation
    *  of current coordinate.
    *
    *  The FieldIterator itself cannot change the position it points to. The ForwardFieldIterator
    *  and ReverseFieldIterator are derived from this class and provide the according operators.
    *  When writing a function that calculates something in a specific cell should
    *  take a FieldIterator as argument ( not Forward- or ReverseFieldIterators ) so that it can be
    *  called with both Forward and Reverse iterators.
    *  For this case there is also a more generic way: see walberla::field::FieldPointer
    *
    * \image html field/doc/FieldIterator.png "Implementation"
    */
   //*******************************************************************************************************************
   template <typename T, uint_t fieldFSize>
   class FieldIterator
   {
   public:
      using iterator_category = std::forward_iterator_tag;
      using value_type = T;
      using difference_type = std::ptrdiff_t;
      using pointer = T*;
      using reference = T&;

      using FieldType = Field<typename std::remove_const<T>::type, fieldFSize>;

      static const uint_t F_SIZE = fieldFSize;

      //** Copy Operations *********************************************************************************************
      /*!\name Copy Operations */
      //@{
      FieldIterator                         ( const FieldIterator<T,fieldFSize> & other );
      FieldIterator<T,fieldFSize>& operator=( const FieldIterator<T,fieldFSize> & other );
      //@}
      //****************************************************************************************************************


      //**Operators*****************************************************************************************************
      /*!\name Operators */
      //@{
      inline bool           operator==( const FieldIterator& it ) const;
      inline bool           operator!=( const FieldIterator& it ) const;

      operator const FieldIterator<const T, fieldFSize> & () const {
         const FieldIterator<const T, fieldFSize> * ptr;
         ptr = reinterpret_cast< const FieldIterator<const T, fieldFSize>* > ( this );
         return *ptr;
      }

      //@}
      //****************************************************************************************************************


      //**Access Functions**********************************************************************************************
      /*!\name Access Functions */
      //@{
      inline T & operator*()    const { return *linePtr_; }
      inline T * operator->()   const { return  linePtr_; }

      inline T & getF( cell_idx_t cf ) const;
      inline T & getF( uint_t     cf ) const;

      inline T & operator[] ( cell_idx_t cf ) const { return getF(cf); }
      inline T & operator[] ( uint_t     cf ) const { return getF(cf); }

      inline T & neighbor( stencil::Direction d, cell_idx_t cf = 0 ) const;
      inline T & neighbor( stencil::Direction d, uint_t cf ) const;
      inline T & neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, cell_idx_t cf = 0 ) const;
      inline T & neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, uint_t cf ) const;

      //@}
      //****************************************************************************************************************


      //** Coordinates of current position ****************************************************************************
      /*! \name Coordinates of current position */
      //@{
      inline cell_idx_t x() const;
      inline cell_idx_t y() const;
      inline cell_idx_t z() const;
      inline cell_idx_t f() const;

      inline Cell cell()    const;
      //@}
      //****************************************************************************************************************


      //**Utility Functions ********************************************************************************************
      /*!\name Utility Functions */
      //@{
      void print(std::ostream & str) const;
      const FieldType * getField() const { return f_; }
      //@}
      //****************************************************************************************************************

   protected:

      //**Constructor/Destructor****************************************************************************************
      /*!\name Constructor/Destructor */
      //@{

      explicit FieldIterator();

      explicit FieldIterator(const FieldType * f,
                             cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                             uint_t xSize, uint_t ySize, uint_t zSize, uint_t fSize, bool forward );

      //@}
      //****************************************************************************************************************


      friend class Field<typename std::remove_const<T>::type, fieldFSize>;

      void incrementLine();
      void decrementLine();

      T * lineBegin_;  ///< Points to begin of fastest coordinate line
      T * linePtr_;    ///< Point to current element
      T * lineEnd_;    ///< Points to end of current line

      /// In the following vectors [0] is the slowest and [3] the fastest coordinate
      /// Current values of the coordinates, forth coordinate implicitly stored in linePtr_
      /// and if needed written to fastestCoord_
      cell_idx_t cur_[3];

      /// Number of elements to skip when coordinate wraps around
      uint_t  skips_[4];

      /// Size of each coordinate
      uint_t  sizes_[4];

      /// Field where iterator belongs to
      const FieldType * f_;

      /// Following offset values are only used in functions x(),y(),z(),f()
      cell_idx_t xBegin_;
      cell_idx_t yBegin_;
      cell_idx_t zBegin_;
      cell_idx_t fBegin_;


      //** Speeding up x(), y(), z(), f() functions *******************************************************************
      /**
      \name Speeding up x(), y(), z(), f() functions
      The following members are only needed to speed up the x(), y(), z(), f() functions.
      In order to get x(), y(), z(), f() function as fast a possible, no if clause for the layout was introduced.
      Instead there are the cur[XYZF]_ members, that point to the cur_ array. The cur_ array does not store the
      fastest coordinate, because it is implicitly stored in (linePtr_ - lineBegin_).  If it would be stored explicitly
      there would have to be an extra update operation in operator++() , which should be as fast as possible.
      The curX_ or curF_ pointer points to the fastestCoord_ member, which always has to be updated before curX_ or
      curF_ is dereferenced.
      @{ */
      void initCoordinateAccessOptimizationPointers();
      mutable cell_idx_t   fastestCoord_;
      cell_idx_t * curX_;
      cell_idx_t * curY_;
      cell_idx_t * curZ_;
      cell_idx_t * curF_;
      //@}
      //****************************************************************************************************************

   };


   template <typename T, uint_t fieldFSize>
   class ForwardFieldIterator : public FieldIterator<T,fieldFSize>
   {
   public:
      using Parent = FieldIterator<T, fieldFSize>;
      using FieldType = Field<typename std::remove_const<T>::type, fieldFSize>;

      //**Constructor/Destructor****************************************************************************************
      /*!\name Constructor/Destructor */
      //@{
      explicit ForwardFieldIterator(const FieldType * field,
                                    cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                                    uint_t xs, uint_t ys, uint_t zs, uint_t fs )
         : FieldIterator<T,fieldFSize>( field, xBeg, yBeg,zBeg, fBeg,xs ,ys, zs, fs, true ) {}


      explicit ForwardFieldIterator()
         : FieldIterator<T,fieldFSize> ()  {}
      //@}
      //****************************************************************************************************************


      //**Operators*****************************************************************************************************
      /*!\name Operators */
      //@{
      inline ForwardFieldIterator& operator++();
      inline ForwardFieldIterator& operator--();
      //@}
      //****************************************************************************************************************

      //**Fast Iteration ***********************************************************************************************
      /*!\name Fast Iteration */
      //@{
      inline void incrOuter();
      inline void incrInner()        { ++Parent::linePtr_; }
      inline bool testInner() const  { return Parent::linePtr_ != Parent::lineEnd_; }
      //@}
      //****************************************************************************************************************

   };



   template <typename T, uint_t fieldFSize>
   class ReverseFieldIterator : public FieldIterator<T,fieldFSize>
   {
   public:
       using Parent = FieldIterator<T, fieldFSize>;
       using FieldType = Field<typename std::remove_const<T>::type, fieldFSize>;

      //**Constructor/Destructor****************************************************************************************
      /*!\name Constructor/Destructor */
      //@{
      explicit ReverseFieldIterator(const FieldType * field,
                                    cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                                    uint_t xs, uint_t ys, uint_t zs, uint_t fs )
          : FieldIterator<T,fieldFSize>( field, xBeg, yBeg,zBeg, fBeg,xs ,ys, zs, fs, false )
      { }

      explicit ReverseFieldIterator()
         : FieldIterator<T,fieldFSize> ()  {}
      //@}
      //****************************************************************************************************************


      //**Operators*****************************************************************************************************
      /*!\name Operators */
      //@{

      inline ReverseFieldIterator& operator++();
      inline ReverseFieldIterator& operator--();
      //@}
      //****************************************************************************************************************


      //** Fast Iteration **********************************************************************************************
      /*!\name Fast Iteration */
      //@{
      inline void incrOuter()        { Parent::decrementLine();               }
      inline void incrInner()        { --Parent::linePtr_;                    }
      inline bool testInner() const  { return Parent::linePtr_ != Parent::lineBegin_; }
      //@}
      //****************************************************************************************************************
   };



} // namespace field
} // namespace walberla


#include "FieldIterator.impl.h"

