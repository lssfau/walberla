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
//! \file Field.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Definitions of Field members
//
//======================================================================================================================

#include "field/iterators/IteratorMacros.h"

#include "core/math/Utility.h" // for equal() in operator==

#include <algorithm>  // std::copy


namespace walberla {
namespace field {


   //===================================================================================================================
   //
   //  CONSTRUCTION
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*!Creates an uninitialized field of size zero (no allocated memory)
    *
    * This field has to be initialized before it can be used using the init() method
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_>::Field( )
       : values_( nullptr ), valuesWithOffset_( nullptr ),
         xSize_(0), ySize_(0), zSize_(0),
         xAllocSize_(0), yAllocSize_(0), zAllocSize_(0), fAllocSize_(0)
   {
   }

   //*******************************************************************************************************************
   /*!Creates an uninitialized field of given size
    *
    * \param _xSize  size of x dimension
    * \param _ySize  size of y dimension
    * \param _zSize  size of z dimension
    * \param l       memory layout of the field (see Field::Layout)
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_>::Field( uint_t _xSize, uint_t _ySize, uint_t _zSize, const Layout & l,
                           const shared_ptr<FieldAllocator<T> > &alloc )
       : values_( nullptr ), valuesWithOffset_( nullptr )
   {
      init(_xSize,_ySize,_zSize,l,alloc);
   }


   //*******************************************************************************************************************
   /*! Creates a field and initializes it with constant
    *
    * \param _xSize  size of x dimension
    * \param _ySize  size of y dimension
    * \param _zSize  size of z dimension
    * \param initVal every element of the field is set to initVal
    * \param l       memory layout of the field (see Field::Layout)
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_>::Field( uint_t _xSize, uint_t _ySize, uint_t _zSize, const T & initVal, const Layout & l,
                           const shared_ptr<FieldAllocator<T> > &alloc )
          : values_( nullptr ), valuesWithOffset_( nullptr )
   {
      init(_xSize,_ySize,_zSize,l,alloc);
      set(initVal);
   }


   //*******************************************************************************************************************
   /*! Creates a field and initializes f coordinate with vector values
    *
    * \param _xSize  size of x dimension
    * \param _ySize  size of y dimension
    * \param _zSize  size of z dimension
    * \param fValues initializes f coordinate with values from vector (see set(std::vector&) )
    * \param l       memory layout of the field (see Field::Layout)
    * \param alloc  class that describes how to allocate memory for the field, see FieldAllocator
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_>::Field( uint_t _xSize, uint_t _ySize, uint_t _zSize,
                           const std::vector<T> & fValues, const Layout & l,
                           const shared_ptr<FieldAllocator<T> > &alloc)
        : values_( NULL ), valuesWithOffset_( NULL )
   {
      init(_xSize,_ySize,_zSize,l,alloc);
      set(fValues);
   }


   //*******************************************************************************************************************
   /*! Deletes all stored data, and resizes the field
    *
    *  The resized field is uninitialized.
    * \param _xSize  size of x dimension
    * \param _ySize  size of y dimension
    * \param _zSize  size of z dimension
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::resize( uint_t _xSize, uint_t _ySize, uint_t _zSize )
   {
      if ( _xSize == xSize_ &&  _ySize == ySize_ && _zSize == zSize_ )
         return;

      allocator_->decrementReferenceCount( values_ );
      values_ = nullptr;
      valuesWithOffset_ = nullptr;
      init( _xSize, _ySize, _zSize, layout_, allocator_ );
   }



   //*******************************************************************************************************************
   /*! Returns a shallow copy of the current field.
    *
    *  Shallow copy means, that the new field internally uses the same memory as this field.
    *
    * \return a new field, that has to be freed by caller
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_> * Field<T,fSize_>::cloneShallowCopy() const
   {
      return cloneShallowCopyInternal() ;
   }


   //*******************************************************************************************************************
   /*! Returns a flattened shallow copy of the current field.
    *
    *  Shallow copy means, that the new field internally uses the same memory as this field.
    *  Flattened means that any VectorTrait-compatible containers are absorbed into the fSize.
    *
    * \return a new field, that has to be freed by caller
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   typename Field<T,fSize_>::FlattenedField * Field<T,fSize_>::flattenedShallowCopy() const
   {
      return flattenedShallowCopyInternal();
   }


   //*******************************************************************************************************************
   /*!\brief Does the same as cloneShallowCopy (but is virtual)
    *
    * This version has to be implemented by derived classes. The cloneShallowCopy() itself cannot be
    * virtual, since the implementation of cloneShallowCopy() of derived classes has a different signature.
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_> * Field<T,fSize_>::cloneShallowCopyInternal() const
   {
      return new Field<T,fSize_>(*this) ;
   }

   //*******************************************************************************************************************
   /*!\brief Does the same as flattenedShallowCopy (but is virtual)
    *
    * This version has to be implemented by derived classes. The flattenedShallowCopy() itself cannot be
    * virtual, since the implementation of flattenedShallowCopy() of derived classes has a different signature.
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   typename Field<T,fSize_>::FlattenedField * Field<T,fSize_>::flattenedShallowCopyInternal() const
   {
      return new FlattenedField(*this) ;
   }

   //*******************************************************************************************************************
   /*! Creates a new field that has equal size and layout as this field. The memory of the new
    *        field is uninitialized.
    *
    * \return a new field, that has to be freed by caller
    *******************************************************************************************************************/
   template <typename T, uint_t fSize_>
   Field<T,fSize_> * Field<T,fSize_>::cloneUninitialized() const
   {
      Field<T,fSize_> * res = cloneShallowCopy();
      res->allocator_->decrementReferenceCount( res->values_ );
      res->values_ = res->allocator_->allocate ( res->allocSize() );

      const auto offset = res->xOff_*res->xfact_+ res->yOff_*res->yfact_+ res->zOff_*res->zfact_;
      res->valuesWithOffset_ = res->values_ + offset;

      WALBERLA_ASSERT ( hasSameSize     ( *res ) );
      WALBERLA_ASSERT ( hasSameAllocSize( *res ) );    
      
      return res;
   }


   //*******************************************************************************************************************
   /*! Returns a deep copy of the current field. The data is copied over.
    *
    * \return a new field, that has to be freed by caller
    *******************************************************************************************************************/
   template <typename T, uint_t fSize_>
   Field<T,fSize_> * Field<T,fSize_>::clone() const
   {
      Field<T,fSize_> * res = cloneUninitialized();
      WALBERLA_ASSERT_EQUAL ( allocSize_, res->allocSize_ );

      std::copy( values_, values_ + allocSize_, res->values_ );

      return res;
   }


   //*******************************************************************************************************************
   /*! Private copy constructor that creates a shallow copy
    *        i.e. reuses the memory of the copied field
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_>::Field( const Field<T,fSize_> & other )
      : values_           ( other.values_ ),
        valuesWithOffset_ ( other.valuesWithOffset_ ),
        xOff_             ( other.xOff_),
        yOff_             ( other.yOff_),
        zOff_             ( other.zOff_),
        xSize_            ( other.xSize_ ),
        ySize_            ( other.ySize_ ),
        zSize_            ( other.zSize_ ),
        xAllocSize_       ( other.xAllocSize_ ),
        yAllocSize_       ( other.yAllocSize_ ),
        zAllocSize_       ( other.zAllocSize_ ),
        fAllocSize_       ( other.fAllocSize_ ),
        layout_           ( other.layout_ ),
        allocSize_        ( other.allocSize_ ),
        ffact_            ( other.ffact_ ),
        xfact_            ( other.xfact_ ),
        yfact_            ( other.yfact_ ),
        zfact_            ( other.zfact_ ),
        allocator_        ( other.allocator_ )
   {
      allocator_->incrementReferenceCount ( values_ );
   }


   //*******************************************************************************************************************
   /*! Private copy constructor that creates a flattened shallow copy
    *        i.e. reuses the memory of the copied field
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   template <typename T2, uint_t fSize2>
   Field<T,fSize_>::Field( const Field<T2,fSize2> & other )
      : values_           ( other.values_[0].data() ),
        valuesWithOffset_ ( other.valuesWithOffset_[0].data() ),
        xOff_             ( other.xOff_),
        yOff_             ( other.yOff_),
        zOff_             ( other.zOff_),
        xSize_            ( other.xSize_ ),
        ySize_            ( other.ySize_ ),
        zSize_            ( other.zSize_ ),
        xAllocSize_       ( other.xAllocSize_ ),
        yAllocSize_       ( other.yAllocSize_ ),
        zAllocSize_       ( other.zAllocSize_ ),
        fAllocSize_       ( other.fAllocSize_*fSize_/fSize2 ),
        layout_           ( other.layout_ ),
        allocSize_        ( other.allocSize_*fSize_/fSize2 ),
        ffact_            ( other.ffact_ ),
        xfact_            ( other.xfact_*cell_idx_t(fSize_/fSize2) ),
        yfact_            ( other.yfact_*cell_idx_t(fSize_/fSize2) ),
        zfact_            ( other.zfact_*cell_idx_t(fSize_/fSize2) ),
        allocator_        ( std::shared_ptr<FieldAllocator<T>>(other.allocator_, reinterpret_cast<FieldAllocator<T>*>(other.allocator_.get())) )
   {
      WALBERLA_CHECK_EQUAL(layout_, Layout::zyxf)
      static_assert(fSize_ % fSize2 == 0, "number of field components do not match");
      static_assert(std::is_same<typename Field<T2,fSize2>::FlattenedField, Field<T,fSize_>>::value, "field types are incompatible for flattening");
      allocator_->incrementReferenceCount ( values_ );
   }


   //*******************************************************************************************************************
   /*! Initializes the field with a given size, in a given layout
    *
    * Must be called exactly once!  This is automatically called by all constructors
    * that take at least one argument
    *
    * \param _xSize  size of x dimension
    * \param _ySize  size of y dimension
    * \param _zSize  size of z dimension
    * \param l       memory layout of the field (see Field::Layout)
    * \param alloc   the allocator to use. If a NULL shared pointer is given, a sensible default is selected,
    *                depending on layout
    * \param innerGhostLayerSizeForAlignedAlloc
    *                This parameter should be set to zero for field that have no ghost layers.
    *                This parameter is passed to the allocator and can there be used to ensure
    *                alignment of the first INNER cell in each line
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void Field<T, fSize_>::init( uint_t _xSize, uint_t _ySize, uint_t _zSize,
                                const Layout & l, shared_ptr<FieldAllocator<T> > alloc,
                                uint_t innerGhostLayerSizeForAlignedAlloc )
   {
      WALBERLA_ASSERT_NULLPTR( values_ );
      WALBERLA_ASSERT_NULLPTR( valuesWithOffset_ );

      // Automatically select allocator if none was given
      if ( alloc == nullptr )
      {
#if defined(__ARM_FEATURE_SVE) && defined(__ARM_FEATURE_SVE_BITS) && __ARM_FEATURE_SVE_BITS > 0
         const uint_t alignment = __ARM_FEATURE_SVE_BITS/8;
#elif defined(__ARM_FEATURE_SVE)
         const uint_t alignment = 64;
#elif defined(__ARM_NEON)
         const uint_t alignment = 16;
#elif defined(__AVX512F__)
         const uint_t alignment = 64;
#elif defined(__AVX__)
         const uint_t alignment = 32;
#elif defined(__SSE__) || defined(_MSC_VER)
         const uint_t alignment = 16;
#elif defined(__BIGGEST_ALIGNMENT__)
         const uint_t alignment = __BIGGEST_ALIGNMENT__;
#else
         const uint_t alignment = 64;
#endif

         // aligned allocator only used (by default) if ...
         if ( l == fzyx                      && // ... we use a structure of arrays layout
              _xSize * sizeof(T) > alignment && // ... the inner coordinate is sufficiently large
              sizeof(T) < alignment          && // ... the stored data type is smaller than the alignment
              alignment % sizeof(T) == 0 )      // ... there is an integer number of elements fitting in one aligned line
            alloc = make_shared<AllocateAligned<T,alignment> >();
         else
            alloc = make_shared<StdFieldAlloc<T> > ();
      }

      allocator_ = alloc;
      allocator_->setInnerGhostLayerSize( innerGhostLayerSizeForAlignedAlloc );
      values_ = nullptr;
      xSize_ = _xSize;
      ySize_ = _ySize;
      zSize_ = _zSize;
      xAllocSize_ = yAllocSize_ = zAllocSize_ = fAllocSize_ = 0; // is set in alloc->allocate()

      layout_ = l;

      WALBERLA_ASSERT(layout_ == zyxf || layout_ == fzyx)

      if (layout_ == fzyx ) {
         values_ = allocator_->allocate(fSize_, zSize_, ySize_, xSize_, zAllocSize_, yAllocSize_, xAllocSize_);
         fAllocSize_ = fSize_;

         WALBERLA_CHECK_LESS_EQUAL( fSize_ * xAllocSize_ * yAllocSize_ * zAllocSize_ + xSize_ + ySize_ * xAllocSize_ + zSize_ * xAllocSize_ * yAllocSize_,
                                    std::numeric_limits< int64_t >::max(),
                                    "The data type 'int64_t' is too small for your field size! Your field is too large." );

         ffact_ = int64_t(xAllocSize_) * int64_t(yAllocSize_) * int64_t(zAllocSize_);
         zfact_ = int64_t(xAllocSize_) * int64_t(yAllocSize_);
         yfact_ = int64_t(xAllocSize_);
         xfact_ = 1;
      } else {
         values_ = allocator_->allocate(zSize_, ySize_, xSize_, fSize_, yAllocSize_, xAllocSize_, fAllocSize_);
         zAllocSize_ = zSize_;

         WALBERLA_CHECK_LESS_EQUAL( fSize_ + xSize_ * fAllocSize_ + ySize_ * fAllocSize_ * xAllocSize_ + zSize_ * fAllocSize_ * xAllocSize_ * yAllocSize_,
                                    std::numeric_limits< int64_t >::max(),
                                    "The data type 'int64_t' is too small for your field size! Your field is too large." );

         zfact_ = int64_t (fAllocSize_) * int64_t(xAllocSize_) * int64_t(yAllocSize_);
         yfact_ = int64_t(fAllocSize_) * int64_t(xAllocSize_);
         xfact_ = int64_t (fAllocSize_);
         ffact_ = 1;
      }

      WALBERLA_ASSERT(xAllocSize_ >= xSize_);
      WALBERLA_ASSERT(yAllocSize_ >= ySize_);
      WALBERLA_ASSERT(zAllocSize_ >= zSize_);

      allocSize_ = fAllocSize_ * xAllocSize_ * yAllocSize_ * zAllocSize_;

      xOff_ = yOff_ = zOff_ = 0;
      valuesWithOffset_ = values_;
   }


   //*******************************************************************************************************************
   /*! Destructor, using Allocator template parameter
    *******************************************************************************************************************/
   template<typename T, uint_t fSize>
   Field<T,fSize>::~Field()
   {
      allocator_->decrementReferenceCount( values_ );
   }




   //===================================================================================================================
   //
   //  ITERATORS
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*! Returns iterator, which can iterate over complete field in a suitable order depending on layout
    *
    * Use this when iterating over a complete field, faster than 4 nested loops and operator() calls in the
    * innermost loop because no index calculations have to be done.
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::iterator Field<T,fSize_>::begin()
   {
      return iterator( this,0,0,0,0, xSize(), ySize(), zSize(), fSize() );
   }

   //*******************************************************************************************************************
   /*! Returns const_iterator, see begin()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::const_iterator Field<T,fSize_>::begin() const
   {
       return const_iterator ( const_cast< Field<T,fSize_> * >(this), 0,0,0,0,
                               xSize(), ySize(), zSize(), fSize() );
   }


   //*******************************************************************************************************************
   /*! Returns iterator which iterates over a sub-block of the field
    *
    * - Iterator over block defined by ( xBeg <= x < xEnd,  yBeg <= y < yEnd, ....)
    * - layout aware
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::iterator
   Field<T,fSize_>::beginSlice( cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                                cell_idx_t xEnd, cell_idx_t yEnd, cell_idx_t zEnd, cell_idx_t fEnd )
   {
      WALBERLA_ASSERT_LESS( xBeg, xEnd );
      WALBERLA_ASSERT_LESS( yBeg, yEnd );
      WALBERLA_ASSERT_LESS( zBeg, zEnd );
      assertValidCoordinates( xBeg  , yBeg  , zBeg  , fBeg   );
      assertValidCoordinates( xEnd-1, yEnd-1, zEnd-1, fEnd-1 ); // -1 since end points behind valid coordinates

      return iterator( this, xBeg, yBeg, zBeg, fBeg, uint_c(xEnd-xBeg), uint_c(yEnd-yBeg), uint_c(zEnd-zBeg), uint_c(fEnd-fBeg) );
   }

   //*******************************************************************************************************************
   /*! Const variant of beginSlice()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::const_iterator
   Field<T,fSize_>::beginSlice( cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                                cell_idx_t xEnd, cell_idx_t yEnd, cell_idx_t zEnd, cell_idx_t fEnd ) const
   {
      WALBERLA_ASSERT_LESS( xBeg, xEnd);
      WALBERLA_ASSERT_LESS( yBeg, yEnd);
      WALBERLA_ASSERT_LESS( zBeg, zEnd);
      assertValidCoordinates( xBeg  , yBeg  , zBeg  , fBeg   );
      assertValidCoordinates( xEnd-1, yEnd-1, zEnd-1, fEnd-1 ); // -1 since end points behind valid coordinates

      return const_iterator( this, xBeg, yBeg, zBeg, fBeg, uint_c( xEnd-xBeg ), uint_c( yEnd-yBeg ), uint_c( zEnd-zBeg ), uint_c( fEnd-fBeg ) );
   }




   //*******************************************************************************************************************
   /*! Returns iterator which iterates over a slice, but only in x,y,z coordinates.
    * \param f fixed value of f coordinate, where iterator points to in each cell
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::iterator
   Field<T,fSize_>::beginSliceXYZ( const CellInterval & ci, cell_idx_t f )
   {
      return ci.empty() ? end() : iterator( this, ci.xMin(),  ci.yMin(),  ci.zMin(),  f,
                                                  ci.xSize(), ci.ySize(), ci.zSize(), 1 );
   }

   //*******************************************************************************************************************
   /*! Const variant of beginSliceXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::const_iterator
   Field<T,fSize_>::beginSliceXYZ ( const CellInterval & ci, cell_idx_t f ) const
   {
      return ci.empty() ? end() : const_iterator( this, ci.xMin(),  ci.yMin(),  ci.zMin(),  f,
                                                        ci.xSize(), ci.ySize(), ci.zSize(), 1 );
   }


   //*******************************************************************************************************************
   /*! Iterates only over XYZ coordinate, f is always 0
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::iterator
   Field<T,fSize_>::beginXYZ()
   {
      return iterator( this, 0,0,0,0, xSize(), ySize(), zSize(), 1 );
   }

   //*******************************************************************************************************************
   /*! Const version of beginXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::const_iterator
   Field<T,fSize_>::beginXYZ() const
   {
      return const_iterator( this, 0,0,0,0, xSize(), ySize(), zSize(), 1 );
   }

   //*******************************************************************************************************************
   /*! End iterator, can be used with begin() and beginBlock()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   const ForwardFieldIterator<T,fSize_> Field<T,fSize_>::staticEnd = ForwardFieldIterator<T,fSize_>();

   template<typename T, uint_t fSize_>
   inline const typename Field<T,fSize_>::iterator & Field<T,fSize_>::end()
   {
      return staticEnd;
   }

   //*******************************************************************************************************************
   /*! Const end iterator, see end()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   const ForwardFieldIterator<const T,fSize_> Field<T,fSize_>::staticConstEnd = ForwardFieldIterator<const T,fSize_>();

   template<typename T, uint_t fSize_>
   inline const typename Field<T,fSize_>::const_iterator & Field<T,fSize_>::end() const
   {
      return staticConstEnd;
   }






   //===================================================================================================================
   //
   //  REVERSE ITERATION
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*! Returns reverse iterator, which can iterate over complete field in a suitable order depending on layout
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::reverse_iterator Field<T,fSize_>::rbegin()
   {
      return reverse_iterator( this,0,0,0,0, xSize(),   ySize(),   zSize(),  fSize() );
   }

   //*******************************************************************************************************************
   /*! Returns const_reverse_iterator, see begin()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::const_reverse_iterator Field<T,fSize_>::rbegin() const
   {
       return const_reverse_iterator ( this, 0,0,0,0, xSize(),   ySize(),   zSize(),  fSize() );
   }


   //*******************************************************************************************************************
   /*! Iterates only over XYZ coordinate, f is always 0
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::reverse_iterator
   Field<T,fSize_>::rbeginXYZ()
   {
      return reverse_iterator( this, 0,0,0, 0, xSize()  , ySize()  , zSize()  , 1 );
   }

   //*******************************************************************************************************************
   /*! Const version of beginXYZ()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline typename Field<T,fSize_>::const_reverse_iterator
   Field<T,fSize_>::rbeginXYZ() const
   {
      return const_reverse_iterator( this, 0,0,0, 0, xSize()  , ySize()  , zSize()  , 1 );
   }

   //*******************************************************************************************************************
   /*! End iterator, can be used with begin() and beginBlock()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   const ReverseFieldIterator<T,fSize_> Field<T,fSize_>::staticREnd = ReverseFieldIterator<T,fSize_>();

   template<typename T, uint_t fSize_>
   inline const typename Field<T,fSize_>::reverse_iterator & Field<T,fSize_>::rend()
   {
      return staticREnd;
   }

   //*******************************************************************************************************************
   /*! Const end iterator, see end()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   const ReverseFieldIterator<const T,fSize_> Field<T,fSize_>::staticConstREnd = ReverseFieldIterator<const T,fSize_>();

   template<typename T, uint_t fSize_>
   inline const typename Field<T,fSize_>::const_reverse_iterator & Field<T,fSize_>::rend() const
   {
      return staticConstREnd;
   }


   //===================================================================================================================
   //
   //  SIZE INFORMATION
   //
   //===================================================================================================================

   template<typename T, uint_t fSize_>
   inline CellInterval Field<T,fSize_>::xyzSize() const
   {
      return CellInterval (0,0,0, cell_idx_c( xSize() )-1,
                                  cell_idx_c( ySize() )-1,
                                  cell_idx_c( zSize() )-1 );
   }


   template<typename T, uint_t fSize_>
   inline CellInterval Field<T,fSize_>::xyzAllocSize() const
   {
      return CellInterval( -xOff_,
                           -yOff_,
                           -zOff_,
                            cell_idx_c( xAllocSize() ) - xOff_-1,
                            cell_idx_c( yAllocSize() ) - yOff_-1,
                            cell_idx_c( zAllocSize() ) - zOff_-1 );
   }

   template<typename T, uint_t fSize_>
   inline uint_t  Field<T,fSize_>::size( uint_t coord )  const
   {
      switch (coord) {
         case 0: return this->xSize();
         case 1: return this->ySize();
         case 2: return this->zSize();
         case 3: return this->fSize();
         default: WALBERLA_ASSERT(false); return 0;
      }
   }



   //===================================================================================================================
   //
   //  ELEMENT ACCESS
   //
   //===================================================================================================================
#ifndef NDEBUG
   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::assertValidCoordinates( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const
   {
      //Bounds checks
      const cell_idx_t xEff = xOff_ + x;
      const cell_idx_t yEff = yOff_ + y;
      const cell_idx_t zEff = zOff_ + z;

      WALBERLA_ASSERT_GREATER_EQUAL( xEff, 0, "Field access out of bounds: x too small: " << x << " < " << - xOff_ );
      WALBERLA_ASSERT_GREATER_EQUAL( yEff, 0, "Field access out of bounds: y too small: " << y << " < " << - yOff_ );
      WALBERLA_ASSERT_GREATER_EQUAL( zEff, 0, "Field access out of bounds: z too small: " << z << " < " << - zOff_ );
      WALBERLA_ASSERT_GREATER_EQUAL( f,    0, "Field access out of bounds: f too small: " << f << " < " << 0 );

      WALBERLA_ASSERT_LESS( xEff, cell_idx_c( xAllocSize() ), "Field access out of bounds: x too big: " << x << " >= " << (xAllocSize() - uint_c(xOff_)) );
      WALBERLA_ASSERT_LESS( yEff, cell_idx_c( yAllocSize() ), "Field access out of bounds: y too big: " << y << " >= " << (yAllocSize() - uint_c(yOff_)) );
      WALBERLA_ASSERT_LESS( zEff, cell_idx_c( zAllocSize() ), "Field access out of bounds: z too big: " << z << " >= " << (zAllocSize() - uint_c(zOff_)) );
      WALBERLA_ASSERT_LESS( f,    cell_idx_c( fSize_)       , "Field access out of bounds: f too big: " << f << " >= " << fSize_ );
   }
#else
   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::assertValidCoordinates( cell_idx_t,  cell_idx_t,  cell_idx_t,  cell_idx_t ) const
   { }
#endif


   template<typename T, uint_t fSize_>
   bool Field<T,fSize_>::coordinatesValid( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const
   {
      //Bounds checks
      const cell_idx_t xEff = xOff_ + x;
      const cell_idx_t yEff = yOff_ + y;
      const cell_idx_t zEff = zOff_ + z;

      if (xEff < 0 || yEff < 0 || zEff < 0 || f < 0 )
         return false;
      if ( xEff >= cell_idx_c( xAllocSize() ) ||
           yEff >= cell_idx_c( yAllocSize() ) ||
           zEff >= cell_idx_c( zAllocSize() ) ||
           f    >= cell_idx_c( fSize_ ) )
         return false;

      return true;
   }

   //*******************************************************************************************************************
   /*! Accesses the value at given coordinate
    *
    * When WALBERLA_FIELD_MONITORED_ACCESS is defined, all registered monitor functions are executed when get is called
    *
    * \note operator() is equivalent to this function
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::get( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f ) const
   {
      assertValidCoordinates( x, y, z, f );

      const int64_t index = f*int64_t(ffact_) + int64_t(x)*int64_t(xfact_) + int64_t(y)*int64_t(yfact_) + int64_t(z)*int64_t(zfact_);

      WALBERLA_ASSERT_LESS( int64_c(index) + int64_c(valuesWithOffset_ - values_), int64_c(allocSize_) );
      WALBERLA_ASSERT_GREATER_EQUAL( int64_c(index) + int64_c(valuesWithOffset_ - values_), int64_c(0) );

#     ifdef WALBERLA_FIELD_MONITORED_ACCESS
      for(uint_t i=0; i< monitorFuncs_.size(); ++i )
         monitorFuncs_[i] (x,y,z,f, *(valuesWithOffset_ + index) );
#     endif

      return *( valuesWithOffset_ + index );
   }

   //*******************************************************************************************************************
   /*! Non-Const variant of get()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T&  Field<T,fSize_>::get(cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f)
   {
      const Field<T,fSize_>& const_this = *this;
      return const_cast<T&>( const_this.get(x,y,z,f) );
   }


   //*******************************************************************************************************************
   /*! get() variant which takes a uint_t as last coordinate, as for example Stencil::toIdx() returns
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::get( cell_idx_t x, cell_idx_t y, cell_idx_t z, uint_t f ) const
   {
      return get(x,y,z,cell_idx_c(f));
   }

   //*******************************************************************************************************************
   /*! get() variant which takes a uint_t as last coordinate, as for example Stencil::toIdx() returns
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::get( cell_idx_t x, cell_idx_t y, cell_idx_t z, uint_t f )
   {
      return get(x,y,z,cell_idx_c(f));
   }


   //*******************************************************************************************************************
   /*! get function with only (x,y,z) coordinates, assumes fSize=1
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::get( cell_idx_t x, cell_idx_t y, cell_idx_t z)
   {
      static_assert(fSize_ == 1, "f coordinate omitted for field with fSize > 1 ");
      return get(x,y,z,0);
   }

   //*******************************************************************************************************************
   /*! get function with only (x,y,z) coordinates, assumes fSize=1
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::get( cell_idx_t x, cell_idx_t y, cell_idx_t z) const
   {
      static_assert(fSize_ == 1, "f coordinate omitted for field with fSize > 1 ");
      return get(x,y,z,0);
   }

   //*******************************************************************************************************************
   /*! get overload using a cell as input, only possible if fSize=1
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::get( const Cell & cell )
   {
      return get( cell.x(), cell.y(), cell.z() );
   }

   //*******************************************************************************************************************
   /*! get overload using a cell as input, only possible if fSize=1
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::get( const Cell & cell ) const
   {
      return get( cell.x(), cell.y(), cell.z() );
   }



   //*******************************************************************************************************************
   /*! get overload, where position is specified using an iterator of another field with equal size
    *
    * Do not use this for iterator that belong to this field, here *iterator is more convenient
    *
    * \param iter Iterator that belongs to another field that has equal size
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::get( const base_iterator & iter )
   {
      WALBERLA_ASSERT( hasSameAllocSize( *iter.getField() ) );
      WALBERLA_ASSERT( hasSameSize     ( *iter.getField() ) );
      WALBERLA_ASSERT( layout() == iter.getField()->layout() );
      WALBERLA_ASSERT( this != iter.getField() ); // use *iterator instead!
      return  *( valuesWithOffset_ + (&(*iter) - iter.getField()->valuesWithOffset_) );
   }

   //*******************************************************************************************************************
   /*! get overload, where position is specified using an iterator of another field with equal size
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::get( const base_iterator & iter ) const
   {
      WALBERLA_ASSERT( hasSameAllocSize( *iter.getField() ) );
      WALBERLA_ASSERT( hasSameSize     ( *iter.getField() ) );
      WALBERLA_ASSERT( layout() == iter.getField()->layout() );
      WALBERLA_ASSERT( this != iter.getField() ); // use *iterator instead!
      return *(valuesWithOffset_ + (&(*iter) - iter.getField()->valuesWithOffset_) );
   }

   //*******************************************************************************************************************
   /*! returns neighboring value of cell in the given direction
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::getNeighbor( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d )
   {
      return getNeighbor(x,y,z,cell_idx_t(0),d);
   }

   //*******************************************************************************************************************
   /*! returns neighboring value of cell in the given direction
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::getNeighbor( cell_idx_t x, cell_idx_t y, cell_idx_t z, stencil::Direction d ) const
   {
      return getNeighbor(x,y,z,cell_idx_t(0),d);
   }

   //*******************************************************************************************************************
   /*! returns neighboring value of cell in the given direction
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::getNeighbor( cell_idx_t x, cell_idx_t y, cell_idx_t z, uint_t f, stencil::Direction d )
   {
      return get( x + stencil::cx[d],
                  y + stencil::cy[d],
                  z + stencil::cz[d],
                  f);
   }

   //*******************************************************************************************************************
   /*! returns neighboring value of cell in the given direction
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::getNeighbor( cell_idx_t x, cell_idx_t y, cell_idx_t z, uint_t f, stencil::Direction d ) const
   {
      return get( x + stencil::cx[d],
                  y + stencil::cy[d],
                  z + stencil::cz[d],
                  f);
   }

   //*******************************************************************************************************************
   /*! get overload using a cell as input, only possible if fSize=1, with neighbor access
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::getNeighbor( const Cell & cell, stencil::Direction d )
   {
      return get( cell.x() + stencil::cx[d],
                  cell.y() + stencil::cy[d],
                  cell.z() + stencil::cz[d] );
   }

   //*******************************************************************************************************************
   /*! get overload using a cell as input, only possible if fSize=1, with neighbor access
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::getNeighbor( const Cell & cell, stencil::Direction d ) const
   {
      return get( cell.x() + stencil::cx[d],
                  cell.y() + stencil::cy[d],
                  cell.z() + stencil::cz[d] );
   }

   template<typename T, uint_t fSize_>
   inline T & Field<T,fSize_>::getF( T * const xyz0, const cell_idx_t f )
   {
      WALBERLA_ASSERT_LESS( f, cell_idx_c(fSize_) );
      WALBERLA_ASSERT( addressInsideAllocedSpace( xyz0 ) );
      WALBERLA_ASSERT( addressInsideAllocedSpace( xyz0 + f * ffact_ ) );
      return *( xyz0 + f * ffact_ );
   }

   template<typename T, uint_t fSize_>
   inline const T & Field<T,fSize_>::getF( const T * const xyz0, const cell_idx_t f ) const
   {
      WALBERLA_ASSERT_LESS( f, cell_idx_c(fSize_) );
      WALBERLA_ASSERT( addressInsideAllocedSpace( xyz0 ) );
      WALBERLA_ASSERT( addressInsideAllocedSpace( xyz0 + f * ffact_ ) );
      return *( xyz0 + f * ffact_ );
   }

   //*******************************************************************************************************************
   /*! Sets all entries of the field to given value
    *
    * Works only in the regions specified by size(), not in the complete allocated region as specified by allocSize()
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::set (const T & value)
   {
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wundefined-bool-conversion"
#endif
      // take care of proper thread<->memory assignment (first-touch allocation policy !)
      WALBERLA_FOR_ALL_CELLS_XYZ( this,

         for( uint_t f = uint_t(0); f < fSize_; ++f )
            get(x,y,z,f) = value;
      
      ) // WALBERLA_FOR_ALL_CELLS_XYZ
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
   }

   //*******************************************************************************************************************
   /*! Initializes the f coordinate to values from vector
    * Sets the entry (x,y,z,f) to fValues[f]
    *
    * Works only in the regions specified by size(), not in the complete allocated region as specified by allocSize()
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::set (const std::vector<T> & fValues)
   {
      WALBERLA_ASSERT(fValues.size() == fSize_);
      
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wundefined-bool-conversion"
#endif
      // take care of proper thread<->memory assignment (first-touch allocation policy !)
      WALBERLA_FOR_ALL_CELLS_XYZ( this,

         for( uint_t f = uint_t(0); f < fSize_; ++f )
            get(x,y,z,f) = fValues[f];
      
      ) // WALBERLA_FOR_ALL_CELLS_XYZ
#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
   }

   //*******************************************************************************************************************
   /*! Copies all entries of the given field to this field
    *
    * Only works when xSize(), ySize(), zSize() and allocSize() are equal
    * Copies complete allocated region as specified by allocSize()
    *
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline void Field<T,fSize_>::set (const Field<T,fSize_> & other )
   {
      WALBERLA_ASSERT_EQUAL( xyzSize(), other.xyzSize() );
      WALBERLA_ASSERT_EQUAL( allocSize(), other.allocSize() );
      std::copy( other.values_, other.values_ + allocSize_, values_ );
   }

   //*******************************************************************************************************************
   /*! Swap two fields efficiently by exchanging only values_ pointer
    * The two fields have to have identical sizes and same layout.
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline void Field<T,fSize_>::swapDataPointers( Field<T,fSize_> & other)
   {
      WALBERLA_ASSERT( hasSameAllocSize(other) );
      WALBERLA_ASSERT( hasSameSize(other) );
      WALBERLA_ASSERT( layout() == other.layout() );
      std::swap( values_, other.values_ );
      std::swap( valuesWithOffset_, other.valuesWithOffset_ );
   }




   //===================================================================================================================
   //
   //  EQUALITY CHECKS
   //
   //===================================================================================================================


   //*********************************************************************************************************************
   /*! Equality operator compares element-wise
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline bool Field<T,fSize_>::operator==(const Field<T,fSize_> & other) const
   {
      if ( !hasSameSize(other) )
         return false;

      const_iterator lhsIt = this->begin();
      const_iterator rhsIt = other.begin();
      while( lhsIt != this->end() )
      {
         if( !math::equal( *lhsIt, *rhsIt ) )
            return false;

         ++lhsIt;
         ++rhsIt;
      }

      WALBERLA_ASSERT_EQUAL( rhsIt, other.end() );

      return true;
   }

   //*********************************************************************************************************************
   /*! Inequality operator compares element-wise
   *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline bool Field<T, fSize_>::operator!=( const Field<T, fSize_> & other ) const
   {
      return !( *this == other );
   }

   //*******************************************************************************************************************
   /*! True if allocation sizes of all dimensions match
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline bool Field<T,fSize_>::hasSameAllocSize( const Field<T,fSize_> & other ) const
   {
      return xAllocSize_ == other.xAllocSize_ &&
             yAllocSize_ == other.yAllocSize_ &&
             zAllocSize_ == other.zAllocSize_ &&
             fAllocSize_ == other.fAllocSize_;
   }

   //*******************************************************************************************************************
   /*! True if sizes of all dimensions match
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   inline bool Field<T,fSize_>::hasSameSize( const Field<T,fSize_> & other ) const
   {
      return xSize_ == other.xSize_ &&
             ySize_ == other.ySize_ &&
             zSize_ == other.zSize_;
   }


   //===================================================================================================================
   //
   //  CHANGING OFFSETS
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*! Moves the coordinate system of the field
    *
    * Can be used by derived classes, to use only a subset of the field.
    * The complete field remains accessible, but now has coordinates smaller than 0 or
    * bigger than [xyzf]Size()
    * This is used for example in the constructor of the GhostLayerField
    *
    * Internally this is implementing by adding an offset to the values_ pointer, and by adapting the size_ members.
    * \param xOffset The x coordinate that is afterwards mapped to zero
    * \param xs      The new size of the x coordinate. Has to be smaller than (old xSize())-xOffset
    * \param yOffset The y coordinate that is afterwards mapped to zero
    * \param ys      The new size of the y coordinate. Has to be smaller than (old ySize())-yOffset
    * \param zOffset The z coordinate that is afterwards mapped to zero
    * \param zs      The new size of the z coordinate. Has to be smaller than (old zSize())-zOffset
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::setOffsets(uint_t xOffset, uint_t xs,
                                    uint_t yOffset, uint_t ys,
                                    uint_t zOffset, uint_t zs)
   {
      xOff_ = cell_idx_c( xOffset );
      yOff_ = cell_idx_c( yOffset );
      zOff_ = cell_idx_c( zOffset );

      WALBERLA_ASSERT_LESS_EQUAL( uint_c(xOff_) + xs, xAllocSize() );
      WALBERLA_ASSERT_LESS_EQUAL( uint_c(yOff_) + ys, yAllocSize() );
      WALBERLA_ASSERT_LESS_EQUAL( uint_c(zOff_) + zs, zAllocSize() );

      valuesWithOffset_ = values_;
      xSize_ = xs;
      ySize_ = ys;
      zSize_ = zs;
      const int64_t offset = int64_t(xOff_)*int64_t(xfact_) + int64_t(yOff_)*int64_t(yfact_) + int64_t(zOff_)*int64_t(zfact_);
      valuesWithOffset_ = values_ + offset;
   }



   template<typename T, uint_t fSize_>
   inline bool Field<T,fSize_>::addressInsideAllocedSpace(const T * const value) const
   {
      return ( value >= values_ && value < values_ + allocSize_ );
   }

   //===================================================================================================================
   //
   //  SLICING
   //
   //===================================================================================================================


   template<typename T, uint_t fSize_>
   void Field<T,fSize_>::shiftCoordinates( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz )
   {
      WALBERLA_ASSERT_LESS_EQUAL ( uint_c(xOff_ + cx) + xSize(), xAllocSize() );
      WALBERLA_ASSERT_LESS_EQUAL ( uint_c(yOff_ + cy) + ySize(), yAllocSize() );
      WALBERLA_ASSERT_LESS_EQUAL ( uint_c(zOff_ + cz) + zSize(), zAllocSize() );

      setOffsets( uint_c( xOff_ + cx ),  xSize(),  uint_c( yOff_ + cy ), ySize(), uint_c( zOff_ + cz ), zSize() );
   }


   //*******************************************************************************************************************
   /*! Create a different "view" of the field, where the created field has the size of the given sliceInterval
    *
    * The ownership of the returned pointer is transfered to the caller i.e. the caller is responsible for
    * deleting the returned object.
    *
    * The returned field uses the same data as the original field. However the returned field has a
    * different coordinate system defined by the given slice-interval. Modifying the returned slice field
    * modifies also the original field!
    * The sliced field has the size as given by the slice-interval.
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   Field<T,fSize_> * Field<T,fSize_>::getSlicedField( const CellInterval & interval ) const
   {
      auto slicedField = cloneShallowCopy();
      slicedField->slice ( interval );
      return slicedField;
   }


   //*******************************************************************************************************************
   /*! Changes the coordinate system of the field.
    *
    * The origin of the new coordinates is at the cell given by min() of the CellInterval.
    * The new size of the field, is the size of the the CellInterval, however the alloc size does not change.
    * Cells that are not in this cell interval can still be accessed
    * ( by coordinates smaller 0, or bigger than [xyz]Size)
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_ >
   void Field<T,fSize_>::slice( const CellInterval & interval )
   {
      setOffsets ( uint_c( interval.xMin() + xOff_ ), interval.xSize(),
                   uint_c( interval.yMin() + yOff_ ), interval.ySize(),
                   uint_c( interval.zMin() + zOff_ ), interval.zSize() );
   }


   //*******************************************************************************************************************
   /*! Returns the number of objects that internally use the same data.
    *******************************************************************************************************************/
   template<typename T, uint_t fSize_>
   uint_t Field<T,fSize_>::referenceCount( ) const
   {
      return allocator_->referenceCount( values_ );
   }

   //===================================================================================================================
   //
   //  MONITORING
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*! Registers a monitoring function
    *
    * Monitoring works only if compiled with WALBERLA_FIELD_MONITORED_ACCESS.
    * When a field is accessed either by get() or operator() the monitoring function is called
    *******************************************************************************************************************/
#  ifdef WALBERLA_FIELD_MONITORED_ACCESS
   template<typename T, uint_t fSize>
   void Field<T,fSize>::addMonitoringFunction(const MonitorFunction & func)
   {
      monitorFuncs_.push_back(func);
   }
#  else
   template<typename T, uint_t fSize>
   void Field<T,fSize>::addMonitoringFunction(const MonitorFunction & )
   {
   }
#  endif




   //===================================================================================================================
   //
   //  Low Level Functions - do not use if not absolutely necessary
   //
   //===================================================================================================================


   //*******************************************************************************************************************
   /*! Returns internal data allocator
    *
    * The allocator can for example be used to prevent free() on the field data when class is deleted.
    * This is useful if you keep a pointer to the internal data, to make sure it remains valid
    * ( also after possible swapDataPointers() etc. )
    *
    *
    * field->getAllocator()->incrementReferenceCount( field->data() );
    *    Use the data pointer
    * field->getAllocator()->decrementReferenceCount( field->data() );
    *
   */
   //*******************************************************************************************************************
   template<typename T, uint_t fSize_ >
   shared_ptr< FieldAllocator<T> > Field<T,fSize_>::getAllocator() const
   {
      return this->allocator_;
   }


}
}
