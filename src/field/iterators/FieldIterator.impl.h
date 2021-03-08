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
//! \file FieldIterator.impl.h
//! \ingroup field
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================



namespace walberla {
namespace field {



//**********************************************************************************************************************
/*!\brief Constructs an iterator for the given slice of a field
 *
 * The constructed iterator goes over the specified slice of the field.
 * Iterator coordinates (f(),z(),y(),x() ) return coordinates inside that slice (they do not start at 0)
 *
 **********************************************************************************************************************/
template <typename T, uint_t fs>
FieldIterator<T,fs>::FieldIterator( const typename FieldIterator<T,fs>::FieldType * field,
                                    cell_idx_t xBeg, cell_idx_t yBeg, cell_idx_t zBeg, cell_idx_t fBeg,
                                    uint_t sx, uint_t sy, uint_t sz, uint_t sf, bool forward )
   :  f_(field), xBegin_(xBeg), yBegin_(yBeg), zBegin_(zBeg), fBegin_(fBeg)
{
   if ( field->xyzSize().empty() )
   {
      linePtr_ = nullptr;
      lineEnd_ = nullptr;
      f_       = nullptr;
      return;
   }

   using NonConstT = typename std::remove_const<T>::type;

   cur_[0] = cur_[1] = cur_[2] = 0;
   if( f_->layout() == fzyx )
   {
      skips_[0] = ( f_->fAllocSize() - sf ) * uint_c( f_->ffact_ );
      skips_[1] = ( f_->zAllocSize() - sz ) * uint_c( f_->zfact_ );
      skips_[2] = ( f_->yAllocSize() - sy ) * uint_c( f_->yfact_ );
      skips_[3] = ( f_->xAllocSize() - sx ) * uint_c( f_->xfact_ );
      sizes_[0] = sf;
      sizes_[1] = sz;
      sizes_[2] = sy;
      sizes_[3] = sx;

      if ( !forward ) {
         cur_[0] = cell_idx_c( sf - 1 );
         cur_[1] = cell_idx_c( sz - 1 );
         cur_[2] = cell_idx_c( sy - 1 );
      }
   }
   else
   {
      skips_[0] = (f_->zAllocSize() - sz) * uint_c( f_->zfact_ );
      skips_[1] = (f_->yAllocSize() - sy) * uint_c( f_->yfact_ );
      skips_[2] = (f_->xAllocSize() - sx) * uint_c( f_->xfact_ );
      skips_[3] = (f_->fAllocSize() - sf) * uint_c( f_->ffact_ );
      sizes_[0] = sz;
      sizes_[1] = sy;
      sizes_[2] = sx;
      sizes_[3] = sf;

      if ( !forward ) {
         cur_[0] = cell_idx_c( sz - 1 );
         cur_[1] = cell_idx_c( sy - 1 );
         cur_[2] = cell_idx_c( sx - 1 );
      }
   }


   if ( forward )
   {
      lineBegin_ = const_cast<NonConstT *>(& f_->get(xBeg,yBeg,zBeg,fBeg) );
      linePtr_   = lineBegin_;
      lineEnd_   = linePtr_ + sizes_[3];
   }
   else
   {
      linePtr_   = const_cast<NonConstT *>(& f_->get(xBeg + cell_idx_c(sx) - 1,
                                                     yBeg + cell_idx_c(sy) - 1,
                                                     zBeg + cell_idx_c(sz) - 1,
                                                     fBeg + cell_idx_c(sf) - 1) );
      lineEnd_   = linePtr_ + 1;
      lineBegin_ = linePtr_ - sizes_[3] + 1;
   }

   initCoordinateAccessOptimizationPointers();
}





//**********************************************************************************************************************
/*!\brief Constructs an end iterator, which is represented by NULL pointers
 **********************************************************************************************************************/
template <typename T, uint_t fs>
FieldIterator<T,fs>::FieldIterator()
   : linePtr_(nullptr), lineEnd_(nullptr), f_(nullptr)
{
}


//**********************************************************************************************************************
/*!\brief Copy Constructor. Required for pointer member cur*_
 **********************************************************************************************************************/
template <typename T, uint_t fs>
FieldIterator<T,fs>::FieldIterator( const FieldIterator<T,fs> & o )
   : lineBegin_     ( o.lineBegin_ ),
     linePtr_       ( o.linePtr_      ),
     lineEnd_       ( o.lineEnd_      ),
     f_             ( o.f_            ),
     xBegin_        ( o.xBegin_       ),
     yBegin_        ( o.yBegin_       ),
     zBegin_        ( o.zBegin_       ),
     fBegin_        ( o.fBegin_       )
{
   // no need to copy fastestCoord_, since it is updated before read
   for(int i=0; i<3; ++i)
      cur_[i] = o.cur_[i];

   for( int i=0; i<4; ++i ) {
      skips_[i] = o.skips_[i];
      sizes_[i] = o.sizes_[i];
   }

   if( f_ )
      initCoordinateAccessOptimizationPointers();
}

//**********************************************************************************************************************
/*!\brief Assignment operator. Required for pointer member cur*_
 **********************************************************************************************************************/
template <typename T, uint_t fs>
FieldIterator<T,fs> & FieldIterator<T,fs>::operator= ( const FieldIterator<T,fs> & o )
{
   if ( &o == this)
      return *this;

   lineBegin_ = o.lineBegin_;
   linePtr_   = o.linePtr_  ;
   lineEnd_   = o.lineEnd_  ;
   f_         = o.f_        ;
   xBegin_    = o.xBegin_   ;
   yBegin_    = o.yBegin_   ;
   zBegin_    = o.zBegin_   ;
   fBegin_    = o.fBegin_   ;

   for(int i=0; i<3; ++i)
      cur_[i] = o.cur_[i];

   for( int i=0; i<4; ++i ) {
      skips_[i] = o.skips_[i];
      sizes_[i] = o.sizes_[i];
   }

   if( f_ )
      initCoordinateAccessOptimizationPointers();

   return *this;
}


//**********************************************************************************************************************
/*!\brief Initializes pointers required for the optimized x(),y(),z(),f() functions
 *        See documentation of fastestCoord_, curX_, curY_, curZ_ and curF_
 **********************************************************************************************************************/
template <typename T, uint_t fs>
void FieldIterator<T,fs>::initCoordinateAccessOptimizationPointers( )
{
   if( f_->layout() == fzyx )
   {
      curF_ = &( cur_[0] );
      curZ_ = &( cur_[1] );
      curY_ = &( cur_[2] );
      curX_ = &( fastestCoord_ );
   }
   else
   {
      curZ_ = &( cur_[0] );
      curY_ = &( cur_[1] );
      curX_ = &( cur_[2] );
      curF_ = &( fastestCoord_ );
   }
}


//**********************************************************************************************************************
/*!\brief Increments the slower 3 coordinates, if innermost coordinate is at end
 *
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline void FieldIterator<T,fs>::incrementLine()
{
   WALBERLA_ASSERT_EQUAL( linePtr_, lineEnd_ );

   linePtr_ += skips_[3];
   cur_[2]++;

   if(cur_[2] == cell_idx_c(sizes_[2]) )
   {
      linePtr_ += skips_[2];
      cur_[2] = 0;
      cur_[1]++;
      if(cur_[1] == cell_idx_c(sizes_[1])  )
      {
         linePtr_ += skips_[1];
         cur_[1] = 0;
         cur_[0]++;
         if(cur_[0] == cell_idx_c(sizes_[0]) )
         {
            // iterator at end
            linePtr_ = nullptr;
            return;
         }
      }
   }

   lineEnd_   = linePtr_ + sizes_[3];
   lineBegin_ = linePtr_;
}


//**********************************************************************************************************************
/*!\brief Decrements the slower 3 coordinates, if innermost coordinate is at beginning
 *
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline void FieldIterator<T,fs>::decrementLine()
{
   WALBERLA_ASSERT_EQUAL( linePtr_, lineBegin_-1 );

   linePtr_ -= skips_[3];
   cur_[2]--;

   if(cur_[2] < 0 )
   {
      linePtr_ -= skips_[2];
      cur_[2] = cell_idx_c(sizes_[2])-1;
      cur_[1]--;
      if(cur_[1] < 0  )
      {
         linePtr_ -= skips_[1];
         cur_[1] = cell_idx_c(sizes_[1])-1;
         cur_[0]--;
         if(cur_[0] < 0 )
         {
            // iterator at end
            linePtr_ = nullptr;
            return;
         }
      }
   }

   lineEnd_   = linePtr_+1;
   lineBegin_ = linePtr_ - sizes_[3]+1;
}


//**********************************************************************************************************************
/*!\brief Equal operator.
 *
 * Test equality only by comparing the internal pointer
 *
 * \return true if both iterators are equal
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline bool FieldIterator<T,fs>::operator==( const FieldIterator<T,fs>& it ) const
{
   return it.linePtr_ == this->linePtr_;
}



//**********************************************************************************************************************
/*!\brief Unequal operator.
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline bool FieldIterator<T,fs>::operator!=( const FieldIterator<T,fs>& it ) const
{
   return it.linePtr_ != this->linePtr_;
}



//**********************************************************************************************************************
/*!\brief Neighbor access relative to current position
 * \param d Direction enumeration which defines deltas for x,y,z
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline T & FieldIterator<T,fs>::neighbor( stencil::Direction d, cell_idx_t cf ) const
{
   using namespace stencil;
   return neighbor(cx[d],cy[d],cz[d],cf);
}


//**********************************************************************************************************************
/*!\brief uint_t variant of above function
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline T & FieldIterator<T,fs>::neighbor( stencil::Direction d, uint_t cf ) const
{
   return neighbor( d, cell_idx_c (cf) );
}


//**********************************************************************************************************************
/*!\brief Neighbor access relative to current position
 * \param d Direction enumeration which defines deltas for x,y,z,f
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline T & FieldIterator<T,fs>::neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, cell_idx_t cf ) const
{
   T * res = linePtr_;

   res += cx * f_->xfact_ +
          cy * f_->yfact_ +
          cz * f_->zfact_ +
          cf * f_->ffact_;

   WALBERLA_ASSERT ( f_->addressInsideAllocedSpace( res ) );

   return *res;
}


//**********************************************************************************************************************
/*!\brief Neighbor variant that takes unsigned int as f parameter,
 *        needed since the stencil toIdx() is an unsigned int
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline T & FieldIterator<T,fs>::neighbor( cell_idx_t cx, cell_idx_t cy, cell_idx_t cz, uint_t cf ) const
{
   return neighbor ( cx, cy, cz, cell_idx_c( cf ) );
}



//**********************************************************************************************************************
/*!\brief For beginXYZ iterators, one often needs a specific f
 * Assumes that iterator stands at f==0
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline  T & FieldIterator<T,fs>::getF( cell_idx_t cf ) const
{
   WALBERLA_ASSERT_EQUAL( f(), 0 );
   WALBERLA_ASSERT_LESS( cf, cell_idx_t ( f_->fSize() ) );
   T * res = linePtr_;
   res += cf * f_->ffact_;
   return *res;
}


//**********************************************************************************************************************
/*!\brief Equivalent to neighbor(cell_idx_t) see above.
 *        Takes an uint_t instead a cell_idx_t, since stencil::toIndex() returns uint_t
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline  T & FieldIterator<T,fs>::getF( uint_t cf ) const
{
   return getF ( cell_idx_c ( cf ) );
}


//======================================================================================================================
//
//  PRINTING
//
//======================================================================================================================

template <typename T, uint_t fs>
inline void FieldIterator<T,fs>::print( std::ostream & os ) const
{
   os << "(" << x() << "," << y() << "," << z() << "/" << f() << ")";
}

template< typename T, uint_t fs>
std::ostream & operator<< ( std::ostream & os, const FieldIterator<T,fs> & it ) {
   it.print(os);
   return os;
}


//======================================================================================================================
//
//  COORDINATES OF CURRENT POSITIONS
//
//======================================================================================================================

/**
 * In order to get x(), y(), z(), f() function as fast a possible, no if clause for the layout was introduced.
 * Instead there are the cur[XYZF]_ members, that point to the cur_ array. The cur_ array does not store the
 * fastest coordinate, because it is implicitly stored in (linePtr_ - lineBegin_).  If it would be stored explicitly
 * there would have to be an extra update operation in operator++() , which should be as fast as possible.
 * The curX_ or curF_ pointer points to the fastestCoord_ member, which always has to be updated before curX_ or
 * curF_ is dereferenced.
 */

template <typename T, uint_t fs>
inline cell_idx_t FieldIterator<T,fs>::x() const
{
   fastestCoord_ = cell_idx_c(linePtr_ - lineBegin_ );
   return xBegin_ + *curX_;
}

template <typename T, uint_t fs>
inline cell_idx_t FieldIterator<T,fs>::y() const
{
   // no fastestCoord_ update required here, since y is never fastest coordinate
   return yBegin_ + *curY_;
}

template <typename T, uint_t fs>
inline cell_idx_t FieldIterator<T,fs>::z() const
{
   // no fastestCoord_ update required here, since z is never fastest coordinate
   return zBegin_ + *curZ_;
}

template <typename T, uint_t fs>
inline cell_idx_t FieldIterator<T,fs>::f() const
{
   fastestCoord_ = cell_idx_c(linePtr_ - lineBegin_ );
   return fBegin_ + *curF_;
}

template <typename T, uint_t fs>
inline Cell FieldIterator<T,fs>::cell() const
{
   fastestCoord_ = cell_idx_c( linePtr_ - lineBegin_ );
   return Cell ( xBegin_ + *curX_,
                 yBegin_ + *curY_,
                 zBegin_ + *curZ_ );
}



//======================================================================================================================
//
//  FORWARD ITERATOR
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\brief Pre-increment operator.
 *
 * \return Reference to the incremented pointer iterator.
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline ForwardFieldIterator<T,fs>& ForwardFieldIterator<T,fs>::operator++()
{

   ++Parent::linePtr_;
   if( Parent::linePtr_ != Parent::lineEnd_)
      return *this;

   // Iteration through line has finished - switch to next
   Parent::incrementLine();

   return *this;
}

//**********************************************************************************************************************
/*!\brief Pre-decrement operator.
 *
 * \return Reference to the decremented pointer iterator.
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline ForwardFieldIterator<T,fs>& ForwardFieldIterator<T,fs>::operator--()
{

   --Parent::linePtr_;
   if( Parent::linePtr_ >= Parent::lineBegin_ )
      return *this;

   // Iteration through line has finished - switch to next
   Parent::decrementLine();

   return *this;
}

//**********************************************************************************************************************
/*!\brief Increments the second inner coordinate c2
 * Use if innermost loop is self written, which is slightly faster
 * \code
   for( const_iterator i = field.begin(); i != field.end(); i.incrOuter() )
      for (; i.testInner(); i.incrInner() )
      {}
   Instead of
  \code
   for(DoubleField::const_iterator i = field.begin(); i != field.end(); ++i) {
   }
  \endcode
 **********************************************************************************************************************/
template <typename T,  uint_t fs>
inline void ForwardFieldIterator<T,fs>::incrOuter()
{
   // incrementing line pointer was done in "inner" iterator
   Parent::incrementLine();
}


/*
template <typename T, uint_t fs>
void ForwardFieldIterator<T,fs>::setCoordinates( cell_idx_t x, cell_idx_t y, cell_idx_t z, cell_idx_t f )
{
   typedef typename Parent::NonConstT NonConstT;
   Parent::linePtr_ = const_cast< NonConstT *>(& Parent::f_->get(x,y,z,f) );;

   if ( Parent::f_->layout() == fzyx )
      Parent::lineBegin_ = const_cast<NonConstT *>(& Parent::f_->get( Parent::xBegin_,y,z,f) );
   else
      Parent::lineBegin_ = const_cast<NonConstT *>(& Parent::f_->get(x,y,z,Parent::fBegin_) );

   Parent::lineEnd_   = Parent::lineBegin_ + Parent::sizes_[3];

   *Parent::curX_ = x - Parent::xBegin_;
   *Parent::curY_ = y - Parent::yBegin_;
   *Parent::curZ_ = z - Parent::zBegin_;
   *Parent::curF_ = f - Parent::fBegin_;

   WALBERLA_ASSERT_GREATER_EQUAL( *Parent::curX_, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( *Parent::curY_, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( *Parent::curZ_, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( *Parent::curF_, 0 );
   WALBERLA_ASSERT_GREATER_EQUAL( Parent::linePtr_, Parent::lineBegin_ );

   WALBERLA_ASSERT_LESS( Parent::cur_[0], Parent::sizes_[0] );
   WALBERLA_ASSERT_LESS( Parent::cur_[1], Parent::sizes_[1] );
   WALBERLA_ASSERT_LESS( Parent::cur_[2], Parent::sizes_[2] );
   WALBERLA_ASSERT_LESS( Parent::lineEnd_ - Parent::linePtr_, Parent::sizes_[3] );
}*/


//======================================================================================================================
//
//  REVERSE ITERATOR
//
//======================================================================================================================


//**********************************************************************************************************************
/*!\brief Pre-increment operator.
 *
 * \return Reference to the incremented pointer iterator.
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline ReverseFieldIterator<T,fs>& ReverseFieldIterator<T,fs>::operator--()
{

   ++Parent::linePtr_;
   if( Parent::linePtr_ != Parent::lineEnd_)
      return *this;

   // Iteration through line has finished - switch to next
   Parent::incrementLine();

   return *this;
}


//**********************************************************************************************************************
/*!\brief Pre-decrement operator.
 *
 * \return Reference to the decremented pointer iterator.
 **********************************************************************************************************************/
template <typename T, uint_t fs>
inline ReverseFieldIterator<T,fs>& ReverseFieldIterator<T,fs>::operator++()
{
   --Parent::linePtr_;

   if( Parent::linePtr_ >= Parent::lineBegin_ )
      return *this;


   // Iteration through line has finished - switch to next
   Parent::decrementLine();

   return *this;
}



} // namespace field
} // namespace walberla
