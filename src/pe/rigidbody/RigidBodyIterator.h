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
//! \file RigidBodyIterator.h
//! \author Klaus Iglberger
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "RigidBody.h"

#include <memory>
#include <type_traits>
#include <vector>

namespace walberla {
namespace pe {

//*************************************************************************************************
/*!\brief Implementation of an iterator for pointer vectors.
 *
 * This iterator implementation will automatically dereference twice at dereference calls!
 * Therefore the return types are not unique_ptrs but RigidBodys.
 */
class RigidBodyIterator
{
   friend inline bool operator==( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs );
   friend inline bool operator!=( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs );
   friend inline bool operator>( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs );
   friend inline bool operator<( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs );
   friend inline bool operator>=( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs );
   friend inline bool operator<=( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs );
public:
   using ContainerType         = std::vector< std::unique_ptr<RigidBody> >;
   //**Type definitions****************************************************************************
   // STL iterator requirements
   using iterator_category     = std::random_access_iterator_tag;
   using value_type            = RigidBody;
   using pointer               = RigidBody*;
   using reference             = RigidBody&;
   using difference_type       = std::ptrdiff_t;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline RigidBodyIterator() = default;
   explicit inline RigidBodyIterator( const typename ContainerType::iterator& it ) : it_(it) {}

   RigidBodyIterator( const RigidBodyIterator& it) = default;
   RigidBodyIterator( RigidBodyIterator&& it) = default;

   RigidBodyIterator& operator=(const RigidBodyIterator& it) = default;
   RigidBodyIterator& operator=(RigidBodyIterator&& it) = default;

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   inline RigidBodyIterator&   operator++()                    {++it_; return *this;}
   inline RigidBodyIterator    operator++( int )               {RigidBodyIterator tmp(*this); ++(*this); return tmp;}
   inline RigidBodyIterator&   operator--()                    {--it_; return *this;}
   inline RigidBodyIterator    operator--( int )               {RigidBodyIterator tmp(*this); --(*this); return tmp;}
   inline RigidBodyIterator&   operator+=( difference_type n ) {it_ += n; return *this;}
   inline RigidBodyIterator    operator+ ( difference_type n ) const {return RigidBodyIterator( it_ + n );}
   inline RigidBodyIterator&   operator-=( difference_type n ) {it_ -= n; return *this;}
   inline RigidBodyIterator    operator- ( difference_type n ) const {return RigidBodyIterator( it_ - n );}
   inline difference_type      operator- ( const RigidBodyIterator& it ) const {return it_ - it.it_;}
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline reference operator[]( difference_type n ) const {return *it_[n];}
   inline reference operator*()                     const {return **it_;}
   inline pointer   operator->()                    const {return it_->get();}
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline pointer getBodyID() {return it_->get();}
   inline ContainerType::iterator get() const {return it_;}
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   typename ContainerType::iterator it_;  //!< wrapped iterator
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the iterators point to the same element, \a false if not.
 */
inline bool operator==( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs )
{
   return lhs.it_ == rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the iterators don't point to the same element, \a false if they do.
 */
inline bool operator!=( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs )
{
   return !operator==(lhs, rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a lower element, \a false if not.
 */
inline bool operator<( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs )
{
   return lhs.it_ < rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a higher element, \a false if not.
 */
inline bool operator>( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs )
{
   return lhs.it_ > rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a lower or the same element, \a false if not.
 */
inline bool operator<=( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs )
{
   return lhs.it_ <= rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a higher or the same element, \a false if not.
 */
inline bool operator>=( const RigidBodyIterator& lhs, const RigidBodyIterator& rhs )
{
   return lhs.it_ >= rhs.it_;
}
//*************************************************************************************************












class ConstRigidBodyIterator
{
   friend inline bool operator==( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs );
   friend inline bool operator!=( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs );
   friend inline bool operator>( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs );
   friend inline bool operator<( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs );
   friend inline bool operator>=( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs );
   friend inline bool operator<=( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs );
public:
   using ContainerType         = std::vector< std::unique_ptr<RigidBody> >;
   //**Type definitions****************************************************************************
   // STL iterator requirements
   using iterator_category     = std::random_access_iterator_tag;
   using value_type            = RigidBody const;
   using pointer               = RigidBody const *;
   using reference             = RigidBody const &;
   using difference_type       = std::ptrdiff_t;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline ConstRigidBodyIterator() = default;
   inline ConstRigidBodyIterator( const RigidBodyIterator& it ) : it_(it.get()) {}
   explicit inline ConstRigidBodyIterator( const typename ContainerType::iterator& it ) : it_(it) {}
   explicit inline ConstRigidBodyIterator( const typename ContainerType::const_iterator& it ) : it_(it) {}

   ConstRigidBodyIterator( const ConstRigidBodyIterator& it) = default;
   ConstRigidBodyIterator( ConstRigidBodyIterator&& it) = default;

   ConstRigidBodyIterator& operator=(const ConstRigidBodyIterator& it) = default;
   ConstRigidBodyIterator& operator=(ConstRigidBodyIterator&& it) = default;

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   inline ConstRigidBodyIterator&   operator++()                    {++it_; return *this;}
   inline ConstRigidBodyIterator    operator++( int )               {ConstRigidBodyIterator tmp(*this); ++(*this); return tmp;}
   inline ConstRigidBodyIterator&   operator--()                    {--it_; return *this;}
   inline ConstRigidBodyIterator    operator--( int )               {ConstRigidBodyIterator tmp(*this); --(*this); return tmp;}
   inline ConstRigidBodyIterator&   operator+=( difference_type n ) {it_ += n; return *this;}
   inline ConstRigidBodyIterator    operator+ ( difference_type n ) const {return ConstRigidBodyIterator( it_ + n );}
   inline ConstRigidBodyIterator&   operator-=( difference_type n ) {it_ -= n; return *this;}
   inline ConstRigidBodyIterator    operator- ( difference_type n ) const {return ConstRigidBodyIterator( it_ - n );}
   inline difference_type operator- ( const ConstRigidBodyIterator& it ) const {return it_ - it.it_;}
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline reference operator[]( difference_type n ) const {return *it_[n];}
   inline reference operator*()                     const {return **it_;}
   inline pointer   operator->()                    const {return it_->get();}
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline pointer getBodyID() const {return it_->get();}
   inline ContainerType::const_iterator get() const {return it_;}
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   typename ContainerType::const_iterator it_;  //!< wrapped iterator
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the iterators point to the same element, \a false if not.
 */
inline bool operator==( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs )
{
   return lhs.it_ == rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the iterators don't point to the same element, \a false if they do.
 */
inline bool operator!=( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs )
{
   return !operator==(lhs, rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a lower element, \a false if not.
 */
inline bool operator<( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs )
{
   return lhs.it_ < rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a higher element, \a false if not.
 */
inline bool operator>( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs )
{
   return lhs.it_ > rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a lower or the same element, \a false if not.
 */
inline bool operator<=( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs )
{
   return lhs.it_ <= rhs.it_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two PtrIterator objects.
 *
 * \param lhs The left-hand side pointer iterator.
 * \param rhs The right-hand side pointer iterator.
 * \return \a true if the left-hand side iterator points to a higher or the same element, \a false if not.
 */
inline bool operator>=( const ConstRigidBodyIterator& lhs, const ConstRigidBodyIterator& rhs )
{
   return lhs.it_ >= rhs.it_;
}
//*************************************************************************************************

} // namespace pe
} // namespace walberla
