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
//! \file RigidBodyCastIterator.h
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
/*!\brief Dynamic cast iterator for polymorphic pointer vectors.
 *
 * The RigidBodyCastIterator is a forward iterator which only selects elements of type C.
 * Dereferencing this iterator will implicitly cast to C.
 */
template <typename C>
class RigidBodyCastIterator
{
   static_assert(std::is_base_of<RigidBody, C>::value && !std::is_base_of<C, RigidBody>::value,
                 "C has to be strictly derived from RigidBody");

   template <typename C2>
   friend inline bool operator==( const RigidBodyCastIterator<C2>& lhs, const RigidBodyCastIterator<C2>& rhs );
   template <typename C2>
   friend inline bool operator!=( const RigidBodyCastIterator<C2>& lhs, const RigidBodyCastIterator<C2>& rhs );
public:
   using ContainerType         = std::vector< std::unique_ptr<RigidBody> >;
   //**Type definitions****************************************************************************
   // STL iterator requirements
   using iterator_category     = std::forward_iterator_tag;
   using value_type            = C;
   using pointer               = C*;
   using reference             = C&;
   using difference_type       = std::ptrdiff_t;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline RigidBodyCastIterator() = default;
   explicit inline RigidBodyCastIterator( const typename ContainerType::iterator& begin, const typename ContainerType::iterator& end );

   RigidBodyCastIterator( const RigidBodyCastIterator<C>& it) = default;
   RigidBodyCastIterator( RigidBodyCastIterator<C>&& it) = default;

   RigidBodyCastIterator& operator=(const RigidBodyCastIterator<C>& it) = default;
   RigidBodyCastIterator& operator=(RigidBodyCastIterator<C>&& it) = default;

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   inline RigidBodyCastIterator<C>&   operator++();
   inline RigidBodyCastIterator<C>    operator++( int ) {RigidBodyCastIterator<C> tmp(*this); ++(*this); return tmp;}
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline reference operator*()   {return static_cast<reference>( *(cur_->get()) );}
   inline pointer   operator->()  {return static_cast<pointer>(     cur_->get()  );}
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline pointer getBodyID() {return static_cast<pointer>(cur_->get());}
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   typename ContainerType::iterator cur_;  //!< Pointer to the current memory location.
   typename ContainerType::iterator end_;  //!< Pointer to the element one past the last element in the element range.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Standard constructor for CastIterator.
 *
 * \param begin The beginning of the element range.
 * \param end The end of the element range.
 */
template< typename C >  // Cast type
inline RigidBodyCastIterator<C>::RigidBodyCastIterator( const typename ContainerType::iterator& begin,
                                                        const typename ContainerType::iterator& end )
   : cur_(begin)  // Pointer to the current memory location
   , end_(end)    // Pointer to the element one past the last element in the element range
{
   cur_ = std::find_if(cur_, end_, [](const ContainerType::iterator::value_type& p) {return p->getTypeID() == C::getStaticTypeID();});
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Pre-increment operator.
 *
 * \return Reference to the incremented cast iterator.
 */
template< typename C >
inline RigidBodyCastIterator<C>& RigidBodyCastIterator<C>::operator++()
{
   cur_ = std::find_if(++cur_, end_, [](const ContainerType::iterator::value_type& p) {return p->getTypeID() == C::getStaticTypeID();});

   return *this;
}
//*************************************************************************************************


//**********************************************************************************************
/*!\brief Equality comparison between two CastIterator objects.
//
// \param lhs The left hand side cast iterator.
// \param rhs The right hand side cast iterator.
// \return \a true if the iterators point to the same element, \a false if not.
*/
template< typename C >
inline bool operator==( const RigidBodyCastIterator<C>& lhs, const RigidBodyCastIterator<C>& rhs )
{
   return lhs.cur_ == rhs.cur_;
}
//**********************************************************************************************

//**********************************************************************************************
/*!\brief Inequality comparison between two CastIterator objects.
//
// \param lhs The left hand side cast iterator.
// \param rhs The right hand side cast iterator.
// \return \a true if the iterators don't point to the same element, \a false if they do.
*/
template< typename C >
inline bool operator!=( const RigidBodyCastIterator<C>& lhs, const RigidBodyCastIterator<C>& rhs )
{
   return !operator==(lhs, rhs);
}
//**********************************************************************************************





















template <typename C>
class ConstRigidBodyCastIterator
{
   static_assert(std::is_base_of<RigidBody, C>::value && !std::is_base_of<C, RigidBody>::value,
                 "C has to be strictly derived from RigidBody");

   template <typename C2>
   friend bool operator==( const ConstRigidBodyCastIterator<C2>& lhs, const ConstRigidBodyCastIterator<C2>& rhs );
   template <typename C2>
   friend bool operator!=( const ConstRigidBodyCastIterator<C2>& lhs, const ConstRigidBodyCastIterator<C2>& rhs );
public:
   using ContainerType         = std::vector< std::unique_ptr<RigidBody> >;
   //**Type definitions****************************************************************************
   // STL iterator requirements
   using iterator_category     = std::forward_iterator_tag;
   using value_type            = C const;
   using pointer               = C const *;
   using reference             = C const &;
   using difference_type       = std::ptrdiff_t;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline ConstRigidBodyCastIterator() = default;
   explicit inline ConstRigidBodyCastIterator( const typename ContainerType::const_iterator& begin,
                                               const typename ContainerType::const_iterator& end );

   ConstRigidBodyCastIterator( const ConstRigidBodyCastIterator<C>& it) = default;
   ConstRigidBodyCastIterator( ConstRigidBodyCastIterator<C>&& it) = default;

   ConstRigidBodyCastIterator& operator=(const ConstRigidBodyCastIterator<C>& it) = default;
   ConstRigidBodyCastIterator& operator=(ConstRigidBodyCastIterator<C>&& it) = default;

   //**Operators***********************************************************************************
   /*!\name Operators */
   //@{
   inline ConstRigidBodyCastIterator<C>&   operator++();
   inline ConstRigidBodyCastIterator<C>    operator++( int ) {ConstRigidBodyCastIterator<C> tmp(*this); ++(*this); return tmp;}
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline reference operator*()   {return static_cast<reference>( *(cur_->get()) );}
   inline pointer   operator->()  {return static_cast<pointer>(     cur_->get()  );}
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline pointer getBodyID() {return static_cast<pointer>(cur_->get());}
   //@}
   //**********************************************************************************************

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   typename ContainerType::const_iterator cur_;  //!< Pointer to the current memory location.
   typename ContainerType::const_iterator end_;  //!< Pointer to the element one past the last element in the element range.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Standard constructor for CastIterator.
 *
 * \param begin The beginning of the element range.
 * \param end The end of the element range.
 */
template< typename C >  // Cast type
inline ConstRigidBodyCastIterator<C>::ConstRigidBodyCastIterator( const typename ContainerType::const_iterator& begin,
                                                                  const typename ContainerType::const_iterator& end )
   : cur_(begin)  // Pointer to the current memory location
   , end_(end)    // Pointer to the element one past the last element in the element range
{
   cur_ = std::find_if(cur_, end_, [](const ContainerType::const_iterator::value_type& p) {return p->getTypeID() == C::getStaticTypeID();});
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Pre-increment operator.
 *
 * \return Reference to the incremented cast iterator.
 */
template< typename C >
inline ConstRigidBodyCastIterator<C>& ConstRigidBodyCastIterator<C>::operator++()
{
   cur_ = std::find_if(++cur_, end_, [](const ContainerType::const_iterator::value_type& p) {return p->getTypeID() == C::getStaticTypeID();});

   return *this;
}
//*************************************************************************************************


//**********************************************************************************************
/*!\brief Equality comparison between two CastIterator objects.
//
// \param lhs The left hand side cast iterator.
// \param rhs The right hand side cast iterator.
// \return \a true if the iterators point to the same element, \a false if not.
*/
template< typename C >
inline bool operator==( const ConstRigidBodyCastIterator<C>& lhs, const ConstRigidBodyCastIterator<C>& rhs )
{
   return lhs.cur_ == rhs.cur_;
}
//**********************************************************************************************

//**********************************************************************************************
/*!\brief Inequality comparison between two CastIterator objects.
//
// \param lhs The left hand side cast iterator.
// \param rhs The right hand side cast iterator.
// \return \a true if the iterators don't point to the same element, \a false if they do.
*/
template< typename C >
inline bool operator!=( const ConstRigidBodyCastIterator<C>& lhs, const ConstRigidBodyCastIterator<C>& rhs )
{
   return !operator==(lhs, rhs);
}
//**********************************************************************************************

} // namespace pe
} // namespace walberla
