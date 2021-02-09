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
//! \file SelectableObject.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/debug/Debug.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>


namespace walberla {
namespace selectable {



//**********************************************************************************************************************
/*!
*   \brief Abstract base class for storing and selecting objects of type T which are attached with selection attributes
*          of type A
*
*   The sole purpose of this class is to store selectable objects of type T, meaning objects which are attached with
*   additional selection attributes of type A, in order to later retrieve only those objects which match certain
*   requested attributes. For the retrieval (see functions "get") a selection object of type S is used (type S may or
*   may not be equal to type A).
*
*   Requirements for type A:
*
*      - objects of type A must be comparable with operator "!="
*      - A must implement operator "<<" for output via a std::ostream object
*
*   Additionally, it must be possible to store objects of both types T and A in a std::vector.
*
*   Objects can be added by calling the function "add" (for every object that is added an identifier string can be
*   provided) and be retrieved by calling "get". For additional information refer to the documentation of these
*   functions.
*
*   Every class that derives from SelectableObject must implement the function
*      \code void select( std::vector< size_t >& index, const S& selector ) const \endcode
*   After returning from the call to this function, the vector "index" must contain the index of every object in
*   "attributes_" that matches with the provided selector "selector" (for an actual implementation of SelectableObject
*   see "SetSelectableObject.h").
*/
//**********************************************************************************************************************

template< typename T, typename A, typename S >
class SelectableObject {

public:

   //**Iterators********************************************************************************************************
   /*! \name Iterators */
   //@{
   class const_iterator;
   friend class const_iterator;

   class iterator;
   friend class iterator;
   class iterator {

      friend class const_iterator;
      friend class SelectableObject;

   public:

      iterator( const iterator& it ) : selectable_( it.selectable_ ), index_( it.index_ ) {}

      iterator& operator++() { ++index_; return *this; }                            // prefix ++X
      iterator  operator++(int) { iterator it( *this ); operator++(); return it; }; // postfix X++

      bool operator==( const iterator& rhs ) const { return selectable_ == rhs.selectable_ && index_ == rhs.index_; }
      bool operator!=( const iterator& rhs ) const { return !operator==(rhs); }

      T *                        get() { WALBERLA_ASSERT_LESS( index_, selectable_->object_.size() );     return &(selectable_->object_[index_]); }
      const std::string & identifier() { WALBERLA_ASSERT_LESS( index_, selectable_->identifier_.size() ); return selectable_->identifier_[index_]; }

      T & operator*()  { WALBERLA_ASSERT_LESS( index_, selectable_->object_.size() ); return   selectable_->object_[index_]; }
      T * operator->() { WALBERLA_ASSERT_LESS( index_, selectable_->object_.size() ); return &(selectable_->object_[index_]); }

   private:

      iterator( SelectableObject* const selectable, const size_t index ) : selectable_( selectable ), index_( index ) {}

      SelectableObject * const selectable_;
      size_t index_;
   };

   class const_iterator {

      friend class SelectableObject;

   public:

      const_iterator( const       iterator& it ) : selectable_( it.selectable_ ), index_( it.index_ ) {}
      const_iterator( const const_iterator& it ) : selectable_( it.selectable_ ), index_( it.index_ ) {}

      const_iterator& operator++() { ++index_; return *this; }                                  // prefix ++X
      const_iterator  operator++(int) { const_iterator it( *this ); operator++(); return it; }; // postfix X++

      bool operator==( const const_iterator& rhs ) const { return selectable_ == rhs.selectable_ && index_ == rhs.index_; }
      bool operator!=( const const_iterator& rhs ) const { return !operator==(rhs); }

      const T *                  get() const { WALBERLA_ASSERT_LESS( index_, selectable_->object_.size() );     return &(selectable_->object_[index_]); }
      const std::string & identifier() const { WALBERLA_ASSERT_LESS( index_, selectable_->identifier_.size() ); return selectable_->identifier_[index_]; }

      const T & operator*()  const { WALBERLA_ASSERT_LESS( index_, selectable_->object_.size() ); return   selectable_->object_[index_]; }
      const T * operator->() const { WALBERLA_ASSERT_LESS( index_, selectable_->object_.size() ); return &(selectable_->object_[index_]); }

   private:

      const_iterator( const SelectableObject* const selectable, const size_t index ) : selectable_( selectable ), index_( index ) {}

      const SelectableObject * const selectable_;
      size_t index_;
   };
   //@}
   //*******************************************************************************************************************



   virtual ~SelectableObject() = default;

   void add( const T& object, const A& attributes, const std::string& identifier = std::string() );

   iterator       begin()       { return       iterator( this, 0 ); }
   const_iterator begin() const { return const_iterator( this, 0 ); }

   iterator       end()       { return       iterator( this, object_.size() ); }
   const_iterator end() const { return const_iterator( this, object_.size() ); }

   size_t get(              T& object, const S& selector ) const;
   void   get( std::vector<T>& object, const S& selector ) const;

   size_t get(              T& object,                std::string& identifier, const S& selector ) const;
   void   get( std::vector<T>& object, std::vector< std::string >& identifier, const S& selector ) const;

         T* getUnique( const S& selector );
   const T* getUnique( const S& selector ) const;
         T* getUnique( const S& selector, std::string & identifierOut );
   const T* getUnique( const S& selector, std::string & identifierOut ) const;

          void        toStream( std::ostream& os ) const;
   inline std::string toString() const;

   size_t size() const { WALBERLA_ASSERT_EQUAL( object_.size(), identifier_.size() );   return object_.size();  }
   bool  empty() const { WALBERLA_ASSERT_EQUAL( object_.empty(), identifier_.empty() ); return object_.empty(); }

private:

   //*******************************************************************************************************************
   /*!
   *   Every class that derives from SelectableObject must implement this function. After returning from the call to
   *   this function, the vector "index" must contain the index of every object in "attributes_" that matches with the
   *   provided selector "selector".
   */
   //*******************************************************************************************************************
   virtual void select( std::vector< size_t >& index, const S& selector ) const = 0;

   std::vector< T >           object_;
   std::vector< std::string > identifier_;

protected:

   std::vector< A >           attributes_;

}; // class SelectableObject



//**********************************************************************************************************************
/*!
*   This function is used to add an object together with its selection attributes stored in "attributes". Optionally,
*   an identifier string can be provided which is used during output (see "toStream" and "toString"). In debug mode,
*   this function triggers an assertion and fails if another object with the same selection attributes already exists.
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
void SelectableObject<T,A,S>::add( const T& object, const A& attributes, const std::string& identifier ) {

   object_.push_back( object );
   identifier_.push_back( identifier.empty() ? std::string( "[anonymous]" ) : identifier );
   attributes_.push_back( attributes );
}



//**********************************************************************************************************************
/*!
*   This function can be used to retrieve the one object whose attributes match with "selector". Depending on "selector"
*   and the actual implementation of the function "select", no object, one object, or multiple objects may be found. If
*   only one object is found, it is returned via the parameter "object". If multiple objects are found, the first object
*   whose attributes match is returned via "object". In any case, the number of objects that match with "selector" is
*   returned by this function - and only if the return value is equal to '1' exactly one object was found and stored in
*   "object".
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
size_t SelectableObject<T,A,S>::get( T& object, const S& selector ) const {

   std::vector< size_t > index;

   select( index, selector );

   if( !index.empty() ) {
      WALBERLA_ASSERT_LESS( index[0], object_.size() );
      object = object_[ index[0] ];
   }

   return index.size();
}



//**********************************************************************************************************************
/*!
*   This function can be used to retrieve all objects whose attributes match with the selector "selector". The objects
*   are returned using the provided vector "objects" (this function may return none, exactly one, ore multiple objects
*   depending on "selector" and the actual implementation of the function "select").
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
void SelectableObject<T,A,S>::get( std::vector<T>& object, const S& selector ) const {

   std::vector< size_t > index;

   select( index, selector );

   for( size_t i = 0; i != index.size(); ++i ) {
      WALBERLA_ASSERT_LESS( index[i], object_.size() );
      object.push_back( object_[ index[i] ] );
   }
}



//**********************************************************************************************************************
/*!
*   This function can be used to retrieve the one object whose attributes match with "selector". Depending on "selector"
*   and the actual implementation of the function "select", no object, one object, or multiple objects may be found. If
*   only one object is found, it is returned via the parameter "object". If multiple objects are found, the first object
*   whose attributes match is returned via "object". In any case, the number of objects that match with "selector" is
*   returned by this function - and only if the return value is equal to '1' exactly one object was found and stored in
*   "object". Additionally, the corresponding identifier is also returned via the parameter "identifier".
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
size_t SelectableObject<T,A,S>::get( T& object, std::string& identifier, const S& selector ) const {

   std::vector< size_t > index;

   select( index, selector );

   if( !index.empty() ) {
      WALBERLA_ASSERT_LESS( index[0], object_.size() );
      object = object_[ index[0] ];
      identifier = identifier_[ index[0] ];
   }

   return index.size();
}



//**********************************************************************************************************************
/*!
*   This function can be used to retrieve all objects whose attributes match with the selector "selector". The objects
*   are returned using the provided vector "objects" (this function may return none, exactly one, ore multiple objects
*   depending on "selector" and the actual implementation of the function "select"). Additionally, for every object that
*   is found, its corresponding identifier is also returned via "identifier".
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
void SelectableObject<T,A,S>::get( std::vector<T>& object, std::vector< std::string >& identifier, const S& selector ) const {

   std::vector< size_t > index;

   select( index, selector );

   for( size_t i = 0; i != index.size(); ++i ) {
      WALBERLA_ASSERT_LESS( index[i], object_.size() );
      object.push_back( object_[ index[i] ] );
      identifier.push_back( identifier_[ index[i] ] );
   }
}



//**********************************************************************************************************************
/*!
*   This function can be used to retrieve a pointer (!) to the object whose attributes match with "selector". Depending
*   on "selector" and the actual implementation of the function "select", no object, one object, or multiple objects may
*   be found. If no objects or multiple objects are found, NULL is returned.
*   Attention: The pointer may get invalidated by subsequent calls to the member function "add".
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
const T* SelectableObject<T,A,S>::getUnique( const S& selector ) const {

   std::vector< size_t > index;

   select( index, selector );

   if( index.size() == 1 ) {
      WALBERLA_ASSERT_LESS( index[0], object_.size() );
      return &(object_[ index[0] ]);
   }

   return nullptr;
}


//**********************************************************************************************************************
/*!
*   Non-const version of getUnique().
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
T* SelectableObject<T,A,S>::getUnique( const S& selector ) {

   const SelectableObject<T,A,S>& const_this = *this;
   return const_cast<T*>( const_this.getUnique( selector ) );
}



//**********************************************************************************************************************
/*!
*   This function is similar to getUnique(const S &selector) but takes an additional output parameter
*   identifier, which is set to the string identifier of the selected object, if a non-zero pointer is
*   returned, otherwise the identifier string is not changed
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
const T* SelectableObject<T,A,S>::getUnique( const S& selector, std::string & identifier ) const {

   std::vector< size_t > index;

   select( index, selector );

   if( index.size() == 1 ) {
      WALBERLA_ASSERT_LESS( index[0], object_.size() );
      identifier = identifier_[ index[0] ];
      return &(object_[ index[0] ]);
   }

   return nullptr;
}


//**********************************************************************************************************************
/*!
*   Non-const version of getUnique(const S&, std::string& ).
*/
//**********************************************************************************************************************
template< typename T, typename A, typename S >
T* SelectableObject<T,A,S>::getUnique( const S& selector, std::string & identifier ) {

   const SelectableObject<T,A,S>& const_this = *this;
   return const_cast<T*>( const_this.getUnique( selector,identifier ) );
}





template< typename T, typename A, typename S >
void SelectableObject<T,A,S>::toStream( std::ostream& os ) const {

   WALBERLA_ASSERT_EQUAL( identifier_.size(), attributes_.size() );

   os << "[ ";

   for( size_t i = 0; i != identifier_.size(); ++i )
      os << "( " << identifier_[i] << " : " << attributes_[i] << ( ( i+1 == identifier_.size() ) ? " ) " : " ), " );

   os << "]";
}



template< typename T, typename A, typename S >
inline std::string SelectableObject<T,A,S>::toString() const {

   std::ostringstream oss;
   toStream( oss );

   return oss.str();
}



//////////////////////
// Global Functions //
//////////////////////



template< typename T, typename A, typename S >
inline std::ostream& operator<<( std::ostream& os, const SelectableObject<T,A,S>& selectableObject ) {

   selectableObject.toStream( os );
   return os;
}



} // namespace selectable
} // namespace walberla
