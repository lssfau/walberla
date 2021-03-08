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
//! \file UID.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Feichtinger
//
//======================================================================================================================

#pragma once

#include "waLBerlaDefinitions.h"

#include "core/debug/Debug.h"
#include "core/mpi/BufferDataTypeExtensions.h"

#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include <mutex>
#include <shared_mutex>


namespace walberla {
namespace uid {



//**********************************************************************************************************************
/*!
*   \brief Class for creating, storing, and managing unique identifiers
*
*   Unique identifiers (UID) can be compared and follow a strict total order (==, !=, >, >=, <, <=). Additionally,
*   every UID is associated with a unique identifier string, meaning that two UIDs that share the same numerical
*   representation (uid_) are also associated with the same identifier string. Furthermore, it is not possible to create
*   two UIDs that have different numerical representations but are both associated with the same identifier string.
*   In other words: two UIDs / UID objects either have the same numerical representation (uid_) and the same identifier
*   string or they have a different numerical representation and a different identifier string.
*
*   In addition to being constructed via constructor, UIDs can also be constructed via two static member functions:
*      \code getOrConstruct( const std::string& identifier ) and get( const std::string& identifier ) \endcode
*   For more information on how to create UIDs see the documentation of constructors and these two functions (UIDs can
*   also be copied / copy constructed).
*
*   The template T specifies a class that is used to generate the numerical representations (uid_) for all UID<T>s.
*   T must be an implementation of a UID generator class (see "UIDGenerators.h"). Every UID generator implements a
*   static member function "uint_type generatUID()" which returns a new UID every time it is called. Additionally, UID
*   generators implement a static member function "uint_type toIndex( const uint_type uid )" which translates the
*   numerical representations of a given UID into an index which, for example, can be used for array access/storage.
*   Moreover, UID generators implement a static member function "uint_type toBitMask( const uint_type uid )" which
*   translates the numerical UID representations into a bit mask. For further information see the documentation of the
*   two member functions "toIndex()" and "toBitMask()" of this class.
*
*   The class UID is thread safe and follows a multiple-readers / single-writer pattern.
*/
//**********************************************************************************************************************

template< typename T > class UID {

public:

   using uint_type = typename T::uint_type;


            inline UID( const char* const  identifier, const bool newUid = false, const bool appendUIDtoIdentifier = false );
            inline UID( const std::string& identifier, const bool newUid = false, const bool appendUIDtoIdentifier = false );
   explicit        UID( const bool createAnonymousUID = false );

   static inline bool exists( const std::string& identifier ) { return stringToUid().find( identifier ) != stringToUid().end(); }

   static        UID<T> getOrConstruct( const std::string& identifier );
   static inline UID<T> get           ( const std::string& identifier );

   inline       uint_type    getUid()        const { WALBERLA_ASSERT( valid_ ); return uid_; }
   inline const std::string& getIdentifier() const;

   //*******************************************************************************************************************
   /*!
   *   Returns '0' for the first UID created by T::generateUID(), '1' for the second, '2' for the third, etc. The value
   *   returned by this function may or may not be equal to "uid_" - depending on the underlying UID generator T.
   */
   //*******************************************************************************************************************
   uint_type toIndex() const { WALBERLA_ASSERT( valid_ ); return T::toIndex( uid_ ); }

   //*******************************************************************************************************************
   /*!
   *   Returns '[...] 0001' for the first UID created by T::generateUID(), '[...] 0010' for the second, '[...] 0100' for
   *   the third, etc. The value returned by this function may or may not be equal to "uid_" - depending on the
   *   underlying UID generator T.
   */
   //*******************************************************************************************************************
   uint_type toBitMask() const { WALBERLA_ASSERT( valid_ ); return T::toBitMask( uid_ ); }

   inline void        toStream( std::ostream& os ) const;
   inline std::string toString() const;

   inline bool operator==( const UID<T>& uid ) const { WALBERLA_ASSERT( valid_ ); return uid_ == uid.getUid(); }
   inline bool operator!=( const UID<T>& uid ) const { WALBERLA_ASSERT( valid_ ); return uid_ != uid.getUid(); }
   inline bool operator> ( const UID<T>& uid ) const { WALBERLA_ASSERT( valid_ ); return uid_ >  uid.getUid(); }
   inline bool operator>=( const UID<T>& uid ) const { WALBERLA_ASSERT( valid_ ); return uid_ >= uid.getUid(); }
   inline bool operator< ( const UID<T>& uid ) const { WALBERLA_ASSERT( valid_ ); return uid_ <  uid.getUid(); }
   inline bool operator<=( const UID<T>& uid ) const { WALBERLA_ASSERT( valid_ ); return uid_ <= uid.getUid(); }

private:

   explicit inline UID( const uint_type uid ); // must not be public! (used internally by UID::getOrConstruct() & UID::get())

   void init( const std::string& identifier, const bool newUid, const bool appendUIDtoIdentifier );

   static std::map< uint_type, std::string >& uidToString() { static std::map< uint_type, std::string > map; return map; }
   static std::map< std::string, uint_type >& stringToUid() { static std::map< std::string, uint_type > map; return map; }

   static std::shared_timed_mutex& uidToStringMutex() { static std::shared_timed_mutex mutex; return mutex; }
   static std::shared_timed_mutex& stringToUidMutex() { static std::shared_timed_mutex mutex; return mutex; }

   uint_type uid_;

#ifndef NDEBUG
   bool   valid_;
#endif

}; // class UID



//**********************************************************************************************************************
/*!
*   For documentation see member function 'init'.
*
*   Do not delete this constructor! This constructor is required, so that writing "UID<TYPE> uid("id")" does not lead to
*   calling the constructor "UID( const bool createAnonymousUID )".
*/
//**********************************************************************************************************************
template< typename T >
inline UID<T>::UID( const char* const identifier, const bool newUid, const bool appendUIDtoIdentifier )
#ifndef NDEBUG
   : valid_( true )
#endif
{
   init( std::string(identifier), newUid, appendUIDtoIdentifier );
}



//**********************************************************************************************************************
/*!
*   For documentation see member function 'init'.
*/
//**********************************************************************************************************************
template< typename T >
inline UID<T>::UID( const std::string& identifier, const bool newUid, const bool appendUIDtoIdentifier )
#ifndef NDEBUG
   : valid_( true )
#endif
{
   init( identifier, newUid, appendUIDtoIdentifier );
}



//**********************************************************************************************************************
/*!
*   This function is called by two constructors in order to create new UIDs or to retrieve copies of existing UIDs.
*   If newUid is true, a completely new UID is created and associated with the string "identifier". In debug mode,
*   creation of a new UID will trigger an assertion and fail if another UID with the same identifier string already
*   exists.
*   If newUid is false (-> default value, can be omitted), a copy of an existing UID is returned if the identifier
*   string already exists, otherwise a new UID is created and associated with "identifier" (this behavior is identical
*   to the static member function "getOrConstruct").
*   If newUid is true AND appendUIDtoIdentifier is also true, the newly created UID is appended to the provided
*   identifier string.
*
*   Must not be made public! Do not use this function. This function is intended for internal use only.
*/
//**********************************************************************************************************************
template< typename T >
void UID<T>::init( const std::string& identifier, const bool newUid, const bool appendUIDtoIdentifier ) {

   std::unique_lock< std::shared_timed_mutex > stringToUid_lock( stringToUidMutex() );

   WALBERLA_ASSERT( !identifier.empty() );

   if( newUid ) { // another UID with the same identifier must not exist

      uid_ = T::generateUID();

      std::string idString = ( appendUIDtoIdentifier ? ( identifier + " [" + std::to_string( uid_ ) + "]" ) : identifier );

      WALBERLA_ASSERT( stringToUid().find( idString ) == stringToUid().end() ); // 'idString' must not exist

      std::unique_lock< std::shared_timed_mutex > uidToString_lock( uidToStringMutex() );

      uidToString()[ uid_ ] = idString;
      stringToUid()[ idString ] = uid_;
   }
   else { // if another UID with the same identifier already exists this UID is chosen, otherwise a new UID is created

      typename std::map< std::string, uint_type >::iterator it = stringToUid().find( identifier );

      if( it == stringToUid().end() ) {

         uid_ = T::generateUID();

         std::unique_lock< std::shared_timed_mutex > uidToString_lock( uidToStringMutex() );
         uidToString()[ uid_ ] = identifier;
         stringToUid()[ identifier ] = uid_;
      }
      else
         uid_ = it->second;
   }
}



//**********************************************************************************************************************
/*!
*   The default constructor (-> createAnonymousUID is false by default!) is used to create uninitialized, invalid UID
*   objects (by just declaring a variable without initializing it at declaration: "UID<T> uid; [...]; uid = [...];").
*   In debug mode, every operation involving such an uninitialized UID will trigger an assertion and fail.
*   If createAnonymousUID is true, a valid UID with an unspecified identifier string is created.
*/
//**********************************************************************************************************************
template< typename T >
UID<T>::UID( const bool createAnonymousUID ) : uid_( 0 ) {

   if( !createAnonymousUID ) {
#ifndef NDEBUG
      valid_ = false;
#endif
   }
   else {
      std::unique_lock< std::shared_timed_mutex > stringToUid_lock( stringToUidMutex() );

      uid_ = T::generateUID();

      std::string identifier = std::string("[anonymous #") + std::to_string( T::toIndex( uid_ ) ) + std::string("]");

      WALBERLA_ASSERT( stringToUid().find( identifier ) == stringToUid().end() ); // 'identifier' must not exist

      std::unique_lock< std::shared_timed_mutex > uidToString_lock( uidToStringMutex() );

      uidToString()[ uid_ ] = identifier;
      stringToUid()[ identifier ] = uid_;

#ifndef NDEBUG
      valid_ = true;
#endif
   }
}



//**********************************************************************************************************************
/*!
*   Must not be made public! Do not use this constructor. This constructor is intended for internal use only.
*/
//**********************************************************************************************************************
template< typename T >
inline UID<T>::UID( const uint_type uid ) : uid_( uid )
#ifndef NDEBUG
   , valid_( true )
#endif
   {}



//**********************************************************************************************************************
/*!
*   Either returns a copy of an existing UID or creates and returns a completely new UID depending on whether a UID
*   that is associated with the identifier string "identifier" already exists or not.
*/
//**********************************************************************************************************************
template< typename T >
UID<T> UID<T>::getOrConstruct( const std::string& identifier ) {
   std::unique_lock< std::shared_timed_mutex > stringToUid_lock( stringToUidMutex() );

   WALBERLA_ASSERT( !identifier.empty() );

   typename std::map< std::string, uint_type >::iterator it = stringToUid().find( identifier );

   if( it == stringToUid().end() ) {

      uint_type uid = T::generateUID();

      std::unique_lock< std::shared_timed_mutex > uidToString_lock( uidToStringMutex() );

      uidToString()[ uid ] = identifier;
      stringToUid()[ identifier ] = uid;

      return UID( uid );
   }

   return UID( it->second );
}



//**********************************************************************************************************************
/*!
*   Returns the UID that is associated with the identifier string "identifier". In debug mode, this function will
*   trigger an assertion and fail if no UID exists that is associated with this identifier string.
*/
//**********************************************************************************************************************
template< typename T >
inline UID<T> UID<T>::get( const std::string& identifier ) {
   std::shared_lock< std::shared_timed_mutex > lock( stringToUidMutex() );

   WALBERLA_ASSERT( stringToUid().find( identifier ) != stringToUid().end() ); // 'identifier' must exist

   return UID( stringToUid().find( identifier )->second );
}



template< typename T >
inline const std::string& UID<T>::getIdentifier() const {
   std::shared_lock< std::shared_timed_mutex > lock( uidToStringMutex() );

   WALBERLA_ASSERT( valid_ );
   WALBERLA_ASSERT( uidToString().find( uid_ ) != uidToString().end() );

   return uidToString().find( uid_ )->second;
}



template< typename T >
inline void UID<T>::toStream( std::ostream& os ) const {

   //os << "( " << getIdentifier() << " = " << uid_ << " | UID generator: " << T::getType() << " )";
   os << getIdentifier();
}



template< typename T >
inline std::string UID<T>::toString() const {

   std::ostringstream oss;
   toStream( oss );
   return oss.str();
}



//////////////////////
// Global Functions //
//////////////////////



template< typename T >
inline std::ostream& operator<<( std::ostream& os, const UID<T>& uid ) {

   uid.toStream( os );
   return os;
}



} // namespace uid

using uid::UID;

} // namespace walberla



//======================================================================================================================
//
//  Send/Recv Buffer Serialization Specialization
//
//======================================================================================================================

namespace walberla {
namespace mpi {

template< typename T, // Element type of SendBuffer
          typename G, // Growth policy of SendBuffer
          typename GE >
inline mpi::GenericSendBuffer<T,G> & operator<<( mpi::GenericSendBuffer<T,G> & buffer, const walberla::uid::UID<GE> & uid )
{
   buffer.addDebugMarker( "ui" );
   buffer << uid.getIdentifier();
   return buffer;
}

template< typename T, // Element type  of RecvBuffer
          typename GE >
inline mpi::GenericRecvBuffer<T>& operator>>( mpi::GenericRecvBuffer<T> & buffer, walberla::uid::UID<GE> & uid )
{
   buffer.readDebugMarker( "ui" );
   std::string identifier;
   buffer >> identifier;
   uid = walberla::uid::UID<GE>::getOrConstruct( identifier );
   return buffer;
}

template< typename GE >
struct BufferSizeTrait< walberla::uid::UID<GE> > { static const bool constantSize = false;  };

}
}
