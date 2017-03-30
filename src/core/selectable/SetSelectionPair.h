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
//! \file SetSelectionPair.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Set.h"


namespace walberla {
namespace selectable {



//**********************************************************************************************************************
/*!
*   \brief Class for storing two sets: one marked as "included" and one marked as "excluded"
*
*   Objects of type SetSelectionPair are used as attributes for objects stored in a "wrapper" object of type
*   SetSelectableObject. For further information see class SetSelectableObject.
*/
//**********************************************************************************************************************

template< typename T >
class SetSelectionPair {

public:

   friend inline bool operator==( const SetSelectionPair& a, const SetSelectionPair& b ) { return  a.isEqual(b); }
   friend inline bool operator!=( const SetSelectionPair& a, const SetSelectionPair& b ) { return !a.isEqual(b); }

   SetSelectionPair( const Set<T>& include, const Set<T>& exclude = Set<T>::emptySet() ) : include_( include ), exclude_( exclude ) {}

   inline const Set<T>& included() const { return include_; }
   inline const Set<T>& excluded() const { return exclude_; }

   inline bool isEqual( const SetSelectionPair<T>& s ) const { return include_ == s.included() && exclude_ == s.excluded(); }

   inline void        toStream( std::ostream& os ) const { os << "{ included = " << include_ << ", excluded = " << exclude_ << " }"; }
   inline std::string toString() const { std::ostringstream oss; toStream( oss ); return oss.str(); }

private:

   Set<T> include_;
   Set<T> exclude_;

}; // SetSelectionPair



template< typename T >
inline std::ostream& operator<<( std::ostream& os, const SetSelectionPair<T>& setSelectionPair ) {

   setSelectionPair.toStream( os );
   return os;
}



} // namespace selectable
} // namespace walberla
