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
//! \file SetSelectableObject.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "IsSetSelected.h"
#include "SelectableObject.h"
#include "SetSelectionPair.h"
#include "core/Set.h"


namespace walberla {
namespace selectable {



//**********************************************************************************************************************
/*!
*   \brief Implementation of SelectableObject using objects of type SetSelectionPair (see "SetSelectionPair.h") as
*          selection attributes, and objects of type Set (see "Set.h") as selector
*
*   SetSelectableObject is an implementation of SelectableObject that stores objects of type T which are attached with
*   selection attributes of type SetSelectionPair<U>. Every object of type SetSelectionPair<U> contains two sets of
*   type Set<U> - one "include" and one "exclude" set. For information on which objects are selected given a certain
*   selection request see the documentation of the function "select".
*/
//**********************************************************************************************************************

template< typename T, typename U >
class SetSelectableObject : public SelectableObject< T, SetSelectionPair<U>, Set<U> >
{
public:

   SetSelectableObject() = default;

   SetSelectableObject( const T& object, const Set<U>& include, const Set<U>& exclude, const std::string& identifier = std::string() ) {

      add( object, include, exclude, identifier );
   }

   ~SetSelectableObject() override = default;

   void add( const T& object, const Set<U>& include, const Set<U>& exclude, const std::string& identifier = std::string() ) {

      SelectableObject< T, SetSelectionPair<U>, Set<U> >::add( object, SetSelectionPair<U>( include, exclude ), identifier );
   }

private:

   struct Compare {
      Compare( const std::vector< SetSelectionPair<U> >& attributes ) : attributes_( attributes ) {}
      bool operator()( const size_t i, const size_t j ) const
         { return attributes_[i].included() > attributes_[j].included(); }
      const std::vector< SetSelectionPair<U> >& attributes_;
   };

   // added inline qualifier to suppress unjustified MSVC warning C4505
   inline void select( std::vector< size_t >& index, const Set<U>& selector ) const override;

}; // SetSelectableObject



//**********************************************************************************************************************
/*!
*   Given a certain selection request 'R' (= "selector" [set of objects of type U]), this function returns the index 'i'
*   of every attribute 'A_i' (A_0 <=> attributes_[0], A_1 <=> attributes_[1], etc.) that matches the request 'R'.
*   Possible candidates are all A_i ...
*
*      - ... whose "include" sets are completely contained within the set R and ...
*      - ... whose "exclude" sets do not intersect with R.
*
*   Out of these possible candidates, the attribute A_x whose "include" set contains most elements is chosen.
*   Additionally, out of all candidates, every A_i whose "include" set contains as many elements as A_x is also chosen.
*/
//**********************************************************************************************************************

template< typename T, typename U >
void SetSelectableObject<T,U>::select( std::vector< size_t >& index, const Set<U>& selector ) const {

   std::vector< size_t > candidates;

   for( size_t i = 0; i != this->attributes_.size(); ++i ) {

      if( isSetSelected( selector, this->attributes_[i].included(), this->attributes_[i].excluded() ) )
         candidates.push_back( i );
   }

   if( !candidates.empty() ) {

      Compare compare( this->attributes_ );

      std::sort( candidates.begin(), candidates.end(), compare );

      index.push_back( candidates[0] );
      for( size_t i = 1; i != candidates.size() &&
                         this->attributes_[ candidates[i] ].included().equalSize( this->attributes_[ candidates[0] ].included() ); ++i )
         index.push_back( candidates[i] );
   }
}



} // namespace selectable
} // namespace walberla
