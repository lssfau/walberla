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
//! \file AllSetSelectableObject.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "SelectableObject.h"
#include "core/AllSet.h"

#include <algorithm>


namespace walberla {
namespace selectable {



//**********************************************************************************************************************
/*!
*   \brief Implementation of SelectableObject using sets (see "AllSet.h") as selection attributes
*
*   AllSetSelectableObject is an implementation of SelectableObject that stores objects of type T which are attached
*   with selection attributes of type `AllSet<U>`. For information on which objects are selected given a certain set as
*   selection request see the documentation of the function "select".
*/
//**********************************************************************************************************************

template< typename T, typename U >
class AllSetSelectableObject : public SelectableObject< T, AllSet<U>, AllSet<U> > {

private:

   struct Compare {
      Compare( const std::vector< AllSet<U> >& attributes ) : attributes_( attributes ) {}
      bool operator()( const size_t i, const size_t j ) const { return attributes_[i] > attributes_[j]; }
      const std::vector< AllSet<U> >& attributes_;
   };

   virtual void select( std::vector< size_t >& index, const AllSet<U>& selector ) const;

}; // AllSetSelectableObject



//**********************************************************************************************************************
/*!
*   Given a certain selection request 'R' (= "selector" [set of objects of type U]), this function returns the index 'i'
*   of every attribute set 'A_i' (A_0 <=> attributes_[0], A_1 <=> attributes_[1], etc.) that matches the request 'R'.
*   Possible candidates are:
*      - for "normal" sets: sets 'A_i' that are completely contained within the set 'R'
*      - for "all" sets (see "AllSet.h"): sets 'A_i' that completely contain the set 'R'
*   Out of these possible candidates, the set A_x that contains most elements is chosen. Additionally, out of these
*   candidates, every set that contains as many elements as A_x and every set that contains as many or even more
*   elements than the selection request 'R' is also chosen.
*/
//**********************************************************************************************************************

template< typename T, typename U >
void AllSetSelectableObject<T,U>::select( std::vector< size_t >& index, const AllSet<U>& selector ) const {

   std::vector< size_t > candidates;

   for( size_t i = 0; i != this->attributes_.size(); ++i ) {
      if( this->attributes_[i].isAll() ) {
         if( this->attributes_[i].contains( selector ) )
            candidates.push_back( i );
      }
      else if( selector.contains( this->attributes_[i] ) )
         candidates.push_back( i );
   }

   if( !candidates.empty() ) {

      Compare compare( this->attributes_ );

      std::sort( candidates.begin(), candidates.end(), compare );

      index.push_back( candidates[0] );
      for( size_t i = 1; i != candidates.size() && ( this->attributes_[ candidates[i] ] >= selector ||
                                                     this->attributes_[ candidates[i] ].equalSize( this->attributes_[ candidates[0] ] ) ); ++i )
         index.push_back( candidates[i] );
   }
}



} // namespace selectable
} // namespace walberla
