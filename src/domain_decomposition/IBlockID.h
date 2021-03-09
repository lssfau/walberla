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
//! \file IBlockID.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include <ostream>

#include "core/DataTypes.h"

namespace walberla {
namespace domain_decomposition {



//**********************************************************************************************************************
/*!
*   \brief Interface/base class for block IDs (unique identifiers for blocks)
*
*   Just as 'IBlock' serves as a base class for all blocks and 'BlockStorage' serves as a base class for all data
*   structures that manage these blocks, 'IBlockID' serves as a base class for all block IDs. The single purpose of
*   block IDs is to uniquely identify blocks within their corresponding block structure (see class 'BlockStorage').
*/
//**********************************************************************************************************************

class IBlockID {

public:
   /// ID type which can be used as a key in a map
   using IDType = uint64_t;

   virtual ~IBlockID() = default;

   virtual bool operator< ( const IBlockID& rhs ) const = 0;
   virtual bool operator==( const IBlockID& rhs ) const = 0;
   virtual bool operator!=( const IBlockID& rhs ) const = 0;

   /// returns a somehow simplified ID which can be used in maps
   virtual IDType getID() const = 0;

   virtual std::ostream& toStream( std::ostream& os ) const = 0;
};



inline std::ostream& operator<<( std::ostream& os, const IBlockID& id ) {

   return id.toStream( os );
}



} // namespace domain_decomposition

using domain_decomposition::IBlockID;

} // namespace walberla
