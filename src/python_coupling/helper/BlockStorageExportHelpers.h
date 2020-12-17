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
//! \file ExportHelpers.h
//! \ingroup python_coupling
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "python_coupling/PythonWrapper.h"
#include "core/Set.h"
#include "core/uid/UID.h"

#include "domain_decomposition/StructuredBlockStorage.h"
#include <stdexcept>
#include <vector>
#include <string>




namespace walberla {
namespace python_coupling {


   struct NoSuchBlockData : public std::out_of_range {
      explicit NoSuchBlockData ( ) : std::out_of_range( "No blockdata with the given name found" ) {}
      explicit NoSuchBlockData ( const std::string & w ) : std::out_of_range(w) {}
      static void translate( const NoSuchBlockData & e );
   };
   struct BlockDataNotConvertible : public std::runtime_error {
      explicit BlockDataNotConvertible (  ) : std::runtime_error( "This blockdata is not accessible from Python" ) {}
      explicit BlockDataNotConvertible ( const std::string & w ) : std::runtime_error(w) {}
      static void translate( const BlockDataNotConvertible &  e );
   };

   BlockDataID blockDataIDFromString( IBlock & block, const std::string & stringID );
   BlockDataID blockDataIDFromString( BlockStorage & bs, const std::string & stringID );
   BlockDataID blockDataIDFromString( StructuredBlockStorage & bs, const std::string & stringID );



#ifdef WALBERLA_BUILD_WITH_PYTHON

   template< typename FField>
   typename FField::value_type maskFromFlagList(  const shared_ptr<StructuredBlockStorage> & bs,
                                                  ConstBlockDataID flagFieldID,
                                                  const std::vector< std::string > & flagList )
   {
      if ( bs->begin() == bs->end() )
         return 0;

      IBlock & firstBlock = *(  bs->begin() );
      const FField * flagField = firstBlock.getData< const FField > ( flagFieldID );

      typedef typename FField::flag_t flag_t;
      flag_t mask = 0;
      for( auto it = flagList.begin(); it != flagList.end(); ++it )
      {
         if ( ! flagField->flagExists( *it ) )
            throw python_coupling::BlockDataNotConvertible( "Unknown FlagID" );

         mask = flag_t( mask | flagField->getFlag( *it ) );
      }

      return mask;
   }

#else
   template<typename UID, typename StringContainer>
   Set< UID > uidSetFromStringContainer( const StringContainer &  )
   {
      return Set< UID >();
   }

   template< typename FField>
   typename FField::value_type maskFromFlagList(  const shared_ptr<StructuredBlockStorage> & ,
                                                  ConstBlockDataID ,
                                                  const std::vector< std::string > &  )
   {
      return typename FField::value_type();
   }

#endif


} // namespace python_coupling
} // namespace walberla


