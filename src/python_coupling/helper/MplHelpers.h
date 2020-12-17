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
//! \file MplHelpers.h
//! \ingroup python_export
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/IBlock.h"


#include <functional>
#include <map>


namespace walberla {
namespace python_coupling {


template <typename T>
struct NonCopyableWrap {};


#define FunctionExporterClass( functionName, functionSignature ) \
   struct Exporter_##functionName \
   { \
      typedef std::function<  functionSignature > FunctionType;\
      Exporter_##functionName( const IBlock * block, ConstBlockDataID id )\
         : block_( block ), blockDataID_( id )\
      {}\
      template< typename FieldType > \
      void operator()( ::walberla::python_coupling::NonCopyableWrap<FieldType> ) \
      { \
         if ( block_->isDataClassOrSubclassOf< FieldType > ( blockDataID_ ) ) \
             result = static_cast<FunctionType>( functionName< FieldType > ); \
      } \
      FunctionType result; \
      const IBlock * block_; \
      const ConstBlockDataID blockDataID_; \
   }



template< typename F >
void for_each_noncopyable_type( const F & )
{}

template< typename Type, typename... Types, typename F >
void for_each_noncopyable_type( const F & f)
{
   f(NonCopyableWrap<Type>());
   for_each_noncopyable_type<Types...>(f);
}



template<typename Exporter, typename... FieldTypes>
class Dispatcher
{
public:
   typedef typename Exporter::FunctionType FunctionType;

   Dispatcher( const IBlock * block )
      : block_( block )
   {}

   FunctionType operator() ( BlockDataID blockDataID )
   {
      if ( map_.find( blockDataID) != map_.end() )
         return map_[ blockDataID ];

      Exporter exporter( block_, blockDataID );
      for_each_noncopyable_type< FieldTypes...>  ( std::ref(exporter) );
      map_[ blockDataID ] = exporter.result;
      return exporter.result;
   }

private:
   const IBlock * block_;
   std::map< BlockDataID, FunctionType > map_;
};



} // namespace python_coupling
} // namespace walberla


