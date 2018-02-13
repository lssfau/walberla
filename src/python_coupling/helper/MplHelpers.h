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
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/IBlock.h"

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/lambda.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/transform.hpp>

#include <boost/bind.hpp>

#include <map>


namespace walberla {
namespace python_coupling {


template <typename V, typename T, typename Result>
struct list_of_pairs
  : boost::mpl::fold<V, Result,
        boost::mpl::push_back<boost::mpl::_1, boost::mpl::pair<T, boost::mpl::_2> > >
{};

template<typename V1, typename V2>
struct combine_vectors
: boost::mpl::fold<
    V1,
    boost::mpl::vector<>,
    boost::mpl::lambda<list_of_pairs<V2,boost::mpl::_2, boost::mpl::_1> >
>::type
{};



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



template< typename Sequence, typename F >
void for_each_noncopyable_type( const F & f)
{
   boost::mpl::for_each< Sequence, NonCopyableWrap< boost::mpl::placeholders::_1> >  ( f );
}



template<typename FieldTypeList, typename Exporter>
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
      for_each_noncopyable_type< FieldTypeList>  ( boost::ref(exporter) );
      map_[ blockDataID ] = exporter.result;
      return exporter.result;
   }

private:
   const IBlock * block_;
   std::map< BlockDataID, FunctionType > map_;
};



} // namespace python_coupling
} // namespace walberla


