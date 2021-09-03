//==============================================================================================================================================================
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
//! \file ExecutionTree.impl.h
//! \ingroup executiontree
//! \author Martin Bauer <martin.bauer@fau.de>
//
//==============================================================================================================================================================

#pragma once

#include "ExecutionTree.h"
#include <utility>

namespace walberla {
namespace executiontree {

namespace internal {

// Helper to handle functors and shared_ptr's to functors the same way
template<typename T>
struct Caller
{
   template<typename ... Args>
   static void call( T &t, Args&&... args )
   {
      t(std::forward<Args>(args)...);
   }
};

template<typename T>
struct Caller< shared_ptr < T > >
{
   template<typename ... Args>
   static void call( shared_ptr <T> &t, Args&&... args )
   {
      ( *t )(std::forward<Args>(args)...);
   }
};


} // namespace internal


template<typename FunctorType>
IFunctionNodePtr functor( FunctorType t, const std::string &name, const TimingTreePtr &timingTree )
{
   return make_shared< Functor< FunctorType > >( t, name, timingTree );
}

inline shared_ptr <Sequence> sequence( std::initializer_list< IFunctionNodePtr > initializerList, const std::string &name,
                                       const TimingTreePtr &timingTree )
{
   return make_shared< Sequence >( initializerList, name, timingTree, false );
}

inline shared_ptr <Sequence> parallelSequence( std::initializer_list< IFunctionNodePtr > initializerList, const std::string &name,
                                               const TimingTreePtr &timingTree )
{
   return make_shared< Sequence >( initializerList, name, timingTree, true );
}


inline shared_ptr< EveryNth > everyNth( const IFunctionNodePtr &node, uint_t interval, bool onFirst, uint_t startValue )
{
   return make_shared< EveryNth >( node, interval, onFirst, startValue );
}


inline shared_ptr< Loop > loop( const IFunctionNodePtr &body, uint_t iterations, bool logTimeStep )
{
   return make_shared< Loop >( body, iterations, logTimeStep );
}


template<typename FunctorType>
Functor< FunctorType >::Functor( const FunctorType &functor, const std::string &name, const TimingTreePtr &timingTree )
        :functor_( functor ), name_( name ), timingTree_( timingTree ) {}


template<typename FunctorType>
void Functor< FunctorType >::operator()()
{
   if ( timingTree_ )
   {
      timingTree_->start( name_ );
      internal::Caller<FunctorType>::call(functor_);
      timingTree_->stop( name_ );
   }
   else
      internal::Caller<FunctorType>::call(functor_);
}



} // namespace executiontree
} // namespace walberla