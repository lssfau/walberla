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
//! \file SharedFunctor.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "DataTypes.h"


namespace walberla {


template< typename F >
class SharedFunctor
{
public:

   SharedFunctor( const shared_ptr<F> & functorPtr ) : functorPtr_( functorPtr ) { }

   void operator()(){ (*functorPtr_)(); }

private:

   shared_ptr<F> functorPtr_;
};


template< typename F >
SharedFunctor<F> makeSharedFunctor( const shared_ptr<F> & functorPtr ) { return SharedFunctor<F>( functorPtr ); }


} // namespace walberla
