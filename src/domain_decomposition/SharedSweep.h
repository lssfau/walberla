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
//! \file SharedSweep.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"


namespace walberla {
namespace domain_decomposition {

class IBlock;

namespace internal {

template< typename T >
class SharedSweep
{
public:

   SharedSweep( const shared_ptr<T> & sweepPtr ) : sweepPtr_( sweepPtr ) { }

   void operator()( IBlock * const block ){ (*sweepPtr_)( block ); }

private:

   shared_ptr<T> sweepPtr_;
};

} // namespace internal

template< typename T >
internal::SharedSweep<T> makeSharedSweep( const shared_ptr<T> & sweepPtr ) { return internal::SharedSweep<T>( sweepPtr ); }

} // namespace domain_decomposition

using domain_decomposition::makeSharedSweep;

} // namespace walberla
