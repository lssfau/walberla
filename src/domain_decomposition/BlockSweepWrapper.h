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
//! \file BlockSweepWrapper.h
//! \ingroup domain_decomposition
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "StructuredBlockStorage.h"
#include "core/debug/CheckFunctions.h"


namespace walberla {
namespace domain_decomposition {



/// Takes an existing block 'sweep' functor and wrapps it into a void-void functor
template< typename BlockSweep_T >
class BlockSweepWrapper
{
public:

   BlockSweepWrapper( const weak_ptr< StructuredBlockStorage > & blocks, const BlockSweep_T & sweep,
                      const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                      const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), sweep_( sweep ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {}

   void operator()()
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'BlockSweepWrapper' for a block storage object that doesn't exist anymore" );

      for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
         sweep_( block.get() );
   }

private:

   weak_ptr< StructuredBlockStorage > blocks_;
   BlockSweep_T sweep_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



} // namespace domain_decomposition

using domain_decomposition::BlockSweepWrapper;

} // namespace walberla
