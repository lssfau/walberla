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
//! \file RBGSIteration.h
//! \ingroup pde
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/Set.h"
#include "core/uid/SUID.h"
#include "domain_decomposition/BlockStorage.h"

#include <functional>



namespace walberla {
namespace pde {



class RBGSIteration
{
public:

   RBGSIteration( BlockStorage & blocks, const uint_t iterations,
                  const std::function< void () > & communication,
                  const std::function< void ( IBlock * ) > & redUpdate,
                  const std::function< void ( IBlock * ) > & blackUpdate,
                  const std::function< real_t () > & residualNorm,
                  const real_t residualNormThreshold = real_t(0), const uint_t residualCheckFrequency = uint_t(1),
                  const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                  const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), iterations_( iterations ),
      residualNormThreshold_( residualNormThreshold ), residualCheckFrequency_( residualCheckFrequency ),
      iterationsPerformed_( uint_t(0) ), thresholdReached_( false ),
      communication_( communication ), redUpdate_( redUpdate ), blackUpdate_( blackUpdate ), residualNorm_( residualNorm ),
      requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
   {}

   void addBoundaryHandling( const std::function< void () > & boundary ) { boundary_ = boundary; }
      
   void operator()();
   
   uint_t iterationsPerformed() const { return iterationsPerformed_; }
   bool   thresholdReached() const { return thresholdReached_; }

protected:

   BlockStorage & blocks_;

   uint_t iterations_;
   real_t residualNormThreshold_;
   uint_t residualCheckFrequency_;
   
   uint_t iterationsPerformed_;
   bool thresholdReached_;

   std::function< void () >            boundary_;
   std::function< void () >            communication_;
   std::function< void ( IBlock * ) >  redUpdate_;
   std::function< void ( IBlock * ) >  blackUpdate_;
   std::function< real_t () >          residualNorm_;
   
   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;
};



} // namespace pde
} // namespace walberla
