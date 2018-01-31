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
//! \file RefinementFunctorWrapper.h
//! \ingroup lbm
//! \author Christoph Rettinger <christoph.rettinger@fau.de>
//
//======================================================================================================================

#pragma once

#include "domain_decomposition/BlockStorage.h"

#include <functional>

namespace walberla {
namespace lbm {
namespace refinement {

class FunctorWrapper {
public:
   FunctorWrapper( const std::function<void()> & fct)
         : fct_(fct) {
   }

   void operator()(const uint_t /*level*/, const uint_t /*executionCounter*/) {
      fct_();
   }

private:
   std::function<void(void)> fct_;
};

class SweepAsFunctorWrapper {
public:
   SweepAsFunctorWrapper( const std::function<void(IBlock * )> & fct,
                          const shared_ptr <StructuredBlockStorage> & blockStorage )
         : fct_(fct), blockStorage_(blockStorage) {
   }

   void operator()(const uint_t level, const uint_t /*executionCounter*/) {
      for (auto blockIt = blockStorage_->begin(); blockIt != blockStorage_->end(); ++blockIt) {
         if (blockStorage_->getLevel(*blockIt) != level) continue;
         fct_(blockIt.get());
      }
   }

private:
   std::function<void(IBlock * )> fct_;
   shared_ptr <StructuredBlockStorage> blockStorage_;
};

} // namespace refinement
} // namespace lbm
} // namespace walberla
