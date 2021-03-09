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
//! \file ExprSystemInitFunction.h
//! \ingroup lbm
//! \author Tobias Schruff <tobias.schruff@gmail.com>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"

#include <vector>
#include <string>
#include <map>

namespace exprtk {
   template <typename T> class expression;
   template <typename T> class symbol_table;
}

namespace walberla {
namespace lbm {
namespace initializer {


class ExprSystemInitFunction
{
public:

   using SymbolTable = exprtk::symbol_table<real_t>;
   using Expression = exprtk::expression<real_t>;
   using ExprSystem = std::map<std::string, Expression>;

   ExprSystemInitFunction( const shared_ptr<StructuredBlockForest> & blocks );
   ~ExprSystemInitFunction();

   bool parse( const Config::BlockHandle & config );

   real_t          getDensity();
   Vector3<real_t> getVelocity();

   void setGlobalCell( const Cell & cell );

   std::vector<real_t> operator()( const Cell & globalCell );

private:

   ExprSystemInitFunction( const ExprSystemInitFunction & other );
   ExprSystemInitFunction & operator=( const ExprSystemInitFunction & other );

   const shared_ptr<StructuredBlockForest> blocks_;

   Vector3<real_t> cell_;
   SymbolTable     * symbolTable_;
   ExprSystem      * expr_;
};


} // namespace initializer
} // namespace lbm
} // namespace walberla