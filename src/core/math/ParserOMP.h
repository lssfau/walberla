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
//! \file ParserOMP.h
//! \ingroup core
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#pragma once

#include "Parser.h"

#include <memory>

namespace walberla {
namespace math {

class FunctionParserOMP
{
public:
   FunctionParserOMP();
   void   parse   (const std::string & equation);
   double evaluate(const std::map<std::string,double> & symbolTable) const;
   inline const bool & isConstant() const { return parser_[0].isConstant(); }
   inline const bool & isZero()     const { return parser_[0].isZero()    ; }
   inline bool symbolExists(const std::string & symbol) const { return parser_[0].symbolExists(symbol); }

private:
   std::unique_ptr< FunctionParser[] > parser_;
#ifndef NDEBUG
   int num_parsers_;
#endif
};


} // namespace math
} // namespace walberla
