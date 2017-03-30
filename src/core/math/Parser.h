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
//! \file Parser.h
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/CMakeDefs.h"

#include <map>
#include <stdexcept>
#include <string>

namespace exprtk {
   template <typename T> class expression;
   template <typename T> class symbol_table;
}


namespace walberla {
namespace math {

#ifndef WALBERLA_DEACTIVATE_MATH_PARSER

/**
 * Function parser, for parsing and evaluating functions
 *
 * Can parse mathematical expressions containing the standard arithmetic
 * operations(+,-,/,*) , brackets, functions (cos,sin,exp,sqrt) and variables
 *
 * After parsing an expression, values can be bound to variables using a
 * symbol table.
 *
 * Example:
 * \code
 *     FunctionParser p;
 *     p.parse("(t*t +5) + exp(Pi/2) * sin(cos(1)) ");
 *     std::map<string,double> symbolTable;
 *     symbolTable["Pi"]=3.141;
 *     symbolTable["t"] = 24;
 *     p.evaluate(symbolTable);
 * \endcode
 */
class FunctionParser
{
public:
   FunctionParser();
   ~FunctionParser();
   
   void   parse   ( const std::string & equation );
   double evaluate( const std::map<std::string,double> & symbolTable ) const;
   inline const bool & isConstant() const { return isConstant_; }
   inline const bool & isZero()     const { return isZero_    ; }
   bool symbolExists(const std::string & symbol) const;

protected:
   FunctionParser( const FunctionParser & other );
   FunctionParser & operator=( const FunctionParser & other );

   exprtk::expression<double>   * expression_;
   exprtk::symbol_table<double> * symbolTable_;
   bool isConstant_;
   bool isZero_;
};

#else // WALBERLA_DEACTIVATE_MATH_PARSER

class FunctionParser
{
public:
   FunctionParser();
   ~FunctionParser();
   
   void   parse   ( const std::string & );
   double evaluate( const std::map<std::string,double> & ) const;
   const bool & isConstant() const;
   const bool & isZero()     const;

protected:
   FunctionParser( const FunctionParser & );
   FunctionParser & operator=( const FunctionParser & );
   
   bool alwaysFalse_;
};

#endif


} // namespace math
} // namespace walberla
