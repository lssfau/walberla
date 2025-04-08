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
//! \file Parser.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Parser.h"

#ifndef WALBERLA_DEACTIVATE_MATH_PARSER

#include "core/debug/Debug.h"
#include "core/Abort.h"

#if defined WALBERLA_CXX_COMPILER_IS_MSVC
#   pragma warning( push, 1 )
#   pragma warning( disable : 4706 )
#elif ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wpragmas"
#   pragma GCC diagnostic ignored "-Wsign-conversion"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#elif defined WALBERLA_CXX_COMPILER_IS_INTEL
#   pragma warning push
#   pragma warning( disable : 187  )
#   pragma warning( disable : 1599 )
#endif


#include "extern/exprtk.h"
#include <algorithm>


#if defined WALBERLA_CXX_COMPILER_IS_MSVC
#   pragma warning( pop )
#elif ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#elif defined WALBERLA_CXX_COMPILER_IS_INTEL
//#   pragma warning pop // disabled because this leads to spilled warnings
#endif


namespace walberla {
namespace math {

FunctionParser::FunctionParser() = default;

FunctionParser::~FunctionParser() 
{
   delete expression_;
   delete symbolTable_;
}

void FunctionParser::parse( const std::string & eq )
{
   delete expression_;
   delete symbolTable_;

   expression_  = new exprtk::expression<double>();
   symbolTable_ = new exprtk::symbol_table<double>();
   expression_->register_symbol_table( *symbolTable_ );

   exprtk::parser<double> parser;
   parser.enable_unknown_symbol_resolver();

   if( !parser.compile( eq, *expression_ ) )
   {
      std::ostringstream oss;
      oss << "Error while parsing expression \"" << eq << "\". Error: " << parser.error() << '\n';

      for( std::size_t i = 0; i < parser.error_count(); ++i )
      {
         exprtk::parser_error::type error = parser.get_error( i );


         oss << "Error: " << i << " Position: " << error.token.position << " Type: [" << exprtk::parser_error::to_str( error.mode ) << "] Message: " << error.diagnostic << " Expression: " << eq << '\n';
      }

      WALBERLA_ABORT( oss.str() );
   }

   // check if this expression evaluates to constant zero
   std::vector< std::string > variables;
   symbolTable_->get_variable_list( variables );
   isConstant_ = (variables.empty());
   if (isConstant_)
   {
      const double value = evaluate(std::map<std::string, double>());
      isZero_ = floatIsEqual(value, 0);
   }
   else
   {
      isZero_ = false;
   }
}

double FunctionParser::evaluate( const std::map<std::string,double> & symbolTable) const
{
   if( expression_ == nullptr )
   {
      WALBERLA_ASSERT_NULLPTR( symbolTable_ );
      WALBERLA_ABORT( "Error: You are trying to evaluate an expression which you never have parsed!" );
   }

   std::vector< std::string > variables;
   symbolTable_->get_variable_list( variables );

   for(auto & variable : variables)
   {
      auto symbolEntryIt = symbolTable.find( variable );

      if( symbolEntryIt == symbolTable.end() )
         WALBERLA_ABORT( "Error evaluating expression. Variable \"" << variable << "\" not specified in symbol table!" );

      symbolTable_->variable_ref( variable ) = symbolEntryIt->second;
   }

   return expression_->value();
}

bool FunctionParser::symbolExists(const std::string & symbol) const
{
   if( expression_ == nullptr )
   {
      WALBERLA_ASSERT_NULLPTR( symbolTable_ );
      WALBERLA_ABORT( "Error: You are trying to evaluate an expression which you never have parsed!" );
   }

   std::vector< std::string > variables;
   symbolTable_->get_variable_list( variables );

   return std::find(variables.begin(), variables.end(), symbol) != variables.end();

}



} // namespace math
} // namespace walberla



#else // WALBERLA_DEACTIVATE_MATH_PARSER

#include "core/Abort.h"

namespace walberla {
namespace math {

FunctionParser::FunctionParser()
   : alwaysFalse_( false )
{}

FunctionParser::~FunctionParser() 
{}

void FunctionParser::parse( const std::string & )
{
   WALBERLA_ABORT( "math::FunctionParser was deactivated at compile time by CMake and cannot be used!" );
}

double FunctionParser::evaluate( const std::map<std::string,double> & ) const
{
   WALBERLA_ABORT( "math::FunctionParser was deactivated at compile time by CMake and cannot be used!" );
   return 0.0;
}

const bool & FunctionParser::isConstant() const
{
   WALBERLA_ABORT( "math::FunctionParser was deactivated at compile time by CMake and cannot be used!" );
   return alwaysFalse_;
}

const bool & FunctionParser::isZero() const
{
   WALBERLA_ABORT( "math::FunctionParser was deactivated at compile time by CMake and cannot be used!" );
   return alwaysFalse_;
}

bool FunctionParser::symbolExists(const std::string &) const
{
   WALBERLA_ABORT( "math::FunctionParser was deactivated at compile time by CMake and cannot be used!" );
   return alwaysFalse_;
}

} // namespace math
} // namespace walberla

#endif
