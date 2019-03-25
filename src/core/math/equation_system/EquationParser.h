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
//! \file EquationParser.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "EquationSystem.h"
#include "core/logging/Logging.h"

#include <sstream>

namespace walberla {
namespace math {

   class EquationParser
   {
   private:
      EquationSystem& es_;

   public:
      EquationParser( EquationSystem& es) : es_(es) { }

   private:
      EquationParser& operator=( const EquationParser& ) { return *this; }

   private:
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // PARSE NUMBER
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      bool    checkNumber( const std::string& str, size_t& index ) const;
      NodePtr parseNumber( const std::string& str, size_t& index ) const;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // PARSE NT_VARIABLE
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      bool    checkVariable( const std::string& str, size_t& index ) const;
      NodePtr parseVariable( const std::string& str, size_t& index ) const;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // PARSE FUNCTION
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      enum OpFunction{
         OP_FUNC_EXP,
         OP_FUNC_LN,
         OP_FUNC_SQRT
      };

      bool    checkFunction( const std::string& str, size_t& index ) const;
      NodePtr parseFunction( const std::string& str, size_t& index ) const;

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // PARSE EXPRESSION
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      NodePtr parseTerm      ( const std::string& str, size_t& index ) const;
      NodePtr parseExpression( const std::string& str, size_t& index ) const;

   public:
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // PARSE EQUATION
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      EquationPtr parseEquation( const std::string& str, size_t& index );
   };

} // namespace math
} // namespace walberla
