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
//! \file Operator.cpp
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "Operator.h"


namespace walberla {
namespace math {

   // no valid operators
   OpNo    OP_NO   ( 'n', "no op",  0u );
   OpNo    OP_EQUAL( '=', "equal",  0u );

   // operators
   OpPlus  OP_PLUS ( '+', "plus",  10u );
   OpMinus OP_MINUS( '-', "minus", 10u );
   OpMult  OP_MULT ( '*', "mult",  30u );
   OpDiv   OP_DIV  ( '/', "div",   30u );
   OpProd  OP_PROD ( '^', "prod",  40u );

   // functions
   OpLog   OP_LOG  ( '$', "log",   50u );
   //OpRoot  OP_ROOT ( '%', "root",  50u );


   int isop( const char c )
   {
      return
         ( OP_PLUS  == c ||
         OP_MINUS == c ||
         OP_MULT  == c ||
         OP_DIV   == c ||
         OP_PROD  == c    ) ? c : 0;
   }

   OpType& getOp ( const char c )
   {
      if (OP_PLUS  == c) return OP_PLUS;
      if (OP_MINUS == c) return OP_MINUS;
      if (OP_MULT  == c) return OP_MULT;
      if (OP_DIV   == c) return OP_DIV;
      if (OP_PROD  == c) return OP_PROD;
      WALBERLA_ABORT( "Found no operator" );
      return OP_NO; // has no effect
   }

   std::ostream& operator<<( std::ostream& os, const OpType & type ){
      if( type == '$' )
         return os << type.name_;
      return os << type.sign_;
   }

} // namespace math
} // namespace walberla
