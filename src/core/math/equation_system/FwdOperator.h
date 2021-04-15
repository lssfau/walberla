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
//! \file FwdOperator.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include <map>


namespace walberla {
namespace math {
   class OpType;
   class OpNo;
   class OpPlus;
   class OpMinus;
   class OpMult;
   class OpDiv;
   class OpProd;
   class OpRoot;
   class OpLog;

   // no valid operators
   extern OpNo    OP_NO;
   extern OpNo    OP_EQUAL;

   // operators
   extern OpPlus  OP_PLUS;
   extern OpMinus OP_MINUS;
   extern OpMult  OP_MULT;
   extern OpDiv   OP_DIV;
   extern OpProd  OP_PROD;

   // functions
   extern OpLog   OP_LOG;
   extern OpRoot  OP_ROOT;

} // namespace math
} // namespace walberla
