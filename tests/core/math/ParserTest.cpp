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
//! \file ParserTest.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "core/debug/TestSubsystem.h"
#include "core/math/Parser.h"

#include <iostream>
#include <map>
#include <string>


using namespace std;

int main(int argc, char**argv)
{
   using walberla::math::FunctionParser;

   walberla::debug::enterTestMode();

   FunctionParser p;
   map<string,double> st;
   st["t"] = 42.0;
   st["someVar"] = 2;

   p.parse("t");
   WALBERLA_CHECK_FLOAT_EQUAL(42.0, p.evaluate(st) );

   p.parse("t/2 * someVar");
   WALBERLA_CHECK_FLOAT_EQUAL(42.0, p.evaluate(st) );

   string expr;
   for(int i=1; i<argc; ++i)
      expr += string(argv[i]);

   if(!expr.empty()) {
      p.parse(expr);
      cout << p.evaluate(st) << endl;
   }

   return 0;
}
