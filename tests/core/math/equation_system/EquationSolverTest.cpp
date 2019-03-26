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
//! \file EquationSolverTest.cpp
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/Environment.h"
#include "core/math/equation_system/EquationParser.h"
#include "core/math/equation_system/EquationSystem.h"

#include <iostream>
#include <string>
#include <vector>


using namespace walberla;
using namespace walberla::math;

/*
int directInput(){
   EquationSystem es;
   EquationParser ep(es);

   std::string str;
   uint_t index = 0;

   std::cout << "\nWrite equations to solve without any blank!\nTo quit enter 'exit'\nTo clear all known variable enter 'clear'"<< std::endl;

   bool run = true;
   do {
      std::cin >> str;
      index = 0;
      if ( strcmp(str.c_str(),"exit") == 0 )
         return 0;
      if ( strcmp(str.c_str(),"clear") == 0 ){
         std::cout << "Clear known variables and equations" << std::endl;
         es.clear();
         continue;
      }
      try {
         EquationPtr eqPtr = ep.parseEquation(str, index);
         es.add(str, eqPtr);
         std::cout << "Convert  '" << str << "'  to  '" << (*eqPtr) << "'" << std::endl;
         if (eqPtr->isComputable()){
            VarPtr varPtr = eqPtr->compute();
            std::cout << "Equation is computable: " << varPtr->name() << "=" << varPtr->value() << std::endl;
         } else if (eqPtr->isEvaluatable()){
            std::cout << "Equation is evaluatable: " << (eqPtr->evaluate() ? "true" : "false") << std::endl;
         } else {
            std::cout << "Equation is neither computable nor evaluatable!" << std::endl;
         }
      } catch (std::runtime_error re) {
         std::cerr << re.what() << std::endl;
      }
   } while (run);
   return 0;
}
 */


int equationInput(){
   std::vector<std::string> eqStringList;

   //// Parameters
   eqStringList.emplace_back("dt = 2e-7");
   eqStringList.emplace_back("dx = 5e-6");
   eqStringList.emplace_back("eta = 0.0001");
   eqStringList.emplace_back("omega = 1.95");
   eqStringList.emplace_back("rho = 1000");

   //// LBM Equations
   eqStringList.emplace_back("'rho_L' = 1.0");
   eqStringList.emplace_back("'dt_L'  = 1.0");
   eqStringList.emplace_back("'dx_L'  = 1.0");
   eqStringList.emplace_back("'c'     = 'dx_L' / 'dt_L'");
   eqStringList.emplace_back("'nu'    = 'eta' / 'rho'");
   eqStringList.emplace_back("'nu_L'  = 'eta_L' / 'rho_L'");
   eqStringList.emplace_back("'dt'    = ( 0.1 * 'dx' ) / 'maxOcurringPhysVel'");
   eqStringList.emplace_back("'cs'    = ( 1.0 / ( 3.0 ^ 0.5 ) ) * 'c'");
   eqStringList.emplace_back("'omega' = 1.0 / 'tau'");
   eqStringList.emplace_back("'nu_L'  = ( 'cs' ^ 2.0 ) * ( 'tau' - ( 0.5 * 'dt_L' ) )");
   /*
   // Unsolvable:
   // Parameters
   eqStringList.push_back( "nu = 3.50E-006");
   eqStringList.push_back( "omega = 1.99");
   eqStringList.push_back( "rho = 1000");
   eqStringList.push_back( "maxOcurringPhysVel = 0.10");


   // LBM Equations
   eqStringList.push_back( "'rho_L' = 1.0");
   eqStringList.push_back( "'dt_L'  = 1.0");
   eqStringList.push_back( "'dx_L'  = 1.0");
   eqStringList.push_back( "'c'     = 'dx_L' / 'dt_L'");
   eqStringList.push_back( "'nu'    = 'eta' / 'rho'");
   eqStringList.push_back( "'nu_L'  = 'eta_L' / 'rho_L'");
   eqStringList.push_back( "'dt'    = ( 0.1 * 'dx' ) / 'maxOcurringPhysVel'");
   eqStringList.push_back( "'cs'    = ( 1.0 / ( 3.0 ^ 0.5 ) ) * 'c'");
   eqStringList.push_back( "'omega' = 1.0 / 'tau'");
   eqStringList.push_back( "'nu_L'  = ( 'cs' ^ 2.0 ) * ( 'tau' - ( 0.5 * 'dt_L' ) )");
   eqStringList.push_back( "'nu_L'  = (nu * dt) / dx^2");
    */

   EquationSystem es;
   EquationParser ep(es);
   size_t index = 0;
   size_t number = 0;

   for (size_t i=0; i<eqStringList.size(); ++i){
      index = 0;
      es.add( std::to_string(++number), ep.parseEquation( eqStringList[i], index ) );
   }

   WALBERLA_CHECK( es.solve() );
   //es.match();
   WALBERLA_LOG_RESULT( es );
   return 0;
}

int unitTest(double v)
{
   std::string s = std::to_string( v );

   std::vector<std::string> eqStringList;
   eqStringList.push_back(       "a = " + s );
   eqStringList.push_back(   "a + 3 =   " + s + " + 3" );
   eqStringList.push_back(   "3 + a =   3 + " + s + "" );
   eqStringList.push_back(   "a - 3 =   " + s + " - 3" );
   eqStringList.push_back(   "3 - a =   3 - " + s + "" );
   eqStringList.push_back(   "a * 3 =   " + s + " * 3" );
   eqStringList.push_back(   "3 * a =   3 * " + s + "" );
   eqStringList.push_back(   "a / 3 =   " + s + " / 3" );
   eqStringList.push_back(   "3 / a =   3 / " + s + "" );
   eqStringList.push_back(   "a ^ 3 =   " + s + " ^ 3" );
   eqStringList.push_back(   "3 ^ a =   3 ^ " + s + "" );
   eqStringList.push_back( "sqrt(a) = sqrt(" + s + ")" );
   eqStringList.push_back(  "exp(a) =  exp(" + s + ")" );
   eqStringList.push_back(   "ln(a) =   ln(" + s + ")" );

   EquationSystem es;
   EquationParser ep(es);
   uint_t index = 0;

   for (uint_t i=0; i<eqStringList.size(); ++i){
      index = 0;
      es.add( eqStringList[i], ep.parseEquation( eqStringList[i], index ) );
   }
   WALBERLA_CHECK( es.solve() );

   return 0;
}

int unitTests(unsigned int count){
   srand( static_cast<unsigned int>(time(nullptr)) );

   double values[] = {0.0, 1.0, 1e-15, 1e+15};
   unsigned int size = 4;

   int test = 0;
   for (unsigned int i=0; i<size && test == 0; ++i){
      test = unitTest( values[i] );
   }

   for (unsigned int i=0; i<count && test == 0; ++i){
      double value = double(rand()) / RAND_MAX;
      int exp = rand() / ( RAND_MAX / 30 ) - 14;
      test = unitTest( pow( value, exp ) );
   }

   return test;
}

int main( int argc, char** argv )
{
   debug::enterTestMode();

   mpi::Environment mpiEnv( argc, argv );

   int value;
   //value = unitTests(100);
   value = equationInput();
   //value = directInput();
   return value;
}
