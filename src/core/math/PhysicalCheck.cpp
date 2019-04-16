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
//! \file PhysicalCheck.cpp
//! \ingroup core
//! \author David Staubach <david.staubach@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"
#ifdef WALBERLA_BUILD_WITH_BOOST

#include "PhysicalCheck.h"
#include "core/Abort.h"
#include "core/logging/Logging.h"
#include "core/math/Parser.h"
#include "equation_system/EquationParser.h"


namespace walberla {
namespace math {

   PhysicalCheck::PhysicalCheck() :
         solved_( false )
   {}

   PhysicalCheck::PhysicalCheck( const std::vector< std::string >& equations, const std::map< std::string, std::string >& unitParameterRelations, const std::vector<std::string>& constraints ) :
         solved_( false )
   {
      addEquations( equations );
      addUnitParameterRelations( unitParameterRelations );
      addConstraints( constraints );
   }

   PhysicalCheck::PhysicalCheck( const Config::BlockHandle& configBlock ) :
      solved_( false )
   {
      addBlock( configBlock );
   }

   void PhysicalCheck::addBlock( const Config::BlockHandle& configBlock )
   {
      auto eqBlock = configBlock.getBlock( "Equations" );

      if( eqBlock )
      {
         std::vector< std::string > equations;
         for( auto i=eqBlock.begin(); i!=eqBlock.end(); ++i )
            equations.push_back(i->second);

         addEquations( equations );
      }

      auto unitBlock = configBlock.getBlock( "Units" );

      if ( unitBlock )
      {
         std::map< std::string, std::string > unitParameterRelations;
         for( auto i=unitBlock.begin(); i!=unitBlock.end(); ++i )
            unitParameterRelations.insert(
                     std::pair<std::string,std::string>(i->first,i->second));


         addUnitParameterRelations( unitParameterRelations );
      }

      auto constraintsBlock = configBlock.getBlock( "Constraints" );

      if( constraintsBlock )
      {
         std::vector<std::string> constraints;
         for( auto i=constraintsBlock.begin(); i!=constraintsBlock.end(); ++i )
            constraints.push_back(i->second);

         addConstraints( constraints );
      }
   }

   void PhysicalCheck::addEquations( const std::vector< std::string >& equations )
   {
      EquationParser ep(es_);
      size_t index = 0;

      solved_ = false;

      for (size_t i=0; i<equations.size(); ++i){
         index = 0;
         es_.add( std::to_string(es_.getNumberOfEquations()+1), ep.parseEquation( equations[i], index ) );
      }
   }

   void PhysicalCheck::addEquation( const std::string& equation )
   {
      EquationParser ep(es_);
      size_t index = 0;

      solved_ = false;

      es_.add( std::to_string(es_.getNumberOfEquations()+1), ep.parseEquation( equation, index ) );
   }

   void PhysicalCheck::addUnitParameterRelations( const std::map< std::string, std::string >& unitParameterRelations )
   {
      // parse units and store unit as int
      for( auto i=unitParameterRelations.begin(); i!=unitParameterRelations.end(); ++i )
      {
//         unitParameterRelations_[i->first]["m"] = 0;
//         unitParameterRelations_[i->first]["s"] = 0;
//         unitParameterRelations_[i->first]["kg"] = 0;
//         unitParameterRelations_[i->first]["A"] = 0;

         std::string unitString = i->second;
         int factor = 1;

         for( size_t j=0; j<unitString.size(); ++j )
         {
            if( unitString[j] == 'm' || unitString[j] == 's' || unitString[j] == 'A' )
            {
               size_t index = j;
               ++j;

               // make this work for exponents larger than 9
               if( j < unitString.size() && unitString[j] == '^' )
               {
                  ++j;
                  int expo = factor * std::atoi( &unitString[j] );
                  if( !setVarUnit( i->first, unitString.substr(index,1), expo ) )
                     WALBERLA_ABORT( "Error in PhysicalCheck::addUnitParameterRelations(). Non-unique description for unit '" << unitString[index] << "' for parameter '" << i->first << "'." );
               }
               else
               {
                  --j;
                  if( !setVarUnit( i->first, unitString.substr(index,1), factor ) )
                     WALBERLA_ABORT( "Error in PhysicalCheck::addUnitParameterRelations(). Non-unique description for unit '" << unitString[index] << "' for parameter '" << i->first << "'." );
               }
            }
            else if( unitString[j] == 'k' )
            {
               size_t index = j;
               ++j;++j;

               if( j < unitString.size() && unitString[j] == '^' )
               {
                  ++j;
                  int expo = factor * std::atoi( &unitString[j] );
                  if( !setVarUnit( i->first, unitString.substr(index,2), expo ) )
                     WALBERLA_ABORT( "Error in PhysicalCheck::addUnitParameterRelations(). Non-unique description for unit 'kg' for parameter '" << i->first << "'.");
               }
               else
               {
                  --j;
                  if( !setVarUnit( i->first, unitString.substr(index,2), factor ) )
                     WALBERLA_ABORT( "Error in PhysicalCheck::addUnitParameterRelations(). Non-unique description for unit 'kg' for parameter '" << i->first << "'.");
               }
            }
            else if( unitString[j] == '/')
               factor = -1;
            else
               continue; // necessary to allow for units of the form: 1/s
         }

         // add equation for the calculation between lattice and physical units
         std::stringstream ss;
         ss << "'" << i->first << "_L' = '" << i->first << "'"<< getParametrizationTerm(i->first);
         addEquation( ss.str() );
      }
   }

   std::string PhysicalCheck::getParametrizationTerm( const std::string& varName )
   {
      std::map<std::string,int> parametrizationTerm;

      parametrizationTerm["dx"]  = - unitParameterRelations_[varName]["m"] - 3*unitParameterRelations_[varName]["kg"] - 5*unitParameterRelations_[varName]["A"];
      parametrizationTerm["dt"]  = - unitParameterRelations_[varName]["s"] + 3*unitParameterRelations_[varName]["A"];
      parametrizationTerm["rho"] = - unitParameterRelations_[varName]["kg"] - unitParameterRelations_[varName]["A"];

      std::stringstream num, denom;
      for( auto i=parametrizationTerm.begin(); i!=parametrizationTerm.end(); ++i )
      {
         if( i->second == 0 )
            continue;

         if( i->second < 0 )
         {
            if( denom.str().size() > 0 )
               denom << " * ";
            denom << i->first;
         }
         else
         {
            num << " * " << i->first;
         }

         if( i->second < -1 )
            denom << " ^ " << std::abs( i->second );
         else if( i->second > 1 )
            num << " ^ " << i->second;
      }

      if( num.str().size() == 0 && denom.str().size() == 0 )
         return std::string();

      if( denom.str().size() == 0 )
      {
         num << " *";
         return num.str();
      }

      if( num.str().size() == 0 )
         num << " * 1";

      num << " / ( ";
      num << denom.str() << " )";
      return num.str();
   }

   void PhysicalCheck::addConstraints( const std::vector<std::string>& constraints )
   {
      for( auto i=constraints.begin(); i!=constraints.end(); ++i )
         constraints_.push_back( *i );
   }

   void PhysicalCheck::completeConfig( const shared_ptr<Config>& config )
   {
      auto globalBlock = config->getWritableGlobalBlock();

      std::map<std::string,double> symbolTable;
      getVarMap(symbolTable);

      completeBlock( globalBlock, symbolTable );
   }

   void PhysicalCheck::completeBlock( Config::Block& configBlock, const std::map<std::string,double>& symbolTable )
   {
      if( configBlock.getKey() == "Physical_Check" )
         return;

      // traverse all parameters in the current block
      for( auto param=configBlock.begin(); param!=configBlock.end(); ++param )
      {
         // check for "'" in the string and erase if present
         std::string expression(param->second);
         if( expression[0] == '\'' )
         {
            expression.erase( std::remove(expression.begin(), expression.end(), '\''), expression.end() );

            FunctionParser funcParser;
            try
            {
               // hand the expression over to the FunctionParser
               funcParser.parse(expression);
            }
            catch( std::exception& )
            {
               WALBERLA_ABORT( "BadSyntaxException when completing Config-File. Block: " << configBlock.getKey()<< ", Parameter: "<< param->first );
            }

            double result=0;
            try
            {
               result = funcParser.evaluate(symbolTable);
            }
            catch( std::exception& )
            {
               WALBERLA_ABORT( "UnknownSymbolException when completing Config-File. Block: " << configBlock.getKey() );
            }

            // set the current parameter to the evaluated result in the current block
            std::ostringstream os;
            os << std::setprecision(16) << result;
            configBlock.setParameter( param->first, os.str() );

            // check for the equality of the parameter and the result
            double val = configBlock.getParameter<double>( param->first );
            if( !floatIsEqual( val, result ) )
               WALBERLA_ABORT( "Error in PhysicalCheck::completeBlock(). Failure when trying to complete Block: " << configBlock.getKey() );
         }
      }

      // recursive call of the inner blocks within the current block
      std::vector<Config::Block*> innerBlocks;
      configBlock.getWritableBlocks( innerBlocks );

      for( auto innerBlock=innerBlocks.begin();innerBlock!=innerBlocks.end();++innerBlock )
         completeBlock( **innerBlock, symbolTable );
   }

   bool PhysicalCheck::isDefined( const std::string& varName )
   {
      if( !solved_ )
         solve();

      return es_.isVarDefined( varName );
   }

   double PhysicalCheck::getVarValue( const std::string& varName )
   {
      if( !isDefined(varName) )
      {
         WALBERLA_ABORT( "Error in PhysicalCheck::getVarValue(). Variable not found: " << varName );
         return 0;
      }

      return es_.getVarValue(varName);
   }

   std::string PhysicalCheck::getVarUnit( const std::string& varName )
   {
      if( !isDefined(varName) )
      {
         WALBERLA_ABORT( "Error in PhysicalCheck::getVarUnit(). Variable not found: " << varName );
         return nullptr;
      }

      std::stringstream num, denom;
      for( auto i=unitParameterRelations_[varName].begin(); i!=unitParameterRelations_[varName].end(); ++i )
      {
         if( i->second == 0 )
            continue;

         if( i->second < 0 )
            denom << i->first;
         else
            num << i->first;

         if( i->second < -1 )
            denom << '^' << std::abs( i->second );
         else if( i->second > 1 )
            num << '^' << i->second;
      }

      if( num.str().size() == 0 && denom.str().size() == 0 )
         return std::string();

      if( denom.str().size() == 0 )
         return num.str();

      if( num.str().size() == 0 )
         num << 1;

      num << '/';
      num << denom.str();
      return num.str();
   }

   bool PhysicalCheck::setVarUnit( const std::string& varName, const std::string& unit, const int expo )
   {
      if( unitParameterRelations_[varName][unit] != 0 )
         return false;

      unitParameterRelations_[varName][unit] = expo;
      return true;
   }

   void PhysicalCheck::getVarMap( std::map<std::string,double>& varMap )
   {
      if( !solved_ )
         solve();

      es_.getVarMap( varMap );
   }

   void PhysicalCheck::solve()
   {
      if( !es_.solve() )
         WALBERLA_ABORT( "Error in PhysicalCheck::solve(). System of equations is not solvable." );

      solved_ = true;
      WALBERLA_LOG_INFO( "System of equations has been solved successfully.\n" << es_.writeVariables() );

      if( checkConstraints() )
         WALBERLA_LOG_INFO( "All constraints are fulfilled." );
   }

   bool PhysicalCheck::checkConstraints()
   {
      enum Operator{ LOWER, LOWER_EQUAL, GREATER, GREATER_EQUAL, UNEQUAL, INVALID };

      std::map<std::string,double> symbolTable;
      getVarMap( symbolTable );

      for( auto i=constraints_.begin(); i!=constraints_.end(); ++i )
      {
         std::string constraintString( *i );

         Operator op = INVALID;
         uint_t lenLHS   = 0;
         uint_t startRHS = 0;

         // parse constraint and search for operator <,>,<=,>=,!=
         for( uint_t j=0; j<constraintString.size(); ++j )
         {
            switch( constraintString[j] )
            {
               case '<':
                  lenLHS = j;
                  if( constraintString[j+1] == '=' )
                  {
                     op = LOWER_EQUAL;
                     startRHS = j+2;
                  }
                  else
                  {
                     op = LOWER;
                     startRHS = j+1;
                  }
                  break;
               case '>':
                  lenLHS = j;
                  if( constraintString[j+1] == '=' )
                  {
                     op = GREATER_EQUAL;
                     startRHS = j+2;
                  }
                  else
                  {
                     op = GREATER;
                     startRHS = j+1;
                  }
                  break;
               case '!':
                  if( constraintString[j+1] == '=' )
                  {
                     op = UNEQUAL;
                     lenLHS = j;
                     startRHS = j+2;
                  }
                  else
                     WALBERLA_ABORT( "Error in PhysicalCheck::checkConstraints(). Invalid operator '!'." );
                  break;
               default:
                  break;
            }
         }

         if ( op == INVALID )
            WALBERLA_ABORT( "Error in PhysicalCheck::checkConstraints(): No Operator found in " << constraintString );

         // use the FunctionParser to parse and solve the lhs and the rhs, respectively
         std::string lhs = constraintString.substr( 0, lenLHS );
         std::string rhs = constraintString.substr( startRHS, constraintString.size()-startRHS );

         FunctionParser funcParser;
         try
         {
            funcParser.parse(lhs);
         }
         catch( std::exception& )
         {
            WALBERLA_ABORT( "BadSyntaxException when checking constraints. Constraint: '" << lhs << "'" );
         }

         double resultLHS=0;
         try
         {
            resultLHS = funcParser.evaluate(symbolTable);
         }
         catch( std::exception& )
         {
            WALBERLA_ABORT( "UnknownSymbolException when checking constraints. Constraint: '" << lhs << "'"  );
         }

         try
         {
            funcParser.parse(rhs);
         }
         catch( std::exception& )
         {
            WALBERLA_ABORT( "BadSyntaxException when checking constraints. Constraint: '" << rhs << "'" );
         }

         double resultRHS=0;
         try
         {
            resultRHS = funcParser.evaluate(symbolTable);
         }
         catch( std::exception& )
         {
            WALBERLA_ABORT( "UnknownSymbolException when checking constraints. Constraint: '" << rhs << "'"  );
         }

         switch( op )
         {
            case LOWER:
               if( !(resultLHS < resultRHS) )
                  WALBERLA_ABORT( "Constraint '" << constraintString << "' failed." );
               break;
            case LOWER_EQUAL:
               if( !(resultLHS <= resultRHS) )
                  WALBERLA_ABORT( "Constraint '" << constraintString << "' failed." );
               break;
            case GREATER:
               if( !(resultLHS > resultRHS) )
                  WALBERLA_ABORT( "Constraint '" << constraintString << "' failed." );
               break;
            case GREATER_EQUAL:
               if( !(resultLHS >= resultRHS) )
                  WALBERLA_ABORT( "Constraint '" << constraintString << "' failed." );
               break;
            case UNEQUAL:
               if( floatIsEqual(resultLHS,resultRHS) )
                  WALBERLA_ABORT( "Constraint '" << constraintString << "' failed." );
               break;
            default:
               WALBERLA_ABORT( "Error in PhysicalCheck::checkConstraints(). Entered unreachable code." );
               break;
         }
      }

      return true;
   }

   void PhysicalCheck::writeUnitParameterRelations()
   {
      std::stringstream ss;
      ss << "Unit-Parameter-Relations:" << std::endl;
      for( auto i=unitParameterRelations_.begin(); i!=unitParameterRelations_.end(); ++i )
      {
         ss << i->first << " (" << getVarUnit(i->first) << ")" << std::endl;
      }
   }

   std::ostream& operator<<( std::ostream& os, PhysicalCheck& pc )
   {
      return os << pc.es_;
   }

} // namespace math
} // namespace walberla

#endif
