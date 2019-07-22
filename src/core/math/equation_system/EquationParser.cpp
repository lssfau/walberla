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
//! \file EquationParser.cpp
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"
#ifdef WALBERLA_BUILD_WITH_BOOST

#include "Equation.h"
#include "EquationParser.h"
#include "Operator.h"
#include "Variable.h"
#include "core/math/Constants.h"
#include "core/StringUtility.h"

#include <memory>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARSE UTIL
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#define THROW(msg, str, index) {\
   std::stringstream ss;\
   ss << (msg) << " -> [" << (str) << "] at [" << (index) << "]";\
   throw std::runtime_error( ss.str() );\
}

namespace walberla {
namespace math {

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARSE NUMBER
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool EquationParser::checkNumber( const std::string& str, size_t& index ) const
{
   if(str[index] == '+' || str[index] == '-')
      return isdigit(str[index+1]) != int(0);
   return isdigit(str[index]) != int(0);
}

NodePtr EquationParser::parseNumber( const std::string& str, size_t& index ) const
{
   size_t start = index;
   double value;

   if(str[index] == '+' || str[index] == '-')
      ++index;

   if( isdigit(str[index]) == int(0) )
      THROW( "No number found", str, index );

   while( isdigit(str[++index]) );

   // numbers are allowed to end with a '.'
   if ( str[index] == '.' )
      while( isdigit(str[++index]) != int(0) );

   if ( str[index] == 'e' || str[index] == 'E' ){
      ++index;
      size_t estart = index;
      if( str[index] == '+' || str[index] == '-' )
         ++index;
      if( isdigit(str[index]) == int(0) )
         THROW( "Number ends with 'e'", str, index );
      while( isdigit(str[++index]) != int(0) );

      value =  std::stod( str.substr(start, estart-start-1) ) *
            pow(10.0, std::stoi( str.substr(estart, index-estart) ) );
   } else {
      value = std::stod( str.substr(start, index-start) );
   }

   return std::make_shared<Node>( value );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARSE NT_VARIABLE
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool EquationParser::checkVariable( const std::string& str, size_t& index ) const
{
   if(str[index] == '+' || str[index] == '-')
      return (isalpha(str[index+1]) != int(0)) || (str[index+1] == '\'');
   return (isalpha(str[index]) != int(0)) || (str[index] == '\'');
}

NodePtr EquationParser::parseVariable( const std::string& str, size_t& index ) const
{
   bool sign = false;
   if(str[index] == '+' || str[index] == '-'){
      sign = str[index] == '-';
      ++index;
   }

   // variables can start with a '
   bool marked = (str[index] == '\'');
   if ( marked )
      ++index;

   if ( isalpha(str[index]) == int(0) )
      THROW( "Variable name has to start with a letter", str, index );

   size_t start = index;
   size_t len;

   for (
         len=1, ++index;
         (isalpha(str[index]) != int(0)) || (isdigit(str[index]) != int(0)) || str[index] == '_';
         ++len, ++index );

   if ( marked )
   {
      if (str[index] == '\'' ){
         ++index;
      } else {
         THROW( "Variable declaration has to end with '", str, index );
      }
   }

   std::string name = str.substr(start, len);

   VarPtr varPtr;
   if ( es_.varMap_.find(name) != es_.varMap_.end() ){
      varPtr = es_.varMap_[name];
   } else {
      varPtr = std::make_shared<Var>( name );
      es_.varMap_[name] = varPtr;
   }

   NodePtr nodePtr;
   if (sign){
      nodePtr.reset( new Node(OP_MULT) );
      nodePtr->left().reset(  new Node(    -1) );
      nodePtr->right().reset( new Node(varPtr) );
   } else {
      nodePtr.reset( new Node(varPtr) );
   }
   return nodePtr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARSE FUNCTION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool EquationParser::checkFunction( const std::string& str, size_t& index ) const
{
   return (str.substr(index, 4) == "exp(") ||
         (str.substr(index, 3) == "ln(") ||
         (str.substr(index, 5) == "sqrt(");
}

NodePtr EquationParser::parseFunction( const std::string& str, size_t& index ) const
{
   OpFunction opFunc;
   if ( str.substr(index, 4) == "exp(" ){
      opFunc = OP_FUNC_EXP;
      index += 4;
   } else if ( str.substr(index, 3) == "ln(" ){
      opFunc = OP_FUNC_LN;
      index += 3;
   } else if ( str.substr(index, 5) == "sqrt(" ){
      opFunc = OP_FUNC_SQRT;
      index += 5;
   } else {
      THROW( "Found no function", str, index );
   }
   NodePtr nodePtr = parseExpression(str, index);
   if ( ! (str[index] == ')') )
      THROW( "Found no enclosing paranthesis", str, index );
   ++index;

   NodePtr funcPtr;

   switch(opFunc)
   {
   case OP_FUNC_EXP:
      funcPtr = std::make_shared<Node>( OP_PROD );
      funcPtr->left()  = std::make_shared<Node>( math::e  );
      funcPtr->right() = nodePtr;
      return funcPtr;
   case OP_FUNC_LN:
      funcPtr = std::make_shared<Node>( OP_LOG );
      funcPtr->right() = std::make_shared<Node>( math::e  );
      funcPtr->left()  = nodePtr;
      return funcPtr;
   case OP_FUNC_SQRT:
      funcPtr = std::make_shared<Node>( OP_PROD );
      funcPtr->left()  = nodePtr;
      funcPtr->right() = std::make_shared<Node>( 0.5 );
      return funcPtr;
   default:
      WALBERLA_ABORT( "Function not yet defined" );
      break;
   }
   return funcPtr; // has no effect
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARSE EXPRESSION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************************************************************
/*!
*   Parses a given term inside the current expression
*
*   Goal: creating a binary tree node for the corresponding term
*
*   \param  str           string representation of the term
*   \param  index         index of the current term within the equation
*/
//**********************************************************************************************************************
NodePtr EquationParser::parseTerm( const std::string& str, size_t& index ) const
{
   NodePtr nodePtr;

   // check for the type of the current term
   if ( str[index] == '(' ){
      nodePtr = parseExpression(str, ++index);
      if ( ! (str[index] == ')') )
         THROW( "Found no enclosing paranthesis", str, index );
      ++index;
   } else if ( checkFunction(str, index) ) {
      nodePtr = parseFunction(str, index);
   } else if ( checkVariable(str, index) ) {
      nodePtr = parseVariable(str, index);
   } else if ( checkNumber(str, index) ) {
      nodePtr = parseNumber(str, index);
   } else {
      THROW( "Found no paranthesis, variable or number", str, index );
   }
   return nodePtr;
}

//**********************************************************************************************************************
/*!
*   Parses a given expression inside the current Equation
*
*   Goal: modeling the current expression in binary tree format
*
*   \param  str           string representation of the expression
*   \param  index         index of the current term within the equation
*/
//**********************************************************************************************************************
NodePtr EquationParser::parseExpression( const std::string& str, size_t& index ) const
{
   NodePtr leftPtr = parseTerm(str, index);

   // index has been shifted to next term by parseTerm function
   size_t indexFstOp = index;
   if ( str[index] == '=' || str.size() == index || str[index] == ')'){
      return leftPtr;
   } else if ( isop(str[index]) ){
      ++index;
   } else {
      THROW( "Found no operator or equal", str, index );
   }

   NodePtr rightPtr;
   bool run = true;
   do {
      rightPtr = parseTerm(str, index);

      size_t indexSndOp = index;
      if ( str[index] == '=' || str.size() == index || str[index] == ')'){
         NodePtr nodePtr ( new Node (getOp(str[indexFstOp])) );
         nodePtr->left()  = leftPtr;
         nodePtr->right() = rightPtr;
         return nodePtr;
      } else if ( isop(str[index]) ){
         ++index;
      } else {
         THROW( "Found no operator or equal", str, index );
      }

      OpType& opFst = getOp(str[indexFstOp]);
      OpType& opSnd = getOp(str[indexSndOp]);

      if (opFst >= opSnd){
         NodePtr nodePtr ( new Node (getOp(str[indexFstOp])) );
         nodePtr->left()  = leftPtr;
         nodePtr->right() = rightPtr;
         leftPtr = nodePtr;
         indexFstOp = indexSndOp;
      } else {
         break;
      }
   } while ( run );

   index = indexFstOp+1;
   rightPtr = parseExpression(str, index);

   NodePtr nodePtr ( new Node (getOp(str[indexFstOp])) );
   nodePtr->left()  = leftPtr;
   nodePtr->right() = rightPtr;
   return nodePtr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PARSE EQUATION
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//**********************************************************************************************************************
/*!
*   Parses a given Equation and constructs a binary tree as representation
*
*   Goal: Equation is modeled in a binary tree format to be solved later on
*
*   \param  str           string representation of the equation
*   \param  index         index of the current term (is always zero here)
*/
//**********************************************************************************************************************
EquationPtr EquationParser::parseEquation( const std::string& str, size_t& index )
{
   // removing leading and trailing spaces of input string
   std::string trimmedStr = string_trim_copy(str);
   // removing spaces inside the trimmed string
   trimmedStr.erase(std::remove(trimmedStr.begin(), trimmedStr.end(), ' '), trimmedStr.end());
   NodePtr leftPtr = parseExpression(trimmedStr, index);
   if ( ! (trimmedStr[index] == '=') )
      THROW( "Found no equal sign in equation", str, index );
   ++index;

   NodePtr rightPtr = parseExpression(trimmedStr, index);

   NodePtr nodePtr ( new Node (OP_EQUAL) );
   nodePtr->left()  = leftPtr;
   nodePtr->right() = rightPtr;

   return std::make_shared<Equation>( nodePtr );
}

} // namespace math
} // namespace walberla

#endif