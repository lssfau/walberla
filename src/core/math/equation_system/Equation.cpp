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
//! \file Equation.cpp
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "Equation.h"
#include "Operator.h"
#include "Variable.h"

#include <algorithm>
#include <cmath>
#include <memory>


namespace walberla {
namespace math {

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // NODE
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   Node::Node( const double  value  ) : nodeType_(NT_CONSTANT), value_(value  ), opType_(OP_NO )           {}
   Node::Node( const VarPtr& var    ) : nodeType_(NT_VARIABLE), value_(FP_NAN), opType_(OP_NO ), var_(var) {}
   Node::Node(       OpType& opType ) : nodeType_(NT_OPERATOR), value_(FP_NAN), opType_(opType)            {}

   void Node::setVar( const VarPtr& var    ){ var_ = var; }

   void Node::collectVariables( VarMap& varMap ){
      switch (nodeType_)
      {
      case NT_CONSTANT:
         break;
      case NT_VARIABLE:
         if ( varMap.find( var_->getName() ) == varMap.end() )
            varMap[var_->getName()] = var_;
         break;
      case NT_OPERATOR:
         left_->collectVariables ( varMap );
         right_->collectVariables( varMap );
         break;
      default:
    	  WALBERLA_ABORT( "No correct node type" );
    	  break;
      }
   }

   uint_t Node::countUnknownVariables(){
      switch (nodeType_)
      {
      case NT_CONSTANT:
         return 0;
      case NT_VARIABLE:
         return var_->valid() ? 0 : 1;
      case NT_OPERATOR:
         return left_->countUnknownVariables() + right_->countUnknownVariables();
      default:
    	 WALBERLA_ABORT( "No correct node type" );
    	 return 0; // has no effect
    	 break;
      }
   }

   double Node::compute(){
      switch (nodeType_)
      {
      case NT_CONSTANT:
         return value_;
      case NT_VARIABLE:
         return var_->getValue();
      case NT_OPERATOR:
         return opType_(left_->compute(), right_->compute());
      default:
    	 WALBERLA_ABORT( "No correct node type" );
    	 return 0; // has no effect
    	 break;
      }
   }

   bool Node::findUnknownVariable(){
      switch (nodeType_)
      {
      case NT_CONSTANT:
         return false;
      case NT_VARIABLE:
         return !var_->valid();
      case NT_OPERATOR:
         if(left_->findUnknownVariable()){
            nodeDir_ = ND_LEFT;
            return true;
         }
         if(right_->findUnknownVariable()){
            nodeDir_ = ND_RIGHT;
            return true;
         }
         return false;
      default:
    	 WALBERLA_ABORT( "No correct node type" );
    	 return false; // has no effect
    	 break;
      }
   }

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // EQUATION
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   Equation::Equation( const NodePtr& root) : root_ (root)
   {
      root_->collectVariables( varMap_ );
   }

   bool Equation::evaluate(){
      if (!isEvaluatable())
    	 WALBERLA_ABORT( "Equation is not evaluatable" );

      double left  = root_->left_->compute();
      double right = root_->right_->compute();

      if ( std::isnan(left) && std::isnan(right) ){
         //WALBERLA_LOG_WARNING( "WARNING: Both values are NAN -> return true" );
         return true;
      } else if ( std::isinf(left) && std::isinf(right) ){
    	 //WALBERLA_LOG_WARNING( "WARNING: Both values are INF -> return true" );
         return true;
      }

      const double border = std::max(
         std::fabs(left/2e12 + right/2e12),
         std::fabs(left/2e12 - right/2e12) );

      return std::fabs( left - right ) < std::max( border, std::numeric_limits<double>::epsilon() );
   }

   VarPtr Equation::compute(){
      if (!isComputable())
    	 WALBERLA_ABORT( "Equation is not computable" );

      sort();

      root_->left_->var_->setValue( root_->right_->compute() );

      return root_->left_->var_;
   }


   void Equation::sort(){
      if ( root_->right_->findUnknownVariable() )
         root_->flip();
      else
         root_->left_->findUnknownVariable();

      while( root_->left_->nodeType_ == NT_OPERATOR ){
         if ( root_->left_->opType_ == OP_PLUS )
         {
            rotate( (root_->left_->nodeDir_ == ND_RIGHT), OP_MINUS, OP_MINUS );
         }
         else if ( root_->left_->opType_ == OP_MINUS )
         {
            rotate( false, OP_PLUS, OP_MINUS );
         }
         else if ( root_->left_->opType_ == OP_MULT )
         {
            rotate( (root_->left_->nodeDir_ == ND_RIGHT), OP_DIV, OP_DIV );
         }
         else if ( root_->left_->opType_ == OP_DIV )
         {
            rotate( false, OP_MULT, OP_DIV );
         }
         else if ( root_->left_->opType_ == OP_PROD )
         {
            //rotate( (root_->left_->nodeDir_ == ND_RIGHT), OP_ROOT, OP_LOG );
            rotate( (root_->left_->nodeDir_ == ND_RIGHT), OP_PROD, OP_LOG );
         }
         else if ( root_->left_->opType_ == OP_LOG )
         {
            //rotate( (root_->left_->nodeDir_ == ND_LEFT), OP_PROD, OP_ROOT );
            rotate( (root_->left_->nodeDir_ == ND_LEFT), OP_PROD, OP_PROD );
         }
         //else if ( root_->left_->opType_ == OP_ROOT )
         //{
         //   rotate( false, OP_PROD, OP_LOG );
         //}
         else
        	WALBERLA_ABORT( "Unknown operator" );
      }
   }

   void Equation::rotate(bool flip, OpType& leftOp, OpType& rightOp){
      NodePtr newNode;
      if ( root_->left_->nodeDir_ == ND_LEFT ){
         newNode = std::make_shared<Node>( leftOp );
         if (flip){
            newNode->left_  = root_->left_->right_;
            newNode->right_ = root_->right_;
         } else {
            newNode->right_ = root_->left_->right_;
            newNode->left_  = root_->right_;
         }
         root_->left_ = root_->left_->left_;
      } else {
         newNode = std::make_shared<Node>( rightOp );
         if (flip){
            newNode->right_ = root_->left_->left_;
            newNode->left_  = root_->right_;
         } else {
            newNode->left_  = root_->left_->left_;
            newNode->right_ = root_->right_;
         }
         root_->left_ = root_->left_->right_;
      }
      root_->right_ = newNode;
   }

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // OUTPUT
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

   std::ostream& operator<<( std::ostream& os, const Node & node ){
      switch (node.nodeType_)
      {
      case NT_CONSTANT:
         os << node.value_;
         break;
      case NT_VARIABLE:
         os << node.var_->getName();
         break;
      case NT_OPERATOR:
         if (node.opType_ == OP_EQUAL)
            os << *node.left_ << node.opType_ << *node.right_;
         else if( node.opType_ == OP_LOG )
            os << "log(" << *node.left_ << ", " << *node.right_ << ")";
         else
            os << "(" << *node.left_ << node.opType_ << *node.right_ << ")";
         break;
      default:
    	 WALBERLA_ABORT( "No correct node type" );
    	 break;
      }
      return os;
   }

   std::ostream& operator<<( std::ostream& os, const Equation & eq ){
      return os << *eq.root_;
   }

} // namespace math
} // namespace walberla
