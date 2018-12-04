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
//! \file Equation.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "FwdEquation.h"
#include "FwdOperator.h"
#include "FwdVariable.h"
#include "core/Abort.h"


namespace walberla {
namespace math {

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // NODE
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   class Node
   {
      friend class Equation;

   private:
      const NodeType    nodeType_;
      const double      value_;
            NodeDir     nodeDir_;
            OpType&     opType_;

      VarPtr    var_;
      NodePtr   left_;
      NodePtr   right_;

   public:
      Node( const double  value  );
      Node( const VarPtr& var   );
      Node(       OpType& opType );
   private:
      Node& operator=( const Node& ){ return *this; }

   public:
      friend std::ostream& operator<<( std::ostream& os, const Node & node );

   private:
      uint_t countUnknownVariables();
      bool         findUnknownVariable();

      void collectVariables( VarMap& varMap );

      void flip(){ left_.swap(right_); }

   public:
      double compute();

      NodePtr& left () { return left_;  }
      NodePtr& right() { return right_; }

      void setVar( const VarPtr& var    );
      void setOp (       OpType& opType );
   };
   // end class Node

   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   // EQUATION
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   class Equation
   {
   public:
      NodePtr root_;
      VarMap  varMap_;

   public:
      Equation( const NodePtr& root );

   private:
      uint_t countUnknownVariables(){ return uint_c( root_->countUnknownVariables() ); }

   public:
      friend std::ostream& operator<<( std::ostream& os, const Equation & eq );

      bool isComputable()  { return countUnknownVariables() == 1; }
      bool isEvaluatable() { return countUnknownVariables() == 0; }

      bool   evaluate();
      VarPtr compute();

   private:
      void sort();
      void rotate(bool flip, OpType& leftOp, OpType& rightOp);

   public:
      const VarMap& getVarMap() { return varMap_; }
   };
   // end class Equation

} // namespace math
} // namespace walberla
