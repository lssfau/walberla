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
//! \file EquationSystem.cpp
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#include "waLBerlaDefinitions.h"
#ifdef WALBERLA_BUILD_WITH_BOOST

#include "Equation.h"
#include "EquationSystem.h"
#include "Variable.h"
#include "core/Abort.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"


#ifdef _MSC_VER
#pragma warning ( push, 1 )
#pragma warning ( disable: 4701 )
#endif

#include <boost/graph/adjacency_list_io.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#ifdef _MSC_VER
#pragma warning ( pop )
#endif

namespace walberla {
namespace math {

void EquationSystem::add(const std::string& key, const EquationPtr& eq)
{
   if ( eqMap_.find(key) != eqMap_.end() )
      WALBERLA_ABORT( "Equation already exists" );

   eqMap_[key]      = eq;
   eqVertices_[key] = boost::add_vertex(eqGraph_);

   for(VarMapIt it = eq->getVarMap().begin(); it != eq->getVarMap().end(); ++it)
   {
      if ( varVertices_.find(it->first) == varVertices_.end() )
      {
         varVertices_[it->first] = boost::add_vertex(eqGraph_);
      }
      boost::add_edge(eqVertices_[key], varVertices_[it->first], eqGraph_);
   }
}

void EquationSystem::clear( )
{
   eqMap_.clear();
   varMap_.clear();
}

void EquationSystem::remove(const std::string& key)
{
   eqMap_.erase(key);
}

void EquationSystem::match()
{
   std::cout << "\nEquation Nodes:\n";
   for(EqMapIt it = eqMap_.begin(); it != eqMap_.end(); ++it)
   {
      //std::cout << *it->second << "\n";
      WALBERLA_LOG_RESULT( *it->second );
   }

   std::cout << "\nVariable Nodes:\n";
   for(VarMapIt it = varMap_.begin(); it != varMap_.end(); ++it)
   {
      //std::cout << *it->second << "\n";
      WALBERLA_LOG_RESULT( *it->second );
   }

   //std::cout << "\nInput Graph:\n" << boost::write( eqGraph_ );
   WALBERLA_LOG_RESULT( "\nInput Graph:\n" << boost::write( eqGraph_ ) );

   std::vector<EqGraph::vertex_descriptor> mate( boost::num_vertices(eqGraph_) );
   WALBERLA_ASSERT( boost::checked_edmonds_maximum_cardinality_matching(eqGraph_, &mate[0]) );

   WALBERLA_LOG_RESULT( "Maximum matching:" );
   EqGraph::vertex_iterator vi;
   EqGraph::vertex_iterator vi_end;
   for(boost::tie(vi,vi_end) = vertices(eqGraph_); vi != vi_end; ++vi)
      if (mate[*vi] != boost::graph_traits<EqGraph>::null_vertex() && *vi < mate[*vi])
         //std::cout << "{" << *vi << ", " << mate[*vi] << "}" << std::endl;
         WALBERLA_LOG_RESULT( "{" << *vi << ", " << mate[*vi] << "}" );
}

bool EquationSystem::solve()
{
   bool change = true;
   while(change){
      change = false;
      EqMapIt it = eqMap_.begin();
      for ( ; it != eqMap_.end(); ++it ){
         EquationPtr eq = it->second;
         if (eq->isComputable()){
            VarPtr varPtr = eq->compute();
            change = true;
         } else if (eq->isEvaluatable()){
            //eq->evaluate();
            if (!eq->evaluate())
               //std::cout << "Equation is not evaluatable! " << *eq << " -> " << eq->root_->left()->compute() << "!=" << eq->root_->right()->compute() << std::endl;
               WALBERLA_ABORT( "Equation is not evaluatable! " << *eq << " -> " << eq->root_->left()->compute() << "!=" << eq->root_->right()->compute() );
         /*} else {
            //std::cout << "Equation '" << *eq << "' is neither computable nor evaluatable!" << std::endl;
            WALBERLA_LOG_RESULT( "Equation '" << *eq << "' is neither computable nor evaluatable!" );*/
         }
      }
   }

   /*for ( VarMapIt it = varMap_.begin(); it != varMap_.end(); ++it ){
      //std::cout << *it->second << std::endl;
      WALBERLA_LOG_RESULT( *it->second );
   }*/
   bool evaluatable = true;
   EqMapIt it = eqMap_.begin();
   for ( ; it != eqMap_.end(); ++it ){
      EquationPtr eq = it->second;
      if ( !eq->isEvaluatable() || !eq->evaluate()){
         evaluatable = false;
         //std::cout << "Equation is not evaluatable! " << *eq << " -> " << eq->root_->left()->compute() << "!=" << eq->root_->right()->compute() << std::endl;
         //WALBERLA_LOG_RESULT( "Equation is not evaluatable! " << *eq << " -> " << eq->root_->left()->compute() << "!=" << eq->root_->right()->compute() );
      }
   }
   /*if (evaluatable)
      //std::cout << "All Equations are evaluatable!" << std::endl;
      WALBERLA_LOG_RESULT( "All Equations are evaluatable!" );*/

   return evaluatable;
}

bool EquationSystem::isVarDefined( const std::string& var ) const
{
   return varMap_.find(var) != varMap_.end();
}

double EquationSystem::getVarValue( const std::string& var ) const
{
   return varMap_.find(var)->second->getValue();
}

void EquationSystem::getVarMap( std::map<std::string,double>& varMap ) const
{
   for( auto it = varMap_.begin(); it != varMap_.end(); ++it )
   {
      varMap.insert( std::pair<std::string,double>(it->first,it->second->getValue()) );
   }
}

size_t EquationSystem::getNumberOfEquations() const
{
   return eqMap_.size();
}

std::string EquationSystem::writeEquations() const
{
   std::stringstream ss;
   ss << "Equations to solve:" << std::endl;
   for( auto it = eqMap_.begin(); it != eqMap_.end(); ++it )
   {
      ss << *it->second << std::endl;
   }
   return ss.str();
}

std::string EquationSystem::writeVariables() const
{
   std::stringstream ss;
   ss << "Solution for each variable:" << std::endl;
   for( auto it = varMap_.begin(); it != varMap_.end(); ++it )
   {
      ss << *it->second << std::endl;
   }
   return ss.str();
}

std::ostream& operator<<( std::ostream& os, EquationSystem& es )
{
   os << es.writeEquations() << es.writeVariables();
   return os;
}

} // namespace math
} // namespace walberla

#endif
