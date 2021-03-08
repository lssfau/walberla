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
//! \file EquationSystem.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "FwdEquation.h"
#include "FwdVariable.h"

#include <map>
#include <string>
#include <vector>


#ifdef _MSC_VER
#pragma warning ( push, 1 )
#endif
#include <boost/graph/adjacency_list.hpp>
#ifdef _MSC_VER
#pragma warning ( pop )
#endif


namespace walberla {
namespace math {

   //===================================================================================================================
   //
   //  CLASS DEFINITION
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*!\brief Wrapper class to handle and solve an equation system, e.g. given by the equations
    * in an input file
    * \author Matthias Markl
    *
    * This class handles and solves the equations given in an input file in terms of a triangular
    * equation system. In order to do so, it employs boost-graphs to organize the equations in data structures.
    * Furthermore, equations and the variables that are solved for are hold in map-structures and can be
    * operated on from outside.
    *
    * The equations need to be given in the following form (e.g.):
    *
    * "'c'     = 'dx_L' / 'dt_L'"
    */
   class EquationSystem
   {
   private:
      // forward declaration of EquationParser class
      friend class EquationParser;

      using EqGraph = boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS>;

      using EqVertexMap = std::map<std::string, EqGraph::vertex_descriptor>;
      using VarVertexMap = std::map<std::string, EqGraph::vertex_descriptor>;

      using EqVertexMapIt = std::map<std::string, EqGraph::vertex_descriptor>::const_iterator;
      using VarVertexMapIt = std::map<std::string, EqGraph::vertex_descriptor>::const_iterator;

      using EqMap = std::map<std::string, EquationPtr>;
      using EqMapIt = std::map<std::string, EquationPtr>::const_iterator;

      EqMap       eqMap_;
      EqGraph     eqGraph_;
      EqVertexMap eqVertices_;

      VarMap       varMap_;
      VarVertexMap varVertices_;
   public:

      //**Get functions****************************************************************************
      /*! \name Get functions */
      //@{
      const EquationPtr& get      (const std::string& key) { return eqMap_[key]; }
      bool   isVarDefined         ( const std::string& var ) const;
      double getVarValue          ( const std::string& var ) const;
      void   getVarMap            ( std::map<std::string,double>& varMap ) const;
      size_t getNumberOfEquations () const;
      //@}
      //****************************************************************************************************************

      //**Output functions*************************************************************************
      /*! \name Output functions */
      //@{
      std::string writeEquations() const;
      std::string writeVariables() const;
      friend std::ostream& operator<<( std::ostream& os, EquationSystem& es );
      //@}
      //****************************************************************************************************************

      //**Utility functions************************************************************************
      /*! \name Utility functions */
      //@{
      void add(const std::string& key, const EquationPtr& eq);
      void remove(const std::string& key);
      bool solve();
      void match();
      void clear();
      //void push();
      //void pop();
      //@}
      //****************************************************************************************************************

   };

   using EquationSystemPtr = shared_ptr<EquationSystem>;

} // namespace math
} // namespace walberla
