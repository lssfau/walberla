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
//! \file PhysicalCheck.h
//! \ingroup core
//! \author David Staubach <david.staubach@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/config/Config.h"
#include "equation_system/EquationSystem.h"

#include <map>
#include <string>
#include <vector>


namespace walberla {
namespace math {

   //===================================================================================================================
   //
   //  CLASS DEFINITION
   //
   //===================================================================================================================

   //*******************************************************************************************************************
   /*!\brief Wrapper class to check for physical properties and correctness of given input parameters
    *
    * This class serves as an interface between the given parameters in an input file and the
    * math module classes that solve the equations and expressions of interest.
    * It can be created either by passing a vector of equations, a map of unit-parameter relations and
    * a vector of constraints, or by passing a BlockHandle to the input file.
    *
    * The blocks for the PhysicalCheck in the input file are bound to the following layout:
    *
    * Physical_Check {
    *    Equations {
    *       eq0 parameter1 = 23;
    *       eq1 parameter2 = 42;
    *
    *       eq2 var1 = parameter1 + cos(parameter2);
    *       eq3 var2 = parameter2 - 23;
    *    }
    *    Units {
    *       parameter1 m;
    *       parameter2 m^2/s;
    *    }
    *    Constraints {
    *       co0 var1 > 20;
    *       co1 var1 <= 30;
    *       co2 var2 != var1;
    *    }
    * }
    *
    * Geometry {
    *    BoundaryConditionXYZ {
    *       velocity 'parameter1';
    *       pressure 'var1 * 29.9';
    *    }
    *
    *    BoundaryConditionABC {
    *       velocity 'parameter1';
    *    }
    * }
    */
   class PhysicalCheck
   {
   public:

      //**Constructors*****************************************************************************
      PhysicalCheck();
      PhysicalCheck( const std::vector<std::string>& equations, const std::map< std::string, std::string >& unitParameterRelations, const std::vector<std::string>& constraints );
      PhysicalCheck( const Config::BlockHandle& configBlock );

      //**Utility functions************************************************************************
      /*! \name Utility functions */
      //@{
      void addBlock                  ( const Config::BlockHandle& configBlock );
      void addEquations              ( const std::vector<std::string>& equations );
      void addUnitParameterRelations ( const std::map< std::string, std::string >& unitParameterRelations );
      void addConstraints            ( const std::vector<std::string>& constraints );
      void completeConfig            ( const shared_ptr<Config>& config );
      //@}
      //****************************************************************************************************************

      //**Get functions****************************************************************************
      /*! \name Get functions */
      //@{
      bool   isDefined  ( const std::string& varName );
      double getVarValue( const std::string& varName );
      //@}
      //****************************************************************************************************************

      //**Output functions*************************************************************************
      /*! \name Output functions */
      //@{
      friend std::ostream& operator<<( std::ostream& os, PhysicalCheck& pc );
      //@}
      //****************************************************************************************************************
   private:
      //**Private functions to setup data layout from the given input file*************************
      /*! \name Private functions */
      //@{
      void addEquation                   ( const std::string& equation );
      std::string getParametrizationTerm ( const std::string& varName );
      void solve                         ();
      bool checkConstraints              ();
      void completeBlock                 ( Config::Block& configBlock, const std::map<std::string,double>& symbolTable );
      std::string getVarUnit             ( const std::string& varName );
      bool setVarUnit                    ( const std::string& varName, const std::string& unit, const int expo );
      void getVarMap                     ( std::map<std::string,double>& varMap );
      void writeUnitParameterRelations   ();
      //@}
      //****************************************************************************************************************

      EquationSystem es_;
      std::map<std::string, std::map<std::string, int> > unitParameterRelations_;
      std::vector<std::string> constraints_;
      bool solved_;
   };

} // namespace math
} // namespace walberla
