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
//! \file ExprSystemInitFunction.cpp
//! \ingroup lbm
//! \author Tobias Schruff <tobias.schruff@gmail.com>
//
//======================================================================================================================

#include "ExprSystemInitFunction.h"


#if defined WALBERLA_CXX_COMPILER_IS_MSVC
#   pragma warning( push, 1 )
#   pragma warning( disable : 4706 )
#elif ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wpragmas"
#   pragma GCC diagnostic ignored "-Wsign-conversion"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wshorten-64-to-32"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#elif defined WALBERLA_CXX_COMPILER_IS_INTEL
#   pragma warning push
#   pragma warning( disable : 187  )
#   pragma warning( disable : 1599 )
#endif

#include "core/math/extern/exprtk.h"

#if defined WALBERLA_CXX_COMPILER_IS_MSVC
#   pragma warning( pop )
#elif ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#elif defined WALBERLA_CXX_COMPILER_IS_INTEL
//#   pragma warning pop // disabled because this leads to spilled warnings
#endif


namespace walberla {
namespace lbm {
namespace initializer {


ExprSystemInitFunction::ExprSystemInitFunction( const shared_ptr<StructuredBlockForest> & blocks )
   : blocks_( blocks ), cell_(), symbolTable_( new exprtk::symbol_table<real_t>() ), expr_( new std::map<std::string, Expression>() )
{
   // add constants
   symbolTable_->add_constant( "n_x", real_c( blocks_->getNumberOfXCells() ) );
   symbolTable_->add_constant( "n_y", real_c( blocks_->getNumberOfYCells() ) );
   symbolTable_->add_constant( "n_z", real_c( blocks_->getNumberOfZCells() ) );
         
   // add variables
   symbolTable_->add_variable( "x" , cell_[0] );
   symbolTable_->add_variable( "y" , cell_[1] );
   symbolTable_->add_variable( "z" , cell_[2] );

   // add global constants (pi, e, etc.)
   symbolTable_->add_constants();

   // register symbols with expressions
   (*expr_)[ "u_x" ].register_symbol_table( *symbolTable_ );
   (*expr_)[ "u_y" ].register_symbol_table( *symbolTable_ );
   (*expr_)[ "u_z" ].register_symbol_table( *symbolTable_ );
   (*expr_)[ "rho" ].register_symbol_table( *symbolTable_ );
}


ExprSystemInitFunction::~ExprSystemInitFunction()
{
   delete expr_;
   delete symbolTable_;
}


bool ExprSystemInitFunction::parse( const Config::BlockHandle & config )
{
   // parse expressions from config and compile them
   exprtk::parser<real_t> parser;

   bool valid = true;
         
   const std::string ux_expr_str = config.getParameter<std::string>( "u_x" , "0.0" );
   if( !parser.compile( ux_expr_str, (*expr_)[ "u_x" ] ) )
   {
      WALBERLA_LOG_WARNING( "Error in expression for u_x\n" <<
                              "Error     : " << parser.error() << "\n" <<
                              "Expression: " << ux_expr_str )

      valid = false;
   }

   const std::string uy_expr_str = config.getParameter<std::string>( "u_y" , "0.0" );
   if( !parser.compile( uy_expr_str, (*expr_)[ "u_y" ] ) )
   {
      WALBERLA_LOG_WARNING( "Error in expression for u_y\n" <<
                              "Error     : " << parser.error() << "\n" <<
                              "Expression: " << uy_expr_str )

      valid = false;
   }

   const std::string uz_expr_str = config.getParameter<std::string>( "u_z" , "0.0" );
   if( !parser.compile( uz_expr_str, (*expr_)[ "u_z" ] ) )
   {
      WALBERLA_LOG_WARNING( "Error in expression for u_z\n" <<
                              "Error     : " << parser.error() << "\n" <<
                              "Expression: " << uz_expr_str )

      valid = false;
   }

   const std::string rho_expr_str = config.getParameter<std::string>( "rho" , "1.0" );
   if( !parser.compile( rho_expr_str, (*expr_)[ "rho" ] ) )
   {
      WALBERLA_LOG_WARNING( "Error in expression for rho\n" <<
                            "Error     : " << parser.error() << "\n" <<
                            "Expression: " << rho_expr_str )

      valid = false;
   }

   return valid;
}


real_t ExprSystemInitFunction::getDensity()
{
   return (*expr_)[ "rho" ].value();
}


Vector3<real_t> ExprSystemInitFunction::getVelocity()
{
   return Vector3<real_t>( (*expr_)[ "u_x" ].value(), (*expr_)[ "u_y" ].value(), (*expr_)[ "u_z" ].value() );
}


void ExprSystemInitFunction::setGlobalCell( const Cell & cell )
{
   cell_[0] = real_c( cell.x() );
   cell_[1] = real_c( cell.y() );
   cell_[2] = real_c( cell.z() );
}


std::vector<real_t> ExprSystemInitFunction::operator()( const Cell & globalCell )
{
   // update variables
   setGlobalCell( globalCell );

   // get results
   const auto velocity = getVelocity();
   const auto density  = getDensity();

   // store in std::vector and return
   std::vector<real_t> results;
   results.push_back( density     );
   results.push_back( velocity[0] );
   results.push_back( velocity[1] );
   results.push_back( velocity[2] );

   return results;
}


} // namespace initializer
} // namespace lbm
} // namespace walberla