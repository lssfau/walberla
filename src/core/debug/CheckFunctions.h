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
//! \file CheckFunctions.h
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "PrintStacktrace.h"
#include "core/DataTypes.h"
#include "core/Macros.h"
#include "core/math/MathTrait.h"
#include "core/math/Utility.h"
#include "core/VectorTrait.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <sstream>
#include <type_traits>



/*******************************************************************************************************************//**
 * Some check functions. The two parameter versions have the advantage of printing the values of their arguments in
 * case of error.
 * Note: The checks will be performed in Debug and in Release mode! Do not use as replacements of
 * assertions in regular code! For assertions, use assert macros!
 *
 * Available macros are:
 * ---------------------
 * - WALBERLA_CHECK( X )                  - X must evaluate to true, otherwise the check fails
 * - WALBERLA_CHECK_NULLPTR( X )          - X must be a null pointer, otherwise the check fails
 * - WALBERLA_CHECK_NOT_NULLPTR( X )      - X must not be a null pointer, otherwise the check fails
 * - WALBERLA_CHECK_EQUAL( X, Y )         - X must be equal to Y, otherwise the check fails
 *                                          (for floating point types, use FLOAT_EQUAL or IDENTICAL)
 * - WALBERLA_CHECK_UNEQUAL( X, Y )       - X must not be equal to Y, otherwise the check fails
 *                                          (for floating point types, use FLOAT_UNEQUAL or NOT_IDENTICAL)
 * - WALBERLA_CHECK_FLOAT_EQUAL( X, Y )   - X must be equal to Y, otherwise the check fails
 *                                          (floating point types are considered equal even if they differ by a small epsilon)
 * - WALBERLA_CHECK_FLOAT_UNEQUAL( X, Y ) - X must not be equal to Y, otherwise the check fails
 *                                          (floating point types are only considered unequal if they differ by more than a small epsilon)
 * - WALBERLA_CHECK_IDENTICAL( X, Y )     - X must be identical to Y, otherwise the check fails
 *                                          (floating point types are only identical if they are bit-identical)
 * - WALBERLA_CHECK_NOT_IDENTICAL( X, Y ) - X must not be identical to Y, otherwise the check fails
 *                                          (floating point types are not identical if they are not bit-identical)
 * - WALBERLA_CHECK_LESS( X, Y )          - check fails if X<Y evaluates to false
 * - WALBERLA_CHECK_GREATER( X, Y )       - check fails if X>Y evaluates to false
 * - WALBERLA_CHECK_LESS_EQUAL( X, Y )    - check fails if X<=Y evaluates to false
 * - WALBERLA_CHECK_GREATER_EQUAL( X, Y ) - check fails if X>=Y evaluates to false
 *
 * For every macro, you can provide a custom message via an additional argument, e.g.:
 *    WALBERLA_CHECK( false, "Will always fail!" )
 * or
 *    WALBERLA_CHECK_LESS( 42, 5, "42 is larger than 5!" )
 *
 **********************************************************************************************************************/

#define WALBERLA_CHECK_1(X)                             { if( !walberla::debug::check_functions_detail::check              ( (X)      ) ) { walberla::debug::check_functions_detail::check              (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_NULLPTR_1(X)                     { if( !walberla::debug::check_functions_detail::check_nullptr      ( (X)      ) ) { walberla::debug::check_functions_detail::check_nullptr      ( (X),      #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_NOT_NULLPTR_1(X)                 { if( !walberla::debug::check_functions_detail::check_not_nullptr  ( (X)      ) ) { walberla::debug::check_functions_detail::check_not_nullptr  (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_EQUAL_2(X,Y)                     { if( !walberla::debug::check_functions_detail::check_equal        ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_equal        ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_UNEQUAL_2(X,Y)                   { if( !walberla::debug::check_functions_detail::check_unequal      ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_unequal      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_FLOAT_EQUAL_2(X,Y)               { if( !walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_FLOAT_UNEQUAL_2(X,Y)             { if( !walberla::debug::check_functions_detail::check_float_unequal( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_float_unequal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_3(X,Y,EPS)   { if( !walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), (EPS) ) ) { walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler(), (EPS) ); } }
#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_3(X,Y,EPS) { if( !walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), (EPS) ) ) { walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler(), (EPS) ); } }
#define WALBERLA_CHECK_IDENTICAL_2(X,Y)                 { if( !walberla::debug::check_functions_detail::check_identical    ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_identical    ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_NOT_IDENTICAL_2(X,Y)             { if( !walberla::debug::check_functions_detail::check_not_identical( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_not_identical( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_LESS_2(X,Y)                      { if( !walberla::debug::check_functions_detail::check_less         ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_less         ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_GREATER_2(X,Y)                   { if( !walberla::debug::check_functions_detail::check_greater      ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_greater      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_LESS_EQUAL_2(X,Y)                { if( !walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_CHECK_GREATER_EQUAL_2(X,Y)             { if( !walberla::debug::check_functions_detail::check_greater_equal( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_greater_equal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }



#define WALBERLA_CHECK_2(X,MSG)                             { if( !walberla::debug::check_functions_detail::check              ( (X)      ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check              (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_NULLPTR_2(X,MSG)                     { if( !walberla::debug::check_functions_detail::check_nullptr      ( (X)      ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_nullptr      ( (X),      #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_NOT_NULLPTR_2(X,MSG)                 { if( !walberla::debug::check_functions_detail::check_not_nullptr  ( (X)      ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_not_nullptr  (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_EQUAL_3(X,Y,MSG)                     { if( !walberla::debug::check_functions_detail::check_equal        ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_equal        ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_UNEQUAL_3(X,Y,MSG)                   { if( !walberla::debug::check_functions_detail::check_unequal      ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_unequal      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_FLOAT_EQUAL_3(X,Y,MSG)               { if( !walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_FLOAT_UNEQUAL_3(X,Y,MSG)             { if( !walberla::debug::check_functions_detail::check_float_unequal( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_unequal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_4(X,Y,EPS,MSG)   { if( !walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), (EPS) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ), (EPS) ); } }
#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_4(X,Y,EPS,MSG) { if( !walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), (EPS) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ), (EPS) ); } }
#define WALBERLA_CHECK_IDENTICAL_3(X,Y,MSG)                 { if( !walberla::debug::check_functions_detail::check_identical    ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_identical    ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_NOT_IDENTICAL_3(X,Y,MSG)             { if( !walberla::debug::check_functions_detail::check_not_identical( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_not_identical( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_LESS_3(X,Y,MSG)                      { if( !walberla::debug::check_functions_detail::check_less         ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_less         ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_GREATER_3(X,Y,MSG)                   { if( !walberla::debug::check_functions_detail::check_greater      ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_greater      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_LESS_EQUAL_3(X,Y,MSG)                { if( !walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_CHECK_GREATER_EQUAL_3(X,Y,MSG)             { if( !walberla::debug::check_functions_detail::check_greater_equal( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_greater_equal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }



#define WALBERLA_CHECK_3(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_NULLPTR_3(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NULLPTR_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NULLPTR_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NULLPTR_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_NOT_NULLPTR_3(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NOT_NULLPTR_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NOT_NULLPTR_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NOT_NULLPTR_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_UNEQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_UNEQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_UNEQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_UNEQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_FLOAT_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_FLOAT_UNEQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_UNEQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_UNEQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_UNEQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_2(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_2(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_IDENTICAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_IDENTICAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_IDENTICAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_IDENTICAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_NOT_IDENTICAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NOT_IDENTICAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NOT_IDENTICAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_NOT_IDENTICAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_LESS_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_LESS_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_LESS_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_LESS_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_GREATER_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_GREATER_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_GREATER_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_GREATER_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_LESS_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_LESS_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_LESS_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_LESS_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK_GREATER_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_GREATER_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_GREATER_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO
#define WALBERLA_CHECK_GREATER_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_A_CHECK_MACRO

#define WALBERLA_CHECK(...)                       WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_, __VA_ARGS__ )
#define WALBERLA_CHECK_NULLPTR(...)               WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_NULLPTR_, __VA_ARGS__ )
#define WALBERLA_CHECK_NOT_NULLPTR(...)           WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_NOT_NULLPTR_, __VA_ARGS__ )
#define WALBERLA_CHECK_EQUAL(...)                 WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_EQUAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_UNEQUAL(...)               WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_UNEQUAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_FLOAT_EQUAL(...)           WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_FLOAT_EQUAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_FLOAT_UNEQUAL(...)         WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_FLOAT_UNEQUAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_FLOAT_EQUAL_EPSILON(...)   WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_FLOAT_EQUAL_EPSILON_, __VA_ARGS__ )
#define WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON(...) WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_FLOAT_UNEQUAL_EPSILON_, __VA_ARGS__ )
#define WALBERLA_CHECK_IDENTICAL(...)             WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_IDENTICAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_NOT_IDENTICAL(...)         WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_NOT_IDENTICAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_LESS(...)                  WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_LESS_, __VA_ARGS__ )
#define WALBERLA_CHECK_GREATER(...)               WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_GREATER_, __VA_ARGS__ )
#define WALBERLA_CHECK_LESS_EQUAL(...)            WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_LESS_EQUAL_, __VA_ARGS__ )
#define WALBERLA_CHECK_GREATER_EQUAL(...)         WALBERLA_MACRO_OVERLOAD( WALBERLA_CHECK_GREATER_EQUAL_, __VA_ARGS__ )



/// \cond internal

namespace walberla {
namespace debug {
namespace check_functions_detail {

struct ExitHandler
{
   ExitHandler() {};
   ExitHandler( const std::string & message ) : message_( message ) {}
   void operator()( const std::string & checkErrorMessage );
private:
   std::string message_;
};

inline bool check( bool b ) { return b; }

template< typename T >
inline bool check_nullptr( T * p );

template< typename T >
inline bool check_nullptr( const shared_ptr<T> & p );

template< typename T >
inline bool check_not_nullptr( T * p );

template< typename T >
inline bool check_not_nullptr( const shared_ptr<T> & p );

template< typename T, typename U >
inline bool check_equal( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_equal( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_equal( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_unequal( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_unequal( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_unequal( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_float_equal( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_float_unequal( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_float_equal_eps( const T & lhs, const U & rhs, const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon );

template< typename T, typename U >
inline bool check_float_unequal_eps( const T & lhs, const U & rhs, const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon );

template< typename T, typename U >
inline bool check_identical( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_identical( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_identical( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_not_identical( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_not_identical( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_not_identical( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_less( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_less( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_less( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_greater( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_greater( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_greater( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_less_equal( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_less_equal( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_less_equal( const T & lhs, const U & rhs, const std::false_type & );

template< typename T, typename U >
inline bool check_greater_equal( const T & lhs, const U & rhs );

template< typename T, typename U >
inline bool check_greater_equal( const T & lhs, const U & rhs, const std::true_type & );

template< typename T, typename U >
inline bool check_greater_equal( const T & lhs, const U & rhs, const std::false_type & );



template< typename T >
void check( const char * const expression, const char * const filename, int line, T failFunc );

template< typename T, typename U >
void check_nullptr( T * p, const char * const ptrExpression, const char * const filename, int line, U failFunc );

template< typename T, typename U >
void check_nullptr( const shared_ptr<T> & p, const char * const ptrExpression, const char * const filename, int line, U failFunc );

template< typename T >
void check_not_nullptr( const char * const ptrExpression, const char * const filename, int line, T failFunc );

template< typename T, typename U, typename V >
void check_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                  const char * const filename, int line, V failFunc);

template< typename T, typename U, typename V >
void check_unequal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                    const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_float_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                        const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_float_unequal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc);

template< typename T, typename U, typename V >
void check_float_equal_eps( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                            const char * const filename, int line, V failFunc,
                            const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon );

template< typename T, typename U, typename V >
void check_float_unequal_eps( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                              const char * const filename, int line, V failFunc,
                              const typename VectorTrait<typename math::MathTrait<T,U>::LowType>::OutputType epsilon );

template< typename T, typename U, typename V >
void check_identical( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                      const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_not_identical( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_less( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                 const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_greater( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                    const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_less_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                       const char * const filename, int line, V failFunc );

template< typename T, typename U, typename V >
void check_greater_equal( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                          const char * const filename, int line, V failFunc );



template< typename T, typename U, typename V >
void printErrorAndExit( const T & lhs, const U & rhs, const char * const lhsExpression, const char * const rhsExpression,
                        const char * const opString, const char * const filename, int line, V failFunc );

template< typename T >
std::ostream & printValue( std::ostream & os, const T & value );

template< typename T >
std::ostream & printValue( std::ostream & os, const T & value, const std::true_type & );

template< typename T >
std::ostream & printValue( std::ostream & os, const T & value, const std::false_type & );

template< typename T >
std::ostream & printValue( std::ostream & os, const T * value );

template< typename T >
std::ostream & printValue( std::ostream & os, T * value );

std::ostream & printValue( std::ostream & os, char value );

std::ostream & printValue( std::ostream & os, unsigned char value );

std::ostream & printValue( std::ostream & os, float value );

std::ostream & printValue( std::ostream & os, double value );

std::ostream & printValue( std::ostream & os, long double value );

} // namespace check_functions_detail
} // namespace debug
} // namespace walberla

#include "CheckFunctions.impl.h"

/// \endcond
