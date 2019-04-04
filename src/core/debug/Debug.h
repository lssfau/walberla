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
//! \file Debug.h
//! \ingroup core
//! \author Christian Feichtinger
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief File for all Debug macros and classes.
//
//======================================================================================================================

#pragma once

#ifndef NDEBUG
#   include "CheckFunctions.h"
#   include <functional>
#   include <string>
#   include "core/Macros.h"
#endif



/**
 * \file Debug.h
 * \details
 * Some advanced assert macros that print a nice error message on failure.
 *
 * Examples:
 * ---------------------
 * - WALBERLA_ASSERT( X )                  - X must evaluate to true, otherwise the assertion fails
 * - WALBERLA_ASSERT_NULLPTR( X )          - X must be a null pointer, otherwise the assertion fails
 * - WALBERLA_ASSERT_NOT_NULLPTR( X )      - X must not be a null pointer, otherwise the assertion fails
 * - WALBERLA_ASSERT_EQUAL( X, Y )         - X must be equal to Y, otherwise the assertion fails
 *                                           (for floating point types, use FLOAT_EQUAL or IDENTICAL)
 * - WALBERLA_ASSERT_UNEQUAL( X, Y )       - X must not be equal to Y, otherwise the assertion fails
 *                                           (for floating point types, use FLOAT_UNEQUAL or NOT_IDENTICAL)
 * - WALBERLA_ASSERT_FLOAT_EQUAL( X, Y )   - X must be equal to Y, otherwise the assertion fails
 *                                           (floating point types are considered equal even if they differ by a small epsilon)
 * - WALBERLA_ASSERT_FLOAT_UNEQUAL( X, Y ) - X must not be equal to Y, otherwise the assertion fails
 *                                           (floating point types are only considered unequal if they differ by more than a small epsilon)
 * - WALBERLA_ASSERT_IDENTICAL( X, Y )     - X must be identical to Y, otherwise the assertion fails
 *                                           (floating point types are only identical if they are bit-identical)
 * - WALBERLA_ASSERT_NOT_IDENTICAL( X, Y ) - X must not be identical to Y, otherwise the assertion fails
 *                                           (floating point types are not identical if they are not bit-identical)
 * - WALBERLA_ASSERT_LESS( X, Y )          - assertion fails if X<Y evaluates to false
 * - WALBERLA_ASSERT_GREATER( X, Y )       - assertion fails if X>Y evaluates to false
 * - WALBERLA_ASSERT_LESS_EQUAL( X, Y )    - assertion fails if X<=Y evaluates to false
 * - WALBERLA_ASSERT_GREATER_EQUAL( X, Y ) - assertion fails if X>=Y evaluates to false
 *
 * For every macro, you can provide a custom message via an additional argument, e.g.:
 *    WALBERLA_ASSERT( false, "Will always fail!" )
 * or
 *    WALBERLA_ASSERT_LESS( 42, 5, "42 is larger than 5!" )
 *
 */
#ifndef NDEBUG

#define WALBERLA_ASSERT_1(X)                 { if( !walberla::debug::check_functions_detail::check              ( (X)      ) ) { walberla::debug::check_functions_detail::check              (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_NULLPTR_1(X)         { if( !walberla::debug::check_functions_detail::check_nullptr      ( (X)      ) ) { walberla::debug::check_functions_detail::check_nullptr      ( (X),      #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_NOT_NULLPTR_1(X)     { if( !walberla::debug::check_functions_detail::check_not_nullptr  ( (X)      ) ) { walberla::debug::check_functions_detail::check_not_nullptr  (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_EQUAL_2(X,Y)         { if( !walberla::debug::check_functions_detail::check_equal        ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_equal        ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_UNEQUAL_2(X,Y)       { if( !walberla::debug::check_functions_detail::check_unequal      ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_unequal      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_FLOAT_EQUAL_2(X,Y)   { if( !walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_2(X,Y) { if( !walberla::debug::check_functions_detail::check_float_unequal( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_float_unequal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_IDENTICAL_2(X,Y)     { if( !walberla::debug::check_functions_detail::check_identical    ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_identical    ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_NOT_IDENTICAL_2(X,Y) { if( !walberla::debug::check_functions_detail::check_not_identical( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_not_identical( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_LESS_2(X,Y)          { if( !walberla::debug::check_functions_detail::check_less         ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_less         ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_GREATER_2(X,Y)       { if( !walberla::debug::check_functions_detail::check_greater      ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_greater      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_LESS_EQUAL_2(X,Y)    { if( !walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }
#define WALBERLA_ASSERT_GREATER_EQUAL_2(X,Y) { if( !walberla::debug::check_functions_detail::check_greater_equal( (X), (Y) ) ) { walberla::debug::check_functions_detail::check_greater_equal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler() ); } }

#define WALBERLA_ASSERT_2(X,MSG)                         { if( !walberla::debug::check_functions_detail::check                  ( (X)      ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check              (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_NULLPTR_2(X,MSG)                 { if( !walberla::debug::check_functions_detail::check_nullptr          ( (X)      ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_nullptr      ( (X),      #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_NOT_NULLPTR_2(X,MSG)             { if( !walberla::debug::check_functions_detail::check_not_nullptr      ( (X)      ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_not_nullptr  (           #X,     __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_EQUAL_3(X,Y,MSG)                 { if( !walberla::debug::check_functions_detail::check_equal            ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_equal        ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_UNEQUAL_3(X,Y,MSG)               { if( !walberla::debug::check_functions_detail::check_unequal          ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_unequal      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_FLOAT_EQUAL_3(X,Y,MSG)           { if( !walberla::debug::check_functions_detail::check_float_equal      ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_equal  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_3(X,Y,MSG)         { if( !walberla::debug::check_functions_detail::check_float_unequal    ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_unequal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_3(X,Y,EPS)   { if( !walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), (EPS) ) ) { walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( ), (EPS) ); } }
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_3(X,Y,EPS) { if( !walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), (EPS) ) ) { walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( ), (EPS) ); } }
#define WALBERLA_ASSERT_IDENTICAL_3(X,Y,MSG)             { if( !walberla::debug::check_functions_detail::check_identical        ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_identical    ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_NOT_IDENTICAL_3(X,Y,MSG)         { if( !walberla::debug::check_functions_detail::check_not_identical    ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_not_identical( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_LESS_3(X,Y,MSG)                  { if( !walberla::debug::check_functions_detail::check_less             ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_less         ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_GREATER_3(X,Y,MSG)               { if( !walberla::debug::check_functions_detail::check_greater          ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_greater      ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_LESS_EQUAL_3(X,Y,MSG)            { if( !walberla::debug::check_functions_detail::check_less_equal        ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_less_equal   ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }
#define WALBERLA_ASSERT_GREATER_EQUAL_3(X,Y,MSG)         { if( !walberla::debug::check_functions_detail::check_greater_equal    ( (X), (Y) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_greater_equal( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ) ); } }

#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_4(X,Y,EPS,MSG)   { if( !walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), (EPS) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_equal_eps  ( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ), (EPS) ); } }
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_4(X,Y,EPS,MSG) { if( !walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), (EPS) ) ) { std::stringstream _ss_; _ss_ << MSG; walberla::debug::check_functions_detail::check_float_unequal_eps( (X), (Y), #X, #Y, __FILE__, __LINE__, walberla::debug::check_functions_detail::ExitHandler( _ss_.str() ), (EPS) ); } }

#define WALBERLA_ASSERT_3(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_NULLPTR_3(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NULLPTR_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NULLPTR_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NULLPTR_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_NOT_NULLPTR_3(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NOT_NULLPTR_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NOT_NULLPTR_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NOT_NULLPTR_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_UNEQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_UNEQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_UNEQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_UNEQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_FLOAT_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_FLOAT_UNEQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_2(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_2(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_IDENTICAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_IDENTICAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_IDENTICAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_IDENTICAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_NOT_IDENTICAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NOT_IDENTICAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NOT_IDENTICAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_NOT_IDENTICAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_LESS_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_LESS_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_LESS_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_LESS_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_GREATER_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_GREATER_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_GREATER_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_GREATER_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_LESS_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_LESS_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_LESS_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_LESS_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT_GREATER_EQUAL_1(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_GREATER_EQUAL_4(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_GREATER_EQUAL_5(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO
#define WALBERLA_ASSERT_GREATER_EQUAL_6(...) THIS_IS_SUPPOSED_TO_FAIL___YOU_MADE_AN_ERROR_WHEN_USING_AN_ASSERT_MACRO

#define WALBERLA_ASSERT(...)                       WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_, __VA_ARGS__ )
#define WALBERLA_ASSERT_NULLPTR(...)               WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_NULLPTR_, __VA_ARGS__ )
#define WALBERLA_ASSERT_NOT_NULLPTR(...)           WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_NOT_NULLPTR_, __VA_ARGS__ )
#define WALBERLA_ASSERT_EQUAL(...)                 WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_EQUAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_UNEQUAL(...)               WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_UNEQUAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_FLOAT_EQUAL(...)           WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_FLOAT_EQUAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_FLOAT_UNEQUAL(...)         WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_FLOAT_UNEQUAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON(...)   WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_FLOAT_EQUAL_EPSILON_, __VA_ARGS__ )
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON(...) WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_FLOAT_UNEQUAL_EPSILON_, __VA_ARGS__ )
#define WALBERLA_ASSERT_IDENTICAL(...)             WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_IDENTICAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_NOT_IDENTICAL(...)         WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_NOT_IDENTICAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_LESS(...)                  WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_LESS_, __VA_ARGS__ )
#define WALBERLA_ASSERT_GREATER(...)               WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_GREATER_, __VA_ARGS__ )
#define WALBERLA_ASSERT_LESS_EQUAL(...)            WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_LESS_EQUAL_, __VA_ARGS__ )
#define WALBERLA_ASSERT_GREATER_EQUAL(...)         WALBERLA_MACRO_OVERLOAD( WALBERLA_ASSERT_GREATER_EQUAL_, __VA_ARGS__ )

#else // NDBEUG

#define WALBERLA_ASSERT_1(X)                 (void(0));
#define WALBERLA_ASSERT_NULLPTR_1(X)         (void(0));
#define WALBERLA_ASSERT_NOT_NULLPTR_1(X)     (void(0));
#define WALBERLA_ASSERT_EQUAL_2(X,Y)         (void(0));
#define WALBERLA_ASSERT_UNEQUAL_2(X,Y)       (void(0));
#define WALBERLA_ASSERT_FLOAT_EQUAL_2(X,Y)   (void(0));
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_2(X,Y) (void(0));
#define WALBERLA_ASSERT_IDENTICAL_2(X,Y)     (void(0));
#define WALBERLA_ASSERT_NOT_IDENTICAL_2(X,Y) (void(0));
#define WALBERLA_ASSERT_LESS_2(X,Y)          (void(0));
#define WALBERLA_ASSERT_GREATER_2(X,Y)       (void(0));
#define WALBERLA_ASSERT_LESS_EQUAL_2(X,Y)    (void(0));
#define WALBERLA_ASSERT_GREATER_EQUAL_2(X,Y) (void(0));

#define WALBERLA_ASSERT_2(X,MSG)                 (void(0));
#define WALBERLA_ASSERT_NULLPTR_2(X,MSG)         (void(0));
#define WALBERLA_ASSERT_NOT_NULLPTR_2(X,MSG)     (void(0));
#define WALBERLA_ASSERT_EQUAL_3(X,Y,MSG)         (void(0));
#define WALBERLA_ASSERT_UNEQUAL_3(X,Y,MSG)       (void(0));
#define WALBERLA_ASSERT_FLOAT_EQUAL_3(X,Y,MSG)   (void(0));
#define WALBERLA_ASSERT_FLOAT_UNEQUAL_3(X,Y,MSG) (void(0));
#define WALBERLA_ASSERT_IDENTICAL_3(X,Y,MSG)     (void(0));
#define WALBERLA_ASSERT_NOT_IDENTICAL_3(X,Y,MSG) (void(0));
#define WALBERLA_ASSERT_LESS_3(X,Y,MSG)          (void(0));
#define WALBERLA_ASSERT_GREATER_3(X,Y,MSG)       (void(0));
#define WALBERLA_ASSERT_LESS_EQUAL_3(X,Y,MSG)    (void(0));
#define WALBERLA_ASSERT_GREATER_EQUAL_3(X,Y,MSG) (void(0));

#define WALBERLA_ASSERT(...)               (void(0));
#define WALBERLA_ASSERT_NULLPTR(...)       (void(0));
#define WALBERLA_ASSERT_NOT_NULLPTR(...)   (void(0));
#define WALBERLA_ASSERT_EQUAL(...)         (void(0));
#define WALBERLA_ASSERT_UNEQUAL(...)       (void(0));
#define WALBERLA_ASSERT_FLOAT_EQUAL(...)   (void(0));
#define WALBERLA_ASSERT_FLOAT_UNEQUAL(...) (void(0));
#define WALBERLA_ASSERT_IDENTICAL(...)     (void(0));
#define WALBERLA_ASSERT_NOT_IDENTICAL(...) (void(0));
#define WALBERLA_ASSERT_LESS(...)          (void(0));
#define WALBERLA_ASSERT_GREATER(...)       (void(0));
#define WALBERLA_ASSERT_LESS_EQUAL(...)    (void(0));
#define WALBERLA_ASSERT_GREATER_EQUAL(...) (void(0));

#endif // NDBEUG



namespace walberla {
namespace debug {

//======================================================================================================================
//
//  ASSERTS SECTION
//
//======================================================================================================================

void myAssert(const char * const file, const int line);

#ifndef NDEBUG
/// \cond internal
class ConditionalExec
{
public:
   ConditionalExec( bool cond, const std::function<void (void)> & function ) : cond_(cond), function_(function) { }
   ~ConditionalExec() { if(cond_) function_(); }
   operator bool() const { return cond_; }
private:
   bool cond_;
   std::function<void (void)> function_;
};
/// \endcond

#define WALBERLA_ASSERT_SECTION(condition) if(true)\
   if(const walberla::debug::ConditionalExec COND_EXEC = walberla::debug::ConditionalExec(!(condition),std::bind(walberla::debug::myAssert,__FILE__,__LINE__)))

#else

#define WALBERLA_ASSERT_SECTION(condition) if(false)

#endif


//======================================================================================================================
//
//  DEBUG Section
//
//======================================================================================================================

#ifndef NDEBUG
#     define WALBERLA_DEBUG_SECTION() if(true)
#else
#     define WALBERLA_DEBUG_SECTION() if(false)
#endif



} // namespace debug
} // namespace walberla
