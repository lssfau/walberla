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
//! \file Operator.h
//! \ingroup core
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "FwdOperator.h"
#include "core/Abort.h"

#include <cmath>
#include <string>


namespace walberla {
namespace math {

   class OpType
   {
   private:
      const char         sign_;
      const std::string  name_;
      const unsigned int strength_;

   public:
      OpType( const char& sign, const std::string& n, const unsigned int strength ) :
         sign_(sign), name_(n), strength_(strength) {}

      virtual ~OpType() = default;

   private:
      OpType& operator=( const OpType& ){ return *this; }

   public:
      bool operator==( const OpType & type ) const { return sign_ == type.sign_; }
      bool operator==( const char   & c    ) const { return sign_ == c;          }

      bool operator<( const OpType & type ) const { return strength_ < type.strength_; }
      bool operator>( const OpType & type ) const { return strength_ > type.strength_; }
      bool operator<=( const OpType & type ) const { return strength_ <= type.strength_; }
      bool operator>=( const OpType & type ) const { return strength_ >= type.strength_; }

      virtual double operator() ( const double&, const double& ) = 0;

      friend std::ostream& operator<<( std::ostream& os, const OpType & type );

   public:
      const std::string & getName() const { return name_; }
   };

   class OpNo : public OpType{
   public:
      OpNo( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double &, const double & ) override { WALBERLA_ABORT( "NO OPERATION" ); return 0; }
   };

   class OpPlus : public OpType{
   public:
      OpPlus( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {};
      double operator() ( const double & a, const double & b ) override { return a + b; }
   };

   class OpMinus : public OpType{
   public:
      OpMinus( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double & a, const double & b ) override { return a - b; }
   };

   class OpMult : public OpType{
   public:
      OpMult( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double & a, const double & b ) override { return a * b; }
   };

   class OpDiv : public OpType{
   public:
      OpDiv( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double & a, const double & b ) override { return a / b; }
   };

   class OpProd : public OpType{
   public:
      OpProd( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double & a, const double & b ) override { return pow( a, b ); }
   };

   class OpRoot : public OpType{
   public:
      OpRoot( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double & a, const double & b ) override { return pow( a, 1/b ); }
   };

   class OpLog : public OpType{
   public:
      OpLog( const char& sign, const std::string& name, const unsigned int strength ) :
         OpType( sign, name, strength ) {}
      double operator() ( const double & a, const double & b ) override { return log10(a) / log10(b); }
   };


   int isop( const char c );

   OpType& getOp ( const char c );

} // namespace math
} // namespace walberla
