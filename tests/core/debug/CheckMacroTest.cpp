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
//! \file CheckMacroTest.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implements Tests for the WALBERLA_CHECK_ macros
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>


// Own data structure with comparison operators but without operator<<
struct MyStructNoStream
{
   MyStructNoStream(const std::string & s) : s_(s) { }
   bool operator==(const MyStructNoStream & other) const { return this->s_ == other.s_; }
   std::string s_;
};

// Own data structure with comparison operators and free function operator<<
struct MyStructFreeStream
{
   MyStructFreeStream(const std::string & s) : s_(s) { }
   bool operator==(const MyStructFreeStream & other) const { return this->s_ == other.s_; }
   std::string s_;
};

inline std::ostream & operator<<(std::ostream & os, const MyStructFreeStream & mystruct) { os << mystruct.s_; return os; }

int main()
{
   using walberla::shared_ptr;
   walberla::debug::enterTestMode();

   {
      // Pointers
      int* nullPointer = nullptr;
      shared_ptr<int> sharedNullPointer;
      WALBERLA_CHECK_NULLPTR(nullPointer);
      WALBERLA_CHECK_NULLPTR(sharedNullPointer);

      int x = 0;
      int* notNullPointer = &x;
      shared_ptr<int> sharedNotNullPointer(new int(0));
      WALBERLA_CHECK_NOT_NULLPTR(notNullPointer);
      WALBERLA_CHECK_NOT_NULLPTR(sharedNotNullPointer);

      int a = 0;
      int b = 0;
      int * pa  = &a;
      int * pac = &a;
      int * pb  = &b;
      WALBERLA_CHECK_EQUAL(pa, pac);
      WALBERLA_CHECK_UNEQUAL(pa, pb);
   }

   {
      // Struct without stream operator
      MyStructNoStream a("Hallo Welt!");
      MyStructNoStream b("Hallo Welt!");
      WALBERLA_CHECK_EQUAL(a, b);
   }

   {
      // Struct with free stream operator
      MyStructFreeStream a("Hallo Welt!");
      MyStructFreeStream b("Hallo Welt!");
      WALBERLA_CHECK_EQUAL(a, b);
   }

   {
      std::vector<int> a( 5, 0 );
      WALBERLA_CHECK_EQUAL( a.begin(), a.begin() );
      WALBERLA_CHECK_EQUAL( std::find( a.begin(), a.end(), 1 ), a.end() );
      WALBERLA_CHECK_UNEQUAL( std::find( a.begin(), a.end(), 0 ), a.end() );
   }

   {
      std::map<int, int> a;
      for( int i = 0; i < 10; ++i )
         a[i] = i * 10;

      WALBERLA_CHECK_EQUAL( a.begin(), a.begin() );
      WALBERLA_CHECK_EQUAL( a.find( 11 ), a.end() );
      WALBERLA_CHECK_UNEQUAL( a.find( 5 ), a.end() );
      WALBERLA_CHECK_EQUAL( a.find( 5 )->second, 50 );
   }

   //WALBERLA_CHECK_(X)
   WALBERLA_CHECK(true);

   //WALBERLA_CHECK__EQUAL(X,Y)
   WALBERLA_CHECK_EQUAL(true,    true);
   WALBERLA_CHECK_EQUAL(false,   false);
   WALBERLA_CHECK_UNEQUAL(false, true);

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<char>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<unsigned char>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned char>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<short>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<unsigned short>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned short>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<int>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<unsigned int>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned int>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<long>(23));

   WALBERLA_CHECK_EQUAL(static_cast<char>(23),           static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<short>(23),          static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<int>(23),            static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<long>(23),           static_cast<unsigned long>(23));
   WALBERLA_CHECK_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned long>(23));


   //WALBERLA_CHECK__UNEQUAL(X,Y)

   WALBERLA_CHECK_UNEQUAL(true,                            false);
   WALBERLA_CHECK_UNEQUAL(false,                           true);

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<char>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<unsigned char>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<unsigned char>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<short>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<unsigned short>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<unsigned short>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<int>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<unsigned int>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<unsigned int>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<long>(11));

   WALBERLA_CHECK_UNEQUAL(static_cast<char>(23),           static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned char>(23),  static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<short>(23),          static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned short>(23), static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<int>(23),            static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned int>(23),   static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<long>(23),           static_cast<unsigned long>(11));
   WALBERLA_CHECK_UNEQUAL(static_cast<unsigned long>(23),  static_cast<unsigned long>(11));

   //WALBERLA_CHECK__GREATER_EQUAL(X,Y)

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<char>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned char>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<short>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned short>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<int>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned int>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<long>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned long>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<float>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<float>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<double>(11));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<double>(11));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<char>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned char>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned char>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<short>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned short>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned short>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<int>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned int>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned int>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<long>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<unsigned long>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<unsigned long>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<float>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<float>(23));

   WALBERLA_CHECK_GREATER_EQUAL(static_cast<char>(23),           static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned char>(23),  static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<short>(23),          static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned short>(23), static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<int>(23),            static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned int>(23),   static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<long>(23),           static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<unsigned long>(23),  static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<float>(23),          static_cast<double>(23));
   WALBERLA_CHECK_GREATER_EQUAL(static_cast<double>(23),         static_cast<double>(23));


   //WALBERLA_CHECK__LESS_EQUAL


   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<char>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned char>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<short>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned short>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<int>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned int>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<long>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned long>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<float>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<float>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<double>(33));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<double>(33));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<char>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned char>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned char>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<short>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned short>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned short>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<int>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned int>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned int>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<long>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<unsigned long>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<unsigned long>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<float>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<float>(23));

   WALBERLA_CHECK_LESS_EQUAL(static_cast<char>(23),           static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned char>(23),  static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<short>(23),          static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned short>(23), static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<int>(23),            static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned int>(23),   static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<long>(23),           static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<unsigned long>(23),  static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<float>(23),          static_cast<double>(23));
   WALBERLA_CHECK_LESS_EQUAL(static_cast<double>(23),         static_cast<double>(23));


   //WALBERLA_CHECK__GREATER(X,Y)


   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<char>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<char>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<unsigned char>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<unsigned char>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<short>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<short>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<unsigned short>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<unsigned short>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<int>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<int>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<unsigned int>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<unsigned int>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<long>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<long>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<unsigned long>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<unsigned long>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<float>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<float>(11));

   WALBERLA_CHECK_GREATER(static_cast<char>(23),           static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned char>(23),  static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<short>(23),          static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned short>(23), static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<int>(23),            static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned int>(23),   static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<long>(23),           static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<unsigned long>(23),  static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<float>(23),          static_cast<double>(11));
   WALBERLA_CHECK_GREATER(static_cast<double>(23),         static_cast<double>(11));


   //WALBERLA_CHECK__LESS


   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<char>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<char>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<unsigned char>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<unsigned char>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<short>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<short>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<unsigned short>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<unsigned short>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<int>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<int>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<unsigned int>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<unsigned int>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<long>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<long>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<unsigned long>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<unsigned long>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<float>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<float>(33));

   WALBERLA_CHECK_LESS(static_cast<char>(23),           static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned char>(23),  static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<short>(23),          static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned short>(23), static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<int>(23),            static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned int>(23),   static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<long>(23),           static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<unsigned long>(23),  static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<float>(23),          static_cast<double>(33));
   WALBERLA_CHECK_LESS(static_cast<double>(23),         static_cast<double>(33));


   // WALBERLA_CHECK__FLOAT_EQUAL


   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< float >( 23 ),       static_cast< float >( 23 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< double >( 23 ),      static_cast< float >( 23 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< long double >( 23 ), static_cast< float >( 23 ) );

   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< float >( 23 ),       static_cast< double >( 23 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< double >( 23 ),      static_cast< double >( 23 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< long double >( 23 ), static_cast< double >( 23 ) );

   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< float >( 23 ),       static_cast< long double >( 23 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< double >( 23 ),      static_cast< long double >( 23 ) );
   WALBERLA_CHECK_FLOAT_EQUAL( static_cast< long double >( 23 ), static_cast< long double >( 23 ) );


   // WALBERLA_CHECK__FLOAT_UNEQUAL


   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< float >( 23 ),       static_cast< float >( 33 ) );
   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< double >( 23 ),      static_cast< float >( 33 ) );
   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< long double >( 23 ), static_cast< float >( 33 ) );

   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< float >( 23 ),       static_cast< double >( 33 ) );
   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< double >( 23 ),      static_cast< double >( 33 ) );
   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< long double >( 23 ), static_cast< double >( 33 ) );

   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< float >( 23 ),       static_cast< long double >( 33 ) );
   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< double >( 23 ),      static_cast< long double >( 33 ) );
   WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< long double >( 23 ), static_cast< long double >( 33 ) );

   bool runtimeErrorThrown = false;
   try
   {
      walberla::Abort::instance()->resetAbortFunction( &walberla::Abort::exceptionAbort );
      WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< float >( 0 ), 0.0 ); // This check is expected to fail
      // Should never be reached
      walberla::Abort::instance()->resetAbortFunction();
      WALBERLA_CHECK(false);
   }
   catch( const std::runtime_error & /*e*/ )
   {
      walberla::Abort::instance()->resetAbortFunction();
      runtimeErrorThrown = true;
   }
   WALBERLA_CHECK( runtimeErrorThrown );

   runtimeErrorThrown = false;
   try
   {
      walberla::Abort::instance()->resetAbortFunction( &walberla::Abort::exceptionAbort );
      WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< double >( 0 ), 0.0 ); // This check is expected to fail
      // Should never be reached
      walberla::Abort::instance()->resetAbortFunction();
      WALBERLA_CHECK(false);
   }
   catch( const std::runtime_error & /*e*/ )
   {
      walberla::Abort::instance()->resetAbortFunction();
      runtimeErrorThrown = true;
   }
   WALBERLA_CHECK( runtimeErrorThrown );

   runtimeErrorThrown = false;
   try
   {
      walberla::Abort::instance()->resetAbortFunction( &walberla::Abort::exceptionAbort );
      WALBERLA_CHECK_FLOAT_UNEQUAL( static_cast< long double >( 0 ), 0.0 ); // This check is expected to fail
      // Should never be reached
      walberla::Abort::instance()->resetAbortFunction();
      WALBERLA_CHECK(false);
   }
   catch( const std::runtime_error & /*e*/ )
   {
      walberla::Abort::instance()->resetAbortFunction();
      runtimeErrorThrown = true;
   }
   WALBERLA_CHECK( runtimeErrorThrown );

   // WALBERLA_CHECK__FLOAT_IDENTICAL


   WALBERLA_CHECK_IDENTICAL( static_cast< float >( 23 ),       static_cast< float >( 23 ) );
   WALBERLA_CHECK_IDENTICAL( static_cast< double >( 23 ),      static_cast< float >( 23 ) );
   WALBERLA_CHECK_IDENTICAL( static_cast< long double >( 23 ), static_cast< float >( 23 ) );

   WALBERLA_CHECK_IDENTICAL( static_cast< float >( 23 ),       static_cast< double >( 23 ) );
   WALBERLA_CHECK_IDENTICAL( static_cast< double >( 23 ),      static_cast< double >( 23 ) );
   WALBERLA_CHECK_IDENTICAL( static_cast< long double >( 23 ), static_cast< double >( 23 ) );

   WALBERLA_CHECK_IDENTICAL( static_cast< float >( 23 ),       static_cast< long double >( 23 ) );
   WALBERLA_CHECK_IDENTICAL( static_cast< double >( 23 ),      static_cast< long double >( 23 ) );
   WALBERLA_CHECK_IDENTICAL( static_cast< long double >( 23 ), static_cast< long double >( 23 ) );


   // WALBERLA_CHECK__FLOAT_NOT_IDENTICAL


   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< float >( 23 ),       static_cast< float >( 33 ) );
   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< double >( 23 ),      static_cast< float >( 33 ) );
   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< long double >( 23 ), static_cast< float >( 33 ) );

   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< float >( 23 ),       static_cast< double >( 33 ) );
   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< double >( 23 ),      static_cast< double >( 33 ) );
   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< long double >( 23 ), static_cast< double >( 33 ) );

   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< float >( 23 ),       static_cast< long double >( 33 ) );
   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< double >( 23 ),      static_cast< long double >( 33 ) );
   WALBERLA_CHECK_NOT_IDENTICAL( static_cast< long double >( 23 ), static_cast< long double >( 33 ) );

}
