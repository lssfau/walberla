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
//! \file SQLiteTest.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================


#include "core/debug/TestSubsystem.h"
#include "core/mpi/Environment.h"

#include "sqlite/SQLite.h"


int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment walberlaEnv( argc, argv );

   std::map<std::string, std::string>       strColumns;
   std::map<std::string, int>               intColumns;
   std::map<std::string, walberla::int64_t> largeColumns;

   for( int i=0; i< 100; ++i )
   {
      strColumns["property1"] = "value1";
      strColumns["property2"] = "value2";
      strColumns["property3"] = "value3";
      intColumns["i"] = int(i);
      walberla::sqlite::storeRunInSqliteDB( "dbFile.sqlite", intColumns, strColumns );
   }

   for( int i=0; i< 100; ++i )
   {
      strColumns["property1"] = "value1";
      strColumns["property2"] = "value2";
      strColumns["property3"] = "value3";
      largeColumns["i"] = 4294967297 + i;
      walberla::sqlite::storeRunInSqliteDB( "dbFile.sqlite", largeColumns, strColumns );
   }

   return 0;
}
