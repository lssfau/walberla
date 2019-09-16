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
//! \file SQLite.h
//! \ingroup postprocessing
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/FPClassify.h"
#include "core/logging/Logging.h"
#include "core/timing/TimingPool.h"
#include "core/timing/TimingTree.h"
#include "core/timing/TimingNode.h"
#include "core/uid/GlobalState.h"

#include <map>
#include <string>

struct sqlite3;

namespace walberla {
namespace sqlite {

using std::string;
using std::map;


// Time in seconds the process waits if database is locked
const int BUSY_TIMEOUT = 30;



class SQLiteDB
{
private:

public:
   SQLiteDB ( const string & dbFile, const int busyTimeout = BUSY_TIMEOUT );
   ~SQLiteDB();

   uint_t storeRun ( const map<string, int>     & integerProperties,
                     const map<string, string > & stringProperties ,
                     const map<string, double > & realProperties);

   uint_t storeRun ( const map<string, int64_t> & integerProperties,
                     const map<string, string > & stringProperties ,
                     const map<string, double > & realProperties);

   void storeAdditionalRunInfo( uint_t runId, const std::string & tableName,
                                const map<string, int>     & integerProperties,
                                const map<string, string > & stringProperties ,
                                const map<string, double > & realProperties );

   void storeAdditionalRunInfo( uint_t runId, const std::string & tableName,
                                const map<string, int64_t> & integerProperties,
                                const map<string, string > & stringProperties ,
                                const map<string, double > & realProperties );

   void storeTimingPool ( uint_t runId, const WcTimingPool & tp, const std::string & name  );
   void storeTimingTree ( uint_t runId, const WcTimingTree & tt, const std::string & timingTreeName  );

private:
   void storeTimingNode ( const uint_t runId,
                          const int    parentId,
                          const WcTimingNode & tn,
                          const std::string & timingTreeName,
                          const std::string & sweep,
                          const double totalTime );

   bool valid_;
   sqlite3 * dbHandle_;
   std::string file_;
};

uint_t storeRunInSqliteDB( const string               & dbFile,
                           const map<string, int>     & integerProperties= map<string,int>(),
                           const map<string, string > & stringProperties = map<string,string>(),
                           const map<string, double > & realProperties   = map<string,double>(),
                           const int                    busyTimeout      = BUSY_TIMEOUT );

uint_t storeRunInSqliteDB( const string               & dbFile,
                           const map<string, int64_t> & integerProperties= map<string,int64_t>(),
                           const map<string, string > & stringProperties = map<string,string>(),
                           const map<string, double > & realProperties   = map<string,double>(),
                           const int                    busyTimeout      = BUSY_TIMEOUT );

void storeAdditionalRunInfoInSqliteDB( const uint_t                 runId,
                                       const string               & dbFile,
                                       const string               & tableName,
                                       const map<string, int>     & integerProperties,
                                       const map<string, string > & stringProperties,
                                       const map<string, double > & realProperties,
                                       const int                    busyTimeout = BUSY_TIMEOUT);

void storeAdditionalRunInfoInSqliteDB( const uint_t                 runId,
                                       const string               & dbFile,
                                       const string               & tableName,
                                       const map<string, int64_t> & integerProperties,
                                       const map<string, string > & stringProperties,
                                       const map<string, double > & realProperties,
                                       const int                    busyTimeout = BUSY_TIMEOUT);

void storeTimingPoolInSqliteDB( const string & dbFile, uint_t runId, const WcTimingPool & tp,
                                const std::string & name, const int busyTimeout = BUSY_TIMEOUT );

void storeTimingTreeInSqliteDB( const string & dbFile, uint_t runId, const WcTimingTree & tt,
                                const std::string & name, const int busyTimeout = BUSY_TIMEOUT );





} // namespace sqlite
} // namespace walberla


