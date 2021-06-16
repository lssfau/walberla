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
//! \file SQLite.cpp
//! \ingroup postprocessing
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "SQLite.h"

#include "sqlite3.h"

#include <core/RandomUUID.h>

#include <sstream>
#include <cassert>


namespace walberla {
namespace sqlite {


SQLiteDB::SQLiteDB( const string & dbFile, const int busyTimeout )
   : valid_(true), dbHandle_(nullptr), file_( dbFile )
{
   static const char * CREATE_RUN_TABLE =
         "CREATE TABLE IF NOT EXISTS runs ("
         " runId     INTEGER PRIMARY KEY, "
         " timestamp DATETIME DEFAULT CURRENT_TIMESTAMP, "
         " uuid      STRING );" ;

   // Create tables if it does not exist
   int retVal = sqlite3_open( file_.c_str(), &dbHandle_ );
   if ( retVal != SQLITE_OK ) {
      WALBERLA_LOG_WARNING( "Failed to open sqlite3 file." << dbFile );
      valid_ = false;
      return;
   }

   sqlite3_busy_timeout( dbHandle_, busyTimeout*1000 );
   sqlite3_exec( dbHandle_, "PRAGMA foreign_keys = ON;",nullptr,nullptr,nullptr );
   sqlite3_exec( dbHandle_, CREATE_RUN_TABLE, nullptr,nullptr,nullptr);

   static const char* UPDATE_RUN_TABLE_CMD = "ALTER TABLE runs ADD COLUMN uuid STRING;";
   sqlite3_exec( dbHandle_, UPDATE_RUN_TABLE_CMD, nullptr, nullptr, nullptr );

}


SQLiteDB::~SQLiteDB ()
{
   if ( valid_ )
      sqlite3_close( dbHandle_ );
}



//*******************************************************************************************************************
/*! Store information about a simulation run into a Sqlite3 Database
   *
   *  The generated database is called "runs". The columns of the table are the keys of the given maps.
   *  If a column does not yet exist, it is appended to the database. The entries that have no value in this column
   *  are set to default value ( zero or empty string ).
   *  Additionally a column is created for each activated global state, and its value is set to 1
   *
   *  \param *Properties  Map of column names to value
   *  \returns            The primary key of the inserted data set.
   */
//*******************************************************************************************************************
template<typename IntType>
uint_t storeRunImpl( sqlite3 * dbHandle, std::string & filename,
                     const map<string, IntType> & integerProperties,
                     const map<string, string > & stringProperties,
                     const map<string, double > & realProperties )
{
   WALBERLA_ASSERT_NOT_NULLPTR( dbHandle );

   sqlite3_exec( dbHandle, "BEGIN;",nullptr,nullptr,nullptr );

   string insertRunCommand = "INSERT INTO runs (timestamp, uuid ";
   std::stringstream values;
   auto uuid = RandomUUID();
   values << " VALUES ( CURRENT_TIMESTAMP, \"" << uuid << "\" ";
   // Add columns for integer properties
   for ( auto i = integerProperties.begin(); i != integerProperties.end(); ++i )
   {
      insertRunCommand += "," + i->first;
      values  << ", " << i->second;
      string command = "ALTER TABLE runs ADD COLUMN " + i->first + " INTEGER ";
      sqlite3_exec ( dbHandle, command.c_str(), nullptr,nullptr,nullptr ); // ignore errors (column can exist already)
   }

   // Add columns for string properties
   for ( auto i = stringProperties.begin(); i != stringProperties.end(); ++i )
   {
      insertRunCommand += "," + i->first;
      values << ", " << "\"" << i->second << "\"";
      string command = "ALTER TABLE runs ADD COLUMN " + i->first + " TEXT ";
      sqlite3_exec ( dbHandle, command.c_str(), nullptr,nullptr,nullptr ); // ignore errors (column can exist already)

   }

   // Add columns for real_t properties
   for ( auto i = realProperties.begin(); i != realProperties.end(); ++i )
   {
      if( math::finite( i->second ) )
      {
         insertRunCommand += "," + i->first;
         values << ", " << i->second;
         string command = "ALTER TABLE runs ADD COLUMN " + i->first + " DOUBLE ";
         sqlite3_exec( dbHandle, command.c_str(), nullptr, nullptr, nullptr ); // ignore errors (column can exist already)
      }
      else
      {
         WALBERLA_LOG_WARNING( "Skipping column \"" << i->first << "\" while inserting a row into the run table of the "
                                                                   "sqlite3 database \"" << filename << "\" due to non-finite value \"" << i->second << "\"." );
      }
   }

   // Add columns for global state selectors
   for( auto i = uid::globalState().begin(); i != uid::globalState().end(); ++i )
   {
      insertRunCommand += "," + i->getIdentifier();
      values << " ,1";
      // no boolean in sqlite3, use integer instead
      string command = "ALTER TABLE runs ADD COLUMN " + i->getIdentifier() + " INTEGER ";
      sqlite3_exec ( dbHandle, command.c_str(), nullptr,nullptr,nullptr ); // ignore errors (column can exist already)
   }

   insertRunCommand += " )  ";
   values << "); ";
   insertRunCommand += values.str();

   int ret = sqlite3_exec ( dbHandle, insertRunCommand.c_str(), nullptr, nullptr, nullptr );
   if ( ret != SQLITE_OK) {
      WALBERLA_LOG_WARNING( "Failed to insert a row into run table of sqlite3 database: " << sqlite3_errmsg(dbHandle) << "\n sql command: " << insertRunCommand.c_str() );
   }
   uint_t generatedPrimaryKey = uint_c ( sqlite3_last_insert_rowid( dbHandle ) );

   sqlite3_exec( dbHandle, "END TRANSACTION;",nullptr,nullptr,nullptr );

   return generatedPrimaryKey;
}

//*******************************************************************************************************************
/*! Stores information in another table, referencing the "run" table
   *
   * \param runId       result of storeRun() member function, primary key of the run to store information for
   * \param tableName   name of the table where the information is stored in
   *                    is created if it does not yet exist
   * \param *Properties Map of column names to value
   */
//*******************************************************************************************************************
template<typename IntType>
void storeAdditionalRunInfoImpl( sqlite3 * dbHandle,
                                 uint_t runId, const std::string & tableName,
                                 const map<string, IntType> & integerProperties,
                                 const map<string, string > & stringProperties ,
                                 const map<string, double > & realProperties )
{
   sqlite3_exec( dbHandle, "BEGIN;",nullptr,nullptr,nullptr );
   std::string CREATE_TABLE =
         "CREATE TABLE IF NOT EXISTS " + tableName +
         " (runId     INTEGER, "
         " FOREIGN KEY (runId) REFERENCES runs(runId) "
         " );" ;

   sqlite3_exec( dbHandle, CREATE_TABLE.c_str(), nullptr,nullptr,nullptr);

   string insertRunCommand = "INSERT INTO " + tableName + "( runId";
   std::stringstream values;
   values << " VALUES (  " << runId;
   // Add columns for integer properties
   for ( auto i = integerProperties.begin(); i != integerProperties.end(); ++i )
   {
      insertRunCommand += "," + i->first;
      values  << ", " << i->second;
      string command = "ALTER TABLE " + tableName + " ADD COLUMN " + i->first + " INTEGER ";
      sqlite3_exec ( dbHandle, command.c_str(), nullptr,nullptr,nullptr ); // ignore errors (column can exist already)
   }

   // Add columns for string properties
   for ( auto i = stringProperties.begin(); i != stringProperties.end(); ++i )
   {
      insertRunCommand += "," + i->first;
      values << ", " << "\"" << i->second << "\"";
      string command = "ALTER TABLE " + tableName + " ADD COLUMN " + i->first + " TEXT ";
      sqlite3_exec ( dbHandle, command.c_str(), nullptr,nullptr,nullptr ); // ignore errors (column can exist already)

   }

   // Add columns for real_t properties
   for ( auto i = realProperties.begin(); i != realProperties.end(); ++i )
   {
      insertRunCommand += "," + i->first;
      values << ", " << i->second ;
      string command = "ALTER TABLE " + tableName + " ADD COLUMN " + i->first + " DOUBLE ";
      sqlite3_exec ( dbHandle, command.c_str(), nullptr,nullptr,nullptr ); // ignore errors (column can exist already)
   }

   insertRunCommand += " )  ";
   values << "); ";
   insertRunCommand += values.str();

   int ret = sqlite3_exec ( dbHandle, insertRunCommand.c_str(), nullptr, nullptr, nullptr );
   if ( ret != SQLITE_OK) {
      WALBERLA_LOG_WARNING( "Failed to insert a row into run table of sqlite3 database: " << sqlite3_errmsg(dbHandle) << "\n sql command: " << insertRunCommand.c_str() );
   }
   sqlite3_exec( dbHandle, "END TRANSACTION;",nullptr,nullptr,nullptr );
}


//*******************************************************************************************************************
/*! Store information about a simulation run into a Sqlite3 Database
   *
   *  The generated database is called "runs". The columns of the table are the keys of the given maps.
   *  If a column does not yet exist, it is appended to the database. The entries that have no value in this column
   *  are set to default value ( zero or empty string ).
   *  Additionally a column is created for each activated global state, and its value is set to 1
   *
   *  \param *Properties  Map of column names to value
   *  \returns            The primary key of the inserted data set.
   */
//*******************************************************************************************************************
uint_t SQLiteDB::storeRun( const map<string, int>     & integerProperties,
                           const map<string, string > & stringProperties,
                           const map<string, double > & realProperties )
{
   return storeRunImpl(dbHandle_, file_, integerProperties, stringProperties, realProperties);
}
/// \see storeRun
uint_t SQLiteDB::storeRun( const map<string, int64_t> & integerProperties,
                           const map<string, string > & stringProperties,
                           const map<string, double > & realProperties )
{
   return storeRunImpl(dbHandle_, file_, integerProperties, stringProperties, realProperties);
}



//*******************************************************************************************************************
/*! Stores information in another table, referencing the "run" table
   *
   * \param runId       result of storeRun() member function, primary key of the run to store information for
   * \param tableName   name of the table where the information is stored in
   *                    is created if it does not yet exist
   * \param *Properties Map of column names to value
   */
//*******************************************************************************************************************
void SQLiteDB::storeAdditionalRunInfo( uint_t runId, const std::string & tableName,
                                       const map<string, int>     & integerProperties,
                                       const map<string, string > & stringProperties ,
                                       const map<string, double > & realProperties )
{
   return storeAdditionalRunInfoImpl( dbHandle_, runId, tableName, integerProperties, stringProperties, realProperties);
}
/// \see storeAdditionalRunInfo
void SQLiteDB::storeAdditionalRunInfo( uint_t runId, const std::string & tableName,
                                       const map<string, int64_t> & integerProperties,
                                       const map<string, string > & stringProperties ,
                                       const map<string, double > & realProperties )
{
   return storeAdditionalRunInfoImpl( dbHandle_, runId, tableName, integerProperties, stringProperties, realProperties);
}



//*******************************************************************************************************************
/*! Stores a TimingPool in a Sqlite3 Database, and links it to a run
   *
   * The generated table is called "timingPools"
   *
   * \param runId   primary key of the run, as returned by storeRun()
   * \param tp      the TimingPool to store
   * \param name    name of the timing pool ( as written to database column )
   */
//*******************************************************************************************************************
void SQLiteDB::storeTimingPool ( uint_t runId,
                                 const WcTimingPool & tp,
                                 const std::string & timingPoolName )
{
   sqlite3_exec( dbHandle_, "BEGIN;",nullptr,nullptr,nullptr );

   assert ( timingPoolName.length() > 0 && timingPoolName.length() < 255 );

   static const char * CREATE_TIMINGPOOL_TABLE =
         "CREATE TABLE IF NOT EXISTS timingPool ("
         " runId       INTEGER, "
         " name        VARCHAR(255),"
         " sweep       VARCHAR(255),"
         " average     DOUBLE, "
         " min         DOUBLE, "
         " max         DOUBLE, "
         " count       INTEGER,"
         " variance    DOUBLE, "
         " percentage  DOUBLE, "
         " FOREIGN KEY (runId) REFERENCES runs(runId)  "
         " );" ;


   static const char * INSERT_STATEMENT =
         " INSERT INTO timingPool (runId,name,sweep,average,min,max,count,variance,percentage) "
         " VALUES( ?, ?, ?, ?, ?, ?, ?, ?, ? )";

   sqlite3_exec( dbHandle_, CREATE_TIMINGPOOL_TABLE, nullptr,nullptr,nullptr );

   sqlite3_stmt *stmt = nullptr;
   auto retVal = sqlite3_prepare_v2(dbHandle_, INSERT_STATEMENT, -1, &stmt, nullptr );
   if ( retVal != SQLITE_OK ) {
      WALBERLA_LOG_WARNING( "Failed to prepare SQL Insert statement." << file_ );
      return;
   }


   double totalTime = 0;
   for ( auto i = tp.begin(); i != tp.end(); ++i )
      totalTime += i->second.total();


   for ( auto i = tp.begin(); i != tp.end(); ++i )
   {
      sqlite3_bind_int64 ( stmt, 1, int64_c(runId) );
      sqlite3_bind_text  ( stmt, 2, timingPoolName.c_str() , -1, SQLITE_STATIC );
      sqlite3_bind_text  ( stmt, 3, i->first.c_str() , -1, SQLITE_STATIC );
      sqlite3_bind_double( stmt, 4, ( ( i->second.getCounter() == uint_t(0) ) ? 0.0 : double_c( i->second.average() ) ) );
      sqlite3_bind_double( stmt, 5, ( ( i->second.getCounter() == uint_t(0) ) ? 0.0 : double_c( i->second.min() ) ) );
      sqlite3_bind_double( stmt, 6, ( ( i->second.getCounter() == uint_t(0) ) ? 0.0 : double_c( i->second.max() ) ) );
      sqlite3_bind_int64 ( stmt, 7, int64_c   ( i->second.getCounter() ));
      sqlite3_bind_double( stmt, 8, ( ( i->second.getCounter() == uint_t(0) ) ? 0.0 : double_c( i->second.variance() ) ) );
      sqlite3_bind_double( stmt, 9, ( ( i->second.getCounter() == uint_t(0) ) ? 0.0 : double_c( i->second.total() / totalTime ) ) );

      sqlite3_step ( stmt );  // execute statement
      sqlite3_reset ( stmt ); // undo binding
   }
   sqlite3_exec( dbHandle_, "END TRANSACTION;",nullptr,nullptr,nullptr );

   sqlite3_finalize( stmt ); // free prepared statement
}

//*******************************************************************************************************************
/*! Stores a TimingTree in a Sqlite3 Database, and links it to a run
   *
   * The generated table is called "timingTree"
   *
   * \param runId   primary key of the run, as returned by storeRun()
   * \param tt      the TimingTree to store
   * \param name    name of the timing tree ( as written to database column )
   */
//*******************************************************************************************************************
void SQLiteDB::storeTimingTree ( uint_t runId,
                                 const WcTimingTree & tt,
                                 const std::string & timingTreeName )
{
   sqlite3_exec( dbHandle_, "BEGIN;",nullptr,nullptr,nullptr );

   assert ( timingTreeName.length() > 0 && timingTreeName.length() < 255 );

   static const char * CREATE_TIMINGTREE_TABLE =
         "CREATE TABLE IF NOT EXISTS timingTree ("
         " id          INTEGER PRIMARY KEY, "
         " runId       INTEGER, "
         " name        VARCHAR(255),"
         " parentId    INTEGER, "
         " sweep       VARCHAR(255),"
         " average     DOUBLE, "
         " min         DOUBLE, "
         " max         DOUBLE, "
         " count       INTEGER,"
         " variance    DOUBLE, "
         " percentage  DOUBLE, "
         " FOREIGN KEY (runId) REFERENCES runs(runId)  "
         " );" ;

   sqlite3_exec( dbHandle_, CREATE_TIMINGTREE_TABLE, nullptr,nullptr,nullptr );

   double totalTime = 0.0;
   for (auto it = tt.getRawData().tree_.begin(); it != tt.getRawData().tree_.end(); ++it)
   {
      totalTime += it->second.timer_.total();
   }

   storeTimingNode(runId, std::numeric_limits<int>::max(), tt.getRawData(), timingTreeName, "Total", totalTime);

   sqlite3_exec( dbHandle_, "END TRANSACTION;",nullptr,nullptr,nullptr );
}

//*******************************************************************************************************************
/*! Stores a TimingNode recursively in a Sqlite3 Database, and links it together
   *
   * \param runId   primary key of the run, as returned by storeRun()
   * \param parentId   parent key of the node
   * \param tn      the TimingNode to store
   * \param name    name of the timing tree ( as written to database column )
   */
//*******************************************************************************************************************
void SQLiteDB::storeTimingNode ( const uint_t runId,
                                 const int    parentId,
                                 const WcTimingNode & tn,
                                 const std::string & timingTreeName,
                                 const std::string & sweep,
                                 const double totalTime )
{
   static const char * INSERT_STATEMENT =
         " INSERT INTO timingTree (runId,name,parentId,sweep,average,min,max,count,variance,percentage) "
         " VALUES( ?, ?, ?, ?, ?, ?, ?, ?, ?, ? )";

   sqlite3_stmt *stmt = nullptr;
   auto retVal = sqlite3_prepare_v2(dbHandle_, INSERT_STATEMENT, -1, &stmt, nullptr );
   if ( retVal != SQLITE_OK ) {
      WALBERLA_LOG_WARNING( "Failed to prepare SQL Insert statement (" << retVal << ")." );
      return;
   }

   sqlite3_bind_int64 ( stmt, 1, int64_c(runId) );
   sqlite3_bind_text  ( stmt, 2, timingTreeName.c_str() , -1, SQLITE_STATIC );
   sqlite3_bind_int64 ( stmt, 3, parentId );
   sqlite3_bind_text  ( stmt, 4, sweep.c_str() , -1, SQLITE_STATIC );
   sqlite3_bind_double( stmt, 5, ( ( tn.timer_.getCounter() == uint_t(0) ) ? 0.0 : double_c( tn.timer_.average() ) ) );
   sqlite3_bind_double( stmt, 6, ( ( tn.timer_.getCounter() == uint_t(0) ) ? 0.0 : double_c( tn.timer_.min() ) ) );
   sqlite3_bind_double( stmt, 7, ( ( tn.timer_.getCounter() == uint_t(0) ) ? 0.0 : double_c( tn.timer_.max() ) ) );
   sqlite3_bind_int64 ( stmt, 8, int64_c   ( tn.timer_.getCounter() ));
   sqlite3_bind_double( stmt, 9, ( ( tn.timer_.getCounter() == uint_t(0) ) ? 0.0 : double_c( tn.timer_.variance() ) ) );
   sqlite3_bind_double( stmt,10, ( ( tn.timer_.getCounter() == uint_t(0) ) ? 0.0 : double_c( tn.timer_.total() / totalTime ) ) );

   sqlite3_step ( stmt );  // execute statement
   sqlite3_reset ( stmt ); // undo binding
   sqlite3_finalize( stmt ); // free prepared statement

   int currentId = int_c( sqlite3_last_insert_rowid( dbHandle_ ) );

   for ( auto i = tn.tree_.begin(); i != tn.tree_.end(); ++i )
   {
      storeTimingNode( runId, currentId, i->second, timingTreeName, i->first, totalTime);
   }
}





uint_t storeRunInSqliteDB( const string               & dbFile,
                           const map<string, int>     & integerProperties,
                           const map<string, string > & stringProperties,
                           const map<string, double > & realProperties,
                           const int                    busyTimeout )
{
   SQLiteDB db ( dbFile, busyTimeout );
   return db.storeRun( integerProperties, stringProperties, realProperties );
}
uint_t storeRunInSqliteDB( const string               & dbFile,
                           const map<string, int64_t> & integerProperties,
                           const map<string, string > & stringProperties,
                           const map<string, double > & realProperties,
                           const int                    busyTimeout )
{
   SQLiteDB db ( dbFile, busyTimeout );
   return db.storeRun( integerProperties, stringProperties, realProperties );
}


void storeAdditionalRunInfoInSqliteDB( const uint_t                 runId,
                                       const string               & dbFile,
                                       const string               & tableName,
                                       const map<string, int>     & integerProperties,
                                       const map<string, string > & stringProperties,
                                       const map<string, double > & realProperties,
                                       const int                    busyTimeout )
{
   SQLiteDB db ( dbFile, busyTimeout );
   return db.storeAdditionalRunInfo( runId, tableName, integerProperties, stringProperties, realProperties );
}
void storeAdditionalRunInfoInSqliteDB( const uint_t                 runId,
                                       const string               & dbFile,
                                       const string               & tableName,
                                       const map<string, int64_t> & integerProperties,
                                       const map<string, string > & stringProperties,
                                       const map<string, double > & realProperties,
                                       const int                    busyTimeout )
{
   SQLiteDB db ( dbFile, busyTimeout );
   return db.storeAdditionalRunInfo( runId, tableName, integerProperties, stringProperties, realProperties );
}


void storeTimingPoolInSqliteDB ( const string & dbFile, uint_t runId,
                                 const WcTimingPool & tp,
                                 const std::string & timingPoolName,
                                 const int           busyTimeout )
{
   SQLiteDB db ( dbFile, busyTimeout );
   db.storeTimingPool( runId, tp, timingPoolName );
}

void storeTimingTreeInSqliteDB ( const string & dbFile, uint_t runId,
                                 const WcTimingTree & tt,
                                 const std::string & timingTreeName,
                                 const int           busyTimeout )
{
   SQLiteDB db ( dbFile, busyTimeout );
   db.storeTimingTree( runId, tt, timingTreeName );
}

} // namespace sqlite
} // namespace walberla
