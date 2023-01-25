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
//! \file BufferTest.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Simple read and write testing of SendBuffer and RecvBuffer
//!
//! This test case writes values of different types to a SendBuffer
//! then simulates the communication with a memcpy and finally
//! extracts the values from the RecvBuffer and compares them
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/cell/CellInterval.h"
#include "core/cell/CellSet.h"
#include "core/cell/CellVector.h"
#include "core/debug/TestSubsystem.h"
#include "core/math/Matrix3.h"
#include "core/math/Vector3.h"
#include "core/mpi/BufferDataTypeExtensions.h"
#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"

#include <random>
#include <cstring>
#include <iostream>


using namespace walberla;
using namespace mpi;

template<typename T>
void initIntegerContainer( T & container )
{
   static std::mt19937 rng;
   std::uniform_int_distribution<typename T::value_type> dist;
   std::uniform_int_distribution<size_t> size_dist(10,1000);

   size_t size = size_dist(rng);
   container.clear();
   for(size_t i = 0; i < size; ++i )
      container.push_back( dist(rng) );
}

template<typename T>
void initIntegerAssocContainer( T & container )
{
   static std::mt19937 rng;
   std::uniform_int_distribution<typename T::value_type> dist;
   std::uniform_int_distribution<size_t> size_dist(10,1000);

   size_t size = size_dist(rng);
   container.clear();
   for(size_t i = 0; i < size; ++i )
      container.insert( container.end(), dist(rng) );
}


void initVecBool( std::vector<bool> & vecBool )
{
   static std::mt19937 rng;
   std::uniform_int_distribution<walberla::uint32_t> dist;
   std::uniform_int_distribution<size_t> size_dist(10,1000);

   size_t size = size_dist(rng);
   vecBool.clear();
   for(size_t i = 0; i < size; ++i )
      vecBool.push_back( dist(rng) % 2 != 0);
}


template<typename T>
void initIntegerMap( T & container )
{
   static std::mt19937 rng;
   std::uniform_int_distribution<typename T::mapped_type> mapped_dist;
   std::uniform_int_distribution<typename T::key_type> key_dist;
   std::uniform_int_distribution<size_t> size_dist(10,1000);

   size_t size = size_dist(rng);
   container.clear();
   for(size_t i = 0; i < size; ++i )
      container.insert( container.end(), std::make_pair(key_dist(rng), mapped_dist(rng)) );
}


template<typename T>
void initCellContainer( T & container )
{
   static std::mt19937 rng;
   std::uniform_int_distribution<cell_idx_t> dist;
   std::uniform_int_distribution<size_t> size_dist(10,1000);

   size_t size = size_dist(rng);
   container.clear();
   for(size_t i = 0; i < size; ++i )
      container.insert( container.end(), Cell( dist(rng), dist(rng), dist(rng) ) );
}


template<typename T, std::size_t N>
void initStdArray( std::array< T, N > & array )
{
   static std::mt19937 rng;
   std::uniform_int_distribution<T> dist;

   for( auto it = array.begin(); it != array.end(); ++it )
      *it = dist( rng );
}


/// Simulates one send and receive operation
/// data is put into send buffer and copied to RecvBuffer
/// Then the values are compared to sent values
template<typename T>
void bufferTest()
{
   const double          testDouble = 42.242;
   const int             testInt = 10;
   const Vector3<int>    vec(1,2,3);
   const Matrix3<double> mat(1,2,3,4,5,6,7,8,9);
   const Cell            cell(1,2,3);
   const CellInterval    cellInterval( Cell(1,2,3), Cell(4,5,6) );

   CellVector            cellVector;
   CellSet               cellSet;

   initCellContainer(cellVector);
   initCellContainer(cellSet);

   std::vector       <bool>         boolStdVec,      boolStdVecEmpty;
   std::vector       <unsigned int> stdVec,          stdVecEmpty;
   std::deque        <unsigned int> stdDeque,        stdDequeEmpty;
   std::list         <unsigned int> stdList,         stdListEmpty;
   std::set          <unsigned int> stdSet,          stdSetEmpty;
   std::multiset     <unsigned int> stdMultiSet,     stdMultiSetEmpty;
   std::unordered_set<unsigned int> stdUnorderedSet, stdUnorderedSetEmpty;

   std::map          <unsigned int, walberla::int64_t> stdMap,          stdMapEmpty;
   std::multimap     <unsigned int, walberla::int64_t> stdMultiMap,     stdMultiMapEmpty;
   std::unordered_map<unsigned int, walberla::int64_t> stdUnorderedMap, stdUnorderedMapEmpty;

   std::array  < unsigned int, 19 > stdArray;

   initVecBool(boolStdVec);
   initIntegerContainer(stdVec);
   initIntegerContainer(stdDeque);
   initIntegerContainer(stdList);
   initIntegerAssocContainer(stdSet);
   initIntegerAssocContainer(stdMultiSet);
   initIntegerAssocContainer(stdUnorderedSet);
   initIntegerMap(stdMap);
   initIntegerMap(stdMultiMap);
   initIntegerMap(stdUnorderedMap);
   initStdArray(stdArray);

   // Create send buffer and put two values in it
   GenericSendBuffer<T> sb;
   sb << testDouble       << testInt;
   sb << vec              << mat;
   sb << cell             << cellInterval;
   sb << cellVector       << cellSet;
   sb << boolStdVec       << boolStdVecEmpty;
   sb << stdVec           << stdVecEmpty;
   sb << stdDeque         << stdDequeEmpty;
   sb << stdList          << stdListEmpty;
   sb << stdSet           << stdSetEmpty;
   sb << stdMultiSet      << stdMultiSetEmpty;
   sb << stdUnorderedSet  << stdUnorderedSetEmpty;
   sb << stdMap           << stdMapEmpty;
   sb << stdMultiMap      << stdMultiMapEmpty;
   sb << stdUnorderedMap  << stdUnorderedMapEmpty;
   sb << stdArray;

   // Copying
   //RecvBuffer<T> rb;
   //rb.resize( sb.size() );
   //memcpy(rb.ptr(),sb.ptr(), sb.size()* sizeof(T));

   // Transferring data to receive buffer
   GenericRecvBuffer<T> rb( sb );

   // Extracting Values
   double          recvD;
   int             recvI;
   Vector3<int>    recvVec;
   Matrix3<double> recvMat;
   Cell            recvCell;
   CellInterval    recvCellInterval;
   CellVector      recvCellVector;
   CellSet         recvCellSet;

   std::vector        <bool>         recvBoolStdVec,      recvBoolStdVecEmpty;
   std::vector        <unsigned int> recvStdVec,          recvStdVecEmpty;
   std::deque         <unsigned int> recvStdDeque,        recvStdDequeEmpty;
   std::list          <unsigned int> recvStdList,         recvStdListEmpty;
   std::set           <unsigned int> recvStdSet,          recvStdSetEmpty;
   std::multiset      <unsigned int> recvStdMultiSet,     recvStdMultiSetEmpty;
   std::unordered_set <unsigned int> recvStdUnorderedSet, recvStdUnorderedSetEmpty;

   std::map           <unsigned int, walberla::int64_t> recvStdMap,          recvStdMapEmpty;
   std::multimap      <unsigned int, walberla::int64_t> recvStdMultiMap,     recvStdMultiMapEmpty;
   std::unordered_map <unsigned int, walberla::int64_t> recvStdUnorderedMap, recvStdUnorderedMapEmpty;

   std::array  <unsigned int, 19> recvStdArray;

   rb >> recvD                >> recvI;
   rb >> recvVec              >> recvMat;
   rb >> recvCell             >> recvCellInterval;
   rb >> recvCellVector       >> recvCellSet;
   rb >> recvBoolStdVec       >> recvBoolStdVecEmpty;
   rb >> recvStdVec           >> recvStdVecEmpty;
   rb >> recvStdDeque         >> recvStdDequeEmpty;
   rb >> recvStdList          >> recvStdListEmpty;
   rb >> recvStdSet           >> recvStdSetEmpty;
   rb >> recvStdMultiSet      >> recvStdMultiSetEmpty;
   rb >> recvStdUnorderedSet  >> recvStdUnorderedSetEmpty;
   rb >> recvStdMap           >> recvStdMapEmpty;
   rb >> recvStdMultiMap      >> recvStdMultiMapEmpty;
   rb >> recvStdUnorderedMap  >> recvStdUnorderedMapEmpty;
   rb >> recvStdArray;

   // Validate
   WALBERLA_CHECK_FLOAT_EQUAL(recvD,testDouble);

   WALBERLA_CHECK_EQUAL(recvI,    testInt);
   WALBERLA_CHECK_EQUAL(recvVec,  vec);
   WALBERLA_CHECK_EQUAL(recvMat,  mat);
   WALBERLA_CHECK_EQUAL(recvCell, cell);
   WALBERLA_CHECK_EQUAL(recvCellInterval, cellInterval);
   WALBERLA_CHECK_EQUAL(recvCellVector,   cellVector);
   WALBERLA_CHECK(recvCellSet == cellSet);

   WALBERLA_CHECK_EQUAL(recvBoolStdVec,      boolStdVec);
   WALBERLA_CHECK_EQUAL(recvBoolStdVecEmpty, boolStdVecEmpty);
   WALBERLA_CHECK_EQUAL(recvStdVec,          stdVec);
   WALBERLA_CHECK_EQUAL(recvStdVecEmpty,     stdVecEmpty);
   WALBERLA_CHECK_EQUAL(recvStdDeque,        stdDeque);
   WALBERLA_CHECK_EQUAL(recvStdDequeEmpty,   stdDequeEmpty);
   WALBERLA_CHECK_EQUAL(recvStdList,         stdList);
   WALBERLA_CHECK_EQUAL(recvStdListEmpty,    stdListEmpty);

   WALBERLA_CHECK_EQUAL(recvStdSet,               stdSet);
   WALBERLA_CHECK_EQUAL(recvStdSetEmpty,          stdSetEmpty);
   WALBERLA_CHECK_EQUAL(recvStdMultiSet,          stdMultiSet);
   WALBERLA_CHECK_EQUAL(recvStdMultiSetEmpty,     stdMultiSetEmpty);
   WALBERLA_CHECK_EQUAL(recvStdUnorderedSet,      stdUnorderedSet);
   WALBERLA_CHECK_EQUAL(recvStdUnorderedSetEmpty, stdUnorderedSetEmpty);

   WALBERLA_CHECK_EQUAL(recvStdMap,               stdMap);
   WALBERLA_CHECK_EQUAL(recvStdMapEmpty,          stdMapEmpty);
   WALBERLA_CHECK_EQUAL(recvStdMultiMap,          stdMultiMap);
   WALBERLA_CHECK_EQUAL(recvStdMultiMapEmpty,     stdMultiMapEmpty);
   WALBERLA_CHECK_EQUAL(recvStdUnorderedMap,      stdUnorderedMap);
   WALBERLA_CHECK_EQUAL(recvStdUnorderedMapEmpty, stdUnorderedMapEmpty);

   WALBERLA_CHECK_EQUAL(recvStdArray,   stdArray);
}


// Special test function to check serializations that only work with buffers that utilize
// uint8_t as underlying data type
void bufferTestUInt8()
{
   std::string stdString("Hello World!"), stdStringEmpty;

   walberla::optional<int> optional0, optional1, optional2, optional3;
   optional2 = 23;
   optional3 = 42;

   // Create send buffer and put two values in it
   GenericSendBuffer<walberla::uint8_t> sb;
   sb << stdString << stdStringEmpty;
   sb << optional0 << optional1 << optional2 << optional3;

   // Copying
   GenericRecvBuffer<walberla::uint8_t> rb( sb );

   std::string recvStdString, recvStdStringEmpty;

   walberla::optional<int> recvOptional0 = 123;
   walberla::optional<int> recvOptional1 = 123;
   walberla::optional<int> recvOptional2 = 456;
   walberla::optional<int> recvOptional3 = 456;

   rb >> recvStdString >> recvStdStringEmpty;
   rb >> recvOptional0 >> recvOptional1;
   rb >> recvOptional2 >> recvOptional3;

   WALBERLA_CHECK_EQUAL(recvStdString,       stdString);
   WALBERLA_CHECK_EQUAL(recvStdStringEmpty,  stdStringEmpty);

   WALBERLA_CHECK( optional0 == recvOptional0 );
   WALBERLA_CHECK( optional1 == recvOptional1 );
   WALBERLA_CHECK( optional2 == recvOptional2 );
   WALBERLA_CHECK( optional3 == recvOptional3 );
}


void bufferOverwriteTest()
{
   //! [SendBuffer Overwrite Test]
   int a = 1;
   int b = 2;
   int c = 3;
   SendBuffer sb;
   auto ptr = sb.allocate<int>(2);
   sb << b << c;
   *ptr   = a;
   ptr[1] = c;
   //ptr[2] = c; // will fail in debug: out of bounds
   //! [SendBuffer Overwrite Test]

   // Copying
   RecvBuffer rb( sb );

   int recv;

   rb >> recv;
   WALBERLA_CHECK_EQUAL(a, recv);
   rb >> recv;
   WALBERLA_CHECK_EQUAL(c, recv);
   rb >> recv;
   WALBERLA_CHECK_EQUAL(b, recv);
   rb >> recv;
   WALBERLA_CHECK_EQUAL(c, recv);
}


int main()
{
   debug::enterTestMode();

   bufferTest<unsigned char>() ;
   bufferTest<unsigned short>();
   bufferTest<unsigned int>()  ;

   bufferTestUInt8();

   // This should not compile, since the containing type has be be smaller
   // than the type that is stored
   //bufferTest<double>();

   bufferOverwriteTest();

   return 0;
}

