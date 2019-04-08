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
//! \file VoxelFileTest.cpp
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================


#include "geometry/structured/VoxelFileReader.h"
#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include "core/Filesystem.h"

#include <random>

typedef std::mersenne_twister_engine< walberla::uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 0xa37d3c92 > mt11213b;

/// randomize the memory underlying the vector up the maximal size (== capacity)
template<typename T>
void randomizeVector( std::vector<T> & v )
{
   static_assert(sizeof(T) > sizeof(char), "cannot use char");

   mt11213b rng;
   std::uniform_int_distribution<T> dist( std::numeric_limits<T>::min(), std::numeric_limits<T>::max() );

   size_t oldSize = v.size();
   v.resize( v.capacity() );
   for(typename std::vector<T>::iterator it = v.begin(); it != v.end(); ++it)
      *it = dist(rng);
   v.resize(oldSize);
}

template<>
void randomizeVector( std::vector<unsigned char> & v )
{
   mt11213b rng;
   std::uniform_int_distribution<walberla::int16_t> dist( std::numeric_limits<unsigned char>::min(), std::numeric_limits<unsigned char>::max() );

   size_t oldSize = v.size();
   v.resize( v.capacity() );
   for(typename std::vector<unsigned char>::iterator it = v.begin(); it != v.end(); ++it)
      *it = static_cast<unsigned char>( dist(rng) );
   v.resize(oldSize);
}

template<>
void randomizeVector( std::vector<char> & v )
{
   mt11213b rng;
   std::uniform_int_distribution<int16_t> dist( std::numeric_limits<char>::min(), std::numeric_limits<char>::max() );

   size_t oldSize = v.size();
   v.resize( v.capacity() );
   for(typename std::vector<char>::iterator it = v.begin(); it != v.end(); ++it)
      *it = static_cast<char>( dist(rng) );
   v.resize(oldSize);
}

template<typename T>
void makeRandomMultiArray( std::vector<T> & ma)
{
   static_assert(sizeof(T) > sizeof(char), "cannot use char");

   mt11213b rng;
   std::uniform_int_distribution<T> dist( std::numeric_limits<T>::min(), std::numeric_limits<T>::max() );

   for(auto it = ma.begin(); it != ma.end(); ++it)
      *it = dist(rng);
}

template<>
void makeRandomMultiArray( std::vector<unsigned char> & ma) {
   mt11213b rng;
   std::uniform_int_distribution<walberla::int16_t> dist(std::numeric_limits<unsigned char>::min(), std::numeric_limits<unsigned char>::max());

   for (auto it = ma.begin(); it != ma.end(); ++it)
      *it = static_cast<unsigned char>( dist(rng));
}

template<>
void makeRandomMultiArray( std::vector<char> & ma)
{
   mt11213b rng;
   std::uniform_int_distribution<int16_t> dist( std::numeric_limits<char>::min(), std::numeric_limits<char>::max() );
   for (auto it = ma.begin(); it != ma.end(); ++it)
      *it = static_cast<char>( dist(rng) );
}

template<typename T>
void runTests(const std::string & filename, size_t xSize, size_t ySize, size_t zSize);

void modifyHeader(std::string inputFilename, std::string outputFilename,
                  size_t xSize, size_t ySize, size_t zSize);


int main(int argc, char** argv)
{
   walberla::debug::enterTestMode();
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );

   WALBERLA_LOG_INFO("Starting test!");

   bool longrun = std::find(argv, argv + argc, std::string("--longrun")) != argv + argc;

   std::vector<size_t> sizes;
   if(longrun)
   {
      sizes.push_back(1);
      sizes.push_back(2);
      sizes.push_back(9);
      sizes.push_back(10);
   }
   sizes.push_back(111);

   for(std::vector<size_t>::const_iterator xSize = sizes.begin(); xSize != sizes.end(); ++xSize)
      for(std::vector<size_t>::const_iterator ySize = sizes.begin(); ySize != sizes.end(); ++ySize)
         for(std::vector<size_t>::const_iterator zSize = sizes.begin(); zSize != sizes.end(); ++zSize)
         {
            std::stringstream ss;
            ss << "geometry_testfile_" << *xSize << "_" << *ySize << "_" << *zSize << ".dat";

            std::string filename = ss.str();

            runTests<unsigned char>(filename, *xSize, *ySize, *zSize);

            if(longrun)
            {
               runTests<char>(filename, *xSize, *ySize, *zSize);

               runTests<short>(filename, *xSize, *ySize, *zSize);
               runTests<unsigned short>(filename, *xSize, *ySize, *zSize);

               runTests<int>(filename, *xSize, *ySize, *zSize);
               runTests<unsigned int>(filename, *xSize, *ySize, *zSize);

               runTests<long>(filename, *xSize, *ySize, *zSize);
               runTests<unsigned long>(filename, *xSize, *ySize, *zSize);
            }
         }
}

template<typename T>
void runTests(const std::string & filename, size_t xSize, size_t ySize, size_t zSize)
{
   using namespace walberla;
   using geometry::VoxelFileReader;

   WALBERLA_LOG_INFO( "Running Test with size " << xSize << "x" << ySize << "x" << zSize << " T = " << typeid(T).name() );

   filesystem::path path(filename);

   if( filesystem::exists( path ) )
      filesystem::remove( path );

   CellInterval aabb(0, 0, 0, cell_idx_c(xSize - 1), cell_idx_c(ySize - 1), cell_idx_c(zSize - 1));
   uint_t numCells = aabb.numCells();

   VoxelFileReader<T> geometryFile(filename, xSize, ySize, zSize);

   WALBERLA_CHECK( geometryFile.isOpen() );
   WALBERLA_CHECK_EQUAL( geometryFile.filename(), filename );
   WALBERLA_CHECK_EQUAL( geometryFile.xSize(), xSize );
   WALBERLA_CHECK_EQUAL( geometryFile.ySize(), ySize );
   WALBERLA_CHECK_EQUAL( geometryFile.zSize(), zSize );
   WALBERLA_CHECK_EQUAL( geometryFile.numCells(), xSize * ySize* zSize );


   std::vector<T> data;
   randomizeVector(data);
   geometryFile.read(aabb, data);
   WALBERLA_CHECK_EQUAL( numCells, data.size() );

   WALBERLA_CHECK_EQUAL( xSize, geometryFile.xSize() );
   WALBERLA_CHECK_EQUAL( ySize, geometryFile.ySize() );
   WALBERLA_CHECK_EQUAL( zSize, geometryFile.zSize() );
   WALBERLA_CHECK_EQUAL( geometryFile.numCells(), xSize * ySize* zSize );

   ptrdiff_t numDefaultConstructedElements = std::count( data.begin(), data.end(), T() );
   WALBERLA_CHECK_EQUAL( numCells, numDefaultConstructedElements );

   geometryFile.close();
   WALBERLA_CHECK( !geometryFile.isOpen() );

   data.clear();

   WALBERLA_CHECK( filesystem::exists( path ) );
   WALBERLA_CHECK( filesystem::is_regular_file( path ) );

   geometryFile.open(filename);
   WALBERLA_CHECK( geometryFile.isOpen() );
   WALBERLA_CHECK_EQUAL( geometryFile.filename(), filename );
   WALBERLA_CHECK_EQUAL( geometryFile.xSize(), xSize );
   WALBERLA_CHECK_EQUAL( geometryFile.ySize(), ySize );
   WALBERLA_CHECK_EQUAL( geometryFile.zSize(), zSize );
   WALBERLA_CHECK_EQUAL( geometryFile.numCells(), xSize * ySize* zSize );

   randomizeVector(data);
   geometryFile.read(aabb, data);
   WALBERLA_CHECK_EQUAL( numCells, data.size() );
   WALBERLA_CHECK_EQUAL( xSize, geometryFile.xSize() );
   WALBERLA_CHECK_EQUAL( ySize, geometryFile.ySize() );
   WALBERLA_CHECK_EQUAL( zSize, geometryFile.zSize() );

   numDefaultConstructedElements = std::count( data.begin(), data.end(), T() );
   WALBERLA_CHECK_EQUAL( numCells, numDefaultConstructedElements );

   geometryFile.close();
   WALBERLA_CHECK( !geometryFile.isOpen() );

   geometryFile.create(filename, xSize, ySize, zSize, 7);
   WALBERLA_CHECK( geometryFile.isOpen() );
   WALBERLA_CHECK_EQUAL( geometryFile.filename(), filename );
   WALBERLA_CHECK_EQUAL( geometryFile.xSize(), xSize );
   WALBERLA_CHECK_EQUAL( geometryFile.ySize(), ySize );
   WALBERLA_CHECK_EQUAL( geometryFile.zSize(), zSize );
   WALBERLA_CHECK_EQUAL( geometryFile.numCells(), xSize * ySize* zSize );

   data.clear();
   randomizeVector(data);
   geometryFile.read(aabb, data);

   WALBERLA_CHECK_EQUAL( numCells, data.size() );
   WALBERLA_CHECK_EQUAL( xSize, geometryFile.xSize() );
   WALBERLA_CHECK_EQUAL( ySize, geometryFile.ySize() );
   WALBERLA_CHECK_EQUAL( zSize, geometryFile.zSize() );

   numDefaultConstructedElements = std::count( data.begin(), data.end(), numeric_cast<T>(7) );
   WALBERLA_CHECK_EQUAL( numCells, numDefaultConstructedElements );

   data.clear();

   std::vector<T> reference(zSize * ySize * xSize);

   makeRandomMultiArray(reference);

   data.resize(reference.size());
   std::copy( reference.begin(), reference.end(), data.begin() );

   geometryFile.open(filename);
   WALBERLA_CHECK( geometryFile.isOpen() );
   WALBERLA_CHECK_EQUAL( geometryFile.filename(), filename );
   WALBERLA_CHECK_EQUAL( geometryFile.xSize(), xSize );
   WALBERLA_CHECK_EQUAL( geometryFile.ySize(), ySize );
   WALBERLA_CHECK_EQUAL( geometryFile.zSize(), zSize );
   WALBERLA_CHECK_EQUAL( geometryFile.numCells(), xSize * ySize* zSize );

   geometryFile.write(aabb, data);
   data.clear();
   randomizeVector(data);
   geometryFile.read(aabb, data);

   WALBERLA_CHECK_EQUAL(reference.size(), data.size());
   WALBERLA_CHECK( std::equal(reference.begin(), reference.end(), data.begin()) );

   randomizeVector(data);

   geometryFile.create(filename, xSize, ySize, zSize, reference.data() );
   WALBERLA_CHECK( geometryFile.isOpen() );
   WALBERLA_CHECK_EQUAL( geometryFile.filename(), filename );
   WALBERLA_CHECK_EQUAL( geometryFile.xSize(), xSize );
   WALBERLA_CHECK_EQUAL( geometryFile.ySize(), ySize );
   WALBERLA_CHECK_EQUAL( geometryFile.zSize(), zSize );
   WALBERLA_CHECK_EQUAL( geometryFile.numCells(), xSize * ySize* zSize );

   geometryFile.read(aabb, data);

   WALBERLA_CHECK_EQUAL(reference.size(), data.size());
   WALBERLA_CHECK( std::equal(reference.begin(), reference.end(), data.begin()) );

   std::vector<size_t> blockSizesX;
   blockSizesX.push_back(std::max(xSize / 2, size_t(1)));
   if( xSize > size_t(1) )
      blockSizesX.push_back(std::max(xSize / 2 - 1, size_t(1)));

   std::vector<size_t> blockSizesY;
   blockSizesY.push_back(std::max(ySize / 2, size_t(1)));
   if( ySize > size_t( 1 ) )
      blockSizesY.push_back(std::max(ySize / 2 - 1, size_t(1)));

   std::vector<size_t> blockSizesZ;
   blockSizesZ.push_back(std::max(zSize / 2, size_t(1)));
   if( zSize > size_t( 1 ) )
      blockSizesZ.push_back(std::max(zSize / 2 - 1, size_t(1)));

   for( auto xit = blockSizesX.begin(); xit != blockSizesX.end(); ++xit ) { size_t blockSizeX = *xit;
      for( auto yit = blockSizesY.begin(); yit != blockSizesY.end(); ++yit ) { size_t blockSizeY = *yit;
         for( auto zit = blockSizesZ.begin(); zit != blockSizesZ.end(); ++zit ) { size_t blockSizeZ = *zit;
            for(size_t zz = 0; zz <= (zSize - 1) / blockSizeZ; ++zz) {
               for(size_t yy = 0; yy <= (ySize - 1) / blockSizeY; ++yy) {
                  for(size_t xx = 0; xx <= (xSize - 1) / blockSizeX; ++xx)
                  {
                     CellInterval blockAABB;
                     blockAABB.xMin() = cell_idx_c(xx * blockSizeX);
                     blockAABB.yMin() = cell_idx_c(yy * blockSizeY);
                     blockAABB.zMin() = cell_idx_c(zz * blockSizeZ);
                     blockAABB.xMax() = cell_idx_c(std::min(blockAABB.xMin() + cell_idx_c(blockSizeX), cell_idx_c(geometryFile.xSize())) - 1);
                     blockAABB.yMax() = cell_idx_c(std::min(blockAABB.yMin() + cell_idx_c(blockSizeY), cell_idx_c(geometryFile.ySize())) - 1);
                     blockAABB.zMax() = cell_idx_c(std::min(blockAABB.zMin() + cell_idx_c(blockSizeZ), cell_idx_c(geometryFile.zSize())) - 1);
                     WALBERLA_CHECK_LESS(blockAABB.xMin(), xSize);
                     WALBERLA_CHECK_LESS(blockAABB.yMin(), ySize);
                     WALBERLA_CHECK_LESS(blockAABB.zMin(), zSize);
                     WALBERLA_CHECK_LESS(blockAABB.xMax(), xSize);
                     WALBERLA_CHECK_LESS(blockAABB.yMax(), ySize);
                     WALBERLA_CHECK_LESS(blockAABB.zMax(), zSize);

                     geometryFile.read(blockAABB, data);
                     WALBERLA_CHECK_EQUAL( data.size(), blockAABB.numCells() );

                     auto zIndexOffset = walberla::numeric_cast<size_t>(blockAABB.zMin());
                     auto yIndexOffset = walberla::numeric_cast<size_t>(blockAABB.yMin());
                     auto xIndexOffset = walberla::numeric_cast<size_t>(blockAABB.xMin());

                     size_t vectorIdx = 0;
                     for(size_t z = 0; z < blockAABB.zSize(); ++z)
                        for(size_t y = 0; y < blockAABB.ySize(); ++y)
                           for(size_t x = 0; x < blockAABB.xSize(); ++x)
                           {
                              WALBERLA_CHECK_EQUAL(data[vectorIdx], reference[(zIndexOffset+z)*ySize*xSize + (yIndexOffset+y)*xSize + (xIndexOffset+x)]);
                              ++vectorIdx;
                           }
                  }
               }
            }
         }
      }
   }

   geometryFile.close();
   WALBERLA_CHECK( !geometryFile.isOpen() );

   if( filesystem::exists( path ) )
      filesystem::remove( path );

   geometryFile.create(filename, xSize, ySize, zSize);

   for( auto xit = blockSizesX.begin(); xit != blockSizesX.end(); ++xit ) { size_t blockSizeX = *xit;
      for( auto yit = blockSizesY.begin(); yit != blockSizesY.end(); ++yit ) { size_t blockSizeY = *yit;
         for( auto zit = blockSizesZ.begin(); zit != blockSizesZ.end(); ++zit ) { size_t blockSizeZ = *zit;
            for(size_t zz = 0; zz <= (zSize - 1) / blockSizeZ; ++zz) {
               for(size_t yy = 0; yy <= (ySize - 1) / blockSizeY; ++yy) {
                  for(size_t xx = 0; xx <= (xSize - 1) / blockSizeX; ++xx)
                  {
                     CellInterval blockAABB;
                     blockAABB.xMin() = cell_idx_c(xx * blockSizeX);
                     blockAABB.yMin() = cell_idx_c(yy * blockSizeY);
                     blockAABB.zMin() = cell_idx_c(zz * blockSizeZ);
                     blockAABB.xMax() = cell_idx_c(std::min(blockAABB.xMin() + cell_idx_c(blockSizeX), cell_idx_c(geometryFile.xSize())) - 1);
                     blockAABB.yMax() = cell_idx_c(std::min(blockAABB.yMin() + cell_idx_c(blockSizeY), cell_idx_c(geometryFile.ySize())) - 1);
                     blockAABB.zMax() = cell_idx_c(std::min(blockAABB.zMin() + cell_idx_c(blockSizeZ), cell_idx_c(geometryFile.zSize())) - 1);
                     WALBERLA_CHECK_LESS(blockAABB.xMin(), xSize);
                     WALBERLA_CHECK_LESS(blockAABB.yMin(), ySize);
                     WALBERLA_CHECK_LESS(blockAABB.zMin(), zSize);
                     WALBERLA_CHECK_LESS(blockAABB.xMax(), xSize);
                     WALBERLA_CHECK_LESS(blockAABB.yMax(), ySize);
                     WALBERLA_CHECK_LESS(blockAABB.zMax(), zSize);

                     auto zIndexOffset = walberla::numeric_cast<size_t>(blockAABB.zMin());
                     auto yIndexOffset = walberla::numeric_cast<size_t>(blockAABB.yMin());
                     auto xIndexOffset = walberla::numeric_cast<size_t>(blockAABB.xMin());

                     data.resize(blockAABB.numCells());
                     size_t vectorIdx = 0;
                     for(size_t z = 0; z < blockAABB.zSize(); ++z)
                        for(size_t y = 0; y < blockAABB.ySize(); ++y)
                           for(size_t x = 0; x < blockAABB.xSize(); ++x)
                           {
                              data[vectorIdx] = reference[(zIndexOffset+z)*ySize*xSize + (yIndexOffset+y)*xSize + (xIndexOffset+x)];
                              ++vectorIdx;
                           }

                     geometryFile.write(blockAABB, data);
                  }
               }
            }
            geometryFile.read(aabb, data);
            WALBERLA_CHECK_EQUAL(reference.size(), data.size());
            WALBERLA_CHECK( std::equal(reference.begin(), reference.end(), data.begin()) );
         }
      }
   }

   modifyHeader(filename, filename + "0", xSize + 1, ySize, zSize);
   bool runtimeErrorThrown = false;
   try
   {
      WALBERLA_LOG_INFO("The following Error is expected!");
      Abort::instance()->resetAbortFunction( &Abort::exceptionAbort );
      geometryFile.open(filename + "0");
      Abort::instance()->resetAbortFunction();
      WALBERLA_CHECK(false);
   }
   catch( const std::runtime_error & /*e*/ )
   {
      Abort::instance()->resetAbortFunction();
      runtimeErrorThrown = true;
   }
   WALBERLA_CHECK( runtimeErrorThrown );
   geometryFile.close();
   if( filesystem::exists( filesystem::path(filename + "0") ) )
      filesystem::remove( filesystem::path(filename + "0") );

   if(xSize > 0)
   {
      modifyHeader(filename, filename + "1", xSize - 1, ySize, zSize);
      runtimeErrorThrown = false;
      try
      {
         WALBERLA_LOG_INFO("The following Error is expected!");
         Abort::instance()->resetAbortFunction( &Abort::exceptionAbort );
         geometryFile.open(filename + "1");
         Abort::instance()->resetAbortFunction();
         WALBERLA_CHECK(false);
      }
      catch( const std::runtime_error & /*e*/ )
      {
         Abort::instance()->resetAbortFunction();
         runtimeErrorThrown = true;
      }
      WALBERLA_CHECK( runtimeErrorThrown );
      geometryFile.close();
      if( filesystem::exists( filesystem::path(filename + "1") ) )
         filesystem::remove( filesystem::path(filename + "1") );
   }

   if( filesystem::exists( path ) )
      filesystem::remove( path );

}

void modifyHeader(std::string inputFilename, std::string outputFilename,
                  size_t xSize, size_t ySize, size_t zSize)
{
   std::ifstream is(inputFilename.c_str(), std::fstream::in | std::fstream::binary);
   WALBERLA_CHECK( is.is_open() );
   std::string oldHeader;
   std::getline(is, oldHeader);

   std::ofstream os(outputFilename.c_str(), std::fstream::out | std::fstream::binary);
   WALBERLA_CHECK( os.is_open() );
   os << xSize << " " << ySize << " " << zSize << "\n";

   while(!is.eof())
   {
      char buffer[1024];
      is.read(buffer, 1024);
      std::streamsize bytesRead = is.gcount();
      os.write( buffer, bytesRead );
      WALBERLA_CHECK( is.eof() || !is.fail() );
      WALBERLA_CHECK( !os.fail() );
   }

   is.close();
   os.close();
}


