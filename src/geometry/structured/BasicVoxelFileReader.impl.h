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
//! \file BasicVoxelFileReader.impl.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Defines class StructuredGeometryFileBasicReader that provides a low-level reader for waLBerla geometry files.
//
//======================================================================================================================

#include <cassert>
#include <sstream>
#include <stdexcept>
#include <typeinfo>


namespace walberla {
namespace geometry {

/*******************************************************************************************************************//**
 * \brief Constructs on empty geometry file
 *
 * \post isOpen() == false
 **********************************************************************************************************************/
template<typename T>
BasicVoxelFileReader<T>::BasicVoxelFileReader() : xSize_(0), ySize_(0), zSize_(0)
{
   assert( !isOpen() );
   assert( filename_.empty() );
   assert( dataBegin_ == std::streampos() );
   assert( xSize_ == 0 );
   assert( ySize_ == 0 );
   assert( zSize_ == 0 );
}

/*******************************************************************************************************************//**
 * \brief Opens an existing geometry file.
 *
 * \param filename Name (path) of the file.
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template<typename T>
BasicVoxelFileReader<T>::BasicVoxelFileReader( const std::string & _filename)
   :  xSize_(0), ySize_(0), zSize_(0)
{
   open(_filename);

   assert( dataBegin_ != std::streampos(-1) );
   assert( dataBegin_ != std::streampos() );
   assert( isOpen() );
   assert( xSize() != 0 );
   assert( ySize() != 0 );
   assert( zSize() != 0 );
   assert( filename() == _filename );
}

/*******************************************************************************************************************//**
 * \brief Creates a new geometry file with extends xSize x ySize x zSize.
 *
 * \param filename Name (path) of the file.
 * \param xSize    Extend of the geometry file in x direction.
 * \param ySize    Extend of the geometry file in y direction.
 * \param zSize    Extend of the geometry file in z direction.
 * \param value    The value the cells are initialized with. Defaults to T().
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template<typename T>
BasicVoxelFileReader<T>::BasicVoxelFileReader( const std::string & _filename, size_t _xSize, size_t _ySize, size_t _zSize, T value /*= T()*/ )
   : xSize_(0), ySize_(0), zSize_(0)
{
   create(_filename, _xSize, _ySize, _zSize, value);

   assert( dataBegin_ != std::streampos(-1) );
   assert( dataBegin_ != std::streampos() );
   assert( isOpen() );
   assert( xSize_ == _xSize );
   assert( ySize_ == _ySize );
   assert( zSize_ == _zSize );
   assert( filename_ == _filename );
}

/*******************************************************************************************************************//**
 * \brief Creates a new geometry file with extends xSize x ySize x zSize.
 *
 * \param filename Name (path) of the file.
 * \param xSize    Extend of the geometry file in x direction.
 * \param ySize    Extend of the geometry file in y direction.
 * \param zSize    Extend of the geometry file in z direction.
 * \param values   An array of size xSize * ySize * zSize with the values to initialize the
 * 					 geometry file with.
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \pre values != nullptr
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template<typename T>
BasicVoxelFileReader<T>::BasicVoxelFileReader( const std::string & _filename, size_t _xSize, size_t _ySize, size_t _zSize, const T * values )
   : xSize_(0), ySize_(0), zSize_(0)
{
   assert(values != 0);

   create(_filename, _xSize, _ySize, _zSize, values);

   assert( dataBegin_ != std::streampos(-1) );
   assert( dataBegin_ != std::streampos() );
   assert( isOpen() );
   assert( xSize_ == _xSize );
   assert( ySize_ == _ySize );
   assert( zSize_ == _zSize );
   assert( filename_ == _filename );
}

/*******************************************************************************************************************//**
 * \brief Destructor that closes the file if necessary.
 **********************************************************************************************************************/
template<typename T>
BasicVoxelFileReader<T>::~BasicVoxelFileReader()
{
   //filestream_.exceptions( std::ios_base::iostate(0) );
   if(filestream_.is_open())
      filestream_.close();
}

/*******************************************************************************************************************//**
* \brief Opens an existing geometry file.
*
* 	An already opened file gets closed beforehand.
*
* \param filename Name (path) of the file.
*
* \throws std::runtime_error on I/O errors.
* \throws std::runtime_error if the loaded geometry file's format is corrupt
*
* \post isOpen() == true
 **********************************************************************************************************************/
template<typename T>
void BasicVoxelFileReader<T>::open( const std::string & _filename )
{
   if( isOpen() )
      close();

   filename_ = _filename;

   filestream_.open(_filename.c_str(),  std::fstream::in | std::fstream::out | std::fstream::binary);
   if( filestream_.fail() || filestream_.bad())
      throw std::runtime_error("Error opening file \"" + _filename + "\"!");

   std::string firstLine;
   std::getline(filestream_, firstLine);
   if( filestream_.fail() || filestream_.bad() )
      throw std::runtime_error("Error reading header from file \"" + _filename + "\"!");
   {
      std::stringstream ss(firstLine);
      ss >> xSize_;
      ss >> ySize_;
      ss >> zSize_;
   }

   dataBegin_ = filestream_.tellg();
   if( filestream_.fail() || filestream_.bad() || dataBegin_ == std::streampos(-1) )
      throw std::runtime_error("I/O Error while reading file \"" + _filename + "\"!");

   filestream_.seekg(0, std::ios::end);
   if( filestream_.fail() || filestream_.bad() )
      throw std::runtime_error("I/O Error while reading file \"" + _filename + "\"!");

   std::streampos dataEnd = filestream_.tellg();
   if( filestream_.fail() || filestream_.bad() || dataEnd == std::streampos(-1) )
      throw std::runtime_error("I/O Error while reading file \"" + _filename + "\"!");


   assert(dataBegin_ <= dataEnd);
   size_t rawDataLengthBytes = static_cast<size_t>( dataEnd - dataBegin_ );
   if( rawDataLengthBytes % sizeof(T) != 0 )
   {
      std::stringstream ss;
      ss << "Raw data part of file " << _filename << " has size " << rawDataLengthBytes << " bytes. Data is of Type T = "
         << typeid(T).name() << " and sizeof(T) == " << sizeof(T) << "\nError: size % sizeof(T) != 0!";
      throw std::runtime_error(ss.str());
   }
   assert(rawDataLengthBytes % sizeof(T) == 0);
   size_t rawDataLength = rawDataLengthBytes / sizeof(T);
   if( rawDataLength != numCells() )
   {
      std::stringstream ss;
      ss << "Raw data part of file " << _filename << " has size " << rawDataLengthBytes << " bytes. Data is of Type T = "
         << typeid(T).name() << " and sizeof(T) == " << sizeof(T) << ". Error: Number of cells in data section of file = size / sizeof(T) = "
         << rawDataLength << " != Number of cells specified file header = " << numCells();
      throw std::runtime_error(ss.str());
   }

   assert(dataBegin_ != std::streampos(-1));
   assert(dataBegin_ != std::streampos());
   assert(filename_ == _filename);
   assert(xSize_ != 0);
   assert(ySize_ != 0);
   assert(zSize_ != 0);
}

/*******************************************************************************************************************//**
 * \brief Creates a new geometry file with extends xSize x ySize x zSize.
 *
 * An already opened file gets closed beforehand.
 *
 * \param filename Name (path) of the file.
 * \param xSize    Extend of the geometry file in x direction.
 * \param ySize    Extend of the geometry file in y direction.
 * \param zSize    Extend of the geometry file in z direction.
 * \param value    The value the cells are initialized with. Defaults to T().
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template<typename T>
void BasicVoxelFileReader<T>::create( const std::string & _filename, size_t _xSize, size_t _ySize, size_t _zSize, T value /*= T()*/ )
{
   if( isOpen() )
      close();

   filename_ = _filename;

   xSize_ = _xSize;
   ySize_ = _ySize;
   zSize_ = _zSize;

   filestream_.open(_filename.c_str(), std::fstream::out | std::fstream::binary);
   if( filestream_.fail() || filestream_.bad())
      throw std::runtime_error("Error opening file \"" + _filename + "\"!");

   filestream_ << _xSize << " " << _ySize << " " << _zSize << "\n";
   if( filestream_.fail() || filestream_.bad())
      throw std::runtime_error("Error writing header to file \"" + _filename + "\"!");

   dataBegin_ = filestream_.tellp();
   if( filestream_.fail() || filestream_.bad() || dataBegin_ == std::streampos(-1) )
      throw std::runtime_error("I/O Error while writing file \"" + _filename + "\"!");

   const size_t maxChunkSizeBytes = 1024 * 1024; // 1 MegaByte
   const size_t numCellsPerChunk = maxChunkSizeBytes / sizeof(T);
   const std::streamsize bytesPerChunk = static_cast< std::streamsize >( numCellsPerChunk * sizeof(T) );

   const std::vector<T> chunkData(numCellsPerChunk, value);
   assert(!chunkData.empty());
   const char * rawData = reinterpret_cast<const char*>(&chunkData[0]);

   size_t cellsLeft = numCells();
   while( cellsLeft >= numCellsPerChunk )
   {
      filestream_.write(rawData, bytesPerChunk);
      if( filestream_.fail() || filestream_.bad() )
         throw std::runtime_error("I/O Error while writing file \"" + _filename + "\"!");
      cellsLeft -= numCellsPerChunk;
   }

   const std::streamsize bytesLastChunk = static_cast< std::streamsize >( cellsLeft * sizeof(T) );
   assert(bytesLastChunk < bytesPerChunk);
   filestream_.write(rawData, bytesLastChunk);
   if( filestream_.fail() || filestream_.bad() )
      throw std::runtime_error("I/O Error while writing file \"" + _filename + "\"!");

   close();
   open(_filename);

   assert(dataBegin_ != std::streampos(-1));
   assert(dataBegin_ != std::streampos());
   assert(filename_ == _filename);
   assert(xSize_ == _xSize);
   assert(ySize_ == _ySize);
   assert(zSize_ == _zSize);
}

/*******************************************************************************************************************//**
 * \brief Creates a new geometry file with extends xSize x ySize x zSize.
 *
 * An already opened file gets closed beforehand.
 *
 * \param filename Name (path) of the file.
 * \param xSize    Extend of the geometry file in x direction.
 * \param ySize    Extend of the geometry file in y direction.
 * \param zSize    Extend of the geometry file in z direction.
 * \param values   An array of size xSize * ySize * zSize with the values to initialize the
 * 					 geometry file with.
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \pre values != nullptr
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template<typename T>
void BasicVoxelFileReader<T>::create( const std::string & _filename, size_t _xSize, size_t _ySize, size_t _zSize, const T * values )
{
   if( isOpen() )
      close();

   filename_ = _filename;

   xSize_ = _xSize;
   ySize_ = _ySize;
   zSize_ = _zSize;

   filestream_.open(_filename.c_str(), std::fstream::out | std::fstream::binary);
   if( filestream_.fail() || filestream_.bad())
      throw std::runtime_error("Error opening file \"" + _filename + "\"!");

   filestream_ << _xSize << " " << _ySize << " " << _zSize << "\n";
   if( filestream_.fail() || filestream_.bad())
      throw std::runtime_error("Error writing header to file \"" + _filename + "\"!");

   dataBegin_ = filestream_.tellp();
   if( filestream_.fail() || filestream_.bad() || dataBegin_ == std::streampos(-1) )
      throw std::runtime_error("I/O Error while writing file \"" + _filename + "\"!");

   const char * rawData = reinterpret_cast<const char*>(values);
   size_t dataSizeBytes = numCells() * sizeof(T);
   filestream_.write(rawData, static_cast< std::streamsize >( dataSizeBytes ) );
   if( filestream_.fail() || filestream_.bad() )
      throw std::runtime_error("I/O Error while writing file \"" + _filename + "\"!");

   close();
   open(_filename);

   assert(dataBegin_ != std::streampos(-1));
   assert(dataBegin_ != std::streampos());
   assert(filename_ == _filename);
   assert(xSize_ == _xSize);
   assert(ySize_ == _ySize);
   assert(zSize_ == _zSize);
}

/*******************************************************************************************************************//**
 * \brief Closes an opened geometry file.
 *
 * If no geometry file is open, this does nothing.
 *
 * \post isOpen() == false
 **********************************************************************************************************************/
template<typename T>
void BasicVoxelFileReader<T>::close()
{
   if( isOpen() )
   {
      filestream_.close();
      dataBegin_ = std::streampos();
      filename_.clear();
      xSize_ = size_t(0);
      ySize_ = size_t(0);
      zSize_ = size_t(0);
   }

   assert(!filestream_.is_open());
   assert(dataBegin_ == std::streampos());
   assert(filename_.empty());
   assert(xSize_ == 0);
   assert(ySize_ == 0);
   assert(zSize_ == 0);
   assert( !isOpen() );
}

/*******************************************************************************************************************//**
 * \brief Query if a geometry file is open.
 *
 * \return true if a file is opened, false if not.
 **********************************************************************************************************************/
template<typename T>
bool BasicVoxelFileReader<T>::isOpen() const
{
   return filestream_.is_open();
}

/*******************************************************************************************************************//**
 * \brief Gets the filename of the opened geometry file.
 *
 * \pre isOpen() == true
 *
 * \return The filename of the opened geometry file.
 **********************************************************************************************************************/
template<typename T>
const std::string & BasicVoxelFileReader<T>::filename() const
{
   assert( isOpen() );
   return filename_;
}

/*******************************************************************************************************************//**
 * \brief Gets the extend of the geometry file in x direction.
 *
 * \pre isOpen() == true
 *
 * \return The extend of the geometry file in x direction.
 **********************************************************************************************************************/
template<typename T>
size_t BasicVoxelFileReader<T>::xSize() const
{
   assert( isOpen() );
   return xSize_;
}

/*******************************************************************************************************************//**
 * \brief Gets the extend of the geometry file in y direction.
 *
 * \pre isOpen() == true
 *
 * \return The extend of the geometry file in y direction.
 **********************************************************************************************************************/
template<typename T>
size_t BasicVoxelFileReader<T>::ySize() const
{
   assert( isOpen() );
   return ySize_;
}

/*******************************************************************************************************************//**
 * \brief Gets the extend of the geometry file in z direction.
 *
 * \pre isOpen() == true
 *
 * \return The extend of the geometry file in z direction.
 **********************************************************************************************************************/
template<typename T>
size_t BasicVoxelFileReader<T>::zSize() const
{
   assert( isOpen() );
   return zSize_;
}

/*******************************************************************************************************************//**
 * \brief Gets the number of cells of the currently loaded geometry file.
 *
 * \pre isOpen() == true
 *
 * \return xSize() * ySize() * zSize().
 **********************************************************************************************************************/
template<typename T>
size_t BasicVoxelFileReader<T>::numCells() const
{
   assert( isOpen() );
   return xSize() * ySize() * zSize();
}

/*******************************************************************************************************************//**
 * \brief Reads a block of data from the opened geometry file.
 *
 * \param cellAABB   The axis-aligned bounding box of the block of data to be read.
 * \param [out] data The vector the read data is stored to. The Storage order is zyx. (Meaning
 * 				      your innermost loop should iterate over x)
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \pre isOpen() == true
 * \pre cellAABB.xEnd < xSize()
 * \pre cellAABB.yEnd < ySize()
 * \pre cellAABB.zEnd < zSize()
 *
 * \post data.size() == cellAABB.numCells()
 **********************************************************************************************************************/
template<typename T>
void BasicVoxelFileReader<T>::read( const CellAABB & cellAABB, std::vector<T> & data ) const
{
   assert( isOpen() );
   assert( cellAABB.xEnd < xSize() );
   assert( cellAABB.yEnd < ySize() );
   assert( cellAABB.zEnd < zSize() );

   data.clear();
   data.resize( cellAABB.numCells() );
   assert(!data.empty());
   char * buffer = reinterpret_cast<char*>( &data[0] );

   const size_t zFactor = xSize() * ySize();
   const size_t yFactor = xSize();
   std::streamsize lineLength = static_cast< std::streamsize >( cellAABB.xSize() * sizeof(T) );

   for(size_t z = cellAABB.zBegin; z <= cellAABB.zEnd; ++z)
      for(size_t y = cellAABB.yBegin; y <= cellAABB.yEnd; ++y)
      {
         assert( buffer < reinterpret_cast<char*>(&data[0] + data.size()) );
         assert( buffer + lineLength - 1 < reinterpret_cast<char*>(&data[0] + data.size()) );
         std::streamoff offset = static_cast< std::streamoff >( (z * zFactor + y * yFactor + cellAABB.xBegin) * sizeof(T) );

         filestream_.seekg(dataBegin_ + offset);
         if( filestream_.fail() || filestream_.bad() )
            throw std::runtime_error("I/O Error while reading file \"" + filename() + "\"!");

         filestream_.read(buffer, lineLength);
         if( filestream_.fail() || filestream_.bad() )
            throw std::runtime_error("I/O Error while reading file \"" + filename() + "\"!");

         assert(filestream_.gcount() == lineLength);

         buffer += lineLength;
      }

   assert( buffer == reinterpret_cast<char*>(&data[0] + data.size()) );
   assert( data.size() == cellAABB.numCells() );
}

/*******************************************************************************************************************//**
 * \brief Writes a block of data to the geometry file.
 *
 * \param cellAABB The axis-aligned bounding box of the block of data to be written.
 * \param data     The vector holding the data to bw written to the geometry file. The Storage
 * 					 order is zyx. (Meaning your innermost loop should iterate over x)
 *
 * \throws std::runtime_error on I/O errors.
 *
 * \pre isOpen() == true
 * \pre cellAABB.xEnd < xSize()
 * \pre cellAABB.yEnd < ySize()
 * \pre cellAABB.zEnd < zSize()
 * \pre data.size() >= cellAABB.numCells()
 **********************************************************************************************************************/
template<typename T>
void BasicVoxelFileReader<T>::write( const CellAABB & cellAABB, const std::vector<T> & data )
{
   assert( isOpen() );
   assert( cellAABB.xEnd < xSize_ );
   assert( cellAABB.yEnd < ySize_ );
   assert( cellAABB.zEnd < zSize_ );
   assert( !data.empty() );
   assert( data.size() >= cellAABB.numCells() );

   const char * buffer = reinterpret_cast<const char*>( &data[0] );

   const size_t zFactor = xSize() * ySize();
   const size_t yFactor = xSize();
   std::streamsize lineLength = static_cast< std::streamsize >( cellAABB.xSize() * sizeof(T) );

   for(size_t z = cellAABB.zBegin; z <= cellAABB.zEnd; ++z)
      for(size_t y = cellAABB.yBegin; y <= cellAABB.yEnd; ++y)
      {
         assert( buffer < reinterpret_cast<const char*>(&data[0] + data.size()) );
         assert( buffer + lineLength - 1 < reinterpret_cast<const char*>(&data[0] + data.size()) );
         std::streamoff offset = static_cast< std::streamoff >( (z * zFactor + y * yFactor + cellAABB.xBegin) * sizeof(T) );

         filestream_.seekp(dataBegin_ + offset);
         if( filestream_.fail() || filestream_.bad() )
            throw std::runtime_error("I/O Error while writing file \"" + filename() + "\"!");

         filestream_.write(buffer, lineLength);
         if( filestream_.fail() || filestream_.bad() )
            throw std::runtime_error("I/O Error while writing file \"" + filename() + "\"!");

         buffer += lineLength;
      }

   assert( buffer == reinterpret_cast<const char*>(&data[0] + cellAABB.numCells()) );
}

/*******************************************************************************************************************//**
 * \brief Constructor to initialize all members to 0.
 *
 * \post #xBegin = 0.
 * \post #yBegin = 0.
 * \post #zBegin = 0.
 * \post #xEnd = 0.
 * \post #yEnd = 0.
 * \post #zEnd = 0.
 * \post xSize() = 1
 * \post ySize() = 1
 * \post zSize() = 1
 * \post numCells() = 1
 **********************************************************************************************************************/
inline CellAABB::CellAABB()
   : xBegin(0), yBegin(0), zBegin(0), xEnd(0), yEnd(0), zEnd(0)
{ }

/*******************************************************************************************************************//**
 * \brief Constructor to initialize all members.
 *
 * \pre _xBegin <= _xEnd.
 * \pre _yBegin <= _yEnd.
 * \pre _zBegin <= _zEnd.
 *
 * \post #xBegin = _xBegin.
 * \post #yBegin = _yBegin.
 * \post #zBegin = _zBegin.
 * \post #xEnd = _xEnd.
 * \post #yEnd = _yEnd.
 * \post #zEnd = _zEnd.
 *
 * \param _xBegin The minimal x coordinate of all cells included in the AABB.
 * \param _yBegin The minimal y coordinate of all cells included in the AABB.
 * \param _zBegin The minimal z coordinate of all cells included in the AABB.
 * \param _xEnd   The maximal x coordinate of all cells included in the AABB.
 * \param _yEnd   The maximal y coordinate of all cells included in the AABB.
 * \param _zEnd   The maximal z coordinate of all cells included in the AABB.
 **********************************************************************************************************************/
CellAABB::CellAABB( size_t _xBegin, size_t _yBegin, size_t _zBegin,
                    size_t _xEnd,   size_t _yEnd,   size_t _zEnd )
   : xBegin(_xBegin), yBegin(_yBegin), zBegin(_zBegin),
     xEnd(_xEnd), yEnd(_yEnd), zEnd(_zEnd)
{
   assert( xBegin <= xEnd );
   assert( yBegin <= yEnd );
   assert( zBegin <= zEnd );
}

/*******************************************************************************************************************//**
 * \brief Gets the number cells included in the AABB.
 *
 * \pre #xBegin <= #xEnd.
 * \pre #yBegin <= #yEnd.
 * \pre #zBegin <= #zEnd.
 *
 * \return xSize() * ySize() * zSize().
 **********************************************************************************************************************/
size_t CellAABB::numCells() const
{
   assert( xBegin <= xEnd );
   assert( yBegin <= yEnd );
   assert( zBegin <= zEnd );

   return xSize() * ySize() * zSize();
}

/*******************************************************************************************************************//**
 * \brief Gets the extent of the AABB in z direction.
 *
 * \pre #xBegin <= #xEnd.
 *
 * \return #xEnd - #xBegin + 1.
 **********************************************************************************************************************/
size_t CellAABB::xSize() const
{
   assert( xBegin <= xEnd  );
   return xEnd - xBegin + 1;
}

/*******************************************************************************************************************//**
 * \brief Gets the extent of the AABB in y direction.
 *
 * \pre #yBegin <= #yEnd.
 *
 * \return #yEnd - #yBegin + 1.
 **********************************************************************************************************************/
size_t CellAABB::ySize() const
{
   assert( yBegin <= yEnd );
   return yEnd - yBegin + 1;
}

/*******************************************************************************************************************//**
 * \brief Gets the extent of the AABB in z direction.
 *
 * \pre #zBegin <= #zEnd.
 *
 * \return #zEnd - #zBegin + 1.
 **********************************************************************************************************************/
size_t CellAABB::zSize() const
{
   assert( zBegin <= zEnd );
   return zEnd - zBegin + 1;
}

} // namespace geometry
} // namespace walberla
