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
//! \file VoxelFileReader.impl.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Defines class StructuredGeometryFileReader that provides a reader for waLBerla geometry files.
//!
//! Note that the data is stored in binary form. There is no correction for other binary data
//! representations on different architectures (e.g. different endianness)!
//
//======================================================================================================================



namespace walberla {
namespace geometry {

/*******************************************************************************************************************//**
 * \brief Constructs on empty geometry file
 *
 * \post isOpen() == false
 **********************************************************************************************************************/
template <typename T>
VoxelFileReader<T>::VoxelFileReader()
try {  WALBERLA_ASSERT( !isOpen() ); }
catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }

/*******************************************************************************************************************//**
 * \brief Opens an existing geometry file.
 *
 * \param filename Name (path) of the file.
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template <typename T>
VoxelFileReader<T>::VoxelFileReader( const std::string & _filename )
   try : geometryFile_( _filename ) {  WALBERLA_ASSERT( isOpen() ); }
catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }

/*******************************************************************************************************************//**
 * \brief Creates a new geometry file with extends xSize x ySize x zSize.
 *
 * \param filename Name (path) of the file.
 * \param xSize    Extend of the geometry file in x direction.
 * \param ySize    Extend of the geometry file in y direction.
 * \param zSize    Extend of the geometry file in z direction.
 * \param value    The value the cells are initialized with. Defaults to T().
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template <typename T>
VoxelFileReader<T>::VoxelFileReader( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, T value /*= T() */ )
   try : geometryFile_( _filename, numeric_cast<size_t>(_xSize), numeric_cast<size_t>(_ySize), numeric_cast<size_t>(_zSize), value )
{
   WALBERLA_ASSERT( isOpen() );
}
catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }


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
* \pre values != nullptr
*
* \post isOpen() == true
***********************************************************************************************************************/
template <typename T>
VoxelFileReader<T>::VoxelFileReader( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, const T * values )
try : geometryFile_( _filename, numeric_cast<size_t>(_xSize), numeric_cast<size_t>(_ySize), numeric_cast<size_t>(_zSize), values )
{
   WALBERLA_ASSERT( isOpen() );
}
catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }


/*******************************************************************************************************************//**
* \brief Opens an existing geometry file.
*
* 	An already opened file gets closed beforehand.
*
* \param filename Name (path) of the file.
*
* \post isOpen() == true
 **********************************************************************************************************************/
template <typename T>
void VoxelFileReader<T>::open( const std::string & _filename )
{
   try {
      geometryFile_.open(_filename);
   }
   catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }
   WALBERLA_ASSERT( isOpen() );
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
 * \post isOpen() == true
 **********************************************************************************************************************/
template <typename T>
void VoxelFileReader<T>::create( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, T value /*= T() */ )
{
   try {
      geometryFile_.create(_filename, numeric_cast<size_t>(_xSize), numeric_cast<size_t>(_ySize), numeric_cast<size_t>(_zSize), value);
   }
   catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }
   WALBERLA_ASSERT( isOpen() );
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
 * \pre values != nullptr
 *
 * \post isOpen() == true
 **********************************************************************************************************************/
template <typename T>
void VoxelFileReader<T>::create( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, const T * values )
{
   try {
      geometryFile_.create(_filename, numeric_cast<size_t>(_xSize), numeric_cast<size_t>(_ySize), numeric_cast<size_t>(_zSize), values);
      WALBERLA_ASSERT( isOpen() );
   }
   catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }
}

/*******************************************************************************************************************//**
 * \brief Query if a geometry file is open.
 *
 * \return true if a file is opened, false if not.
 **********************************************************************************************************************/
template <typename T>
bool VoxelFileReader<T>::isOpen() const
{
   return geometryFile_.isOpen();
}

/*******************************************************************************************************************//**
 * \brief Closes an opened geometry file.
 *
 * If no geometry file is open, this does nothing.
 *
 * \post isOpen() == false
 **********************************************************************************************************************/
template <typename T>
void VoxelFileReader<T>::close()
{
   try {
      geometryFile_.close();
      WALBERLA_ASSERT( !isOpen() );
   }
   catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }
}

/*******************************************************************************************************************//**
 * \brief Gets the filename of the opened geometry file.
 *
 * \pre isOpen() == true
 *
 * \return The filename of the opened geometry file.
 **********************************************************************************************************************/
template <typename T>
const std::string & VoxelFileReader<T>::filename() const
{
   return geometryFile_.filename();
}

/*******************************************************************************************************************//**
 * \brief Gets the number of cells of the currently loaded geometry file.
 *
 * \pre isOpen() == true
 *
 * \return xSize() * ySize() * zSize().
 **********************************************************************************************************************/
template <typename T>
uint_t VoxelFileReader<T>::numCells() const
{
   return uint_c( geometryFile_.numCells() );
}

/*******************************************************************************************************************//**
 * \brief Gets the extend of the geometry file in x direction.
 *
 * \pre isOpen() == true
 *
 * \return The extend of the geometry file in x direction.
 **********************************************************************************************************************/
template <typename T>
uint_t VoxelFileReader<T>::xSize() const
{
   return uint_c( geometryFile_.xSize() );
}

/*******************************************************************************************************************//**
 * \brief Gets the extend of the geometry file in y direction.
 *
 * \pre isOpen() == true
 *
 * \return The extend of the geometry file in y direction.
 **********************************************************************************************************************/
template <typename T>
uint_t VoxelFileReader<T>::ySize() const
{
   return uint_c( geometryFile_.ySize() );
}

/*******************************************************************************************************************//**
 * \brief Gets the extend of the geometry file in z direction.
 *
 * \pre isOpen() == true
 *
 * \return The extend of the geometry file in z direction.
 **********************************************************************************************************************/
template <typename T>
uint_t VoxelFileReader<T>::zSize() const
{
   return uint_c( geometryFile_.zSize() );
}

/*******************************************************************************************************************//**
 * \brief Reads a block of data from the opened geometry file.
 *
 * \param cellInterval The axis-aligned bounding box of the block of data to be read.
 * \param [out] data The vector the read data is stored to. The Storage order is zyx. (Meaning
 * 				      your innermost loop should iterate over x)
 *
 * \pre isOpen() == true
 *
 * \post data.size() == cellAABB.numCells()
 **********************************************************************************************************************/
template <typename T>
void VoxelFileReader<T>::read( const CellInterval & cellInterval, std::vector<T> & data ) const
{
   validateCellInterval( cellInterval );
   WALBERLA_ASSERT( isOpen() );
   try { geometryFile_.read( toCellAABB(cellInterval), data ); }
   catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }
}

/*******************************************************************************************************************//**
 * \brief Writes a block of data to the geometry file.
 *
 * \param cellInterval The axis-aligned bounding box of the block of data to be written.
 * \param data     The vector holding the data to bw written to the geometry file. The Storage
 * 					 order is zyx. (Meaning your innermost loop should iterate over x)
 *
 * \pre isOpen() == true
 * \pre !cellInterval.empty()
 * \pre cellInterval.positiveIndicesOnly()
 * \pre data.size() >= cellAABB.numCells()
 **********************************************************************************************************************/
template <typename T>
void VoxelFileReader<T>::write( const CellInterval & cellInterval, const std::vector<T> & data )
{
   validateCellInterval( cellInterval );

   try { geometryFile_.write( toCellAABB(cellInterval), data ); }
   catch( std::exception & e) { WALBERLA_ABORT( e.what() ); }
}

/*******************************************************************************************************************//**
 * \brief   Validate if the given CellInterval can be converted to a CellAABB.
 *
 * Logs an Error and aborts if either cellInterval.empty() or !cellInterval.positiveIndicesOnly()
 * since a CellAABB can neither be empty nor contain negative bounds.
 *
 * \param   cellInterval   The validated CellInterval.
 **********************************************************************************************************************/
void validateCellInterval( const CellInterval & cellInterval )
{
   if( cellInterval.empty() )
      WALBERLA_ABORT("CellInterval may not be empty! Given CellInterval is: " << cellInterval);
   if( !cellInterval.positiveIndicesOnly() )
      WALBERLA_ABORT("CellInterval may not have negative indices! Given CellInterval is: " << cellInterval);
}

/*******************************************************************************************************************//**
 * \brief   Converts a CellInterval to a CellAABB.
 *
 * A CellAABB can not hold negative bounds and therefore not be empty. You have to check yourself if
 * this is the case!
 *
 * \param   cellInterval   The CellInterval to be converted.
 *
 * \pre !cellInterval.empty()
 * \pre cellInterval.positiveIndicesOnly()
 *
 * \return  The CellAABB.
 **********************************************************************************************************************/
CellAABB toCellAABB( const CellInterval & cellInterval )
{
    WALBERLA_ASSERT( !cellInterval.empty(), cellInterval );
    WALBERLA_ASSERT( cellInterval.positiveIndicesOnly(), cellInterval );

    return CellAABB ( numeric_cast<size_t>( cellInterval.xMin() ),
                      numeric_cast<size_t>( cellInterval.yMin() ),
                      numeric_cast<size_t>( cellInterval.zMin() ),
                      numeric_cast<size_t>( cellInterval.xMax() ),
                      numeric_cast<size_t>( cellInterval.yMax() ),
                      numeric_cast<size_t>( cellInterval.zMax() ) );
}

} // namespace geometry
} // namespace walberla
