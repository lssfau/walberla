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
//! \file VoxelFileReader.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Declares class StructuredGeometryFileReader that provides a reader for waLBerla geometry files.
//!
//! This reader is basically a wrapper around StructuredGeometryFileBasicReader utilizing the waLBerla data
//! structures
//
//======================================================================================================================

#pragma once

#include "BasicVoxelFileReader.h"
#include "core/Abort.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"


namespace walberla {
namespace geometry {

/*******************************************************************************************************************//**
 *
 * \brief Provides a reader for waLBerla geometry files.
 *
 * Note that the data is stored in binary form. There is no correction for other binary data
 * representations on different architectures (e.g. different endianness)!
 * An opened file is automatically closed upon destruction.
 *
 * \tparam T The underlying datatype that is stored in binary form in the geometry file
 *
 * \ingroup geometry
 *
 **********************************************************************************************************************/
template <typename T>
class VoxelFileReader
{
public:
	VoxelFileReader();
	VoxelFileReader( const std::string & _filename);
	VoxelFileReader( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, T value = T() );
	VoxelFileReader( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, const T * values );

	void open  ( const std::string & _filename );
	void create( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, T value = T() );
	void create( const std::string & _filename, uint_t _xSize, uint_t _ySize, uint_t _zSize, const T * values );
	void close ();

	bool isOpen() const;
	const std::string & filename() const;
	uint_t numCells() const;

	uint_t xSize() const;
	uint_t ySize() const;
	uint_t zSize() const;

	void read ( const CellInterval & cellInterval,       std::vector<T> & data ) const;
	void write( const CellInterval & cellInterval, const std::vector<T> & data );
private:
	BasicVoxelFileReader<T> geometryFile_; ///< The low-level geometry file reader
};

inline CellAABB toCellAABB( const CellInterval & cellInterval );
inline void validateCellInterval( const CellInterval & cellInterval );

} // namespace geometry
} // namespace walberla

#include "VoxelFileReader.impl.h"
