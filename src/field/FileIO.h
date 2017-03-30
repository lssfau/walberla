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
//! \file FileIO.h
//! \ingroup field
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include <core/mpi/MPIWrapper.h>
#include <core/mpi/Reduce.h>

#include <domain_decomposition/BlockStorage.h>

#include <fstream>
#include <string>
#include <vector>

namespace walberla {
namespace field {



//======================================================================================================================
/*!
 *  \brief Writes a field from a BlockStorage to file
 *
 *  Only the inner cells of a Field are written to file, ghost layer cells are ignored.
 *
 *  The generated file is binary with no regard to endianness, so it is only save to read it on the same computer
 *  architecture. 
 *
 *  Blocks are processed in the order they appear in the given BlockStorage. This requires that the BlockStorage is
 *  structured the same way when the file is read again. It is advised that you save the SetupBlockForest that was used
 *  to create the BlockStorage during writing and use it to create the BlockStorage when reading the file, since load
 *  balancers may employ randomness for their decisions.
 *
 *  If the file specified by filename already exists it is overwritten.
 *
 *  This is a collective function, it has to be called by all MPI processes simultaneously.
 *
 *  \param filename     The name of the file to be created
 *  \param blockStorage The BlockStorage the field is registered at
 *  \param fieldID      The ID of the field as returned by the BlockStorage at its registration
 */
//======================================================================================================================
template< typename FieldT >
void writeToFile( const std::string & filename, const BlockStorage & blockStorage, const BlockDataID & fieldID,
                  const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(), const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );



//======================================================================================================================
/*!
*  \brief Reads a field from a file
*
*  Only the inner cells of a Field are read from file, ghost layer cells are ignored.
*
*  The file is binary with no regard to endianness, so it is only save to read it on the same computer architecture as 
*  it was written on.
*
*  Blocks are processed in the order they appear in the given BlockStorage. This requires that the BlockStorage is
*  structured the same way as when the file was written. It is advised that you save the SetupBlockForest that was used
*  to create the BlockStorage during writing and use it to create the BlockStorage when reading the file, since load
*  balancers may employ randomness for their decisions.
*
*  This is a collective function, it has to be called by all MPI processes simultaneously.
*
*  \param filename     The name of the file to be read
*  \param blockStorage The BlockStorage the field is registered at
*  \param fieldID      The ID of the field as returned by the BlockStorage at its registration
*/
//======================================================================================================================
template< typename FieldT >
void readFromFile( const std::string & filename, BlockStorage & blockStorage, const BlockDataID & fieldID,
                   const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(), const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() );



} // namespace walberla
} // namespace field

#include "FileIO.impl.h"

