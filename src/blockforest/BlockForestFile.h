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
//! \file BlockForestFile.h
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

namespace walberla {
namespace blockforest {
namespace internal {

//**********************************************************************************************************************
/*!
 *   \file BlockForestFile.h
 *   \brief Description of the BlockForest save file format.
 *
 *   \section FileFormat File Format
 *
 *   \subsection HEADER HEADER
 *
 *   BYTES                      | DESCRIPTION
 *   ---------------------------|-----------------
 *   6 x ( 3 + sizeof(real_t) ) | domain AABB
 *   3 x 4                      | number of coarse/root blocks in each direction ( max 2^32 = 4 294 967 296 )
 *   3 x 1                      | domain periodicity
 *   1                          | block forest depth (= number of levels - 1)
 *   1                          | treeIdDigits (= number of bits used for storing the tree ID [tree ID marker + tree index])
 *   1                          | processIdBytes (= number of bytes required for storing process IDs)
 *   1                          | insertBuffersIntoProcessNetwork? ( 0=no, 1=yes )
 *   4                          | number of processes ( max 2^32 = 4 294 967 296 )
 *
 *   --> 23 + 6 x ( 3 + sizeof(real_t) ) BYTES
 *
 *   \subsection SUID SUID MAPPING:
 *
 *   1 | number of SUIDs (= \#SUIDs)
 *
 *   \code{.unparsed}
 *   for each SUID:
 *      1                | length of the UID identifier string
 *      length-of-string | UID identifier string
 *   \endcode
 *
 *   --> 1 + \#SUIDs + number-of-characters-of-all-identifiers-combined BYTES
 *
 *   How the mapping works:\n
 *   SUID #1 is assigned bit #1 ( -> [...]0 0000 0001 )\n
 *   SUID #2 is assigned bit #2 ( -> [...]0 0000 0010 )\n
 *   SUID #3 is assigned bit #3 ( -> [...]0 0000 0100 )\n
 *   ...\n
 *   For every block a bit mask containing information about all SUIDs (i.e., is the corresponding SUID set at this block?) is saved.
 *   -> The number of available SUIDs determines the size that is needed to store this bit mask (= SUID-mask-bytes).
 *   One byte is enough to hold 8 SUIDs, two bytes are enough to hold 16 SUIDs, ...
 *
 *   \subsection BLOCKDATA BLOCK DATA
 *
 *   \code{.unparsed}
 *   for each process:
 *      2 | number of blocks (can be '0' -> buffer process! - 2^16 = 65 536 )
 *      if( number-of-blocks > 0 ):
 *         for each block:
 *            block-ID-bytes  | ID of the block (the number of bytes required for storing the block ID largely depends on the size
 *                              of the simulation, the total number of blocks, and the number of refinement levels)
 *            SUID-mask-bytes | state of the block = bit mask containing information about all SUIDs (see "How the mapping works" for SUIDs,
 *                              SUID-mask-bytes can be equal to 0 bytes if no SUIDs exist!)
 *      2 | number of neighbor processes
 *      for each neighbor process:
 *         process-ID-bytes | process ID / rank of the neighbor process (one byte if there are less than 257 processes,
 *                                                                          two bytes if there are less than 65 537 processes, ...)
 *   \endcode
 */
//**********************************************************************************************************************

static const uint_t FILE_HEADER_SIZE = 6 * sizeof( real_t ) + 6 + 12 + 3 * 4 + 3 + 1 + 1 + 1 + 1 + 4;

}
}
}
