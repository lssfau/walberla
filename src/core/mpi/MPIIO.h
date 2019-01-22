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
//! \file MPIIO.h
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/mpi/RecvBuffer.h"
#include "core/mpi/SendBuffer.h"


namespace walberla {
namespace mpi  {

    /*! Writes contents of local buffer to a single binary file via MPIIO
     *
     * Has to be called by all processes
     */
   void writeMPIIO(const std::string & file, SendBuffer & buffer);


   /*! Counterpart to writeMPIIO - has to be called with exactly the same process distribution
    *
    * Reads local part of the data into a buffer
    */
   void readMPIIO(const std::string & file, RecvBuffer & buffer);


} // namespace mpi
} // namespace walberla


