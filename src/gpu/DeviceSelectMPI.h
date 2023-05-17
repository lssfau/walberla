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
//! \file DeviceSelectMPI.h
//! \ingroup gpu
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once


namespace walberla {
namespace gpu
{


/**
 * Selects active GPU device based on MPI rank
 *
 * assumes that on each node there are as many MPI processes started as there are GPUs
 * - if there are more GPUs than processes on a node, a warning is printed and not all GPUs are utilized
 * - if there are more processes than GPUs, also a warning is printed and multiple processes may access the same GPU.
 *   Processes are assigned to GPUs in a round-robin fashion
 */
void selectDeviceBasedOnMpiRank();

} // namespace gpu
} // namespace walberla