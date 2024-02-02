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
//! \file Operation.h
//! \ingroup core
//! \author Michael Zikeli <michael.zikeli@fau.de>
//
//======================================================================================================================
#pragma once

namespace walberla::mpi
{
// Note: I don't like at all that this is an enum and not an enum class, but changing this would be a major change in the framework.
enum Operation { MIN, MAX, SUM, PRODUCT, LOGICAL_AND, BITWISE_AND, LOGICAL_OR, BITWISE_OR, LOGICAL_XOR, BITWISE_XOR };
} // namespace walberla::mpi