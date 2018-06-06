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
//! \file DataProcessor.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief DataProcessor Interface
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"

#include <string>
#include <vector>


namespace walberla {
namespace gather {

/**
 * Interface for classes that operate on "Graph" like data
 * Can be implemented in classes that do output (print a GnuPlot file)
 * or for post-processing classes (calculation of mean value, std deviation..)
 *
 * One could build a chain: first the data are collected, then processed by
 * a statistical reductor, and finally written by a GnuPlot output class
 */
class DataProcessor
{
   public:

      virtual ~DataProcessor() = default;

      /*
       * Process "graph like" data
       * every entry in the outer vector is a data-point.
       */
      virtual void process(const std::vector<std::vector<real_t> > & data) = 0;
};


} // namespace gather
} // namespace walberla

