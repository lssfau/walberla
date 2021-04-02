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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <blockforest/BlockForest.h>
#include <core/DataTypes.h>

#include <map>
#include <string>

namespace walberla {
namespace mesa_pd {

void addBuildInfoToSQL(std::map<std::string, int64_t> &integerProperties,
                       std::map<std::string, double> &realProperties,
                       std::map<std::string, std::string> &stringProperties);

void addDomainPropertiesToSQL(const walberla::blockforest::BlockForest &forest,
                              std::map<std::string, int64_t> &integerProperties,
                              std::map<std::string, double> &realProperties,
                              std::map<std::string, std::string> & /*stringProperties*/ );

void addSlurmPropertiesToSQL(std::map<std::string, int64_t> &integerProperties,
                             std::map<std::string, double> &realProperties,
                             std::map<std::string, std::string> &stringProperties);

} //namespace mesa_pd
} //namespace walberla
