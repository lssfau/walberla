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
//! \file   NodeTimings.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "NodeTimings.h"

namespace walberla {
namespace mesa_pd {

void storeNodeTimings( const uint_t                 runId,
                       const std::string          & dbFile,
                       const std::string          & tableName,
                       const WcTimingPool         & tp )
{
   std::map< std::string, walberla::int64_t > integerProperties;
   std::map< std::string, double >            realProperties;
   std::map< std::string, std::string >       stringProperties;

   walberla::mpi::SendBuffer sb;
   walberla::mpi::RecvBuffer rb;

   sb << walberla::getHostName();
   sb << int64_t(walberla::mpi::MPIManager::instance()->rank());
   sb << tp;

   walberla::mpi::gathervBuffer(sb, rb);

   WALBERLA_ROOT_SECTION()
   {
      while (!rb.isEmpty())
      {
         integerProperties.clear();
         realProperties.clear();
         stringProperties.clear();

         std::string  hostname;
         int64_t      rank;
         WcTimingPool cTP;
         rb >> hostname;
         rb >> rank;
         rb >> cTP;

         stringProperties["hostname"] = hostname;
         integerProperties["rank"]    = rank;
         for (auto& v : cTP)
         {
            realProperties[v.first] = v.second.average();
         }

         sqlite::storeAdditionalRunInfoInSqliteDB( runId,
                                                           dbFile,
                                                           tableName,
                                                           integerProperties,
                                                           stringProperties,
                                                           realProperties );
      }
   }
}

} // namespace mesa_pd
} // namespace walberla
