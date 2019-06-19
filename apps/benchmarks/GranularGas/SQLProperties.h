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
//! \file SQLProperties.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include <blockforest/loadbalancing/DynamicDiffusive.h>
#include <core/DataTypes.h>

#include <map>
#include <string>

namespace walberla { namespace blockforest{
class DynamicParMetis;
} }

namespace walberla {
namespace mesa_pd {

void addBuildInfoToSQL( std::map< std::string, int64_t > &       integerProperties,
                        std::map< std::string, double > &        realProperties,
                        std::map< std::string, std::string > & stringProperties );

void addDomainPropertiesToSQL( const walberla::blockforest::BlockForest& forest,
                               std::map< std::string, int64_t > &       integerProperties,
                               std::map< std::string, double > &        realProperties,
                               std::map< std::string, std::string > & /*stringProperties*/ );

void addLoadBalancingPropertiesToSQL( const walberla::blockforest::BlockForest& forest,
                                      std::map< std::string, int64_t > &       integerProperties,
                                      std::map< std::string, double > &        /*realProperties*/,
                                      std::map< std::string, std::string > & /*stringProperties*/ );

template <typename PhantomBlockWeight_T>
void addDynamicDiffusivePropertiesToSQL( const walberla::blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>& ddb,
                                         std::map< std::string, int64_t > &       integerProperties,
                                         std::map< std::string, double > &        /*realProperties*/,
                                         std::map< std::string, std::string > &   stringProperties )
{
   integerProperties[ "diffMaxIterations" ]                = walberla::int64_c(ddb.getMaxIterations());
   integerProperties[ "diffFlowIterations" ]               = walberla::int64_c(ddb.getFlowIterations());
   integerProperties[ "bDiffAbortEarly" ]                  = ( ddb.checkForEarlyAbort() ? 1 : 0 );
   integerProperties[ "bDiffAdaptInflow" ]                 = ( ddb.adaptInflowWithGlobalInformation() ? 1 : 0 );
   integerProperties[ "bDiffAdaptOutflow" ]                = ( ddb.adaptOutflowWithGlobalInformation() ? 1 : 0 );
   std::string diffModeStr = "unknown";
   if (ddb.getMode() == walberla::blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PUSH) diffModeStr = "push";
   if (ddb.getMode() == walberla::blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PULL) diffModeStr = "pull";
   if (ddb.getMode() == walberla::blockforest::DynamicDiffusionBalance<PhantomBlockWeight_T>::DIFFUSION_PUSHPULL) diffModeStr = "pushpull";
   stringProperties[ "diffMode" ]                          = diffModeStr;
}

void addParMetisPropertiesToSQL( const ::walberla::blockforest::DynamicParMetis&         dpm,
                                 std::map< std::string, int64_t > &        integerProperties,
                                 std::map< std::string, double > &         realProperties,
                                 std::map< std::string, std::string > &    stringProperties );

void addSlurmPropertiesToSQL( std::map< std::string, int64_t > &        integerProperties,
                              std::map< std::string, double > &         realProperties,
                              std::map< std::string, std::string > &    stringProperties );

} //namespace mesa_pd
} //namespace walberla
