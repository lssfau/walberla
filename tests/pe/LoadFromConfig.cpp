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
//! \file LoadFromConfig.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "domain_decomposition/all.h"

#include "pe/basic.h"
#include "pe/utility/CreateWorld.h"

#include "core/debug/TestSubsystem.h"

#include <boost/tuple/tuple.hpp>

using namespace walberla;
using namespace walberla::pe;

int main( int argc, char ** argv )
{
   walberla::debug::enterTestMode();

   Environment env(argc, argv);
   const Config::BlockHandle configBlock  = env.config()->getBlock( "LoadFromConfig" );

   shared_ptr<BodyStorage> globalBodyStorage = make_shared<BodyStorage>();

   // create blocks
   shared_ptr<BlockForest> forest = createBlockForestFromConfig( configBlock );
   WALBERLA_CHECK_EQUAL( forest->getXSize(), 3 );
   WALBERLA_CHECK_EQUAL( forest->getYSize(), 4 );
   WALBERLA_CHECK_EQUAL( forest->getZSize(), 5 );
   WALBERLA_CHECK( !forest->isXPeriodic() );
   WALBERLA_CHECK(  forest->isYPeriodic() );
   WALBERLA_CHECK( !forest->isZPeriodic() );
   WALBERLA_CHECK_FLOAT_EQUAL( forest->getDomain().minCorner(), Vec3(-15, -15, 0) );
   WALBERLA_CHECK_FLOAT_EQUAL( forest->getDomain().maxCorner(), Vec3(-3, 8, 34) );

   BlockDataID blockDataID;
   cr::HCSITS hcsits( globalBodyStorage, forest, blockDataID, blockDataID, blockDataID);
   configure(configBlock, hcsits);
   WALBERLA_CHECK_EQUAL( hcsits.getRelaxationModel(), cr::HCSITS::RelaxationModel::ApproximateInelasticCoulombContactByDecoupling );
   WALBERLA_CHECK_EQUAL( hcsits.getMaxIterations(), 123 );
   WALBERLA_CHECK_FLOAT_EQUAL( hcsits.getRelaxationParameter(), real_t(0.123) );
   WALBERLA_CHECK_FLOAT_EQUAL( hcsits.getErrorReductionParameter(), real_t(0.123) );
   WALBERLA_CHECK_FLOAT_EQUAL( hcsits.getGlobalLinearAcceleration(), Vec3(1,-2,3) );

   return EXIT_SUCCESS;
}
