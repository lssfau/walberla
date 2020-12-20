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
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/ParticleStorage.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <algorithm>
#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

void basic_test()
{
   //init data structures
   data::ParticleStorage ps(100);

   ps.create();
   ps.create();
   ps.create();
   ps.create();
   ps.create();
   ps.create();
   ps.create();

   WALBERLA_CHECK_EQUAL( ps.size(), 7);
   for (size_t i = 0; i < ps.size(); ++i)
   {
      for (size_t j = 0; j < ps.size(); ++j)
      {
         if (i==j)
         {
            WALBERLA_CHECK_EQUAL(ps.getUid(i), ps.getUid(j));
         } else
         {
            WALBERLA_CHECK_UNEQUAL(ps.getUid(i), ps.getUid(j));
         }
      }

      auto it = ps.find(ps.getUid(i));
      WALBERLA_CHECK_EQUAL( it->getUid(), ps.getUid(i));
   }

   auto it  = data::ParticleStorage::iterator(&ps, 3);
   auto uid = ps.getUid(3);
   ps.erase(it);
   WALBERLA_CHECK_EQUAL( ps.size(), 6);
   WALBERLA_CHECK_EQUAL( ps.find(uid), ps.end());
}

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   basic_test();

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
