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
//! \file PackPerformance.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/math/Vector3.h"
#include "core/mpi/BufferSystem.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingPool.h"

#include <array>
#include <iostream>
#include <sstream>

namespace walberla {

int main( int /*argc*/, char ** /*argv*/ )
{
   const size_t numElements = 100000000;
   mpi::SendBuffer sb0;
   mpi::SendBuffer sb1;
   mpi::SendBuffer sb2;
   Vector3<real_t> v(1,2,3);
   WcTimer timer0;
   WcTimer timer1;
   WcTimer timer2;

   for (size_t i = 0; i < numElements; ++i)
   {
      sb0 << v;
      sb1 << v;
      sb2 << v;
   }

   WALBERLA_LOG_DEVEL_VAR(sb0.size());
   sb0.clear();
   sb1.clear();
   sb2.clear();

   timer0.start();
   for (size_t i = 0; i < numElements; ++i)
   {
      sb0 << v;
   }
   timer0.end();

   WALBERLA_LOG_DEVEL_VAR(sb0.size());
   sb0.clear();

   timer1.start();
   for (size_t i = 0; i < numElements; ++i)
   {
      sb1 << v[0] << v[1] << v[2];
   }
   timer1.end();

   WALBERLA_LOG_DEVEL_VAR(sb0.size());
   sb0.clear();

   timer2.start();
   for (size_t i = 0; i < numElements; ++i)
   {
      auto pos = sb2.forward(sizeof(real_t) * 3);
      memcpy(pos, v.data(), sizeof(real_t) * 3);
   }
   timer2.end();

   WALBERLA_LOG_DEVEL_VAR(sb0.size());
   sb0.clear();

   //auto ptr0 = sb0.ptr();
   //auto ptr1 = sb1.ptr();
   //for (auto i = 0; i < numElements; ++i)
   //{
   //   WALBERLA_ASSERT_EQUAL(*ptr0, *ptr1);
   //   ++ptr0;
   //   ++ptr1;
   //}

   WALBERLA_LOG_DEVEL("native:      " << timer0.total());
   WALBERLA_LOG_DEVEL("elementwise: " << timer1.total());
   WALBERLA_LOG_DEVEL("memcpy:      " << timer2.total());

   return 0;
}

} // namespace walberla

int main( int argc, char* argv[] )
{
   walberla::mpi::Environment mpiEnv( argc, argv );
   WALBERLA_UNUSED(mpiEnv);

   return walberla::main( argc, argv );
}
