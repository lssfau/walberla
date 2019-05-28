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
//! \file   SingleCast.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/data/DataTypes.h>

#include <mesa_pd/kernel/SingleCast.h>

#include <core/Environment.h>
#include <core/logging/Logging.h>

#include <iostream>

namespace walberla {

using namespace walberla::mesa_pd;

class SingleParticleAccessor : public data::IAccessor
{
public:
   data::BaseShape* getShape(const size_t p_idx)
   {
      return p_idx == 0 ?
               static_cast<data::BaseShape*>(&halfspace_) :
               static_cast<data::BaseShape*>(&sphere_);
   }

   data::HalfSpace halfspace_ = data::HalfSpace(Vec3(1,0,0));
   data::Sphere    sphere_    = data::Sphere(real_t(0.5));
};

class ShapeTester
{
public:
   void operator()(const size_t idx, data::HalfSpace& /*hs*/)
   {
      WALBERLA_CHECK_EQUAL(idx, 0);
   }

   void operator()(const size_t idx, data::Sphere& /*s*/)
   {
      WALBERLA_CHECK_EQUAL(idx, 1);
   }
};

int main( int argc, char ** argv )
{
   Environment env(argc, argv);
   WALBERLA_UNUSED(env);
   mpi::MPIManager::instance()->useWorldComm();

   //init data structures
   SingleParticleAccessor accessor;

   //init kernels
   kernel::SingleCast singlecast;

   ShapeTester tester;
   singlecast(0, accessor, tester);
   singlecast(1, accessor, tester);

   return EXIT_SUCCESS;
}

} //namespace walberla

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}
