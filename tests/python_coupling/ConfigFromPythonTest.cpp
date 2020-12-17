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
//! \file ConfigFromPythonTest.cpp
//! \ingroup core
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "core/Abort.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "core/mpi/Environment.h"

#include "python_coupling/CreateConfig.h"

using namespace walberla;

int main(int argc, char** argv)
{
   debug::enterTestMode();

   mpi::Environment env(argc, argv);

   int counter = 0;
   for (auto cfg = python_coupling::configBegin(argc, argv); cfg != python_coupling::configEnd(); ++cfg)
   {
      auto config                  = *cfg;
      auto parameters              = config->getOneBlock("DomainSetup");
      const int test_int           = parameters.getParameter< int >("testInt");
      const std::string testString = parameters.getParameter< std::string >("testString");
      const real_t testDouble      = parameters.getParameter< real_t >("testDouble");
      Vector3< real_t > testVector = parameters.getParameter< Vector3< real_t > >("testVector");
      const bool testBool          = parameters.getParameter< bool >("testBool");

      if (counter == 0)
         WALBERLA_CHECK(test_int == 4)
      else
         WALBERLA_CHECK(test_int == 5)

      counter++;

      WALBERLA_CHECK(testString == "someString")
      WALBERLA_CHECK(testDouble > 42 && testDouble < 43)
      WALBERLA_CHECK(testVector == Vector3< real_t >(0.5, 0.5, 0.7))
      WALBERLA_CHECK(testBool == false)

      WALBERLA_LOG_INFO_ON_ROOT(test_int)
   }

   return EXIT_SUCCESS;
}
