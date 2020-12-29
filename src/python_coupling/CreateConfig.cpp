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
//! \file CreateConfigFromPythonScript.cpp
//! \ingroup python
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

#include "CreateConfig.h"

#include "core/StringUtility.h"
#include "core/config/Config.h"
#include "core/config/Create.h"
#include "core/logging/Logging.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON
#   include "python_coupling/helper/ConfigFromDict.h"

#   include "PythonCallback.h"
#   include <pybind11/pybind11.h>

namespace walberla
{
namespace python_coupling
{
namespace py = pybind11;

shared_ptr< Config > createConfigFromPythonScript(const std::string& scriptFile,
                                                         const std::string& pythonFunctionName,
                                                         const std::vector< std::string >& argv)
{
   importModuleOrFile(scriptFile, argv);

   PythonCallback pythonCallback(pythonFunctionName);

   pythonCallback();

   using py::dict;
   using py::object;

   object returnValue = pythonCallback.data().dict()["returnValue"];
   if (returnValue.is(object())) return shared_ptr< Config >();

   bool isDict = py::isinstance< dict >(returnValue);
   if (!isDict) { WALBERLA_ABORT("Python configuration did not return a dictionary object."); }
   dict returnDict = dict(returnValue);
   return configFromPythonDict(returnDict);
}

//===================================================================================================================
//
//  Config Generators and iterators
//
//===================================================================================================================

class PythonMultipleConfigGenerator : public config::ConfigGenerator
{
 public:
   PythonMultipleConfigGenerator(py::object ScenarioConfigGenerator) // NOLINT
      : ScenarioConfigGenerator_(ScenarioConfigGenerator)            // NOLINT
   {}

   shared_ptr< Config > next() override
   {
      shared_ptr< Config > config = make_shared< Config >();
      try
      {
         py::dict configDict = ScenarioConfigGenerator_.attr("__next__")();
         configFromPythonDict(config->getWritableGlobalBlock(), configDict);
         return config;
      }
      catch (py::error_already_set&)
      {
         return shared_ptr<Config>();
      }
   }

 private:
   py::object ScenarioConfigGenerator_;
};

class PythonSingleConfigGenerator : public config::ConfigGenerator
{
 public:
   PythonSingleConfigGenerator(const shared_ptr< Config >& config) : config_(config) {}

   shared_ptr< Config > next() override
   {
      auto res = config_;
      config_.reset();
      return res;
   }

 private:
   shared_ptr< Config > config_;
};

config::Iterator createConfigIteratorFromPythonScript(const std::string& scriptFile,
                                                             const std::string& pythonFunctionName,
                                                             const std::vector< std::string >& argv)
{
   importModuleOrFile(scriptFile, argv);
   PythonCallback pythonCallback(pythonFunctionName);
   pythonCallback();

   py::object returnValue = pythonCallback.data().dict()["returnValue"];
   bool isDict            = py::isinstance< py::dict >(returnValue);

   shared_ptr< config::ConfigGenerator > generator;
   if (isDict)
   {
      auto config            = make_shared< Config >();
      py::dict extractedDict = py::cast< py::dict >(returnValue);
      configFromPythonDict(config->getWritableGlobalBlock(), extractedDict);
      generator = make_shared< PythonSingleConfigGenerator >(config);
   }
   else
   {
      try
      {
         generator = make_shared< PythonMultipleConfigGenerator >(returnValue);
      } catch (py::error_already_set&)
      {
         std::string message = std::string("Error while running Python function ") + pythonFunctionName;
         WALBERLA_ABORT_NO_DEBUG_INFO(message);
      }
   }

   return config::Iterator(generator);
}

} // namespace python_coupling
} // namespace walberla

#else

namespace walberla
{
namespace python_coupling
{
shared_ptr< Config > createConfigFromPythonScript(const std::string&, const std::string&,
                                                  const std::vector< std::string >&)
{
   WALBERLA_ABORT("Tried to run with Python config but waLBerla was built without Python support.");
   return shared_ptr< Config >();
}

config::Iterator createConfigIteratorFromPythonScript(const std::string&, const std::string&,
                                                      const std::vector< std::string >&)
{
   WALBERLA_ABORT("Tried to run with Python config but waLBerla was built without Python support.");
   return config::Iterator();
}

} // namespace python_coupling
} // namespace walberla

#endif

namespace walberla
{
namespace python_coupling
{
shared_ptr< Config > createConfig(int argc, char** argv)
{
   if (argc < 2) throw std::runtime_error(config::usageString(argv[0]));

   shared_ptr< Config > config;
   std::string filename(argv[1]);

   auto argVec = std::vector< std::string >(argv + 1, argv + argc);

   if (string_ends_with(filename, ".py")) { config = createConfigFromPythonScript(filename, "config", argVec); }
   else
   {
      config = make_shared< Config >();
      config::createFromTextFile(*config, filename);
   }

   config::substituteCommandLineArgs(*config, argc, argv);

   return config;
}

config::Iterator configBegin(int argc, char** argv)
{
   if (argc < 2) throw std::runtime_error(config::usageString(argv[0]));

   std::string filename(argv[1]);
   if (string_ends_with(filename, ".py"))
   {
      auto argVec = std::vector< std::string >(argv + 1, argv + argc);
      return createConfigIteratorFromPythonScript(filename, "config", argVec);
   }
   else
   {
      return config::begin(argc, argv);
   }
}

} // namespace python_coupling
} // namespace walberla
