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
//
//======================================================================================================================

#include "CreateConfig.h"


#include "core/config/Config.h"
#include "core/config/Create.h"
#include "core/logging/Logging.h"
#include "core/StringUtility.h"

#include <exception>


#ifdef WALBERLA_BUILD_WITH_PYTHON

#include "PythonCallback.h"
#include "PythonWrapper.h"
#include "DictWrapper.h"

#include "python_coupling/helper/ConfigFromDict.h"
#include "helper/ExceptionHandling.h"


namespace walberla {
namespace python_coupling {

   namespace bp = boost::python;


   shared_ptr<Config> createConfigFromPythonScript( const std::string & scriptFile,
                                                    const std::string & pythonFunctionName,
                                                    const std::vector<std::string> & argv )
   {
      importModuleOrFile( scriptFile, argv );

      PythonCallback pythonCallback ( pythonFunctionName );

      pythonCallback();

      using boost::python::object;
      using boost::python::dict;
      using boost::python::extract;

      object returnValue = pythonCallback.data().dict()[ "returnValue" ];
      if ( returnValue == object() )
         return shared_ptr<Config>();

      bool isDict = extract< dict >( returnValue ).check();
      if ( ! isDict ) {
         WALBERLA_ABORT("Python configuration did not return a dictionary object.");
      }
      dict returnDict = extract<dict>( returnValue );
      return configFromPythonDict( returnDict );
   }


   //===================================================================================================================
   //
   //  Config Generators and iterators
   //
   //===================================================================================================================


   class PythonMultipleConfigGenerator : public config::ConfigGenerator
   {
   public:
      PythonMultipleConfigGenerator( bp::stl_input_iterator< bp::dict > iterator  )  //NOLINT
         : iter_( iterator ), firstTime_(true) //NOLINT
      {}

      shared_ptr<Config> next() override
      {
         // this seemingly unnecessary complicated firstTime variable is used
         // since in alternative version where (++iter_) is at the end of the function
         // the python generator expression for the next time is already
         // called before the current simulation finished
         if ( !firstTime_ )
            ++iter_;
         else
            firstTime_ = false;

         if ( iter_ == bp::stl_input_iterator< bp::dict >()  )
            return shared_ptr<Config>();

         shared_ptr<Config> config = make_shared<Config>();

         bp::dict configDict = *iter_;
         configFromPythonDict( config->getWritableGlobalBlock(), configDict );

         return config;
      }

   private:
      bp::stl_input_iterator< bp::dict > iter_;
      bool firstTime_;
   };



   class PythonSingleConfigGenerator : public config::ConfigGenerator
   {
   public:
      PythonSingleConfigGenerator( const shared_ptr<Config> & config ): config_ ( config ) {}

      shared_ptr<Config> next() override
      {
         auto res = config_;
         config_.reset();
         return res;
      }

   private:
      shared_ptr<Config> config_;
   };


   config::Iterator createConfigIteratorFromPythonScript( const std::string & scriptFile,
                                                          const std::string & pythonFunctionName,
                                                          const std::vector<std::string> & argv )
   {
      importModuleOrFile( scriptFile, argv );

      PythonCallback pythonCallback ( pythonFunctionName );

      pythonCallback();

      bp::object returnValue = pythonCallback.data().dict()[ "returnValue" ];

      bool isDict = bp::extract< bp::dict >( returnValue ).check();

      shared_ptr< config::ConfigGenerator> generator;
      if ( isDict )
      {
         auto config = make_shared<Config>();
         bp::dict extractedDict = bp::extract<bp::dict> ( returnValue );
         configFromPythonDict( config->getWritableGlobalBlock(), extractedDict );
         generator = make_shared<PythonSingleConfigGenerator>( config );
      }
      else {

         try {
            generator= make_shared<PythonMultipleConfigGenerator>( returnValue );
         }
         catch ( bp::error_already_set & ) {
            python_coupling::terminateOnPythonException("Error while running Python config generator");
         }
      }

      return config::Iterator( generator );
   }


} // namespace python_coupling
} // namespace walberla


#else


namespace walberla {
namespace python_coupling {

   shared_ptr<Config> createConfigFromPythonScript( const std::string &, const std::string &, const std::vector<std::string> &  )
   {
      WALBERLA_ABORT( "Tried to run with Python config but waLBerla was built without Python support." );
      return shared_ptr<Config>();
   }


   config::Iterator createConfigIteratorFromPythonScript( const std::string & , const std::string &, const std::vector<std::string> &   )
   {
      WALBERLA_ABORT( "Tried to run with Python config but waLBerla was built without Python support." );
      return config::Iterator();
   }


} // namespace python_coupling
} // namespace walberla


#endif












namespace walberla {
namespace python_coupling {



   shared_ptr<Config> createConfig( int argc, char ** argv )
   {
      if(argc<2)
         throw std::runtime_error( config::usageString(argv[0]) );

      shared_ptr<Config> config;
      std::string filename( argv[1] );

      auto argVec = std::vector<std::string> (argv+1, argv + argc);

      if ( string_ends_with( filename, ".py")  ) {
         config = createConfigFromPythonScript( filename, "config", argVec );
      }
      else {
         config = make_shared<Config>();
         config::createFromTextFile( *config, filename );
      }

      config::substituteCommandLineArgs( *config, argc, argv );

      return config;
   }

   config::Iterator configBegin( int argc, char ** argv )
   {
      if(argc<2)
         throw std::runtime_error( config::usageString(argv[0]) );

      std::string filename( argv[1] );
      if ( string_ends_with( filename, ".py")  ) {
         auto argVec = std::vector<std::string> (argv+1, argv + argc);
         return createConfigIteratorFromPythonScript( filename, "config", argVec );
      }
      else {
         return config::begin( argc, argv );
      }

   }





} // namespace python_coupling
} // namespace walberla

