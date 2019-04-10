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
//! \file Create.cpp
//! \ingroup core
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "Create.h"
#include "core/logging/Tracing.h"
#include "core/Abort.h"
#include "core/StringUtility.h"

#include <iomanip>

namespace walberla {
namespace config {

   /********************************************************************************************************************
    * Returns usage help, is called when program was started with wrong parameters
    *
    * \param executableName  name of the executable, usually taken from argv[0]
    *******************************************************************************************************************/
   std::string usageString(const std::string & executableName)
   {
      std::stringstream ss;
      ss << "Wrong usage" << std::endl;
      ss << "Usage:" << executableName << " ParameterFile" << std::endl;
      ss << "One parameter file is expected!" << std::endl;
      return ss.str();
   }


   shared_ptr<Config> create(int argc, char ** argv)
   {
      WALBERLA_TRACE_IN;

      shared_ptr<Config> config = make_shared<Config>();

      if(argc<2)
         throw std::runtime_error( usageString(argv[0]) );

      // parse command line parameters and parameter file name
      std::string filename;
      for (int i = 1; i < argc; ++i)
      {
         if ( ( i < argc - 1 ) && ( argv[i][0] == '-' ) ) {
            config->addValueReplacement( &argv[i][1], argv[i+1] );
            ++i;
         }
         else {
            if (! filename.empty() ) {
               WALBERLA_LOG_DETAIL( "Ignoring parameter " << argv[i] );
            }
            else
               filename = std::string( argv[i] );
         }
      }

      if ( filename.empty() )
         throw std::runtime_error( usageString( argv[0] ) );


      createFromTextFile( *config, filename );
      substituteCommandLineArgs( *config, argc, argv );
      return config;
   }


   void createFromTextFile( Config & config, const std::string & pathToTextFile )
   {
      config.readParameterFile( pathToTextFile.c_str() );
      if (config.error() != "" )
         throw std::runtime_error( "FileReader returned an error reading the parameter file: \n" + config.error() );
   }


   void substituteCommandLineArgs( Config & config, int argc, char**argv )
   {
      std::vector< std::string > params;
      params.reserve( uint_c( argc - 1 ) );

      for(int i=1; i < argc; ++i ) {
         std::string curArg ( argv[i] );
         if ( curArg[0] == '-' )
            params.push_back( curArg.substr( 1) );
      }
      substituteCommandLineArgs( config, params );
   }

   void substituteCommandLineArgs( Config & config, const std::vector<std::string> & params )
   {
      using std::vector;
      using std::string;

      for( auto param = params.begin(); param != params.end(); ++param )
      {
         vector< string > equalitySignSplitResult = string_split( *param, "=" );

         std::string value;
         if ( equalitySignSplitResult.size() == 0 )
         {
            WALBERLA_LOG_WARNING( "Ignoring empty parameter");
            continue;
         }
         else if ( equalitySignSplitResult.size() == 1 ) {
            value = "1";
         }
         else if ( equalitySignSplitResult.size() == 2 ) {
            value = equalitySignSplitResult[1];
         }
         else
         {
            WALBERLA_LOG_WARNING( "Ignoring illegally formed command line parameter: '" << *param << "'  (Multiple '='s )");
            continue;
         }

         const std::string & blockDescriptor = equalitySignSplitResult[0];

         vector< string > blocks = string_split( blockDescriptor, "." );

         if ( blocks.empty() ) {
            WALBERLA_LOG_WARNING( "Ignoring Parameter: Missing block descriptor on left hand side: '" << *param <<"'" );
            continue;
         }

         Config::Block * currentBlock = & config.getWritableGlobalBlock();
         for( uint_t i=0; i < blocks.size() -1; ++i )
         {
            std::string & blockName = blocks[i];
            string_trim( blockName );

            if ( blockName.empty() )
            {
               currentBlock = nullptr;
               WALBERLA_LOG_WARNING("Ignoring Parameter '" << *param << "' empty block name");
               break;
            }
            vector< Config::Block * > possibleBlocks;
            currentBlock->getWritableBlocks( blockName, possibleBlocks );
            if ( possibleBlocks.size() > 1 )
            {
               currentBlock = nullptr;
               WALBERLA_LOG_WARNING("Ignoring Parameter '" << *param << "' since block is ambiguous: " << blockName );
               break;
            }
            else if ( possibleBlocks.empty() ) {
               currentBlock = & ( currentBlock->createBlock( blockName ) );
            }
            else {
               WALBERLA_ASSERT_EQUAL( possibleBlocks.size(), 1 );
               currentBlock = possibleBlocks[0];
            }
         }
         if ( ! currentBlock )
            continue;

         string_trim( blocks.back() );
         if ( blocks.back().empty() )
         {
            WALBERLA_LOG_WARNING( "Ignoring Parameter '" << *param << "' since key is empty");
            continue;
         }
         else
         {
            if ( currentBlock->isDefined( blocks.back() ))
               currentBlock->setParameter( blocks.back(), value );
            else
               currentBlock->addParameter( blocks.back(), value );
         }
      }
   }


   //===================================================================================================================
   //
   //  Config Iterators
   //
   //===================================================================================================================

   class SingleConfigGenerator : public config::ConfigGenerator
   {
   public:
      SingleConfigGenerator( const shared_ptr<Config> & config ): config_ ( config ) {}

      shared_ptr<Config> next() override
      {
         auto res = config_;
         config_.reset();
         return res;
      }

   private:
      shared_ptr<Config> config_;
   };


   class MultipleConfigGenerator : public config::ConfigGenerator
   {
   public:
      MultipleConfigGenerator( const std::string & baseName, const std::string & extension, int numberOfDigits )
         : baseName_( baseName ), extension_( extension ), numberOfDigits_( numberOfDigits), counter_(-1) {}

      shared_ptr<Config> next() override
      {
         ++counter_;
         std::stringstream ss;
         ss << baseName_ << std::setfill('0') << std::setw(numberOfDigits_) << counter_ << extension_;

         WALBERLA_LOG_PROGRESS( "Simulating " << ss.str()  );
         auto config = make_shared<Config>();
         createFromTextFile( *config, ss.str() );

         return config;
      }

   private:
      std::string baseName_;
      std::string extension_;
      int numberOfDigits_;
      int counter_;
   };


   Iterator begin( int argc, char ** argv)
   {
      if( argc<2 )
         throw std::runtime_error( usageString(argv[0]) );

      // parse command line parameters and parameter file name
      std::string filename;
      for( int i= argc-1; i >=0; --i ) {
         if ( argv[i][0] != '-' ) {
            filename = std::string (argv[i]);
            break;
         }
      }
      if ( filename.empty() )
          throw std::runtime_error( usageString(argv[0]) );


      auto dotPosition = filename.find_last_of('.');

      int numberOfZeros=0;
      if ( dotPosition != std::string::npos )
      {
         auto searchPosition = dotPosition -1;
         while( searchPosition > 0 && filename[searchPosition] == '0' ) {
            ++numberOfZeros;
            --searchPosition;
         }
      }

      if ( numberOfZeros > 0 )
      {
         std::string basename  = filename.substr( 0, dotPosition - uint_c( numberOfZeros ) );
         std::string extension = filename.substr( dotPosition );
         return createConfigIteratorFromTextFiles( basename, extension, numberOfZeros );
      }
      else
      {
         auto config = make_shared<Config>();
         createFromTextFile( *config, filename );
         substituteCommandLineArgs( *config, argc, argv );
         return config::Iterator( make_shared<SingleConfigGenerator> ( config ) );
      }
   }

   Iterator createConfigIteratorFromTextFiles ( const std::string & basename, const std::string & extension, int nrOfDigits )
   {
      // Intel compiler complains when casting shared_ptr<MultipleConfigGenerator> to shared_ptr<ConfigGenerator>
      ConfigGenerator * cg = new MultipleConfigGenerator( basename, extension, nrOfDigits );
      return Iterator( shared_ptr<config::ConfigGenerator>(cg) );
   }


} // namespace config
} // namespace walberla
