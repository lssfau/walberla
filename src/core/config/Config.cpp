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
//! \file Config.cpp
//! \ingroup core
//! \author Klaus Iglberger
//! \brief Implementation of a parameter file parser.
//
//======================================================================================================================

#include "Config.h"
#include "core/mpi/MPIManager.h"
#include "core/StringUtility.h"

#include <fstream>
#include <iostream>
#include <stdexcept>

namespace walberla {
namespace config {



//======================================================================================================================
//
//  CONSTRUCTORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config::Config()
// \brief Default constructor for the Config class.
 */
Config::Config()
     : stateFlag_(true) // Internal status flag
     , error_()         // Container for all error messages
     , block_()         // Global parameter block
{}
//**********************************************************************************************************************



//======================================================================================================================
//
//  DESTRUCTOR
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config::~Config()
// \brief Destructor for the Config class.
 */
Config::~Config()
= default;
//**********************************************************************************************************************




//======================================================================================================================
//
//  OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config& Config::operator=( const Config& c )
// \brief Copy assignment operator for the Config class.
//
// \param c The Config object to be copied.
 */
Config& Config::operator=( const Config& c )
{
   if( &c == this ) return *this;

   stateFlag_ = c.stateFlag_;
   error_.str() = c.error_.str();
   block_ = c.block_;
   valueReplacements_ = c.valueReplacements_;

   return *this;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn void Config::readParameterFile( const char* const filename )
// \brief Extraction of the parameter file \p filename.
//
// \param filename Name of the parameter file.
// \return void
 */
void Config::readParameterFile( const char* const filename )
{

   WALBERLA_ROOT_SECTION(){
      parseFromFile( filename, block_, 1 );
   }

   WALBERLA_MPI_SECTION() {
      std::string fileContent="";
      int size=0;
      WALBERLA_ROOT_SECTION() {
         fileContent=block_.getString();
         size=int(fileContent.size());
         //std::cout<<"Sending contents of parameter file to other processes... ("<<size<<" bytes)"<<std::endl;
         }
      MPI_Bcast(&size,1,MPITrait<int>::type(),0,MPI_COMM_WORLD);
      char* buffer = new char[uint_c(size)];
      WALBERLA_ROOT_SECTION() {
         fileContent.copy(buffer,uint_c(size));
      }
      MPI_Bcast(buffer,size,MPITrait<char>::type(),0,MPI_COMM_WORLD);
      WALBERLA_NON_ROOT_SECTION() {
         std::string content(buffer,uint_c(size));
         parseFromString(content,block_,1);
      }
      delete[] buffer;
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::parseFromFile( const char* filename, Block& block, unsigned int level )
// \brief Parses the given file and extracts the parameters to the given block.
//
// \param filename Name of the parameter file.
// \param block Parameter block for the found parameters.
// \param level The current inclusion level.
 */
void Config::parseFromFile( const char* filename, Block& block, unsigned int level )
{
   std::stringstream input;
   LineVector lineNumbers;
   std::string line;
   std::string key;
   std::string value;
   std::string::size_type pos1;
   std::string::size_type pos2;
   unsigned int lineCounter(0);
   bool comment(false);


   /////////////////////////////////////////////
   // Preprocessing the parameter file stream

   std::ifstream in( filename, std::ifstream::in );
   if( !in.is_open() )
   {
      error_ << "   Error opening parameter input file '" << filename << "' !\n";
      stateFlag_ = false;
      return;
   }

   line.reserve( 100 );
   key.reserve( 100 );
   value.reserve( 100 );

   while( std::getline( in, line ) )
   {
      ++lineCounter;

      // Filtering comments
      if( comment ) {
         if( ( pos1 = line.find( "*/", 0 ) ) != std::string::npos ) {
            line.erase( line.begin(), line.begin()+ numeric_cast< std::string::difference_type >( pos1+2 ) );
            comment = false;
         }
         else continue;
      }

      if( ( pos1 = line.find( "//", 0 ) ) != std::string::npos ) {
         line.erase( line.begin() + numeric_cast< std::string::difference_type >( pos1 ), line.end() );
      }
      if( ( pos1 = line.find( "/*", 0 ) ) != std::string::npos ) {
         if( ( pos2 = line.find( "*/", pos1+2 ) ) != std::string::npos ) {
            line.replace( line.begin() + numeric_cast< std::string::difference_type >( pos1 ),
                          line.begin() + numeric_cast< std::string::difference_type >( pos2+2 ), " " );
         }
         else {
            line.erase( line.begin() + numeric_cast< std::string::difference_type >( pos1 ), line.end() );
            comment = true;
         }
      }

      //Adding whitespaces
      for( pos1=0; ; ++pos1 )
      {
         if( pos1 >= line.size() ) break;

         if( line[pos1] == '{' ) {
            line.replace( pos1, 1, " { " );
            pos1 += 2;
         }
         else if( line[pos1] == '}' ) {
            line.replace( pos1, 1, " } " );
            pos1 += 2;
         }
         else if( line[pos1] == ';' ) {
            line.replace( pos1, 1, " ; " );
            pos1 += 2;
         }
      }

      //Adding the line to the input string
      lineNumbers.push_back( Pair( input.tellp(), lineCounter ) );
      input << line << "\n";
   }

   in.close();


   //////////////////////////
   // Parameter extraction

   while( input >> key )
   {
      unsigned int commandLine = getLineNumber( lineNumbers, input.tellg() );

      if( input >> std::ws && input.eof() ) break;

      //Extraction of a parameter block
      else if( input.peek() == '{' )
      {
         input.ignore( 1 );
         Block& newBlock = block.createBlock( key );
         extractBlock( filename, input, newBlock, lineNumbers,
                       getLineNumber( lineNumbers, input.tellg() ), level );
      }

      else if( string_icompare( key, "include" ) == 0 )
      {
         if( std::getline( input, value, ';' ) && !input.eof() )
         {
            input.ignore( 1 );

            if( level == maxInclusionLevel ) {
               error_ << "   " << filename << ", line " << commandLine
                      << ": Include limit reached! Include directives are limited to "
                      << maxInclusionLevel << " levels!\n"
                      << "     Inclusion of file '" << value << "' is skipped!\n";
               stateFlag_ = false;
            }
            else {
               removeTrailingWhiteSpaces( value );
               std::string file( (filesystem::path( getDirectory( filename ) ) / value).string() );
               parseFromFile( file.c_str(), block, level+1 );
            }
         }
      }

      //Extraction of a single parameter
      else if( std::getline( input, value, ';' ) && !input.eof() )
      {
         input.ignore( 1 );
         while( (value.find("$(") != std::string::npos) && (value.find(')') != std::string::npos) ) {
            size_t s = value.find("$("); size_t e = value.find(')');
            ValueReplacementMap::iterator mkey = valueReplacements_.find( value.substr( s+2, e-s+1-3 ) );
            if(mkey != valueReplacements_.end()) {
               value.replace( s,e-s+1, mkey->second );
            }
            else {
               if(e!=std::string::npos)
                  value.erase(e,1);
               if(s!=std::string::npos)
                  value.erase(s,2);
            }
         }
         if( !block.addParameter( key, value.substr( 0, value.size()-1 ) ) )
         {
            error_ << "   " << filename << ", line " << commandLine
                   << ": Duplicate parameter '" << key << "' in";
            if( block.getKey().empty() ) error_ << " global section!\n";
            else error_ << " block '" << block.getKey() << "'\n";
            stateFlag_ = false;
         }
      }
   }
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn void Config::parseFromString( const std::string str, Block& block, unsigned int level )
// \brief Parses the given file and extracts the parameters to the given block.
//
// \param str Config-generated input string.
// \param block Parameter block for the found parameters.
// \param level The current inclusion level.
//
// By internal convention, each block first contains all parameter, and then all sub-blocks
 */
void Config::parseFromString( const std::string & str, Block& block, unsigned int level )
{
   std::string input = str;

   while (!input.empty()) {
      std::string::size_type f=input.find_first_of(" {");
      if (input[f]==' ') {
         std::string::size_type posSemicolon = input.substr(f+1).find(';');
         block.addParameter(input.substr(0,f),input.substr(f+1,posSemicolon));
         input = input.substr(f+1+posSemicolon+1);
      } else {
         int lev=1;
         bool inParam=false;
         std::string::size_type curly2=std::string::npos;
         for (std::string::size_type pos = f+1; pos<input.size(); ++pos) {
            if (input[pos]==' ' && !inParam) inParam=true;
            if (input[pos]==';' && inParam) inParam=false;
            if (input[pos]=='{' && !inParam) ++lev;
            if (input[pos]=='}' && !inParam) --lev;
            if (lev==0) {
               curly2=pos;
               break;
            }
         }
         Block& newBlock = block.createBlock( input.substr(0,f) );
         parseFromString(input.substr(f+1,curly2-f-1),newBlock,level+1);
         input=input.substr(curly2+1);
      }
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::extractBlock( const char* filename, std::stringstream& input, Block& block,
//                                const LineVector& lineNumbers, unsigned int lineNumber,
//                                unsigned int level )
// \brief Extracting a parameter block from the preprocessed input stream.
//
// \param filename Name of the parameter file
// \param input The preprocessed input stream.
// \param block The block whose parameters are extracted.
// \param lineNumbers The line numbers belonging to the preprocessed input stream.
// \param lineNumber Starting line of the block in the parameter file.
// \param level The current inclusion level.
 */
void Config::extractBlock( const char* filename, std::stringstream& input, Block& block,
                           const LineVector& lineNumbers, unsigned int lineNumber,
                           unsigned int level )
{
   std::string key;
   std::string value;

   while( input >> key )
   {
      if( string_icompare( key, "}" ) == 0 ) {
         return;
      }
      else if( input >> std::ws && input.eof() ) {
         break;
      }

      unsigned int commandLine = getLineNumber( lineNumbers, input.tellg() );

      //Extraction of a parameter block
      if( input.peek() == '{' )
      {
         input.ignore( 1 );
         Block& newBlock = block.createBlock( key );
         extractBlock( filename, input, newBlock, lineNumbers,
               getLineNumber( lineNumbers, input.tellg() ), level );
      }

      else if( string_icompare( key, "include" ) == 0 )
      {
         if( std::getline( input, value, ';' ) && !input.eof() )
         {
            input.ignore( 1 );

            if( level == maxInclusionLevel ) {
               error_ << "   " << filename << ", line " << commandLine
                  << ": Include limit reached! Include directives are limited to "
                  << maxInclusionLevel << " levels!\n"
                  << "     Inclusion of file '" << value << "' is skipped!\n";
               stateFlag_ = false;
            }
            else {
               removeTrailingWhiteSpaces( value );
               std::string file( (filesystem::path( getDirectory( filename ) ) / value).string() );
               parseFromFile( file.c_str(), block, level+1 );
            }
         }
      }

      //Extraction of a single parameter
      else if( std::getline( input, value, ';' ) && !input.eof() )
      {
         input.ignore( 1 );
         while( (value.find("$(") != std::string::npos) && (value.find(')') != std::string::npos) ) {
            size_t s = value.find("$("); size_t e = value.find(')');
            ValueReplacementMap::iterator mkey = valueReplacements_.find( value.substr( s+2, e-s+1-3 ) );
            if(mkey != valueReplacements_.end()) {
               value.replace( s,e-s+1, mkey->second );
            }
            else {
               if(e!=std::string::npos)
                  value.erase(e,1);
               if(s!=std::string::npos)
                  value.erase(s,2);
            }
         }
         if( !block.addParameter( key, value.substr( 0, value.size()-1 ) ) ) {
            error_ << "   " << filename << ", line " << commandLine
               << ": Duplicate parameter '" << key << "' in block '" << block.getKey() << "'!\n";
            stateFlag_ = false;
         }
      }
   }

   error_ << "   Missing '}' for " << block.getKey()
      << " block starting in line " << lineNumber << "\n";
   stateFlag_ = false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn unsigned int Config::getLineNumber( const LineVector& lineNumbers, sstreamPos pos ) const
// \brief Calculation of the line number of a command.
//
// \param lineNumbers ?
// \param pos Position in the preprocessed input stream.
 */
unsigned int Config::getLineNumber( const LineVector& lineNumbers, sstreamPos pos ) const
{
   LineVector::const_iterator it=lineNumbers.begin();
   for( ; it!=lineNumbers.end()-1; ++it ) {
      if( (it+1)->first > pos )
         return it->second;
   }
   return it->second;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  CLASS FILEREADER::BLOCK
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config::Block::Block( const std::string& key )
// \brief Standard constructor for the Block class.
//
// \param key The key of the block (default argument = "").
 */
Config::Block::Block( const std::string& key )
   :key_(key),blocks_(),params_()
{}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Block& Config::Block::operator=( const Block& b )
// \brief The copy assignment operator for the Block class.
//
// \param b The Block object to be copied.
// \return Reference to the Config object.
 */
Config::Block& Config::Block::operator=( const Block& b )
{
   if( &b == this ) return *this;

   key_ = b.key_;
   blocks_ = b.blocks_;
   params_ = b.params_;

   return *this;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Block::~Block()
// \brief Default destructor of the Block class.
 */
Config::Block::~Block()
= default;
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn Config::size_type Config::Block::getNumBlocks( const std::string& key ) const
// \brief Returns the number of contained parameter blocks corresponding to given key.
//
// \param key The key of the blocks to count.
//
// \return The number of contained parameter blocks corresponding to given key.
 */
Config::size_type Config::Block::getNumBlocks( const std::string& key ) const
{
   size_type c = 0;
   for( List::const_iterator it=blocks_.begin(); it!=blocks_.end(); ++it ) {
      if( string_icompare( key, it->getKey() ) == 0 ) {
         ++c;
      }
   }
   return c;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::BlockHandle Config::Block::getBlock( const std::string& key ) const
// \brief Returns one or none block that matches with \p key.
//
// \param key The key of the extracted block.
// \return One or none block that matches with \p key.
 */
Config::BlockHandle Config::Block::getBlock( const std::string& key ) const {

   Blocks blocks;
   getBlocks( key, blocks, 0, 1 );
   if( blocks.empty() )
      return BlockHandle();
   return blocks[0];
}

//**********************************************************************************************************************
/*!\fn Config::BlockHandle Config::Block::getOneBlock( const std::string& key ) const
// \brief Returns exactly one block that matches with \p key.
//
// \param key The key of the extracted block.
// \return One block that matches with \p key.
 */
Config::BlockHandle Config::Block::getOneBlock( const std::string& key ) const
{
   Blocks blocks;
   getBlocks( key, blocks, 1, 1 );
   return blocks[0];
}

//**********************************************************************************************************************
/*!\fn void Config::Block::getBlocks( const std::string& key, Blocks& blocks,
//                                    size_t min, size_t max ) const
// \brief Adds to the given \p blocks all extracted blocks with key \p key.
//
// \param key The key of the extracted blocks.
// \param[out] blocks Reference to the vector of blocks
// \param min minimum number of blocks the caller expects
// \param max maximum number of blocks the caller can handle
// \return void.
 */
void Config::Block::getBlocks( const std::string& key, Blocks& blocks, size_t min, size_t max ) const
{
   size_t c = 0;
   for( List::const_iterator it=blocks_.begin(); it!=blocks_.end(); ++it ) {
      if( string_icompare( key, it->getKey() ) == 0 ) {
         blocks.push_back( BlockHandle( &*it ) );
         ++c;
      }
   }
   if( c < min || c > max ) {
      std::ostringstream oss;
      oss << "You requested at least " << min << " and at most " << max << " Config::Blocks with key "
          << key << " - but there are " << c << ".";
      throw std::range_error( oss.str() );
   }
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn void Config::Block::getBlocks( Blocks& blocks ) const
// \brief Adds to the given \p blocks all blocks contained in this block.
//
// \param[out] blocks Reference to the vector of blocks
// \return void.
 */
void Config::Block::getBlocks( Blocks& blocks ) const
{
   for( List::const_iterator it=blocks_.begin(); it!=blocks_.end(); ++it ) {
      blocks.push_back( BlockHandle( &*it ) );
   }
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn void Config::Block::getWritableBlocks( std::vector<Block*>& blocks )
// \brief Adds to the given \p blocks all blocks contained in this block.
//
// \param[out] blocks Reference to the vector of blocks
// \return void.
 */
void Config::Block::getWritableBlocks( std::vector<Block*> & blocks )
{
   for( List::iterator it=blocks_.begin(); it!=blocks_.end(); ++it ) {
      blocks.push_back( &*it );
   }
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn void Config::Block::getWritableBlocks( const std::string& key,std::vector<Block*>& blocks,size_t min,size_t max )
// \brief Adds to the given \p blocks all extracted blocks with key \p key.
//
// \param key The key of the extracted blocks.
// \param[out] blocks Reference to the vector of blocks
// \param min minimum number of blocks the caller expects
// \param max maximum number of blocks the caller can handle
// \return void.
 */
void Config::Block::getWritableBlocks( const std::string & key, std::vector<Block*> & blocks, size_t min, size_t max )
{
   size_t c = 0;
   for( List::iterator it=blocks_.begin(); it!=blocks_.end(); ++it ) {
      if( string_icompare( key, it->getKey() ) == 0 ) {
         blocks.push_back(  &*it );
         ++c;
      }
   }
   if( c < min || c > max ) {
      std::ostringstream oss;
      oss << "You requested at least " << min << " and at most " << max << " Config::Blocks with key "
          << key << " - but there are " << c << ".";
      throw std::range_error( oss.str() );
   }
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  GET FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn bool Config::Block::addParameter( std::string key, const std::string& value )
// \brief Adding a parameter to a parameter block.
//
// \param key The key of the new parameter.
// \param value The value of the new parameter.
// \return \a true if a new parameter was added, \a false if duplicate key is present.
 */
bool Config::Block::addParameter( const std::string& key, const std::string& value )
{
   return params_.insert( std::pair<Key,Value>(key,value) ).second;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Block& Config::Block::createBlock( std::string key )
// \brief Adding a new block within the parameter block.
//
// \param key The key of the new block.
// \return Reference to the new block.
 */
Config::Block& Config::Block::createBlock( const std::string& key )
{
   blocks_.push_back( Block( key ) );
   return *blocks_.rbegin();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::Block::listParameters() const
// \brief Output function for the contained parameters.
//
// \return void
 */
void Config::Block::listParameters() const
{
   for( Map::const_iterator it=params_.begin(); it!=params_.end(); ++it ) {
      std::cout << " Key = '" << it->first << "' , Value = '" << it->second << "'\n";
   }
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn std::string Config::Block::getString() const
// \brief Returns Config-internal string to use for communication.
//
// \return void
 */
std::string Config::Block::getString() const
{
   std::string ret="";
   for( Map::const_iterator it=params_.begin(); it!=params_.end(); ++it ) {
      ret += it->first + " " + it->second + ";";
   }
   for( List::const_iterator it=blocks_.begin(); it!=blocks_.end(); ++it ) {
      ret += it->key_ + "{" + it->getString() + "}";
   }
   return ret;
}
//**********************************************************************************************************************





//======================================================================================================================
//
//  OUTPUT
//
//======================================================================================================================
static void printConfig( std::ostream & os, const Config::BlockHandle & block, int depth = 0  )
{
   std::vector< Config::BlockHandle > subBlocks;
   block.getBlocks( subBlocks );

   std::stringstream prefix;
   for( int i=0; i < depth; ++i )
      prefix << "   ";


   if ( depth >=0 )
   {
      std::string blockName = block.getKey();
      string_trim( blockName );
      os << prefix.str() << blockName << "\n";
      os << prefix.str() << "{\n";
   }

   for( auto subBlock = subBlocks.begin(); subBlock != subBlocks.end(); ++subBlock )
      printConfig( os, *subBlock, depth+1);

   for( auto paramIt = block.begin(); paramIt != block.end(); ++paramIt ) {
      std::string key = paramIt->first;
      std::string value = paramIt->second;
      string_trim( key );
      string_trim( value );
      os << prefix.str() << "   " << key << " " << value << ";\n";
   }


   if ( depth >=0 )
      os << prefix.str() << "}\n";
}


std::ostream & operator<< ( std::ostream & os, const Config & config )
{
   printConfig( os, config.getGlobalBlock(), -1 );
   return os;
}

std::ostream & operator<< ( std::ostream & os, const Config::BlockHandle & block )
{
   printConfig( os, block );
   return os;
}



} // namespace config
} // namespace walberla
