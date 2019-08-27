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
//! \file Config.h
//! \ingroup core
//! \author Klaus Iglberger
//
//======================================================================================================================

#pragma once

#include "core/Abort.h"
#include "core/Filesystem.h"
#include "core/debug/CheckFunctions.h"
#include "core/StringUtility.h"

#include <cctype>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>






namespace walberla {
namespace config {



struct CaseInsensitiveCompare
{
   bool operator()( const std::string & lhs, const std::string & rhs) const
   {
      return ( string_icompare(lhs, rhs) < 0 );
   }
};



//======================================================================================================================
//
//  CLASS DEFINITION
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Configuration handler & parameter file parser.
// \author Klaus Iglberger
//
// The Config class parses parameter files and stores value-parameter pairs. These files may
// contain an arbitrary number of parameters in any order of the following form:\n
// key value;\n
// Additionally, the files may contain an arbitrary number of parameter blocks of the following
// form:\n
// key { ... },\n
// which can again contain an arbitrary number of parameters and parameter blocks in any order.\n
// During the parameter extraction, the file parser checks for duplicate parameters within a
// single parameter block (which also includes the global block). Furthermore, the file parser is
// case insensitive with respect to the used keys. Therefore 'key', 'KEY' and 'Key' are handled
// equally and are recognized as duplicate parameters.
 */
class Config
{
public:
   //**Declarations*****************************************************************************************************
   class Block;
   class BlockHandle;
   //*******************************************************************************************************************

private:
   //**Type definitions****************************************************************************
   typedef std::string                                  Key;                  //!< Parameter key.
   typedef std::string                                  Value;                //!< Parameter value.
   typedef std::map<Key,Value, CaseInsensitiveCompare>  Map;                  //!< Parameter map.
   typedef std::list<Block>                             List;                 //!< List for parameter blocks.
   typedef std::map<Value,Value>                        ValueReplacementMap;  //!< Parameter value replacement map

   typedef std::stringstream::pos_type         sstreamPos;  //!< Stream position.
   typedef std::pair<sstreamPos,unsigned int>  Pair;        //!< Pair consisting of a stream position
                                                            //!< and a line number.
   typedef std::vector<Pair>                   LineVector;  //!< Vector for Pair.
   //*******************************************************************************************************************

   //**Member constants*************************************************************************************************
   static const unsigned int maxInclusionLevel = 30;
   //*******************************************************************************************************************

public:
   //**Type definitions****************************************************************************
   typedef std::vector<BlockHandle>  Blocks;          //!< Container for block handles.
   typedef List::size_type           size_type;       //!< Size type for a Block count.

   typedef Map::iterator             iterator;        //!< Iterator over the contained parameters.
   typedef Map::const_iterator       const_iterator;  //!< Constant iterator over the contained parameters.
   //*******************************************************************************************************************

   //**Error codes******************************************************************************************************
   enum ErrorCode {
      noerror   = 0,
      undefined = 1,
      badcast   = 2
   };
   //*******************************************************************************************************************

   //**class Parameter**************************************************************************************************
   template< typename Type >
   class Parameter
   {
   public:
      //**Constructors**************************************************************************************************
      explicit Parameter():parameter_(),error_(noerror) {}
      explicit Parameter(Type parameter, ErrorCode err, std::string key )
         : parameter_(std::move(parameter)), error_(err), key_(std::move(key)) {}
      // No explicitly declared copy constructor.
      //****************************************************************************************************************

      //**Destructor****************************************************************************************************
      // No explicitly declared destructor.
      //****************************************************************************************************************

      //**Operators****************************************************************************************************
      // No explicitly declared copy assignment operator.
      inline operator Type();
      //****************************************************************************************************************

      //**Utility functions*********************************************************************************************
      inline ErrorCode error() const { return error_; }
      //****************************************************************************************************************

   private:
      //**Member variables**********************************************************************************************
      Type parameter_;
      ErrorCode error_;
      std::string key_;
      //****************************************************************************************************************
   };
   //*******************************************************************************************************************

   //**class Block******************************************************************************************************
   /*!\brief Parameter block.
   // \author Klaus Iglberger
   //
   // A parameter block contains all parameters of one section of the input parameter file. A
   // single parameter block can again contain other parameter blocks.
    */
   class Block
   {
   private:
      //**Friend declarations*******************************************************************************************
      friend class Config;
      //****************************************************************************************************************

   public:
      //**Constructors**************************************************************************************************
      Block( const std::string& key="" );
      Block( const Block& b ) = default;
      //****************************************************************************************************************

      //**Destructor****************************************************************************************************
      ~Block();
      //****************************************************************************************************************

      //**Operators*****************************************************************************************************
      /*! \name Operators */
      //@{
      Block&              operator= ( const Block& b );
      inline std::string& operator[]( std::string key );
      //@}
      //****************************************************************************************************************

      //**Get functions*************************************************************************************************
      /*! \name Get functions */
      //@{
      inline const std::string& getKey() const;
      inline bool isDefined( std::string key ) const;

      template< typename T >
      inline Parameter<T> getParameter( std::string key ) const;

      template< typename T >
      inline Parameter<T> getParameter( const std::string & key, const T & defaultValue ) const;

      inline bool setParameter( const std::string & key, const std::string & value );
      
      inline void setOrAddParameter( const std::string & key, const std::string & value );

      inline iterator       begin();
      inline const_iterator begin() const;
      inline iterator       end();
      inline const_iterator end() const;

      inline size_type getNumBlocks() const;
             size_type getNumBlocks( const std::string& key ) const;
      BlockHandle getBlock( const std::string& key ) const;
      BlockHandle getOneBlock( const std::string& key ) const;
      void getBlocks( const std::string& key, Blocks& blocks,
                      size_t min = 0, size_t max = std::numeric_limits< size_t >::max() ) const;
      void getBlocks( Blocks& blocks ) const;
      void getWritableBlocks    ( std::vector<Block*> & blocks );
      void getWritableBlocks    ( const std::string & key, std::vector<Block*> & blocks,
                                  size_t min = 0, size_t max = std::numeric_limits< size_t >::max() );

      //@}
      //****************************************************************************************************************

      //**Utility functions*********************************************************************************************
      /*! \name Utility functions */
      //@{
      bool addParameter( const std::string& key, const std::string& value );
      void listParameters() const;
      //@}
      //****************************************************************************************************************

      Block& createBlock( const std::string& key );

   private:
      //**Utility functions*********************************************************************************************
      /*! \name Utility functions */
      //@{
      std::string getString() const;
      //@}
      //****************************************************************************************************************

      //**Member variables**********************************************************************************************
      /*! \name Member variables */
      //@{
      Key key_;      //!< Key of the block.
      List blocks_;  //!< List of contained blocks.
      Map params_;   //!< The parameters of the block.
      //@}
      //****************************************************************************************************************
   };
   //*******************************************************************************************************************

   //**class BlockHandle************************************************************************************************
   /*!\brief Handle for a Block object.
   // \author Klaus Iglberger
   //
   // The BlockHandle class is a simple wrapper for a pointer to a Block object that facilitates
   // the handling of a Block pointer. The class is used transparently in the Blocks container
   // type.
    */
   class BlockHandle
   {
   public:
      //**Constructors**************************************************************************************************
      inline BlockHandle();
      inline BlockHandle( const Block* block );
      // No explicitly declared copy constructor.
      //****************************************************************************************************************

      //**Get functions*************************************************************************************************
      /*! \name Get functions */
      //@{
      inline const std::string& getKey() const;
      inline bool isDefined( std::string key ) const;

      template< typename T >
      inline Parameter<T> getParameter( const std::string & key ) const;

      template< typename T >
      inline Parameter<T> getParameter( const std::string & key, const T& defaultValue ) const;

      bool isValid() const { return block_ != NULL; }
      operator bool() const { return isValid(); }

      inline const_iterator begin() const;
      inline const_iterator end() const;

      inline size_type getNumBlocks() const;
      inline size_type getNumBlocks( const std::string& key ) const;
      BlockHandle getOneBlock( const std::string& key ) const { return block_->getOneBlock( key ); };
      BlockHandle getBlock( const std::string& key ) const { return block_->getBlock( key ); }
      inline void getBlocks( const std::string& key, Blocks& blocks,
                             size_t min = 0, size_t max = std::numeric_limits< size_t >::max() ) const;
      inline void getBlocks( Blocks& blocks ) const;

      inline Block cloneBlock() const;
      //@}
      //****************************************************************************************************************

      //**Utility functions*********************************************************************************************
      /*! \name Utility functions */
      //@{
      inline void listParameters() const;
      //@}
      //****************************************************************************************************************

   private:
      //**Member variables**********************************************************************************************
      const Block * block_;  //!< Pointer to a const parameter block.
      //****************************************************************************************************************
   };
   //*******************************************************************************************************************

   //**Constructors*****************************************************************************************************
   Config( );
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   ~Config();
   //*******************************************************************************************************************

   //**Operators********************************************************************************************************
   /*! \name Operators */
   //@{
   Config& operator=( const Config& c );
   inline  operator bool() const;
   //@}
   //*******************************************************************************************************************

   //**Error functions**************************************************************************************************
   /*! \name Error functions */
   //@{
   inline void clear();
   inline std::string error() const;
   //@}
   //*******************************************************************************************************************

   //**Get functions****************************************************************************************************
   /*! \name Get functions */
   //@{

   inline bool isDefined( std::string key ) const;

   template< typename T >
   inline Parameter<T> getParameter( const std::string& key ) const;

   template< typename T >
   inline Parameter<T> getParameter( const std::string& key, const T & defaultValue ) const;


   inline iterator       begin();
   inline const_iterator begin() const;
   inline iterator       end();
   inline const_iterator end()   const;

   inline size_type getNumBlocks() const { return block_.getNumBlocks(); }
   inline size_type getNumBlocks( const std::string& key ) const { return block_.getNumBlocks( key ); }
   BlockHandle getBlock( const std::string& key ) const { return block_.getBlock( key ); }
   BlockHandle getOneBlock( const std::string& key ) const { return block_.getOneBlock( key ); }
   inline void getBlocks( const std::string& key, Blocks& blocks,
                          size_t min = 0, size_t max = std::numeric_limits< size_t >::max() ) const;
   inline void getBlocks( Blocks& blocks ) const;
   inline BlockHandle getGlobalBlock() const;
   inline Block&      getWritableGlobalBlock();
   //@}
   //*******************************************************************************************************************


   //**Utility functions************************************************************************************************
   /*! \name Utility functions */
   //@{
   void readParameterFile( const char* const filename );
   inline void listParameters() const;
   static inline void convertToLowerCase( std::string& s );
   void addValueReplacement( const Value& original, const Value& replacement ) { valueReplacements_[original] = replacement; }
   //@}
   //*******************************************************************************************************************

private:
   //**Utility functions************************************************************************************************
   /*! \name Utility functions */
   //@{
   static inline std::string  getDirectory( const std::string& path );
   static inline std::string  getFilename( const std::string& path );
   inline void  removeTrailingWhiteSpaces( std::string& s );
   unsigned int getLineNumber( const LineVector& lineNumbers, sstreamPos pos ) const;
   void         parseFromFile( const char* const filename, Block& block, unsigned int level );
   void         parseFromString( const std::string & str, Block& block, unsigned int level );
   void         extractBlock( const char* filename, std::stringstream& input, Block& block,
                              const LineVector& lineNumbers, unsigned int lineNumber, unsigned int level );
   //@}
   //*******************************************************************************************************************

   //**Member variables*************************************************************************************************
   /*! \name Member variables */
   //@{
   bool stateFlag_;             //!< Internal status of the config object.
                                /*!< An error is indicated by \p false. */
   std::ostringstream error_;   //!< Container for all error messages.
   Block block_;                //!< The global parameter block.

   ValueReplacementMap valueReplacements_; //!< value replacements
                                           /*!< "$(mapElements.first)" in input file is replaced
                                            *   by "mapElement.second" */
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************


//======================================================================================================================
//
//  OUTPUT
//
//======================================================================================================================

std::ostream & operator<< ( std::ostream & os, const Config & config );
std::ostream & operator<< ( std::ostream & os, const Config::BlockHandle & blockHandle );


//======================================================================================================================
//
//  INLINE OPERATORS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config::operator bool() const
// \brief Conversion operator for the Config class.
//
// This conversion operator allows for the check of a Config object. The operator returns
// the internal status of the config object.
 */
inline Config::operator bool() const
{
   return stateFlag_;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  INLINE ERROR FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn void Config::clear()
// \brief Clearing the internal status of the file config object.
//
// \return void
//
// This function resets the internal status of the config object and clears all error messages.
 */
inline void Config::clear()
{
   stateFlag_ = true;
   error_.str().clear();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::string Config::error() const
// \brief Returns the error messages.
//
// \return Error messages.
 */
inline std::string Config::error() const
{
   return error_.str();
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  INLINE GET FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn bool Config::isDefined( std::string key ) const
// \brief Checks if a parameter was defined in the parameter file.
//
// \param key The parameter key to be checked.
// \return \a true if the parameter was defined, \a false if the parameter wasn't defined.
 */
inline bool Config::isDefined( std::string key ) const
{
   return block_.isDefined( key );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Parameter<Type> Config::getParameter( const std::string& key ) const
// \brief Returns an extracted parameter.
//
// \param key The key of the extracted parameter.
// \return The extracted parameter.
 */
template< typename Type >
inline Config::Parameter<Type> Config::getParameter( const std::string & key ) const
{
   return block_.getParameter<Type>( key );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Parameter<T> Config::getParameter( const std::string& key, const T& ) const
// \brief Returns an extracted parameter.
//
// \param key The key of the extracted parameter.
// \param def default value
// \return The extracted parameter if key is defined, otherwise default value
 */
template< typename T >
inline Config::Parameter<T> Config::getParameter( const std::string & key, const T & def ) const
{
   if ( isDefined( key ) )
      return block_.getParameter<T>( key );
   else
      return Config::Parameter<T> (def, noerror, key );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::iterator Config::begin()
// \brief Returns an iterator to the first globally defined parameter.
//
// \return Iterator to the first globally defined parameter.
 */
inline Config::iterator Config::begin()
{
   return block_.begin();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::const_iterator Config::begin() const
// \brief Returns an iterator to the first globally defined parameter.
//
// \return Iterator to the first globally defined parameter.
 */
inline Config::const_iterator Config::begin() const
{
   return block_.begin();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::iterator Config::end()
// \brief Returns an iterator one past the last globally defined parameter.
//
// \return Iterator one past the last globally defined parameter.
 */
inline Config::iterator Config::end()
{
   return block_.end();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::const_iterator Config::end() const
// \brief Returns an iterator one past the last globally defined parameter.
//
// \return Iterator one past the last globally defined parameter.
 */
inline Config::const_iterator Config::end() const
{
   return block_.end();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::getBlocks( const std::string& key, Blocks& blocks,
//                             size_t min, size_t max ) const
// \brief Adds to the given \p blocks all extracted blocks with key \p key.
//
// \param key The key of the extracted blocks.
// \param [out] blocks Reference to the vector of blocks.
// \param min minimum number of blocks the caller expects
// \param max maximum number of blocks the caller can handle
// \return void
 */
inline void Config::getBlocks( const std::string& key, Blocks& blocks, size_t min, size_t max ) const
{
   block_.getBlocks( key, blocks, min, max );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::getBlocks( Blocks& blocks ) const
// \brief Adds to the given \p blocks all globally defined blocks.
//
// \param [out] blocks Reference to the vector of blocks.
// \return void
 */
inline void Config::getBlocks( Blocks& blocks ) const
{
   block_.getBlocks( blocks );
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Config::BlockHandle Config::getGlobalBlock() const
// \brief gives a handle to the global block.
//
// \return Handle to global block
 */
inline Config::BlockHandle Config::getGlobalBlock() const
{
   return BlockHandle( &block_ );
}
//*******************************************************

//**********************************************************************************************************************
/*!\fn Config::BlockHandle Config::getWritableGlobalBlock()
// \brief returns the global block.
//
// \return Global block
 */
inline Config::Block& Config::getWritableGlobalBlock()
{
   return block_;
}
//*******************************************************

//======================================================================================================================
//
//  INLINE UTILITY FUNCTIONS
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn void Config::listParameters() const
// \brief Output function for the contained parameters.
//
// \return void
 */
void Config::listParameters() const
{
   block_.listParameters();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::convertToLowerCase( std::string& s )
// \brief Conversion of the character in \p s to lower case characters.
//
// \param s Reference to the string that gets converted.
// \return void
 */
inline void Config::convertToLowerCase( std::string& s )
{
   for( std::string::size_type i=0; i<s.size(); ++i ) {
      s[i] = (char)std::tolower( s[i] );
   }
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::string Config::getDirectory( const std::string& path )
// \brief Returns the directory part of the given path.
//
// \param path The path containing both directory and filename.
// \return The directory part of the given path.
 */
inline std::string Config::getDirectory( const std::string& path )
{
   filesystem::path p( path );
   return p.parent_path().string();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn std::string Config::getFilename( const std::string& path )
// \brief Returns the filename part of the given path.
//
// \param path The path containing both directory and filename.
// \return The filename part of the given path.
 */
inline std::string Config::getFilename( const std::string& path )
{
   filesystem::path p( path );
   return p.filename().string();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::removeTrailingWhiteSpaces( std::string& s )
// \brief Removes the trailing white spaces from the given string.
//
// \param s The string to be cleaned from trailing white spaces.
// \return void
 */
inline void Config::removeTrailingWhiteSpaces( std::string& s )
{
   if( s.empty() )
      return;

   std::string::iterator it=s.end();

   while( it!=s.begin() ) {
      if( isalnum( *(it-1) ) ) {
         s.erase( it, s.end() );
         break;
      }
      --it;
   }
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  CLASS CONFIG::Parameter
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config::Parameter<Type>::operator Type()
// \brief Converts a Parameter to its value's type
//
// If the parameter's value could not be cast to Type or it has not been found an error is logged
// and the program is aborted.
 */
template< typename Type >
Config::Parameter<Type>::operator Type()
{
   switch( error_ )
   {
   case noerror:
      break;
   case undefined:
      WALBERLA_ABORT( "The parameter \"" << key_ << "\" is undefined!" );
   case badcast:
      WALBERLA_ABORT( "The parameter \"" << key_ << "\" could not be cast to the desired type!" );
   default:
      WALBERLA_ABORT( "Unknown error while getting parameter \"" << key_ << "\"" );
   }

   return parameter_;
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  CLASS FILEREADER::BLOCK
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn std::string& Config::Block::operator[]( std::string key )
// \brief Returns a reference to the value of the specified key.
//
// \param key The key of the parameter.
// \return Reference to the value of the specified key.
//
// If the key has not been defined before, a new parameter is added. The return value is a
// reference to the value of the (new) key.
 */
inline std::string& Config::Block::operator[]( std::string key )
{
   return params_[key];
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const std::string& Config::Block::getKey() const
// \brief Returns the key of the block.
//
// \return The block key.
 */
inline const std::string& Config::Block::getKey() const
{
   return key_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Config::Block::isDefined( std::string key ) const
// \brief Checks if a parameter was defined in the parameter file.
//
// \return \a true if the parameter was defined, \a false if the parameter wasn't defined.
 */
inline bool Config::Block::isDefined( std::string key ) const
{
   return !( params_.find( key ) == params_.end() );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::size_type Config::Block::getNumBlocks() const
// \brief Returns the number of contained parameter blocks.
//
// \return The number of contained parameter blocks.
 */
inline Config::size_type Config::Block::getNumBlocks() const
{
   return blocks_.size();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Parameter<Type> Config::Block::getParameter( std::string key ) const
// \brief Returns an extracted parameter.
//
// \param key The key of the parameter.
// \return The extracted parameter.
 */
template< typename Type >
inline Config::Parameter<Type> Config::Block::getParameter( std::string key ) const
{
   Map::const_iterator it = params_.find( key );

   if( it != params_.end() ) {
      Type tmp;
      std::istringstream iss( it->second );
      if( !(iss >> tmp) ) return Parameter<Type>( Type(), badcast, key );
      else return Parameter<Type>( tmp, noerror, key );
   }
   else return Parameter<Type>( Type(), undefined, key );
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Config::Parameter<T> Config::Block::getParameter( const std::string& key, const T& def ) const
// \brief Returns an extracted parameter.
//
// \param key The key of the extracted parameter.
// \param def The default value, if the key does not exist
// \return The extracted parameter if key is defined, otherwise default value
 */
template< typename T >
inline Config::Parameter<T> Config::Block::getParameter( const std::string& key, const T & def ) const
{
   if( isDefined(key) )
      return getParameter<T>( key );
   else
      return Config::Parameter<T> (def, noerror, key );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn Config::Parameter<std::string> Config::Block::getParameter<std::string>( std::string key ) const
// \brief Returns an extracted string parameter.
//
// \param key The key of the string parameter.
// \return The extracted string parameter.
 */
template<>
inline Config::Parameter<std::string> Config::Block::getParameter<std::string>( std::string key ) const
{
   Map::const_iterator it = params_.find( key );

   if( it != params_.end() )
      return Parameter<std::string>( it->second, noerror, key );
   else
      return Parameter<std::string>( "", undefined, key );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn Config::Parameter<bool> Config::Block::getParameter<bool>( std::string key ) const
// \brief Returns an extracted bool parameter.
//
// \param key The key of the bool parameter.
// \return The extracted bool parameter.
 */
template<>
inline Config::Parameter<bool> Config::Block::getParameter<bool>( std::string key ) const
{
   Map::const_iterator it = params_.find( key );

   if( it != params_.end() )
   {
      if( ( it->second.length() == 1u && it->second[0] == '1' )   ||
          string_icompare( it->second, std::string("true") ) == 0 ||
          string_icompare( it->second, std::string("yes") ) == 0  ||
          string_icompare( it->second, std::string("on") ) == 0 )
      {
         return Parameter<bool>( true, noerror, key );
      }
      else if( ( it->second.length() == 1u && it->second[0] == '0' )    ||
               string_icompare( it->second, std::string("false") ) == 0 ||
               string_icompare( it->second, std::string("no") ) == 0    ||
               string_icompare( it->second, std::string("off") ) ==0 )
      {
         return Parameter<bool>( false, noerror, key );
      }
      else
      {
         return Parameter<bool>( false, badcast, key );
      }
   }
   else
      return Parameter<bool>( false, undefined, key );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn Config::Parameter<unsigned int> Config::Block::getParameter<unsigned int>( std::string key ) const
// \brief Returns an extracted unsigned integer parameter.
//
// \param key The key of the unsigned integer parameter.
// \return The extracted unsigned integer parameter.
 */
template<>
inline Config::Parameter<unsigned int> Config::Block::getParameter<unsigned int>( std::string key ) const
{
   Map::const_iterator it = params_.find( key );

   if( it != params_.end() ) {
      int tmp;
      std::istringstream iss( it->second );
      if( !(iss >> tmp) || tmp < 0 )
         return Parameter<unsigned int>( 0, badcast, key );
      else return Parameter<unsigned int>( static_cast<unsigned int>( tmp ), noerror, key );
   }
   else return Parameter<unsigned int>( 0, undefined, key );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Config::Block::setParameter( const std::string& key, const std::string& value )
// \brief Sets a given parameter to a given value.
//
// \param key The key of the parameter.
// \param value The value of the parameter.
// \return True for success, False for failure.
 */
inline bool Config::Block::setParameter( const std::string& key, const std::string& value )
{
   if( isDefined(key) )
   {
      params_[key] = value;
      return true;
   }
   else
      return false;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::Block::setOrAddParameter( const std::string& key, const std::string& value )
// \brief Sets a given parameter to a given value (adds the parameter if it does not yet exist).
//
// \param key The key of the parameter.
// \param value The value of the parameter.
 */
inline void Config::Block::setOrAddParameter( const std::string& key, const std::string& value )
{
   params_[key] = value;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::iterator Config::Block::begin()
// \brief Returns an iterator to the first parameter of the block.
//
// \return Iterator to the first parameter of the block.
 */
inline Config::iterator Config::Block::begin()
{
   return params_.begin();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::const_iterator Config::Block::begin() const
// \brief Returns an iterator to the first parameter of the block.
//
// \return Iterator to the first parameter of the block.
 */
inline Config::const_iterator Config::Block::begin() const
{
   return params_.begin();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::iterator Config::Block::end()
// \brief Returns an iterator one past the last parameter of the block.
//
// \return Iterator one past the last parameter of the block.
 */
inline Config::iterator Config::Block::end()
{
   return params_.end();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::const_iterator Config::Block::end() const
// \brief Returns an iterator one past the last parameter of the block.
//
// \return Iterator one past the last parameter of the block.
 */
inline Config::const_iterator Config::Block::end() const
{
   return params_.end();
}
//**********************************************************************************************************************




//======================================================================================================================
//
//  CLASS CONFIG::BLOCKHANDLE
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\fn Config::BlockHandle::BlockHandle()
// \brief Default constructor for a block handle.
//
// \b Note: The BlockHandle object is uninitialized. Any access to the bound block will result in
// an error until the BlockHandle object is bound to a block.
 */
inline Config::BlockHandle::BlockHandle()
   :block_(0)
{}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::BlockHandle::BlockHandle( Block* block )
// \brief Standard constructor for a block handle.
//
// \param block Pointer to the block to be wrapped by the BlockHandle object.
 */
inline Config::BlockHandle::BlockHandle( const Block* block )
   :block_(block)
{}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn const std::string& Config::BlockHandle::getKey() const
// \brief Returns the key of the block.
//
// \return The block key.
 */
inline const std::string& Config::BlockHandle::getKey() const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get key of a nullpointer block");
   return block_->getKey();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn bool Config::BlockHandle::isDefined( std::string key ) const
// \brief Checks if a parameter was defined in the parameter file.
//
// \return \a true if the parameter was defined, \a false if the parameter wasn't defined.
 */
inline bool Config::BlockHandle::isDefined( std::string key ) const
{
   return block_ && block_->isDefined( key );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::Parameter<Type> Config::BlockHandle::getParameter( const std::string& key ) const
// \brief Returns an extracted parameter.
//
// \param key The key of the extracted parameter.
// \return The extracted parameter.
 */
template< typename Type >
inline Config::Parameter<Type> Config::BlockHandle::getParameter( const std::string& key ) const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get a parameter out of a nullpointer block");
   return block_->getParameter<Type>( key );
}
//**********************************************************************************************************************

//**********************************************************************************************************************
/*!\fn Config::Parameter<T> Config::BlockHandle::getParameter( const std::string& key, const T & ) const
// \brief Returns an extracted parameter.
//
// \param key The key of the extracted parameter.
// \param def default value, used if key is not defined
// \return The extracted parameter if key is defined, default value otherwise
 */
template< typename T >
inline Config::Parameter<T> Config::BlockHandle::getParameter( const std::string& key, const T& def ) const
{
   if( isDefined(key) )
      return block_->getParameter<T>( key );
   else
      return Config::Parameter<T> (def, noerror, key );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn Config::const_iterator Config::BlockHandle::begin() const
// \brief Returns an iterator to the first parameter of the block.
//
// \return Iterator to the first parameter of the block.
 */
inline Config::const_iterator Config::BlockHandle::begin() const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get the begin iterator of a nullpointer block");
   return block_->begin();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::const_iterator Config::BlockHandle::end() const
// \brief Returns an iterator one past the last parameter of the block.
//
// \return Iterator one past the last parameter of the block.
 */
inline Config::const_iterator Config::BlockHandle::end() const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get the end iterator of a nullpointer block");
   return block_->end();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::size_type Config::BlockHandle::getNumBlocks() const
// \brief Returns the number of contained parameter blocks.
//
// \return The number of contained parameter blocks.
 */
inline Config::size_type Config::BlockHandle::getNumBlocks() const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get the number of blocks out of a nullpointer block");
   return block_->getNumBlocks();
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn Config::size_type Config::BlockHandle::getNumBlocks( const std::string& key ) const
// \brief Returns the number of contained parameter blocks corresponding to given key.
//
// \param key The key of the blocks to count.
//
// \return The number of contained parameter blocks corresponding to given key.
 */
inline Config::size_type Config::BlockHandle::getNumBlocks( const std::string& key ) const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get the number of blocks with key out of a nullpointer block");
   return block_->getNumBlocks( key );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::BlockHandle::getBlocks( const std::string& key, Blocks& blocks,
//                                          size_t min, size_t max ) const
// \brief Adds to the given \p blocks all extracted blocks with key \p key.
//
// \param key The key of the extracted blocks.
// \param[out] blocks Reference to the vector of blocks
// \param min minimum number of blocks the caller expects
// \param max maximum number of blocks the caller can handle
// \return void.
 */
inline void Config::BlockHandle::getBlocks( const std::string& key, Blocks& blocks, size_t min, size_t max ) const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get blocks by string from a nullpointer block");
   block_->getBlocks( key, blocks, min, max );
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::BlockHandle::getBlocks( Blocks& blocks ) const
// \brief Adds to the given \p blocks all blocks inside the handles block
//
// \param[out] blocks Reference to the vector of blocks
// \return void.
 */
inline void Config::BlockHandle::getBlocks( Blocks& blocks ) const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to get all blocks out of a nullpointer block");
   block_->getBlocks( blocks );
}
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\fn Block Config::BlockHandle::cloneBlock() const
// \brief Creates a copy of the bound block.
//
// \return The copy of the bound block.
 */
inline Config::Block Config::BlockHandle::cloneBlock() const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to clone nullpointer block");
   return *block_;
}
//**********************************************************************************************************************


//**********************************************************************************************************************
/*!\fn void Config::BlockHandle::listParameters() const
// \brief Output function for the contained parameters.
//
// \return void
 */
inline void Config::BlockHandle::listParameters() const
{
   WALBERLA_CHECK_NOT_NULLPTR(block_, "Tried to list parameters of a nullpointer block");
   block_->listParameters();
}
//**********************************************************************************************************************



} // namespace config

using config::Config;

} // namespace walberla
