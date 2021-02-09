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
//! \file Singleton.h
//! \ingroup core
//! \author Klaus Iglberger
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/NonCopyable.h"

#include <mutex>


/// \cond internal

namespace walberla {
namespace singleton {



//======================================================================================================================
//
//  BEFRIEND_SINGLETON MACRO
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Friendship declaration for the Singleton class template.
//
// This macro has to be used in order to declare the Singleton functionality as friend of the
// class deriving from Singleton.
*/
#define WALBERLA_BEFRIEND_SINGLETON \
   template< typename > friend class walberla::singleton::Singleton
//**********************************************************************************************************************



//======================================================================================================================
//
//  CLASS SINGLETON
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T >  // Type of the singleton
class Singleton : private NonCopyable
{
protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
   = default;
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   = default;
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
      // this implementation is thread safe
      // https://stackoverflow.com/questions/1661529/is-meyers-implementation-of-the-singleton-pattern-thread-safe
      static shared_ptr<T> object( new T() );
      isInstantiated_ = true;
      return object;
   }
   //@}
   //*******************************************************************************************************************

   static bool isInstantiated() { return isInstantiated_; }

private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static std::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//======================================================================================================================

template< typename T >
std::mutex Singleton<T>::instanceMutex_;

template< typename T >
bool Singleton<T>::isInstantiated_ = false;



} // namespace singleton
} // namespace walberla

/// \endcond
