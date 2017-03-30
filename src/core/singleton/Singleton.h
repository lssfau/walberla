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

#include "NullType.h"
#include "TypeList.h"
#include "core/DataTypes.h"
#include "core/NonCopyable.h"

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
#include <boost/thread/mutex.hpp>
#endif


/// \cond internal

namespace walberla {
namespace singleton {



#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#pragma GCC diagnostic ignored "-Wpragmas"
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wunused-local-typedefs"
#endif



//**********************************************************************************************************************
/*!\brief Helper macro for macro concatenation.
//
// The following code was borrowed from the Boost C++ framework (www.boost.org). This piece of
// macro magic joins the two arguments together, even when one of the arguments is itself a
// macro (see 16.3.1 in C++ standard). The key is that macro expansion of macro arguments does
// not occur in WALBERLA_DO_JOIN2 but does in WALBERLA_DO_JOIN.
*/
#define WALBERLA_JOIN( X, Y ) WALBERLA_DO_JOIN( X, Y )
#define WALBERLA_DO_JOIN( X, Y ) WALBERLA_DO_JOIN2(X,Y)
#define WALBERLA_DO_JOIN2( X, Y ) X##Y
//**********************************************************************************************************************



//======================================================================================================================
//
//  NAMESPACE FORWARD DECLARATIONS
//
//======================================================================================================================

template< typename T, typename TL, bool C > struct HasCyclicDependency;



//======================================================================================================================
//
//  CLASS HASCYCLICDEPENDENCYHELPER
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Auxiliary helper struct for the HasCyclicDependency class template.
//
//
// Helper template class for the HasCyclicDependency template class to resolve all lifetime
// dependencies represented by means of a dependency type list.
 */
template< typename TL                      // Type list of checked lifetime dependencies
        , typename D                       // Type list of lifetime dependencies to check
        , size_t   N = Length<D>::value >  // Length of the dependency type list
struct HasCyclicDependencyHelper;
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the given
// dependency type list is empty. In this case no cyclic lifetime dependency could be detected.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , size_t   N >  // Length of the dependency type list
struct HasCyclicDependencyHelper<TL,NullType,N>
{
   enum { value = 0 };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 1.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,1>
{
   typedef typename TypeAt<D,0>::Result  D1;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 2.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,2>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 3.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,3>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;
   typedef typename TypeAt<D,2>::Result  D3;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value ||
                  HasCyclicDependency<D3,TL,Contains<TL,D3>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 4.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,4>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;
   typedef typename TypeAt<D,2>::Result  D3;
   typedef typename TypeAt<D,3>::Result  D4;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value ||
                  HasCyclicDependency<D3,TL,Contains<TL,D3>::value>::value ||
                  HasCyclicDependency<D4,TL,Contains<TL,D4>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 5.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,5>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;
   typedef typename TypeAt<D,2>::Result  D3;
   typedef typename TypeAt<D,3>::Result  D4;
   typedef typename TypeAt<D,4>::Result  D5;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value ||
                  HasCyclicDependency<D3,TL,Contains<TL,D3>::value>::value ||
                  HasCyclicDependency<D4,TL,Contains<TL,D4>::value>::value ||
                  HasCyclicDependency<D5,TL,Contains<TL,D5>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 6.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,6>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;
   typedef typename TypeAt<D,2>::Result  D3;
   typedef typename TypeAt<D,3>::Result  D4;
   typedef typename TypeAt<D,4>::Result  D5;
   typedef typename TypeAt<D,5>::Result  D6;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value ||
                  HasCyclicDependency<D3,TL,Contains<TL,D3>::value>::value ||
                  HasCyclicDependency<D4,TL,Contains<TL,D4>::value>::value ||
                  HasCyclicDependency<D5,TL,Contains<TL,D5>::value>::value ||
                  HasCyclicDependency<D6,TL,Contains<TL,D6>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 7.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,7>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;
   typedef typename TypeAt<D,2>::Result  D3;
   typedef typename TypeAt<D,3>::Result  D4;
   typedef typename TypeAt<D,4>::Result  D5;
   typedef typename TypeAt<D,5>::Result  D6;
   typedef typename TypeAt<D,6>::Result  D7;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value ||
                  HasCyclicDependency<D3,TL,Contains<TL,D3>::value>::value ||
                  HasCyclicDependency<D4,TL,Contains<TL,D4>::value>::value ||
                  HasCyclicDependency<D5,TL,Contains<TL,D5>::value>::value ||
                  HasCyclicDependency<D6,TL,Contains<TL,D6>::value>::value ||
                  HasCyclicDependency<D7,TL,Contains<TL,D7>::value>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependencyHelper class template.
//
// This specialization of the HasCyclicDependencyHelper class is selected in case the length
// of the given type list is 8.
 */
template< typename TL   // Type list of checked lifetime dependencies
        , typename D >  // Type list of lifetime dependencies to check
struct HasCyclicDependencyHelper<TL,D,8>
{
   typedef typename TypeAt<D,0>::Result  D1;
   typedef typename TypeAt<D,1>::Result  D2;
   typedef typename TypeAt<D,2>::Result  D3;
   typedef typename TypeAt<D,3>::Result  D4;
   typedef typename TypeAt<D,4>::Result  D5;
   typedef typename TypeAt<D,5>::Result  D6;
   typedef typename TypeAt<D,6>::Result  D7;
   typedef typename TypeAt<D,7>::Result  D8;

   enum { value = HasCyclicDependency<D1,TL,Contains<TL,D1>::value>::value ||
                  HasCyclicDependency<D2,TL,Contains<TL,D2>::value>::value ||
                  HasCyclicDependency<D3,TL,Contains<TL,D3>::value>::value ||
                  HasCyclicDependency<D4,TL,Contains<TL,D4>::value>::value ||
                  HasCyclicDependency<D5,TL,Contains<TL,D5>::value>::value ||
                  HasCyclicDependency<D6,TL,Contains<TL,D6>::value>::value ||
                  HasCyclicDependency<D7,TL,Contains<TL,D7>::value>::value ||
                  HasCyclicDependency<D8,TL,Contains<TL,D8>::value>::value };
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  CLASS HASCYCLICDEPENDENCY
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Class template for the detection of cyclic lifetime dependencies.
//
// This class template checks the given type \a T for cyclic lifetime dependencies. In case a
// cyclic lifetime dependency is detected, the \a value member enumeration is set to 1. Otherwise
// it is set to 0.
 */
template< typename T                      // The type to be checked for cyclic lifetime dependencies
        , typename TL                     // Type list of checked lifetime dependencies
        , bool C=Contains<TL,T>::value >  // Flag to indicate whether T is contained in TL
struct HasCyclicDependency
{
   typedef typename Append<TL,T>::Result  ETL;
   enum { value = HasCyclicDependencyHelper<ETL,typename T::Dependencies>::value };
};
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Specialization of the HasCyclicDependency class template.
//
// This specialization of the HasCyclicDependency class is selected in case the given type \a T
// is contained in the given lifetime dependency type list \a TL. In this case a cyclic lifetime
// dependency was detected and the \a value member enumeration is set to 1.
 */
template< typename T     // The type to be checked for cyclic lifetime dependencies
        , typename TL >  // Type list of checked lifetime dependencies
struct HasCyclicDependency<T,TL,true>
{
   enum { value = 1 };
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  CLASS CYCLIC_LIFETIME_DEPENDENCY_TEST
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Compile time constraint wrapper class.
//
// Helper class for the pe::CYCLIC_LIFETIME_DEPENDENCY_DETECTED class template. This class is
// used as a wrapper for the instantiation of the pe::CYCLIC_LIFETIME_DEPENDENCY_DETECTED
// constraint class. It serves the purpose to force the instantiation of either the defined
// specialization or the undefined basic template during the compilation. In case the compile
// time condition is met, the type pe::CYCLIC_LIFETIME_DEPENDENCY_TEST<1> is defined.
 */
template< int > struct CYCLIC_LIFETIME_DEPENDENCY_TEST {};
//**********************************************************************************************************************



//======================================================================================================================
//
//  DETECT_CYCLIC_LIFETIME_DEPENDENCY CONSTRAINT
//
//======================================================================================================================

//**********************************************************************************************************************
/*!\brief Compile time constraint.
//
// Helper template class for the compile time constraint enforcement. Based on the compile time
// constant expression used for the template instantiation, either the undefined basic template
// or the specialization is selected. If the undefined basic template is selected, a compilation
// error is created.
 */
template< bool > struct CYCLIC_LIFETIME_DEPENDENCY_DETECTED;
template<> struct CYCLIC_LIFETIME_DEPENDENCY_DETECTED<false> { enum { value = 1 }; };
//**********************************************************************************************************************



//**********************************************************************************************************************
/*!\brief Constraint on the data type.
//
// In case the given data type \a T is not an integral data type, a compilation error is created.
 */
#define WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY(T) \
      typedef \
      walberla::singleton::CYCLIC_LIFETIME_DEPENDENCY_TEST< \
      walberla::singleton::CYCLIC_LIFETIME_DEPENDENCY_DETECTED< walberla::singleton::HasCyclicDependency<T,walberla::singleton::NullType>::value >::value > \
      WALBERLA_JOIN( DETECT_CYCLIC_LIFETIME_DEPENDENCY_TYPEDEF, __LINE__ )
//**********************************************************************************************************************



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
   template< typename, typename, typename, typename, typename, typename,typename, typename, typename > friend class walberla::singleton::Singleton; \
   template< typename, typename, bool > friend struct walberla::singleton::HasCyclicDependency
//**********************************************************************************************************************


//======================================================================================================================
//
//  CLASS SINGLETON
//
//======================================================================================================================

//**********************************************************************************************************************
/// Base class for all lifetime managed singletons.
template< typename T                // Type of the singleton (CRTP pattern)
        , typename D1 = NullType    // Type of the first lifetime dependency
        , typename D2 = NullType    // Type of the second lifetime dependency
        , typename D3 = NullType    // Type of the third lifetime dependency
        , typename D4 = NullType    // Type of the fourth lifetime dependency
        , typename D5 = NullType    // Type of the fifth lifetime dependency
        , typename D6 = NullType    // Type of the sixth lifetime dependency
        , typename D7 = NullType    // Type of the seventh lifetime dependency
        , typename D8 = NullType >  // Type of the eighth lifetime dependency
class Singleton : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_8( D1, D2, D3, D4, D5, D6, D7, D8 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
        , dependency3_( D3::instance() )  // Handle to the third lifetime dependency
        , dependency4_( D4::instance() )  // Handle to the fourth lifetime dependency
        , dependency5_( D5::instance() )  // Handle to the fifth lifetime dependency
        , dependency6_( D6::instance() )  // Handle to the sixth lifetime dependency
        , dependency7_( D7::instance() )  // Handle to the seventh lifetime dependency
        , dependency8_( D8::instance() )  // Handle to the eighth lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D3 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D4 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D5 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D6 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D7 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D8 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.
   shared_ptr<D3> dependency3_;  //!< Handle to the third lifetime dependency.
   shared_ptr<D4> dependency4_;  //!< Handle to the fourth lifetime dependency.
   shared_ptr<D5> dependency5_;  //!< Handle to the fifth lifetime dependency.
   shared_ptr<D6> dependency6_;  //!< Handle to the sixth lifetime dependency.
   shared_ptr<D7> dependency7_;  //!< Handle to the seventh lifetime dependency.
   shared_ptr<D8> dependency8_;  //!< Handle to the eighth lifetime dependency.
   
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif
   
   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************


//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (7 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1    // Type of the first lifetime dependency
        , typename D2    // Type of the second lifetime dependency
        , typename D3    // Type of the third lifetime dependency
        , typename D4    // Type of the fourth lifetime dependency
        , typename D5    // Type of the fifth lifetime dependency
        , typename D6    // Type of the sixth lifetime dependency
        , typename D7 >  // Type of the eighth lifetime dependency
class Singleton<T,D1,D2,D3,D4,D5,D6,D7,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_7( D1, D2, D3, D4, D5, D6, D7 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
        , dependency3_( D3::instance() )  // Handle to the third lifetime dependency
        , dependency4_( D4::instance() )  // Handle to the fourth lifetime dependency
        , dependency5_( D5::instance() )  // Handle to the fifth lifetime dependency
        , dependency6_( D6::instance() )  // Handle to the sixth lifetime dependency
        , dependency7_( D7::instance() )  // Handle to the seventh lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D3 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D4 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D5 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D6 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D7 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.
   shared_ptr<D3> dependency3_;  //!< Handle to the third lifetime dependency.
   shared_ptr<D4> dependency4_;  //!< Handle to the fourth lifetime dependency.
   shared_ptr<D5> dependency5_;  //!< Handle to the fifth lifetime dependency.
   shared_ptr<D6> dependency6_;  //!< Handle to the sixth lifetime dependency.
   shared_ptr<D7> dependency7_;  //!< Handle to the seventh lifetime dependency.
   
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif
   
   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (6 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1    // Type of the first lifetime dependency
        , typename D2    // Type of the second lifetime dependency
        , typename D3    // Type of the third lifetime dependency
        , typename D4    // Type of the fourth lifetime dependency
        , typename D5    // Type of the fifth lifetime dependency
        , typename D6 >  // Type of the eighth lifetime dependency
class Singleton<T,D1,D2,D3,D4,D5,D6,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_6( D1, D2, D3, D4, D5, D6 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
        , dependency3_( D3::instance() )  // Handle to the third lifetime dependency
        , dependency4_( D4::instance() )  // Handle to the fourth lifetime dependency
        , dependency5_( D5::instance() )  // Handle to the fifth lifetime dependency
        , dependency6_( D6::instance() )  // Handle to the sixth lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D3 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D4 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D5 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D6 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.
   shared_ptr<D3> dependency3_;  //!< Handle to the third lifetime dependency.
   shared_ptr<D4> dependency4_;  //!< Handle to the fourth lifetime dependency.
   shared_ptr<D5> dependency5_;  //!< Handle to the fifth lifetime dependency.
   shared_ptr<D6> dependency6_;  //!< Handle to the sixth lifetime dependency.

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (5 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1    // Type of the first lifetime dependency
        , typename D2    // Type of the second lifetime dependency
        , typename D3    // Type of the third lifetime dependency
        , typename D4    // Type of the fourth lifetime dependency
        , typename D5 >  // Type of the fifth lifetime dependency
class Singleton<T,D1,D2,D3,D4,D5,NullType,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_5( D1, D2, D3, D4, D5 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
        , dependency3_( D3::instance() )  // Handle to the third lifetime dependency
        , dependency4_( D4::instance() )  // Handle to the fourth lifetime dependency
        , dependency5_( D5::instance() )  // Handle to the fifth lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D3 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D4 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D5 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.
   shared_ptr<D3> dependency3_;  //!< Handle to the third lifetime dependency.
   shared_ptr<D4> dependency4_;  //!< Handle to the fourth lifetime dependency.
   shared_ptr<D5> dependency5_;  //!< Handle to the fifth lifetime dependency.

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (4 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1    // Type of the first lifetime dependency
        , typename D2    // Type of the second lifetime dependency
        , typename D3    // Type of the third lifetime dependency
        , typename D4 >  // Type of the fourth lifetime dependency
class Singleton<T,D1,D2,D3,D4,NullType,NullType,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_4( D1, D2, D3, D4 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
        , dependency3_( D3::instance() )  // Handle to the third lifetime dependency
        , dependency4_( D4::instance() )  // Handle to the fourth lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D3 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D4 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.
   shared_ptr<D3> dependency3_;  //!< Handle to the third lifetime dependency.
   shared_ptr<D4> dependency4_;  //!< Handle to the fourth lifetime dependency.

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (3 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1    // Type of the first lifetime dependency
        , typename D2    // Type of the second lifetime dependency
        , typename D3 >  // Type of the third lifetime dependency
class Singleton<T,D1,D2,D3,NullType,NullType,NullType,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_3( D1, D2, D3 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
        , dependency3_( D3::instance() )  // Handle to the third lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D3 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.
   shared_ptr<D3> dependency3_;  //!< Handle to the third lifetime dependency.

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (2 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1    // Type of the first lifetime dependency
        , typename D2 >  // Type of the second lifetime dependency
class Singleton<T,D1,D2,NullType,NullType,NullType,NullType,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_2( D1, D2 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
        : dependency1_( D1::instance() )  // Handle to the first lifetime dependency
        , dependency2_( D2::instance() )  // Handle to the second lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D2 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the first lifetime dependency.
   shared_ptr<D2> dependency2_;  //!< Handle to the second lifetime dependency.

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (1 LIFETIME DEPENDENCY)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T     // Type of the singleton (CRTP pattern)
        , typename D1 >  // Type of the lifetime dependency
class Singleton<T,D1,NullType,NullType,NullType,NullType,NullType,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef WALBERLA_TYPELIST_1( D1 )  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
      : dependency1_( D1::instance() )  // Handle to the lifetime dependency
   {
      WALBERLA_DETECT_CYCLIC_LIFETIME_DEPENDENCY( D1 );
   }
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
   shared_ptr<D1> dependency1_;  //!< Handle to the lifetime dependency.

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

   static bool isInstantiated_;
   //@}
   //*******************************************************************************************************************
};
//**********************************************************************************************************************



//======================================================================================================================
//
//  SINGLETON SPECIALIZATION (0 LIFETIME DEPENDENCIES)
//
//======================================================================================================================

//**********************************************************************************************************************
template< typename T >  // Type of the singleton (CRTP pattern)
class Singleton<T,NullType,NullType,NullType,NullType,NullType,NullType,NullType,NullType> : private NonCopyable
{
public:
   //**Type definitions****************************************************************************
   //! Type list of all lifetime dependencies.
   typedef NullType  Dependencies;
   //*******************************************************************************************************************

protected:
   //**Constructor******************************************************************************************************
   /*!\brief Constructor for the Singleton class.
   //
   // In case a cyclic lifetime dependency is detected, a compilation error is created.
    */
   explicit Singleton()
   {}
   //*******************************************************************************************************************

   //**Destructor*******************************************************************************************************
   /*!\brief Destructor for the Singleton class.
    */
   ~Singleton()
   {}
   //*******************************************************************************************************************

public:
   //**Instance function***************************************************************************
   /*!\name Instance function */
   //@{
   static const shared_ptr<T>& instance()
   {
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
      boost::mutex::scoped_lock lock( instanceMutex_ );
#endif
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
#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
   static boost::mutex instanceMutex_;  //!< Synchronization mutex for access to the singleton.
#endif

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

#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H >
boost::mutex Singleton<T,A,B,C,D,E,F,G,H>::instanceMutex_;
#endif

template< typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H >
bool Singleton<T,A,B,C,D,E,F,G,H>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G >
boost::mutex Singleton<T,A,B,C,D,E,F,G,NullType>::instanceMutex_;
#endif

template< typename T, typename A, typename B, typename C, typename D, typename E, typename F, typename G >
bool Singleton<T,A,B,C,D,E,F,G,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B, typename C, typename D, typename E, typename F >
boost::mutex Singleton<T,A,B,C,D,E,F,NullType,NullType>::instanceMutex_;
#endif

template< typename T, typename A, typename B, typename C, typename D, typename E, typename F >
bool Singleton<T,A,B,C,D,E,F,NullType,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B, typename C, typename D, typename E >
boost::mutex Singleton<T,A,B,C,D,E,NullType,NullType,NullType>::instanceMutex_;
#endif

template< typename T, typename A, typename B, typename C, typename D, typename E >
bool Singleton<T,A,B,C,D,E,NullType,NullType,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B, typename C, typename D >
boost::mutex Singleton<T,A,B,C,D,NullType,NullType,NullType,NullType>::instanceMutex_;
#endif

template< typename T, typename A, typename B, typename C, typename D >
bool Singleton<T,A,B,C,D,NullType,NullType,NullType,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B, typename C >
boost::mutex Singleton<T,A,B,C,NullType,NullType,NullType,NullType,NullType>::instanceMutex_;
#endif

template< typename T, typename A, typename B, typename C >
bool Singleton<T,A,B,C,NullType,NullType,NullType,NullType,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A, typename B >
boost::mutex Singleton<T,A,B,NullType,NullType,NullType,NullType,NullType,NullType>::instanceMutex_;
#endif

template< typename T, typename A, typename B >
bool Singleton<T,A,B,NullType,NullType,NullType,NullType,NullType,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T, typename A >
boost::mutex Singleton<T,A,NullType,NullType,NullType,NullType,NullType,NullType,NullType>::instanceMutex_;
#endif

template< typename T, typename A >
bool Singleton<T,A,NullType,NullType,NullType,NullType,NullType,NullType,NullType>::isInstantiated_ = false;


#ifdef WALBERLA_BUILD_WITH_BOOST_THREAD
template< typename T >
boost::mutex Singleton<T,NullType,NullType,NullType,NullType,NullType,NullType,NullType,NullType>::instanceMutex_;
#endif

template< typename T >
bool Singleton<T,NullType,NullType,NullType,NullType,NullType,NullType,NullType,NullType>::isInstantiated_ = false;



#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

} // namespace singleton
} // namespace walberla

/// \endcond
