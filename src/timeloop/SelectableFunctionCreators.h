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
//! \file SelectableFunctionCreators.h
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/selectable/SetSelectableObject.h"
#include "core/uid/SUID.h"
#include "domain_decomposition/BlockStorage.h"

#include <functional>
#include <string>


namespace walberla {
namespace timeloop {



   template < typename FuncType >
   struct FuncCreator
   {
      FuncCreator( std::function< FuncType > fct,
                  const std::string& identifier            = std::string(),
                  const Set<SUID>&   requiredSelectors     = Set<SUID>::emptySet(),
                  const Set<SUID>&   incompatibleSelectors = Set<SUID>::emptySet() )
      :  function_( fct ), identifier_( identifier ),
         requiredSelectors_( requiredSelectors ), incompatibleSelectors_( incompatibleSelectors )
      {}

      FuncCreator() = default;

      std::function< FuncType > function_;
      std::string                 identifier_;
      Set<SUID>                   requiredSelectors_;
      Set<SUID>                   incompatibleSelectors_;
   };


   template <typename FuncType>
   struct SelectableFunction
   {
      SelectableFunction() = default;

      SelectableFunction ( std::function< FuncType > fct,
                           const std::string& identifier            = std::string(),
                           const Set<SUID>&   requiredSelectors     = Set<SUID>::emptySet(),
                           const Set<SUID>&   incompatibleSelectors = Set<SUID>::emptySet() )
      {
         selectableFunc_.add( fct, requiredSelectors, incompatibleSelectors, identifier );
      }

      SelectableFunction& operator<<( const FuncCreator<FuncType> & fct )
      {
         selectableFunc_.add( fct.function_, fct.requiredSelectors_, fct.incompatibleSelectors_, fct.identifier_ );
         return *this;
      }

      selectable::SetSelectableObject< std::function<FuncType>, SUID > selectableFunc_;
   };


   struct BeforeFunction : public SelectableFunction< void () >
   {
      BeforeFunction() = default;
      BeforeFunction(std::function< void () > fct,
                     const std::string& id      = std::string(),
                     const Set<SUID>&   req    = Set<SUID>::emptySet(),
                     const Set<SUID>&   incomp = Set<SUID>::emptySet())
         : SelectableFunction<void ()> (fct, id, req, incomp )
      {}

      BeforeFunction& operator<<( const FuncCreator< void () > & fct )
      {
         SelectableFunction<void()> ::operator<<( fct );
         return *this;
      }

   };

   struct AfterFunction : public SelectableFunction< void () >
   {
       AfterFunction() = default;
       AfterFunction(std::function< void () > fct,
                     const std::string& id      = std::string(),
                     const Set<SUID>&   req    = Set<SUID>::emptySet(),
                     const Set<SUID>&   incomp = Set<SUID>::emptySet())
         : SelectableFunction<void ()> (fct, id, req, incomp )
      {}

      AfterFunction& operator<<( const FuncCreator< void () > & fct )
      {
         SelectableFunction<void()> ::operator<<( fct );
         return *this;
      }
   };


   template< typename T >
   struct SweepOnBlock
   {
      SweepOnBlock( std::function< T* ( IBlock* const block ) > function,
                     const std::string& identifier            = std::string(),
                     const Set<SUID>&   required     = Set<SUID>::emptySet(),
                     const Set<SUID>&   incompatible = Set<SUID>::emptySet() ) :
      function_( function ), identifier_( identifier ),
      requiredSelectors_( required ), incompatibleSelectors_( incompatible ) {}

      std::function< T* ( IBlock* const block ) > function_;

      std::string                                   identifier_;
      Set<SUID>                                     requiredSelectors_;
      Set<SUID>                                     incompatibleSelectors_;
   };

   using Sweep = FuncCreator<void (IBlock *)>;


   template<typename SweepClass>
   void executeSweepOnBlock ( IBlock * block, BlockDataID bdID )
   {
      SweepClass * sweep = block->getData<SweepClass>( bdID );
      if( sweep ) // not allocated on block
         (*sweep)();
   }


   /*-----------------------------------------------------------------------------------------------------------------*/


   struct SweepAdder
   {
      SweepAdder( BlockStorage & storage, uint_t id )
         : bs_ (storage ), id_(id)
      {}


      SweepAdder& operator<<( const BeforeFunction & bf )
      {
         beforeFuncs.push_back( bf );
         return *this;
      }

      SweepAdder& operator<<( const AfterFunction & af )
      {
         afterFuncs.push_back( af );
         return *this;
      }

      SweepAdder& operator<<( const Sweep & sw )
      {
         sweep.add ( sw.function_, sw.requiredSelectors_, sw.incompatibleSelectors_, sw.identifier_ );
         return *this;
      }


      template < typename SweepClass >
      SweepAdder & operator<< ( const SweepOnBlock<SweepClass> & sw )
      {
         BlockDataCreator<SweepClass> bdCreator (sw.function_, sw.identifier_, sw.requiredSelectors_, sw.incompatibleSelectors_);
         BlockDataID bdId = bs_.addBlockData() << bdCreator;

         // add a sweep function that fetches the block data sweep and executes it
         auto sweepFunc = std::bind ( executeSweepOnBlock<SweepClass>,  std::placeholders::_1, bdId );
         ( *this ) << Sweep ( sweepFunc, sw.identifier_, sw.requiredSelectors_, sw.incompatibleSelectors_  );

         return *this;
      }


   private:
      template < typename TimingPolicy > friend class SweepTimeloop;

      BlockStorage & bs_;

      std::vector<BeforeFunction> beforeFuncs;
      std::vector<AfterFunction>  afterFuncs;

      selectable::SetSelectableObject<Sweep, SUID > sweep;

      uint_t id_;
   };



} // namespace timeloop
} // namespace walberla


