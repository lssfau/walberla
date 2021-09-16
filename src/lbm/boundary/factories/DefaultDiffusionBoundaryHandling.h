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
//! \file DefaultDiffusionBoundaryHandling.h
//! \ingroup lbm
//! \author Matthias Markl <matthias.markl@fau.de>
//
//======================================================================================================================

#pragma once

#include "blockforest/StructuredBlockForest.h"

#include "boundary/BoundaryHandling.h"

#include "core/debug/Debug.h"

#include "lbm/boundary/DiffusionDirichlet.h"
#include "lbm/boundary/SimpleDiffusionDirichlet.h"
#include "lbm/boundary/NoDiffusion.h"
#include "lbm/boundary/FreeDiffusion.h"

#include "field/FlagFunctions.h"

#include <functional>


namespace walberla {
namespace lbm{

template< typename LatticeModel_T, typename FlagField_T >
class DefaultDiffusionBoundaryHandlingFactory
{
private:
   using flag_t = typename FlagField_T::flag_t;

   using Stencil = typename LatticeModel_T::Stencil;

public:
   using DiffusionDirichlet_T = lbm::DiffusionDirichlet<LatticeModel_T, flag_t>;
   using SimpleDiffusionDirichlet_T = lbm::SimpleDiffusionDirichlet<LatticeModel_T, flag_t>;
   using NoDiffusion_T = lbm::NoDiffusion<LatticeModel_T, flag_t>;
   using FreeDiffusion_T = lbm::FreeDiffusion<LatticeModel_T, FlagField_T>;

public:
   using BoundaryHandling_T = walberla::boundary::BoundaryHandling<FlagField_T, Stencil, DiffusionDirichlet_T, SimpleDiffusionDirichlet_T, NoDiffusion_T, FreeDiffusion_T, SimpleDiffusionDirichlet_T, SimpleDiffusionDirichlet_T>;
   using BoundaryHandling = BoundaryHandling_T;

   const static FlagUID& getDiffusionDirichletFlagUID()        { static FlagUID uid( "DiffusionDirichlet"        ); return uid; }
   const static FlagUID& getSimpleDiffusionDirichletFlagUID()  { static FlagUID uid( "SimpleDiffusionDirichlet"  ); return uid; }
   const static FlagUID& getSimpleDiffusionDirichletFlagUID1() { static FlagUID uid( "SimpleDiffusionDirichlet1" ); return uid; }
   const static FlagUID& getSimpleDiffusionDirichletFlagUID2() { static FlagUID uid( "SimpleDiffusionDirichlet2" ); return uid; }
   const static FlagUID& getNoDiffusionFlagUID()               { static FlagUID uid( "NoDiffusion"               ); return uid; }
   const static FlagUID& getFreeDiffusionFlagUID()             { static FlagUID uid( "FreeDiffusion"             ); return uid; }

   const static BoundaryUID& getDiffusionDirichletBoundaryUID()        { static BoundaryUID uid( "DiffusionDirichlet"        ); return uid; }
   const static BoundaryUID& getSimpleDiffusionDirichletBoundaryUID()  { static BoundaryUID uid( "SimpleDiffusionDirichlet"  ); return uid; }
   const static BoundaryUID& getSimpleDiffusionDirichletBoundaryUID1() { static BoundaryUID uid( "SimpleDiffusionDirichlet1" ); return uid; }
   const static BoundaryUID& getSimpleDiffusionDirichletBoundaryUID2() { static BoundaryUID uid( "SimpleDiffusionDirichlet2" ); return uid; }
   const static BoundaryUID& getNoDiffusionBoundaryUID()               { static BoundaryUID uid( "NoDiffusion"               ); return uid; }
   const static BoundaryUID& getFreeDiffusionBoundaryUID()             { static BoundaryUID uid( "FreeDiffusion"             ); return uid; }

   static Set<FlagUID> getFlagUIDs(){
      Set<FlagUID> flagUIDs;
      flagUIDs += getDiffusionDirichletFlagUID();
      flagUIDs += getSimpleDiffusionDirichletFlagUID();
      flagUIDs += getSimpleDiffusionDirichletFlagUID1();
      flagUIDs += getSimpleDiffusionDirichletFlagUID2();
      flagUIDs += getNoDiffusionFlagUID();
      flagUIDs += getFreeDiffusionFlagUID();
      return flagUIDs;
   }
private:

   static BoundaryHandling_T* createDefaultDiffusionBoundaryHandlingFactory( 
      IBlock* const block, const StructuredBlockStorage* const /*bs*/, const BlockDataID& flagFieldID, const Set< FlagUID > & domainFlagUIDs, const BlockDataID& pdfFieldID, const Set< FlagUID > & initFlagUIDs )
   {
      using PDFField_T = lbm::PdfField<LatticeModel_T>;

      WALBERLA_ASSERT_NOT_NULLPTR( block );

      FlagField_T* flagField = block->getData< FlagField_T >( flagFieldID );
      PDFField_T*  pdfField  = block->getData< PDFField_T  >( pdfFieldID );

      WALBERLA_ASSERT_NOT_NULLPTR( flagField );
      WALBERLA_ASSERT_NOT_NULLPTR( pdfField );

      flag_t domainMask(0u);
      for( auto domainFlagUID = domainFlagUIDs.begin(); domainFlagUID != domainFlagUIDs.end(); ++domainFlagUID )
         field::addMask( domainMask, flagField->getOrRegisterFlag( *domainFlagUID ) );

      BoundaryHandling_T * handling = new BoundaryHandling_T( "Diffusion Boundary Handling", flagField, domainMask,
         DiffusionDirichlet_T      ( getDiffusionDirichletBoundaryUID(),        getDiffusionDirichletFlagUID(),        pdfField, flagField ),
         SimpleDiffusionDirichlet_T( getSimpleDiffusionDirichletBoundaryUID(),  getSimpleDiffusionDirichletFlagUID(),  pdfField ),
         NoDiffusion_T             ( getNoDiffusionBoundaryUID(),               getNoDiffusionFlagUID(),               pdfField ),
         FreeDiffusion_T           ( getFreeDiffusionBoundaryUID(),             getFreeDiffusionFlagUID(),             pdfField, flagField, domainMask ),
         SimpleDiffusionDirichlet_T( getSimpleDiffusionDirichletBoundaryUID1(), getSimpleDiffusionDirichletFlagUID1(), pdfField ),
         SimpleDiffusionDirichlet_T( getSimpleDiffusionDirichletBoundaryUID2(), getSimpleDiffusionDirichletFlagUID2(), pdfField ) );

      if( initFlagUIDs.size() > size_t(0u) )
      {
         flag_t initMask(0u);
         for( auto initFlagUID = initFlagUIDs.begin(); initFlagUID != initFlagUIDs.end(); ++initFlagUID )
            if( flagField->flagExists( *initFlagUID ) )
               field::addMask( initMask, flagField->getFlag( *initFlagUID ) );
            else
               WALBERLA_ABORT( "Try to init flag field with a non registered flag: " << *initFlagUID );
         handling->fillWithDomain( initMask, uint_t(0u) );
      }

      return handling;
   }

public:
   static BlockDataID addDefaultDiffusionBoundaryHandlingToStorage(
      const shared_ptr< StructuredBlockStorage >& bs, const std::string & identifier, const BlockDataID& flagFieldID, const Set<FlagUID>& domainFlagUIDs, const BlockDataID& pdfFieldID, const bool fillWithDomain = true )
   {
      return addDefaultDiffusionBoundaryHandlingToStorage( bs, identifier, flagFieldID, domainFlagUIDs, pdfFieldID, fillWithDomain ? domainFlagUIDs : Set<FlagUID>() );
   }

   static BlockDataID addDefaultDiffusionBoundaryHandlingToStorage(
      const shared_ptr< StructuredBlockStorage >& bs, const std::string & identifier, const BlockDataID& flagFieldID, const Set<FlagUID>& domainFlagUIDs, const BlockDataID& pdfFieldID, const Set<FlagUID>& initFlagUIDs )
   {
      auto func = std::bind( createDefaultDiffusionBoundaryHandlingFactory, std::placeholders::_1, std::placeholders::_2, flagFieldID, domainFlagUIDs, pdfFieldID, initFlagUIDs );
      return bs->addStructuredBlockData< BoundaryHandling_T >( func, identifier );
   }

}; // class DefaultDiffusionBoundaryHandlingFactory


} // namespace lbm
} // namespace walberla
