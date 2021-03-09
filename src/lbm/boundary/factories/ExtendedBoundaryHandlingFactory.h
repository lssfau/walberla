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
//! \file ExtendedBoundaryHandlingFactory.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/boundary/FreeSlip.h"
#include "lbm/boundary/NoSlip.h"
#include "lbm/boundary/Pressure.h"
#include "lbm/boundary/ParserUBB.h"
#include "lbm/boundary/Outlet.h"
#include "lbm/boundary/Curved.h"

#include "core/math/Vector3.h"
#include "domain_decomposition/BlockDataID.h"
#include "boundary/BoundaryHandling.h"
#include "lbm/field/PdfField.h"

namespace walberla {
namespace lbm {


//**********************************************************************************************************************
/*!
* \brief Creates a default boundary handling for LBM simulations
*
* \ingroup lbm
*
* This functor is usually used like
\code{.cpp}
typedef lbm::ExtendedBoundaryHandlingFactory< LatticeModel_T, FlagField_T > Factory;
BlockDataID bid = Factory::addBoundaryHandlingToStorage( blocks, "boundary handling",
                                              flagFieldId, pdfFieldId, setOfDomainFlags );
\endcode
*
* Use this boundary handling if the DefaultBoundaryHandling is not flexible enough for you.
* This boundary handling stores boundary velocities/pressure for each cell separately resulting in a memory
* consumption of 4 real values per cell. The Handling itself is also more costly since these values have to be
* loaded from memory.
*
*
* The following boundary conditions are available:
*   - NoSlip
*   - FreeSlip
*   - Pressure
*   - UBB
*   - Outlet
*   - Curved
*
* \tparam LatticeModel  The lattice model used for the simulation
* \tparam FlagFieldT    Type of the used flag field
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT >
class ExtendedBoundaryHandlingFactory
{
public:
   using flag_t = typename FlagFieldT::flag_t;
   using Stencil = typename LatticeModel::Stencil;
   using Velocity = Vector3<real_t>;
   using PdfFieldLM = PdfField<LatticeModel>;

   using BcNoSlip = NoSlip<LatticeModel, flag_t>;
   using BcFreeSlip = FreeSlip<LatticeModel, FlagFieldT>;
   using BcPressure = Pressure<LatticeModel, flag_t>;
   using BcUBB = ParserUBB<LatticeModel, flag_t>;
   using BcOutlet = Outlet<LatticeModel, FlagFieldT>;
   using BcCurved = Curved<LatticeModel, FlagFieldT>;

   using BoundaryHandling = walberla::boundary::BoundaryHandling<FlagFieldT, Stencil, BcNoSlip, BcFreeSlip, BcPressure, BcUBB, BcOutlet, BcCurved>;

   static BlockDataID addBoundaryHandlingToStorage( const shared_ptr< StructuredBlockStorage > & bs, const std::string & identifier,
                                                    BlockDataID flagFieldID, BlockDataID pdfFieldID, const Set< FlagUID > & flagUIDSet)
   {
      return addBoundaryHandlingToStorage(bs, identifier, pdfFieldID, pdfFieldID, flagUIDSet,BoundaryHandling::Mode::OPTIMIZED_SPARSE_TRAVERSAL );
   }

   static BlockDataID addBoundaryHandlingToStorage( const shared_ptr< StructuredBlockStorage > & bs, const std::string & identifier,
                                                    BlockDataID flagFieldID, BlockDataID pdfFieldID, const Set< FlagUID > & flagUIDSet,
                                                    const typename BoundaryHandling::Mode boundaryHandlingTraversalMode )
   {
      ExtendedBoundaryHandlingFactory factory ( flagFieldID, pdfFieldID, flagUIDSet, boundaryHandlingTraversalMode );

      return bs->addStructuredBlockData< BoundaryHandling >( factory, identifier );
   }

   static const walberla::FlagUID & getNoSlip()    { static FlagUID uid( "NoSlip" );    return uid; }
   static const walberla::FlagUID & getFreeSlip()  { static FlagUID uid( "FreeSlip" );  return uid; }
   static const walberla::FlagUID & getPressure()  { static FlagUID uid( "Pressure" );  return uid; }
   static const walberla::FlagUID & getUBB()       { static FlagUID uid( "UBB" );       return uid; }
   static const walberla::FlagUID & getOutlet()    { static FlagUID uid( "Outlet" );    return uid; }
   static const walberla::FlagUID & getCurved()    { static FlagUID uid( "Curved" );    return uid; }

   static const walberla::BoundaryUID & getNoSlipBoundaryUID()    { static BoundaryUID uid( "NoSlip" );    return uid; }
   static const walberla::BoundaryUID & getFreeSlipBoundaryUID()  { static BoundaryUID uid( "FreeSlip" );  return uid; }
   static const walberla::BoundaryUID & getPressureBoundaryUID()  { static BoundaryUID uid( "Pressure" );  return uid; }
   static const walberla::BoundaryUID & getUBBBoundaryUID()       { static BoundaryUID uid( "UBB" );       return uid; }
   static const walberla::BoundaryUID & getOutletBoundaryUID()    { static BoundaryUID uid( "Outlet" );    return uid; }
   static const walberla::BoundaryUID & getCurvedBoundaryUID()    { static BoundaryUID uid( "Curved" );    return uid; }


   ExtendedBoundaryHandlingFactory( const BlockDataID & flagField, const BlockDataID & pdfField, const Set< FlagUID > & flagUIDSet,
                                    const typename BoundaryHandling::Mode boundaryHandlingTraversalMode );

   BoundaryHandling * operator()( walberla::IBlock * const block, const walberla::StructuredBlockStorage * const storage ) const;

private:
   BlockDataID flagField_;
   BlockDataID pdfField_;

   const Set< FlagUID > flagUIDSet_;

   const typename BoundaryHandling::Mode boundaryHandlingTraversalMode_;

}; // class ExtendedBoundaryHandlingFactory


//**********************************************************************************************************************
/*!
* \ingroup lbm
*
* \param flagField  BlockDataID of the flag field used in the simulation
* \param pdfField   BlockDataID of the PDF field used in the simulation
*/
//**********************************************************************************************************************
template <typename LatticeModel, typename FlagFieldT >
ExtendedBoundaryHandlingFactory<LatticeModel, FlagFieldT>::ExtendedBoundaryHandlingFactory(
                                                   const BlockDataID & flagField, const BlockDataID & pdfField, const Set< FlagUID > & flagUIDSet,
                                                   const typename BoundaryHandling::Mode boundaryHandlingTraversalMode ) :
   flagField_( flagField ), pdfField_( pdfField ), flagUIDSet_(flagUIDSet), boundaryHandlingTraversalMode_(boundaryHandlingTraversalMode)
{
}


template <typename LatticeModel, typename FlagFieldT >
typename ExtendedBoundaryHandlingFactory<LatticeModel, FlagFieldT>::BoundaryHandling *
ExtendedBoundaryHandlingFactory<LatticeModel, FlagFieldT>::operator()( IBlock * const block,
                                                                      const walberla::StructuredBlockStorage * const storage ) const
{
   PdfFieldLM * const pdfField  = block->getData< PdfFieldLM >( pdfField_  );
   FlagFieldT * const flagField = block->getData< FlagFieldT >( flagField_ );

   flag_t mask = 0;
   for( auto flag = flagUIDSet_.begin(); flag != flagUIDSet_.end(); ++flag )
      mask = static_cast< flag_t >( mask | flagField->getOrRegisterFlag( *flag ) );



   BoundaryHandling * handling = new BoundaryHandling( "extended lbm boundary handling", flagField, mask,
        BcNoSlip    ( getNoSlipBoundaryUID(),   getNoSlip(),   pdfField ),
        BcFreeSlip  ( getFreeSlipBoundaryUID(), getFreeSlip(), pdfField, flagField, mask ),
        BcPressure  ( getPressureBoundaryUID(), getPressure(), pdfField ),
        BcUBB       ( getUBBBoundaryUID(),      getUBB(),      pdfField, flagField, storage->getLevel(*block), block->getAABB() ),
        BcOutlet    ( getOutletBoundaryUID(),   getOutlet(),   pdfField, flagField, mask ),
        BcCurved    ( getCurvedBoundaryUID(),   getCurved(),   pdfField, flagField, mask ),
        boundaryHandlingTraversalMode_
    );

   return handling;
}


} // namespace lbm
} // namespace walberla
