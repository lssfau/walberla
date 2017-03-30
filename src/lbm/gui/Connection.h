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
//! \file Connection.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "lbm/field/Adaptors.h"
#include "lbm/field/PdfField.h"
#include "gui/Gui.h"


#ifdef WALBERLA_ENABLE_GUI
#include "gui/BlockSliceView/ScalarFieldDisplayAdaptor.h"
#include "gui/BlockSliceView/VectorFieldDisplayAdaptor.h"
#include "PdfFieldDisplayAdaptor.h"
#endif

namespace walberla {
namespace lbm      {


#ifdef WALBERLA_ENABLE_GUI

   namespace internal
   {
      template<typename LatticeModel_T>
      gui::DisplayAdaptor * addAdaptors( const IBlock & block, ConstBlockDataID blockDataID )
      {
         const uint_t stencilSize = LatticeModel_T::Stencil::Size;
         // Pdf fields
         if ( block.isDataOfType< GhostLayerField<real_t, stencilSize > >( blockDataID) )
            return new PdfFieldDisplayAdaptor<GhostLayerField<real_t, stencilSize>, typename LatticeModel_T::Stencil >( blockDataID );

         // Pdf fields
         if ( block.isDataOfType< PdfField<LatticeModel_T> >( blockDataID) )
            return new PdfFieldDisplayAdaptor< PdfField<LatticeModel_T>, typename LatticeModel_T::Stencil >( blockDataID );


         // Field adaptors
         typedef typename lbm::Adaptor<LatticeModel_T>::Density                      DensityAdaptor;
         typedef typename lbm::Adaptor<LatticeModel_T>::VelocityVector               VelocityVectorAdaptor;
         typedef typename lbm::Adaptor<LatticeModel_T>::StreamMomentumDensityVector  MomentumVectorAdaptorStream;
         typedef typename lbm::Adaptor<LatticeModel_T>::CollideMomentumDensityVector MomentumVectorAdaptorCollide;


         if ( block.isDataOfType< DensityAdaptor >( blockDataID) )
            return new gui::ScalarFieldDisplayAdaptor<DensityAdaptor>( blockDataID );

         if ( block.isDataOfType< VelocityVectorAdaptor >( blockDataID) )
            return new gui::VectorFieldDisplayAdaptor<VelocityVectorAdaptor>( blockDataID );

         if ( block.isDataOfType< MomentumVectorAdaptorStream >( blockDataID) )
            return new gui::VectorFieldDisplayAdaptor<MomentumVectorAdaptorStream>( blockDataID );

         if ( block.isDataOfType< MomentumVectorAdaptorCollide >( blockDataID) )
            return new gui::VectorFieldDisplayAdaptor<MomentumVectorAdaptorCollide>( blockDataID );

         return NULL;
      }
   }

   template<typename LatticeModel_T>
   void connectToGui( gui::GUI & guiObject )
   {
      guiObject.registerDisplayAdaptorCreator( internal::addAdaptors<LatticeModel_T> );
   }
#else

   template<typename LatticeModel_T>
   void connectToGui( gui::GUI &  )
   {
   }

#endif



} // namespace lbm
} // namespace walberla






