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
//! \file PdfFieldDisplayAdaptor.h
//! \ingroup lbm
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "gui/BlockSliceView/DisplayPropertiesItem.h"
#include "gui/BlockSliceView/FieldDisplayAdaptor.h"

#include <type_traits>


namespace walberla {
namespace lbm {


   //*******************************************************************************************************************
   /*!
   * Class for drawing a slice of a D3Q19 PDF field.
   *
   *
   * \ingroup gui
   *
   */
   //*******************************************************************************************************************
   template<typename field_t, typename stencil_t >
   class PdfFieldDisplayAdaptor : public gui::FieldDisplayAdaptor< field_t >
   {
   public:

      PdfFieldDisplayAdaptor( ConstBlockDataID pdfFieldID );
      virtual ~PdfFieldDisplayAdaptor();


      virtual void addConfigurationItem( QTreeWidgetItem * parentItem );

      virtual void draw( QVector<QVector< gui::CellView*> > & grid, int nrOfGhostLayers );

      void configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                      Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers );

   private:
      typedef typename field_t::value_type T;
      static_assert( (std::is_same<T,double>::value || std::is_same<T,float>::value),
                     "Only floating point fields are supported" );

      using gui::FieldDisplayAdaptor<field_t>::sliceDim_;
      using gui::FieldDisplayAdaptor<field_t>::field_;
      using gui::FieldDisplayAdaptor<field_t>::sliceInterval_;
      using gui::FieldDisplayAdaptor<field_t>::blockDataId_;

      gui::DisplayPropertiesItem * displayProperties_;
   };




} // namespace lbm
} // namespace walberla


#include "PdfFieldDisplayAdaptor.impl.h"
