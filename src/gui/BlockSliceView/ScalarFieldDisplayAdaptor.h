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
//! \file ScalarFieldDisplayAdaptor.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "DisplayPropertiesItem.h"
#include "FieldDisplayAdaptor.h"


namespace walberla {
namespace gui {


   //*******************************************************************************************************************
   /*! Class for drawing a slice of a scalar field
   *
   * Scalar fields can be displayed either by showing the numeric values as text, or by using a colormap.
   *
   * \ingroup gui
   *
   */
   //*******************************************************************************************************************
   template< typename field_t>
   class ScalarFieldDisplayAdaptor : public FieldDisplayAdaptor<field_t>
   {
   public:
      ScalarFieldDisplayAdaptor( ConstBlockDataID scalarFieldID );
      ~ScalarFieldDisplayAdaptor();

      virtual void addConfigurationItem( QTreeWidgetItem * parentItem );

      virtual void draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers );

      void configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                      Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers );

   private:
      typedef typename field_t::value_type T;

      static const uint_t F_SIZE = field_t::F_SIZE;

      void drawScalarFieldNumeric( CellView * cell, const typename field_t::const_iterator & it );
      void drawScalarFieldColormap( CellView * cell, const typename field_t::const_iterator & it, T min, T max );

      using FieldDisplayAdaptor<field_t>::sliceDim_;
      using FieldDisplayAdaptor<field_t>::field_;
      using FieldDisplayAdaptor<field_t>::sliceInterval_;
      using FieldDisplayAdaptor<field_t>::blockDataId_;
      cell_idx_t f;
      DisplayPropertiesItem * displayProperties_;
   };




} // namespace gui
} // namespace walberla


#include "ScalarFieldDisplayAdaptor.impl.h"
