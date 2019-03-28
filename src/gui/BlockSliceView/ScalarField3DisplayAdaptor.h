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
//! \file ScalarField3DisplayAdaptor.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "DisplayPropertiesItem.h"
#include "FieldDisplayAdaptor.h"
#include "core/math/Vector3.h"

#include <type_traits>


namespace walberla {
namespace gui {



   //*******************************************************************************************************************
   /*!
   * Class for drawing slices of vector fields
   *
   * \ingroup gui
   *
   */
   //*******************************************************************************************************************
   template < typename field_t >
   class ScalarField3DisplayAdaptor : public FieldDisplayAdaptor< field_t >
   {
   public:
      typedef typename field_t::value_type T;
      typedef FieldDisplayAdaptor< field_t > base_t;
      typedef real_t                length_t;

       ScalarField3DisplayAdaptor( ConstBlockDataID vectorFieldID );
      ~ScalarField3DisplayAdaptor();


      virtual void addConfigurationItem( QTreeWidgetItem * parentItem );

      virtual void draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers );

      void configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                      Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers );


   protected:
      // first value_type for field, second for Vector3
      static_assert( field_t::F_SIZE == 3, "Can only display fields with 3 components ( vector )" );

      void drawVectorFieldNumeric ( CellView * cell, const typename field_t::const_iterator & it );
      void drawVectorFieldColormap( CellView * cell, const typename field_t::const_iterator & it, length_t min, length_t max );
      void drawVectorFieldArrow   ( CellView * cell, const typename field_t::const_iterator & it, length_t max );

      using base_t::blockDataId_;
      using base_t::sliceDim_;
      using base_t::sliceInterval_;
      using base_t::field_;


      DisplayPropertiesItem * displayProperties_;

   };




} // namespace gui
} // namespace walberla


#include "ScalarField3DisplayAdaptor.impl.h"


