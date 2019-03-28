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
//! \file FlagFieldDisplayAdaptor.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "DisplayPropertiesItem.h"
#include "FieldDisplayAdaptor.h"
#include "field/FlagField.h"

#include <QVector>
#include <QMap>
#include <type_traits>

namespace walberla {
namespace gui {


   //*******************************************************************************************************************
   /*!
   * Displays a slice of a flag field
   *
   * \ingroup gui
   */
   //*******************************************************************************************************************
   template< typename field_t >
   class FlagFieldDisplayAdaptor : public FieldDisplayAdaptor<field_t>
   {
   public:

      FlagFieldDisplayAdaptor( ConstBlockDataID flagFieldID );
      ~FlagFieldDisplayAdaptor();

      virtual void addConfigurationItem( QTreeWidgetItem * parentItem );

      virtual void draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers );

      void configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                      Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers );


   protected:
      typedef typename field_t::value_type T;

      static_assert( field_t::F_SIZE == 1, "Flag fields have no F_SIZE > 1 " );
      static_assert( std::is_same<T,uint8_t >::value ||
                     std::is_same<T,uint16_t>::value ||
                     std::is_same<T,uint32_t>::value ||
                     std::is_same<T,uint64_t>::value,
                     "Flag fields have unsigned integers as value");

      using FieldDisplayAdaptor<field_t>::sliceDim_;
      using FieldDisplayAdaptor<field_t>::field_;
      using FieldDisplayAdaptor<field_t>::sliceInterval_;
      using FieldDisplayAdaptor<field_t>::blockDataId_;

      DisplayPropertiesItem * displayProperties_;
      QVector<DisplayPropertiesItem*> flagProperties_;

      const FlagField<T> * flagField_;

      const IBlock * lastBlock_;

      static QMap< QString, QColor > flagNameToColor;
      static QMap< QString, QColor > createFlagNameToColorMap();

      static QVector< QColor > defaultFlagColors;
      static QVector< QColor > createDefaultFlagColors();
      static int defaultFlagColorCounter;
   };




} // namespace gui
} // namespace walberla


#include "FlagFieldDisplayAdaptor.impl.h"
