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
//! \file DisplayAdaptor.h
//! \ingroup gui
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#ifndef Q_MOC_RUN
#include "core/math/Vector3.h"
#include "domain_decomposition/IBlock.h"
#include "field/Field.h"
#include "field/GhostLayerField.h"
#endif

#include <QPoint>
#include <QSize>
#include <QTreeWidgetItem>
#include <QVector>


namespace walberla {
namespace gui {

   class CellView;

   //*******************************************************************************************************************
   /*! Abstract base class for display adaptors.
   *
   * A display adaptor describes how to display some structured grid like data structure in a BlockSliceView.
   *
   * The DisplayAdaptor can use a QTreeWidget branch for configuration of display properties.
   * Here the convenience class DisplayPropertiesItem can be used.
   *
   * Using this options the DisplayAdaptor draws the data structure, in a BlockSliceView.
   *
   * \ingroup gui
   *
   */
   //*******************************************************************************************************************
   class DisplayAdaptor : public QObject
   {
   Q_OBJECT

   public:

      /**
       * Configures the DisplayAdaptor by passing in a block and slice information
       *
       * @param [in]  block         the block that should be displayed
       * @param [in]  sliceDim      The coordinate value of this dimension is kept constant. Defines the slice direction
       *                            0 means x, 1 means y, 2 means z
       * @param [in]  sliceValue    The fixed value of the sliceDim coordinate. If the value is not in the allowed range
       *                            nothing is displayed
       * @param [out] innerSize     size of the field, without ghost layer.
       *                            innerSize[0] is the horizontal extent and will be the first coordinate in the grid array of draw()
       *                            innerSize[1] is the vertical extend and will be the second coordinate in the grid array
       *                            innerSize[2] is used for limiting the allowed sliceValues
       * @param [out] ghostLayers   number of ghost layers for all three dimensions of innerSize
       */
      virtual void configure( const IBlock * block, int sliceDim, cell_idx_t sliceValue,
                              Vector3<uint_t> & innerSize, cell_idx_t & ghostLayers ) = 0;

      /**
       * Draws data structure represented by this Adaptor to a grid of cells
       * Block has to be set first!
       *
       * @param [in,out] grid          Grid of cells. Use functions of CellView to draw
       *                               The Grid has dimension of innerSize + ghostLayers
       * @param [in] nrOfGhostLayers   the number of ghost layers in each direction. This value may be bigger
       *                               than the returned value of configure() since other fields might have more ghost layers.
       *                               The cell at grid[0][0] lies at the outermost ghost layer
       */
      virtual void draw( QVector<QVector<CellView*> > & grid, int nrOfGhostLayers ) = 0;



      /**
       * Optionally adds a QTreeWidgetItem to a tree, containing configuration options
       */
      virtual void addConfigurationItem( QTreeWidgetItem *   ) {}


   signals:

      /// This signal requests a redraw, for example after a configuration change
      void requestRedraw();

   };







} // namespace gui
} // namespace walberla


