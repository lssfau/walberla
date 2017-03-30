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
//! \file MarchingCubes.h
//! \ingroup postprocessing
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/Vector3.h"
#include "field/GhostLayerField.h"


namespace walberla {

   namespace geometry {
      class TriangleMesh;
   }

   namespace postprocessing
   {

      /**
       * \brief Creates an iso-surface mesh from a scalar field using the Marching Cubes algorithm
       *
       * The mesh is not cleared, all vertices and faces are added.
       *
       * \param field         Scalar field, that should be triangulated, has to be a GhostLayerField or
       *                      a class with similar interface (adaptor)
       * \param threshold     all areas with values greater than threshold, are inside the generated surface
       * \param mesh          the generated vertices and triangles are added to this mesh
       *                      duplicate vertices are not recognized, but added again
       * \param dx            length of a cell i.e. grid spacing in dx/dy/dz
       * \param fCoord        fixed f coordinate, only scalar fields can be processed, so here a fixed f
       *                      has to be chosen
       * \param offset        this vector is added to all vertices that are inserted into the mesh
       * \param cellInterval  use only the part specified by this CellInterval to create the mesh
       *                      if the cellInterval is empty the complete field is used
       * \tparam Field_T      Has to be a GhostLayerField<real_t,1> or a class with similar interface (adaptor)
       */
      template<typename Field_T>
      void generateIsoSurface( const Field_T & field,
                               real_t threshold,
                               geometry::TriangleMesh & mesh,
                               const Vector3<real_t> & dx = Vector3<real_t>( real_t(1) ),
                               uint_t fCoord = uint_t(0),
                               const Vector3<real_t> & offset = Vector3<real_t>( real_t(0) ),
                               const CellInterval & cellInterval = CellInterval(),
                               bool calcNormals = true );



   } // namespace postprocessing
} // namespace walberla


#include "MarchingCubes.impl.h"
