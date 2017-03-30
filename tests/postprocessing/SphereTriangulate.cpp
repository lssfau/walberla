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
//! \file SphereTriangulate.cpp
//! \ingroup postprocessing
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#include "postprocessing/MarchingCubes.h"

#include "core/DataTypes.h"
#include "core/debug/TestSubsystem.h"

#include "field/GhostLayerField.h"

#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/TriangleMeshIO.h"

#include <fstream>
#include <iostream>


using namespace walberla;
using namespace postprocessing;

using geometry::TriangleMesh;

using std::cout;
using std::endl;
using std::ofstream;

typedef GhostLayerField<real_t,1> ScalarField;


template< typename T, uint_t fs>
shared_ptr<Field<shared_ptr<GhostLayerField<T,fs> >,1 > >
splitGhostLayerField(const GhostLayerField<T,fs> & fieldToSplit,
                     const std::vector<size_t> & cutsX,
                     const std::vector<size_t> & cutsY,
                     const std::vector<size_t> & cutsZ)
{
   typedef Field<shared_ptr<GhostLayerField<T,fs> >,1 > FieldOfFieldPtrs;
   cell_idx_t gl = cell_idx_c ( fieldToSplit.nrOfGhostLayers() );

   cell_idx_t fieldsX = cell_idx_c ( cutsX.size()+1 );
   cell_idx_t fieldsY = cell_idx_c ( cutsY.size()+1 );
   cell_idx_t fieldsZ = cell_idx_c ( cutsZ.size()+1 );
   shared_ptr<FieldOfFieldPtrs> out = make_shared<FieldOfFieldPtrs>(fieldsX,fieldsY,fieldsZ);

   for(cell_idx_t fz = 0; fz < fieldsZ; ++fz)
      for(cell_idx_t fy = 0; fy < fieldsY; ++fy)
         for(cell_idx_t fx = 0; fx < fieldsX; ++fx)
         {
            cell_idx_t glbStartX = (fx == 0) ? 0 : cell_idx_c ( cutsX[uint_c(fx-1)] );
            cell_idx_t glbStartY = (fy == 0) ? 0 : cell_idx_c ( cutsY[uint_c(fy-1)] );
            cell_idx_t glbStartZ = (fz == 0) ? 0 : cell_idx_c ( cutsZ[uint_c(fz-1)] );
            cell_idx_t glbEndX = (fx == fieldsX-1) ? cell_idx_c( fieldToSplit.xSize() ): cell_idx_c (cutsX[uint_c(fx)] );
            cell_idx_t glbEndY = (fy == fieldsY-1) ? cell_idx_c( fieldToSplit.ySize() ): cell_idx_c (cutsY[uint_c(fy)] );
            cell_idx_t glbEndZ = (fz == fieldsZ-1) ? cell_idx_c( fieldToSplit.zSize() ): cell_idx_c (cutsZ[uint_c(fz)] );

            out->get(fx,fy,fz) = make_shared<GhostLayerField<T,fs> >(glbEndX-glbStartX, glbEndY-glbStartY, glbEndZ - glbStartZ, gl );
            GhostLayerField<T,fs> & smallField = * (out->get(fx,fy,fz) );
            for(cell_idx_t z = -gl; z < cell_idx_c(smallField.zSize())+gl; ++z )
               for(cell_idx_t y = -gl; y <  cell_idx_c(smallField.ySize())+gl; ++y )
                  for(cell_idx_t x = -gl; x <  cell_idx_c(smallField.xSize())+gl; ++x )
                     for(uint_t f = 0; f < fs; ++f)
                        smallField.get(x,y,z,f) = fieldToSplit.get(x+cell_idx_c(glbStartX),
                                                                             y+cell_idx_c(glbStartY),
                                                                             z+cell_idx_c(glbStartZ));

         }

   return out;
}


void singleField()
{
   const uint_t radius = 20;
   const uint_t fieldSize = uint_t(2.6 * radius);
   const uint_t midpoint  = uint_t(1.3 * radius);

   GhostLayerField<real_t,1> f (fieldSize,fieldSize,fieldSize,1,0);

   for(ScalarField::iterator i = f.begin(); i != f.end(); ++i)
   {
      uint_t dx = uint_c(i.x()) - midpoint;
      uint_t dy = uint_c(i.y()) - midpoint;
      uint_t dz = uint_c(i.z()) - midpoint;
      if( dx*dx + dy*dy + dz*dz < radius * radius)
         *i = 1;
   }

   TriangleMesh m;
   generateIsoSurface(f,0.5,m);
   //ofstream file ("out.obj");
   //writeMeshObj(file,m);
}

void multipleFields()
{
   const uint_t radius = 20;
   // the expectedDoubles where determined by an external tool (blender) and depend on the radius
   const uint_t expectedDoubles = 474;
   const uint_t fieldSize = uint_t(3 * radius);
   const uint_t midpoint  = uint_t(1.5 * radius);

   typedef GhostLayerField<real_t,1> GlF;
   GlF f (fieldSize,fieldSize,fieldSize,1,0);

   for(ScalarField::iterator i = f.begin(); i != f.end(); ++i)
   {
      uint_t dx = uint_c(i.x()) - midpoint;
      uint_t dy = uint_c(i.y()) - midpoint;
      uint_t dz = uint_c(i.z()) - midpoint;
      if( dx*dx + dy*dy + dz*dz < radius * radius)
         *i = 1;
   }

   std::vector<size_t> cuts (1, fieldSize/2);
   shared_ptr<Field<shared_ptr<GlF>,1 > > arr = splitGhostLayerField(f, cuts, cuts, cuts);

   TriangleMesh m;
   for( Field<shared_ptr<GlF>,1 >::iterator i = arr->begin(); i != arr->end(); ++i) {
      Vector3<real_t> offset(0);
      offset[0] = i.x() ? fieldSize * real_c(0.5) : real_t(0);
      offset[1] = i.y() ? fieldSize * real_c(0.5) : real_t(0);
      offset[2] = i.z() ? fieldSize * real_c(0.5) : real_t(0);

      generateIsoSurface( **i, real_c(0.5), m, Vector3<real_t>(real_t(1)), uint_t(0), offset );
   }

   //ofstream file1("meshWithDoubles.obj");
   //writeMeshObj(file1,m);
   size_t verticesRemoved =  m.removeDuplicateVertices();
   WALBERLA_CHECK_EQUAL(verticesRemoved, expectedDoubles);
   //ofstream file2("meshWithoutDoubles.obj");
   //writeMeshObj(file2,m);
}






int main()
{
   debug::enterTestMode();

   singleField();
   multipleFields();
   return 0;
}
