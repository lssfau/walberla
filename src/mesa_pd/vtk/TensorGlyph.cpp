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
//! \file TensorGlyph.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include <mesa_pd/vtk/TensorGlyph.h>

namespace walberla {
namespace mesa_pd {
namespace vtk {

TensorGlyph createTensorGlyph(const Vec3& semiAxes, const Rot3& rot)
{
   // compute tensor glyph for visualization with ParaView (tensorGlyph)
   const Mat3& rotMat = rot.getMatrix();
   Vector3<real_t> directionVectorX(rotMat[0], rotMat[3], rotMat[6]);
   Vector3<real_t> directionVectorY(rotMat[1], rotMat[4], rotMat[7]);
   Vector3<real_t> directionVectorZ(rotMat[2], rotMat[5], rotMat[8]);
   Mat3 axa = math::dyadicProduct(directionVectorX, directionVectorX);
   Mat3 bxb = math::dyadicProduct(directionVectorY, directionVectorY);
   Mat3 cxc = math::dyadicProduct(directionVectorZ, directionVectorZ);
   Mat3 tensor = axa * semiAxes[0] + bxb * semiAxes[1] + cxc * semiAxes[2];
   // use symmetry to only write 6 of the 9 elements: XX YY ZZ XY YZ XZ
   return {{tensor(0,0), tensor(1,1), tensor(2,2), tensor(0,1), tensor(1,2), tensor(0,2)}};
}

} // namespace vtk
} // namespace pe
} // namespace walberla
