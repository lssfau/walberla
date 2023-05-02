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
//! \file QCriterionTest.cpp
//! \author Lukas Werner <lks.werner@fau.de>
//
//======================================================================================================================


#include "blockforest/all.h"
#include "core/all.h"
#include "field/all.h"
#include "lbm/all.h"

namespace walberla {

class Filter {
public:
   explicit Filter(uint_t numberOfCells) : numberOfCells_(numberOfCells) {}

   void operator()( const IBlock & /*block*/ ){

   }

   bool operator()( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z ) const {
      return x >= -1 && x <= cell_idx_t(numberOfCells_) &&
             y >= -1 && y <= cell_idx_t(numberOfCells_) &&
             z >= -1 && z <= cell_idx_t(numberOfCells_);
   }

private:
   uint_t numberOfCells_;
};

using VelocityField_T = GhostLayerField<Vector3<real_t>, 1>;
using FluidFilter_T = Filter;

int main( int argc, char ** argv )
{
   debug::enterTestMode();
   Environment env( argc, argv );

   auto numberOfCells = uint_t(40);

   VelocityField_T velocityField(numberOfCells, numberOfCells, numberOfCells, uint_t(1), field::zyxf);

   FluidFilter_T filter(numberOfCells);

   //// initialize field

   real_t xv = 0, yv = 0, zv = 0;

   for (auto cellIt = velocityField.beginWithGhostLayerXYZ(); cellIt != velocityField.end(); ++cellIt) {
      xv = math::realRandom<real_t>(real_t(-1), real_t(1));
      yv = math::realRandom<real_t>(real_t(-1), real_t(1));
      zv = math::realRandom<real_t>(real_t(-1), real_t(1));

      *cellIt = Vector3<real_t>(real_t(xv), real_t(yv), real_t(zv));
   }

   //// evaluate field

   const auto one = cell_idx_t(1);
   for (auto cellIt = velocityField.beginXYZ(); cellIt != velocityField.end(); ++cellIt) {
      cell_idx_t x = cellIt.x();
      cell_idx_t y = cellIt.y();
      cell_idx_t z = cellIt.z();

      // ParaView's formula serves as comparison.
      // From: https://github.com/Kitware/VTK/blob/8b08cf7a0523d88a5602a98c0456f7e397560698/Filters/General/vtkGradientFilter.cxx

      const Vector3<real_t> xa = velocityField.get(x+one,y,z);
      const Vector3<real_t> xb = velocityField.get(x-one,y,z);
      const Vector3<real_t> ya = velocityField.get(x,y+one,z);
      const Vector3<real_t> yb = velocityField.get(x,y-one,z);
      const Vector3<real_t> za = velocityField.get(x,y,z+one);
      const Vector3<real_t> zb = velocityField.get(x,y,z-one);

      const real_t duxdx = (xa[0] - xb[0]) * real_t(0.5);
      const real_t duxdy = (ya[0] - yb[0]) * real_t(0.5);
      const real_t duxdz = (za[0] - zb[0]) * real_t(0.5);

      const real_t duydx = (xa[1] - xb[1]) * real_t(0.5);
      const real_t duydy = (ya[1] - yb[1]) * real_t(0.5);
      const real_t duydz = (za[1] - zb[1]) * real_t(0.5);

      const real_t duzdx = (xa[2] - xb[2]) * real_t(0.5);
      const real_t duzdy = (ya[2] - yb[2]) * real_t(0.5);
      const real_t duzdz = (za[2] - zb[2]) * real_t(0.5);

      real_t q_paraview = -real_t(0.5)*(duxdx*duxdx + duydy*duydy + duzdz*duzdz)
                          -(duxdy*duydx + duxdz*duzdx + duydz*duzdy);

      WALBERLA_CHECK_FLOAT_EQUAL(q_paraview, lbm::getQCriterion(velocityField, filter, x, y, z));
   }

   return 0;
}
}

int main( int argc, char ** argv )
{
   return walberla::main(argc, argv);
}