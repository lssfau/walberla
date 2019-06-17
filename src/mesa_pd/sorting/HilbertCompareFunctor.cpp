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
//! \file HilbertCompareFunctor.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "HilbertCompareFunctor.h"

#include <blockforest/HilbertCurveConstruction.h>
#include <core/DataTypes.h>
#include <core/Environment.h>
#include <core/logging/Logging.h>
#include <core/math/Vector3.h>
#include <core/mpi/MPIManager.h>

#include <stack>
#include <vector>

namespace walberla {
namespace mesa_pd {
namespace sorting {

math::Vector3<uint_t> getChildPosition(const uint_t currentSize, const math::Vector3<uint_t>& parent, const uint_t childId)
{
   auto size = currentSize >> 1;
   auto child = parent;

   if (!(childId & 0b001))
   {
      child[0] -= size;
   }

   if (!(childId & 0b010))
   {
      child[1] -= size;
   }

   if (!(childId & 0b100))
   {
      child[2] -= size;
   }

   return child;
}

std::vector< math::Vector3<uint_t> > generateHilbertCurve(const uint_t size)
{
   using Node = std::pair<uint_t, Vector3<uint_t>>; //level, coordinate

   std::vector< math::Vector3<uint_t> > result;

   std::stack< Node > stack;
   std::stack< uint_t > orientation;

   stack.push( std::make_pair(size, math::Vector3<uint_t>(size)) );
   orientation.push( uint_t(0) );

   while( !stack.empty() )
   {
      const Node node = stack.top();
      uint_t hilbertIndex = orientation.top();

      stack.pop();
      orientation.pop();

      if( node.first > 1 )
      {
         for( uint_t c = 8; c-- != 0; )
         {
            stack.push( std::make_pair( node.first >> 1,
                                        getChildPosition(node.first, node.second, blockforest::hilbertOrder[hilbertIndex][c]) ));
            orientation.push( blockforest::hilbertOrientation[hilbertIndex][c] );
         }
      }
      else
      {
         result.push_back( node.second );
      }
   }
   return result;
}


HilbertCompareFunctor::HilbertCompareFunctor(const math::AABB& domain, const uint_t cells)
   : domain_(domain)
   , cells_(cells)
   , hilbertLookup_(cells_ * cells_ * cells_)
{
   WALBERLA_CHECK(math::uintIsPowerOfTwo(cells));
   inverse_dx[0] = real_t(1.0) / (domain_.xSize() / real_c(cells_));
   inverse_dx[1] = real_t(1.0) / (domain_.ySize() / real_c(cells_));
   inverse_dx[2] = real_t(1.0) / (domain_.zSize() / real_c(cells_));
   initializeLookup();
}

bool HilbertCompareFunctor::operator()(const data::Particle p1, const data::Particle p2) const
{
   const auto hash1 = discretize(p1.getPosition() - domain_.minCorner());
   WALBERLA_ASSERT_LESS( hash1, cells_*cells_*cells_);
   const auto hash2 = discretize(p2.getPosition() - domain_.minCorner());
   WALBERLA_ASSERT_LESS( hash2, cells_*cells_*cells_);

   return hilbertLookup_[hash1] < hilbertLookup_[hash2];
}

uint_t HilbertCompareFunctor::discretize(const Vec3& pos) const
{
   const real_t& x = pos[0];
   const real_t& y = pos[1];
   const real_t& z = pos[2];

   const uint_t hashMask = (cells_ - 1);

   math::Vector3<uint_t> hash;

   if( x < 0 ) {
      real_t i = ( -x ) * inverse_dx[0];
      hash[0]  = cells_ - 1 - ( static_cast<uint_t>( i ) & hashMask );
   }
   else {
      real_t i = x * inverse_dx[0];
      hash[0]  = static_cast<uint_t>( i ) & hashMask;
   }

   if( y < 0 ) {
      real_t i = ( -y ) * inverse_dx[1];
      hash[1]  = cells_ - 1 - ( static_cast<uint_t>( i ) & hashMask );
   }
   else {
      real_t i = y * inverse_dx[1];
      hash[1]  = static_cast<uint_t>( i ) & hashMask;
   }

   if( z < 0 ) {
      real_t i = ( -z ) * inverse_dx[2];
      hash[2]  = cells_ - 1 - ( static_cast<uint_t>( i ) & hashMask );
   }
   else {
      real_t i = z * inverse_dx[2];
      hash[2]  = static_cast<uint_t>( i ) & hashMask;
   }

   WALBERLA_ASSERT_LESS(hash[0], cells_);
   WALBERLA_ASSERT_LESS(hash[1], cells_);
   WALBERLA_ASSERT_LESS(hash[2], cells_);

   return hash[2] * cells_ * cells_ + hash[1] * cells_ + hash[0];
}


void HilbertCompareFunctor::initializeLookup()
{
   const auto sfc = generateHilbertCurve(cells_);
   uint_t counter = 0;
   for (auto p : sfc)
   {
      WALBERLA_ASSERT_GREATER(p[0], 0);
      WALBERLA_ASSERT_GREATER(p[1], 0);
      WALBERLA_ASSERT_GREATER(p[2], 0);
      p[0] -= 1;
      p[1] -= 1;
      p[2] -= 1;

      hilbertLookup_[p[2] * cells_ * cells_ + p[1] * cells_ + p[0]] = counter;
      ++counter;
   }
   WALBERLA_ASSERT_EQUAL(counter, sfc.size());
   WALBERLA_ASSERT_EQUAL(counter, hilbertLookup_.size());
}

} //namespace sorting
} //namespace mesa_pd
} //namespace walberla
