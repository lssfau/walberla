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
//! \file SCIterator.h
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include <iterator>

namespace walberla {
namespace grid_generator {

/// Helper class to generate points in a simple cubic structure within a certain domain.
/// The lattice is fixed by a point of reference (x).
/// \code
/// . . . . . . . .
///        +-----+
/// . . . .|. . .|.
///        |     |
/// . . . .|. . .|.
///        +-----+
/// . . . . . . . .
///
/// . x . . . . . .
///
/// . . . . . . . .
/// \endcode
/// Usage:
/// \code for (auto it = SCIterator::begin(...); it != SCIterator::end(); ++it) \endcode
class SCIterator : public std::iterator< std::forward_iterator_tag, Vector3<real_t> >
{
public:
   /**
    * @brief begin iterator
    * @param domain volume were lattice points will be returned
    * @param pointOfReference point somewhere in the world which fixes the lattice
    * @param spacing spacing between grid points in x, y and z direction
    */
   SCIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing);
   /**
    * @brief begin iterator
    * @param domain volume were lattice points will be returned
    * @param pointOfReference point somewhere in the world which fixes the lattice
    * @param spacing spacing between grid points in x, y and z direction
    */
   SCIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const Vector3<real_t>& spacing);
   /**
    * @brief end iterator
    */
   SCIterator();

   SCIterator& operator++();
   SCIterator  operator++(int);
   Vector3<real_t> operator*() const;
   bool operator==(const SCIterator& rhs) const;
   bool operator!=(const SCIterator& rhs) const;

   static inline real_t getUnitCellX(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellY(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellZ(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellX(const Vector3<real_t> spacing) { return spacing[0]; }
   static inline real_t getUnitCellY(const Vector3<real_t> spacing) { return spacing[1]; }
   static inline real_t getUnitCellZ(const Vector3<real_t> spacing) { return spacing[2]; }

private:
   void updatePoint();

   int i_;
   int iReturn_;
   int j_;
   int jReturn_;
   int k_;

   AABB aabb_;
   Vector3<real_t> pointOfReference_;
   Vector3<real_t> spacing_;

   Vector3<real_t> point_;

   bool ended_;
};

} // namespace grid_generator
} // namespace walberla
