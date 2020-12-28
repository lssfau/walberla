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
//! \file HCPIterator.h
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

/// Helper class to generate points in a hexagonal close packing structure within a certain domain
///
/// Usage:
/// \code for (auto it = HCPIterator::begin(...); it != HCPIterator::end(); ++it) \endcode
class HCPIterator
{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = Vector3<real_t>;
    using difference_type = std::ptrdiff_t;
    using pointer = Vector3<real_t>*;
    using reference = Vector3<real_t>&;

   /**
    * @brief begin iterator
    * @param domain volume were lattice points will be returned
    * @param pointOfReference point somewhere in the world which fixes the lattice
    * @param spacing spacing between grid points in x direction
    */
   HCPIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing);
   /**
    * @brief end iterator
    */
   HCPIterator();

   HCPIterator& operator++();
   HCPIterator  operator++(int);
   Vector3<real_t> operator*() const;
   bool operator==(const HCPIterator& rhs) const;
   bool operator!=(const HCPIterator& rhs) const;

   static inline real_t getUnitCellX(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellY(const real_t spacing) { return real_c(sqrt(3.0) * 2) * spacing; }
   static inline real_t getUnitCellZ(const real_t spacing) { return real_c(2 * sqrt(6.0) / 3.0) * spacing; }

private:
   void updatePoint();

   unsigned int i_;
   unsigned int iReturn_;
   unsigned int j_;
   unsigned int jReturn_;
   unsigned int k_;

   AABB aabb_;
   Vector3<real_t> pointOfReference_;
   real_t radius_;

   Vector3<real_t> point_;

   bool ended_;
};

/// Convenience class to enable range based for loops over grid points.
/// Usage:
/// \code for (const auto& pt : HCPGrid(...) ) \endcode
class HCPGrid
{
public:
   using iterator = HCPIterator;
   using value_type = iterator::value_type;

   /**
    * @param domain volume were lattice points will be returned
    * @param pointOfReference point somewhere in the world which fixes the lattice
    * @param spacing spacing between grid points in x direction
    */
   HCPGrid(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing)
      : domain_(domain)
      , pointOfReference_(pointOfReference)
      , spacing_(spacing)
   {}

   HCPIterator begin() {return HCPIterator(domain_, pointOfReference_, spacing_);}
   HCPIterator begin()  const {return HCPIterator(domain_, pointOfReference_, spacing_);}
   HCPIterator cbegin() const {return HCPIterator(domain_, pointOfReference_, spacing_);}

   HCPIterator end() {return HCPIterator();}
   HCPIterator end()  const {return HCPIterator();}
   HCPIterator cend() const {return HCPIterator();}

private:
   AABB domain_;
   Vector3<real_t> pointOfReference_;
   real_t spacing_;
};

} // namespace grid_generator
} // namespace walberla
