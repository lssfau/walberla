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
//! \file SCIterator.cpp
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "core/debug/Debug.h"
#include "core/DataTypes.h"
#include "core/math/AABB.h"
#include "core/math/Vector3.h"

#include <iterator>

namespace walberla {
namespace grid_generator {

/// Helper class to generate points in a simple cubic structure within a certain domain
///
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
   inline SCIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing);
   /**
    * @brief begin iterator
    * @param domain volume were lattice points will be returned
    * @param pointOfReference point somewhere in the world which fixes the lattice
    * @param spacing spacing between grid points in x, y and z direction
    */
   inline SCIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const Vector3<real_t>& spacing);
   /**
    * @brief end iterator
    */
   inline SCIterator();

   inline SCIterator& operator++();
   inline SCIterator  operator++(int);
   inline Vector3<real_t> operator*() const;
   inline bool operator==(const SCIterator& rhs) const;
   inline bool operator!=(const SCIterator& rhs) const;

   static inline real_t getUnitCellX(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellY(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellZ(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellX(const Vector3<real_t> spacing) { return spacing[0]; }
   static inline real_t getUnitCellY(const Vector3<real_t> spacing) { return spacing[1]; }
   static inline real_t getUnitCellZ(const Vector3<real_t> spacing) { return spacing[2]; }

private:
   inline void updatePoint();

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

SCIterator::SCIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing)
   : i_(0)
   , j_(0)
   , k_(0)
   , aabb_( domain )
   , pointOfReference_( pointOfReference )
   , spacing_( spacing, spacing, spacing )
   , ended_(false)
{
   auto min = domain.min() - pointOfReference_;
   iReturn_ = int_c (ceil( min[0] / spacing ));
   i_ = iReturn_;
   jReturn_ = int_c (ceil( min[1] / spacing ));
   j_ = jReturn_;
   k_ = int_c (ceil( min[2] / spacing ));

   updatePoint();

   if (!aabb_.contains(point_))
      ended_ = true;
}

SCIterator::SCIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const Vector3<real_t>& spacing)
   : i_(0)
   , j_(0)
   , k_(0)
   , aabb_( domain )
   , pointOfReference_( pointOfReference )
   , spacing_( spacing )
   , ended_(false)
{
   auto min = domain.min() - pointOfReference_;
   iReturn_ = int_c( ceil( min[0] / spacing[0] ) );
   i_ = iReturn_;
   jReturn_ = int_c(ceil( min[1] / spacing[1] ) + real_c(0.1));
   j_ = jReturn_;
   k_ = int_c (ceil( min[2] / spacing[2] ));

   updatePoint();

   if (!aabb_.contains(point_))
      ended_ = true;
}

SCIterator::SCIterator()
   : ended_(true)
{}

SCIterator& SCIterator::operator++()
{
   WALBERLA_ASSERT( aabb_.contains(point_), "Trying to increment an already out of bounds iterator!");

   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;

   i_ = iReturn_;
   ++j_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;

   j_ = jReturn_;
   ++k_;
   updatePoint();
   if (!aabb_.contains(point_)) ended_ = true;

   return *this;
}

SCIterator SCIterator::operator++(int)
{
   SCIterator temp(*this);
   operator++();
   return temp;
}

Vector3<real_t> SCIterator::operator*() const
{
   return point_;
}

bool SCIterator::operator==(const SCIterator &rhs) const
{
   if (ended_ || rhs.ended_)
   {
      if (ended_ == rhs.ended_)
      {
         return true;
      }
      else
      {
         return false;
      }
   }

//   WALBERLA_ASSERT_FLOAT_EQUAL(aabb_, rhs.aabb_, "Comparing iterators for different starting configurations!");
   WALBERLA_ASSERT_FLOAT_EQUAL(pointOfReference_, rhs.pointOfReference_, "Comparing iterators for different starting configurations!");
   WALBERLA_ASSERT_FLOAT_EQUAL(spacing_, rhs.spacing_, "Comparing iterators for different starting configurations!");

   return (i_==rhs.i_) && (j_==rhs.j_) && (k_==rhs.k_);
}

bool SCIterator::operator!=(const SCIterator &rhs) const
{
   return !(*this == rhs);
}

void SCIterator::updatePoint()
{
   point_ = pointOfReference_ + Vector3<real_t>( real_c(i_) * spacing_[0], real_c(j_) * spacing_[1], real_c(k_) * spacing_[2] );
}

}
}
