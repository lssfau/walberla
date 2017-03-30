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
//! \file HCPIterator.cpp
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

/// Helper class to generate points in a hexagonal close packing structure within a certain domain
///
/// Usage:
/// \code for (auto it = HCPIterator::begin(...); it != HCPIterator::end(); ++it) \endcode
class HCPIterator : public std::iterator< std::forward_iterator_tag, Vector3<real_t> >
{
public:
   /**
    * @brief begin iterator
    * @param domain volume were lattice points will be returned
    * @param pointOfReference point somewhere in the world which fixes the lattice
    * @param spacing spacing between grid points in x direction
    */
   inline HCPIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing);
   /**
    * @brief end iterator
    */
   inline HCPIterator();

   inline HCPIterator& operator++();
   inline HCPIterator  operator++(int);
   inline Vector3<real_t> operator*() const;
   inline bool operator==(const HCPIterator& rhs) const;
   inline bool operator!=(const HCPIterator& rhs) const;

   static inline real_t getUnitCellX(const real_t spacing) { return spacing; }
   static inline real_t getUnitCellY(const real_t spacing) { return real_c(sqrt(3.0) * 2) * spacing; }
   static inline real_t getUnitCellZ(const real_t spacing) { return real_c(2 * sqrt(6.0) / 3.0) * spacing; }

private:
   inline void updatePoint();

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

HCPIterator::HCPIterator(const AABB& domain, const Vector3<real_t>& pointOfReference, const real_t spacing)
   : i_(0)
   , j_(0)
   , k_(0)
   , aabb_( domain )
   , pointOfReference_( pointOfReference )
   , radius_(real_c(0.5) * spacing)
   , ended_(false)
{
   auto min = domain.min() - pointOfReference_;

//   if ( min[0] <= 2 * spacing ) throw std::invalid_argument("Domain needs to have a certain distance to origin!");
//   if ( min[1] <= 2 * spacing ) throw std::invalid_argument("Domain needs to have a certain distance to origin!");
//   if ( min[2] <= 2 * spacing ) throw std::invalid_argument("Domain needs to have a certain distance to origin!");

   while ( min[0] <= 2 * spacing )
   {
      pointOfReference_[0] -= 2 * spacing;
      min[0] = domain.min()[0] - pointOfReference_[0];
   }

   while ( min[1] <= 2 * spacing )
   {
      pointOfReference_[1] -= real_c(sqrt(3.0) * 2) * spacing;
      min[1] = domain.min()[1] - pointOfReference_[1];
   }

   while ( min[2] <= 2 * spacing )
   {
      pointOfReference_[2] -= real_c(2 * sqrt(6.0) / 3.0) * spacing;
      min[2] = domain.min()[2] - pointOfReference_[2];
   }

   k_ = static_cast<unsigned int>(ceil(min[2] / radius_ * real_c(3.0 / (2.0 * sqrt(6.0)))));
   jReturn_ = static_cast<unsigned int>(
            ceil( min[1] / (real_c(sqrt(3.0)) * radius_) - real_c(1.0 / 3.0) * real_c(k_ % 2) )
         );
   j_ = jReturn_;
   iReturn_ = static_cast<unsigned int>(
            ceil( (min[0] / radius_ - real_c( (j_ + k_) % 2)) * 0.5 )
         );
   i_ = iReturn_;

   updatePoint();

   if (!aabb_.contains(point_))
      ended_ = true;
}

HCPIterator::HCPIterator()
   : ended_(true)
{}

HCPIterator& HCPIterator::operator++()
{
   //WALBERLA_ASSERT( aabb_.contains(point_), "Trying to increment an already out of bounds iterator!");

   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;

   i_ = iReturn_ - 1;
   ++j_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   i_ = iReturn_ - 1;
   ++j_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;

   i_ = iReturn_ - 1;
   j_ = jReturn_ - 1;
   ++k_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   i_ = iReturn_ - 1;
   ++j_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   i_ = iReturn_ - 1;
   ++j_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;
   ++i_;
   updatePoint();
   if (aabb_.contains(point_)) return *this;

   ended_ = true;

   return *this;
}

HCPIterator HCPIterator::operator++(int)
{
   HCPIterator temp(*this);
   operator++();
   return temp;
}

Vector3<real_t> HCPIterator::operator*() const
{
   return point_;
}

bool HCPIterator::operator==(const HCPIterator &rhs) const
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
   WALBERLA_ASSERT_FLOAT_EQUAL(radius_, rhs.radius_, "Comparing iterators for different starting configurations!");

   return (i_==rhs.i_) && (j_==rhs.j_) && (k_==rhs.k_);
}

bool HCPIterator::operator!=(const HCPIterator &rhs) const
{
   return !(*this == rhs);
}

void HCPIterator::updatePoint()
{
   point_ = pointOfReference_ + Vector3<real_t>( real_c( 2 * i_ + ( ( j_ + k_ ) % 2) ),
                                       real_c( std::sqrt(3.0) ) * real_c( ( real_c( j_ ) + real_c( 1.0 / 3.0) * real_c( k_ % 2 ) ) ),
                                       real_c(2) * real_c( std::sqrt( 6.0 ) ) / real_c(3.0) * real_c(k_) ) * radius_;
}

}
}
