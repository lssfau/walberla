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
//! \file RandomPointGenerator.h
//! \ingroup geometry
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/math/Vector3.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <vector>


template < typename scalar_t >
struct RandomPointGenerator
{
   typedef walberla::Vector3<scalar_t> vector_t;
   typedef boost::mt11213b RandomNumberEngine;
   typedef boost::normal_distribution<scalar_t> NormalDistribution;
   typedef boost::variate_generator< RandomNumberEngine, NormalDistribution > Generator;

   RandomPointGenerator( const vector_t & mu, const vector_t & sigma )
   {
      for( walberla::uint_t i = 0; i < 3; ++i )
         normalDistributions.push_back( Generator( RandomNumberEngine(), NormalDistribution( mu[i], sigma[i] ) ) );
   }

   vector_t operator()()
   {
      return vector_t( normalDistributions[0](), normalDistributions[1](), normalDistributions[2]() );
   }

private:
   std::vector< Generator > normalDistributions;
};
