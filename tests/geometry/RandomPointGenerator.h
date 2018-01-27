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

#include <random>
#include <vector>


template < typename scalar_t >
struct RandomPointGenerator
{
   typedef walberla::Vector3<scalar_t> vector_t;
   typedef std::mersenne_twister_engine< walberla::uint32_t, 32, 351, 175, 19, 0xccab8ee7, 11, 0xffffffff, 7, 0x31b6ab00, 15, 0xffe50000, 17, 0xa37d3c92 > mt11213b;
   typedef mt11213b RandomNumberEngine;
   typedef std::normal_distribution<scalar_t> NormalDistribution;

   RandomPointGenerator( const vector_t & mu, const vector_t & sigma )
   {
      for( walberla::uint_t i = 0; i < 3; ++i )
      {
         engines.push_back( RandomNumberEngine() );
         normalDistributions.push_back( NormalDistribution( mu[i], sigma[i] ) );
      }
   }

   vector_t operator()()
   {
      return vector_t( normalDistributions[0](engines[0]), normalDistributions[1](engines[1]), normalDistributions[2](engines[2]) );
   }

private:
   std::vector< RandomNumberEngine > engines;
   std::vector< NormalDistribution > normalDistributions;
};
