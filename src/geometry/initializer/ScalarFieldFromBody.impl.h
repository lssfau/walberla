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
//! \file ScalarFieldFromBody.cpp
//! \ingroup geometry
//! \author Michael Kuron <mkuron@icp.uni-stuttgart.de>
//
//======================================================================================================================

#include "geometry/bodies/AABBBody.h"
#include "geometry/bodies/BodyFromConfig.h"
#include "geometry/bodies/Cylinder.h"
#include "geometry/bodies/Ellipsoid.h"
#include "geometry/bodies/Sphere.h"
#include "geometry/bodies/Torus.h"

#include "core/Abort.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/mpi/Reduce.h"
#include "core/stringToNum.h"


namespace walberla {
namespace geometry {
namespace initializer {


   template <typename Field_T>
   ScalarFieldFromBody<Field_T>::ScalarFieldFromBody( StructuredBlockStorage & structuredBlockStorage, std::vector<BlockDataID> scalarFieldID )
      : structuredBlockStorage_( structuredBlockStorage ), scalarFieldID_( scalarFieldID ), addKeyword_("add"), setKeyword_("set")
   {}


   template <typename Field_T>
   void ScalarFieldFromBody<Field_T>::init( BlockStorage & /*blockStorage*/, const Config::BlockHandle & subBlock )
   {
      bool addOrSet = true;

      bool addDefined      = subBlock.isDefined( addKeyword_ );
      bool setDefined      = subBlock.isDefined( setKeyword_ );

      if ( addDefined && setDefined )
         WALBERLA_ABORT( "Specify only one of " << addKeyword_ << " and " << setKeyword_ << "!\n"
                         << "Both are defined.");

      if ( setDefined )
         addOrSet = false;
      
      auto        id         = subBlock.getParameter< std::vector<BlockDataID>::size_type > ( "id", 0 );
      std::string shape      = subBlock.getParameter< std::string >                         ( "shape" );
      std::string expression = subBlock.getParameter< std::string >                         ( "value" );
      
      try
      {
         Value_T value = stringToNum<Value_T>(expression);
         init ( *bodyFromConfig ( subBlock ), value, addOrSet, id );
      }
      
      catch(std::invalid_argument&)
      {
         math::FunctionParser p;
         p.parse(expression);         
         init ( *bodyFromConfig ( subBlock ), p, addOrSet, id );
      }
   }


   template< typename Field_T >
   template< typename Body >
   void ScalarFieldFromBody<Field_T>::init( const Body & body, Value_T value, bool addOrSet, std::vector<BlockDataID>::size_type id )
   {
      real_t V = 0;

      for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
      {
         IBlock * block = &(*blockIt);
         
         const uint_t level = structuredBlockStorage_.getLevel(*block);
         const real_t dx = structuredBlockStorage_.dx() / real_c(1 << level);
         const real_t dy = structuredBlockStorage_.dy() / real_c(1 << level);
         const real_t dz = structuredBlockStorage_.dz() / real_c(1 << level);
         const real_t dV = dx*dy*dz;

         auto ff = block->getData<Field_T>( scalarFieldID_[id] );
         auto gl = cell_idx_c( ff->nrOfGhostLayers() );

         // If Block (extended with ghost layers) does not intersect body - skip the complete block
         AABB blockBB = block->getAABB();
         blockBB.extend( math::Vector3< real_t >( dx * real_c( gl ), dy * real_c( gl ), dz * real_c( gl ) ) );
         if( fastOverlapCheck( body, blockBB ) == geometry::COMPLETELY_OUTSIDE )
            continue;

         AABB firstCellBB;
         structuredBlockStorage_.getBlockLocalCellAABB( *block, ff->beginWithGhostLayer().cell(), firstCellBB );
         Vector3<real_t> firstCellMidpoint;
         for( uint_t i = 0; i < 3; ++i )
            firstCellMidpoint[i] = firstCellBB.min(i) + real_t(0.5) * firstCellBB.size(i);

         Vector3<real_t> currentMidpoint;
         currentMidpoint[2] = firstCellMidpoint[2];
         for ( cell_idx_t z = -gl; z < cell_idx_c(ff->zSize())+gl; ++z, currentMidpoint[2] += dz )
         {
            currentMidpoint[1] = firstCellMidpoint[1];
            for ( cell_idx_t y = -gl; y < cell_idx_c(ff->ySize())+gl; ++y, currentMidpoint[1] += dy )
            {
               currentMidpoint[0] = firstCellMidpoint[0];
               for( cell_idx_t x = -gl; x < cell_idx_c(ff->xSize())+gl; ++x, currentMidpoint[0] += dx )
               {
                  if ( !contains(body, currentMidpoint) )
                     continue;
                  
                  if (ff->isInInnerPart(Cell(x,y,z)))
                     V += dV;
                  real_t & val = ff->get(x,y,z);

                  if ( addOrSet )
                     val += value;
                  else
                     val =  value;

               }
            }
         }
      }
      mpi::allReduceInplace(V, mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("Total body volume: " << V);
   }
   
   template< typename Field_T >
   template< typename Body >
   void ScalarFieldFromBody<Field_T>::init( const Body & body, math::FunctionParser & parser, bool addOrSet, std::vector<BlockDataID>::size_type id )
   {
      real_t V = 0;

      for( auto blockIt = structuredBlockStorage_.begin(); blockIt != structuredBlockStorage_.end(); ++blockIt )
      {
         IBlock * block = &(*blockIt);
         
         const uint_t level = structuredBlockStorage_.getLevel(*block);
         const real_t dx = structuredBlockStorage_.dx() / real_c(1 << level);
         const real_t dy = structuredBlockStorage_.dy() / real_c(1 << level);
         const real_t dz = structuredBlockStorage_.dz() / real_c(1 << level);
         const real_t dV = dx*dy*dz;

         auto ff = block->getData<Field_T>( scalarFieldID_[id] );
         auto gl = cell_idx_c( ff->nrOfGhostLayers() );

         // If Block (extended with ghost layers) does not intersect body - skip the complete block
         AABB blockBB = block->getAABB();
         blockBB.extend( math::Vector3< real_t >( dx * real_c( gl ), dy * real_c( gl ), dz * real_c( gl ) ) );
         if( fastOverlapCheck( body, blockBB ) == geometry::COMPLETELY_OUTSIDE )
            continue;

         AABB firstCellBB;
         structuredBlockStorage_.getBlockLocalCellAABB( *block, ff->beginWithGhostLayer().cell(), firstCellBB );
         Vector3<real_t> firstCellMidpoint;
         for( uint_t i = 0; i < 3; ++i )
            firstCellMidpoint[i] = firstCellBB.min(i) + real_t(0.5) * firstCellBB.size(i);

         Vector3<real_t> currentMidpoint;
         currentMidpoint[2] = firstCellMidpoint[2];
         for ( cell_idx_t z = -gl; z < cell_idx_c(ff->zSize())+gl; ++z, currentMidpoint[2] += dz )
         {
            currentMidpoint[1] = firstCellMidpoint[1];
            for ( cell_idx_t y = -gl; y < cell_idx_c(ff->ySize())+gl; ++y, currentMidpoint[1] += dy )
            {
               currentMidpoint[0] = firstCellMidpoint[0];
               for( cell_idx_t x = -gl; x < cell_idx_c(ff->xSize())+gl; ++x, currentMidpoint[0] += dx )
               {
                  if ( !contains(body, currentMidpoint) )
                     continue;
                  
                  Cell cell(x,y,z);
                  if (ff->isInInnerPart(cell))
                     V += dV;
                  real_t & val = ff->get(x,y,z);
                  
                  std::map<std::string,Value_T> params;
                  structuredBlockStorage_.transformBlockLocalToGlobalCell(cell, *block);
                  params["x"] = cell.x();
                  params["y"] = cell.y();
                  params["z"] = cell.z();
                  Value_T value = parser.evaluate(params);

                  if ( addOrSet )
                     val += value;
                  else
                     val =  value;

               }
            }
         }
      }
      mpi::allReduceInplace(V, mpi::SUM);
      WALBERLA_LOG_INFO_ON_ROOT("Total body volume: " << V);
   }

} // namespace initializer
} // namespace geometry
} // namespace walberla
