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
//! \file CurveGatherPackInfo.impl.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Implementation of CurveGatherPackInfo
//
//======================================================================================================================

#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/Parser.h"
#include "core/mpi/MPIManager.h"

#include <algorithm>
#include <cmath>


namespace walberla {
namespace gather {


//======================================================================================================================
//
//  Constructors
//
//======================================================================================================================


template< typename GlF, typename IP>
CurveGatherPackInfo<GlF,IP>::CurveGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                                  ConstBlockDataID fieldID,
                                                  const std::string & cx, const std::string & cy, const std::string & cz,
                                                  real_t tStart, real_t tEnd, uint_t numSamples,
                                                  const shared_ptr<DataProcessor> & dp )
   : blocks_(bs), fieldID_(fieldID), dataProcessor_(dp)
{
   WALBERLA_ASSERT_LESS( tStart, tEnd );

   real_t tIncr = ( tEnd - tStart ) / real_c(numSamples-1);

   using math::FunctionParser;

   FunctionParser parserX;
   FunctionParser parserY;
   FunctionParser parserZ;
   parserX.parse(cx);
   parserY.parse(cy);
   parserZ.parse(cz);

   std::map<std::string, double> symbolTable;

   for(real_t t=tStart; t <= tEnd; t+=tIncr )
   {
      symbolTable["t"] = double(t);
      RealVec3 p ( real_c( parserX.evaluate(symbolTable) ),
                   real_c( parserY.evaluate(symbolTable) ),
                   real_c( parserZ.evaluate(symbolTable) ) );

      addSample(t,p);
   }
}

template< typename GlF, typename IP>
CurveGatherPackInfo<GlF,IP>::CurveGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                                  ConstBlockDataID fieldID,
                                                  CurveCallback curveCallback,
                                                  real_t tStart, real_t tEnd, uint_t numSamples,
                                                  const shared_ptr<DataProcessor> & dp )
   : blocks_(bs), fieldID_(fieldID), dataProcessor_(dp)
{
   WALBERLA_ASSERT_LESS( tStart, tEnd );
   real_t tIncr = (tEnd-tStart) / real_c(numSamples-1);

   // sample along curve
   for(real_t t=tStart; t <= tEnd; t+=tIncr )
   {
      RealVec3 p = curveCallback(t);
      addSample(t,p);
   }
}


template< typename GlF, typename IP>
CurveGatherPackInfo<GlF,IP>::CurveGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                                                  ConstBlockDataID fieldID,
                                                  const std::vector<Vector3<real_t> > & inpPoints,
                                                  const shared_ptr<DataProcessor> & dp )
   : blocks_(bs), fieldID_(fieldID), dataProcessor_(dp)
{
   for( unsigned int i=0; i<inpPoints.size(); ++i )
      addSample( real_c(i),inpPoints[i] );
}


template< typename GlF, typename IP>
void CurveGatherPackInfo<GlF,IP>::addSample( real_t t, const RealVec3 & p )
{
   if ( ! blocks_->getDomain().contains( p[0], p[1], p[2] ) ) {
      WALBERLA_LOG_WARNING("CurveGatherPackInfo: Skipping " << p << " because does not fit into domain " << blocks_->getDomain() );
      return;
   }

   // Find block for global coordinate and translate sample to local cell coordinates
   IBlock * block = blocks_->getBlock( p[0], p[1], p[2] );
   if ( ! block ) // non-local block
      return;

   const AABB & blockBB = block->getAABB();
   RealVec3 pLocal ( p[0] - blockBB.xMin(),
                     p[1] - blockBB.yMin(),
                     p[2] - blockBB.zMin() );

   pLocal[0] /= blocks_->dx();
   pLocal[1] /= blocks_->dy();
   pLocal[2] /= blocks_->dz();

   Samples & s = localSamplePoints_[ block ];
   s.t.push_back(t);
   s.coord.push_back( pLocal );
}




//======================================================================================================================
//
//  Packing Interface
//
//======================================================================================================================



template< typename GlF, typename Interpolator>
void CurveGatherPackInfo<GlF,Interpolator>::packData( const IBlock * sender, mpi::SendBuffer & outBuffer )
{
   auto i = localSamplePoints_.find( sender );

   if(i == localSamplePoints_.end() ) //nothing to send
      return;

   const std::vector<RealVec3> & points = i->second.coord;
   const std::vector<real_t> & t        = i->second.t;
   WALBERLA_ASSERT_EQUAL(t.size(), points.size() );

   const GlF* field = sender->getData<GlF>( fieldID_ );
   Interpolator ip ( *field );
   const size_t fieldSize = Interpolator::F_SIZE;

   outBuffer << points.size();
   outBuffer << fieldSize;

   // Write data for each component of the field
   for(unsigned int j=0; j< t.size(); ++j)
   {
      outBuffer << t[j];
      for(uint_t f=0; f< fieldSize; ++f )
         outBuffer << ip( points[j][0], points[j][1], points[j][2], cell_idx_c( f ) );
   }

   // Message format ( example with fieldSize=#f=4 and  points.size()=#Points=2  )
   // #Points |  #f  | t0 (f0 f1 f2 f3) | t1 (f0 f1 f2 f3)
}

template< typename GlF, typename IP>
void CurveGatherPackInfo<GlF,IP>::unpackData( mpi::RecvBuffer & buffer )
{
   size_t nrPoints;
   buffer >> nrPoints;
   size_t fieldSize;
   buffer >> fieldSize;

   for( size_t i=0; i< nrPoints; ++i )
   {
      receivedData.push_back(std::vector<real_t>(fieldSize+1)); //+1 because we also store t value as first entry
      std::vector<real_t> & pointVec = receivedData[receivedData.size()-1];

      real_t t;
      real_t val;

      buffer >> t;
      pointVec[0] = t;
      for( size_t f=0; f<fieldSize; ++f )
      {
         buffer >> val;
         pointVec[f+1] = val;
      }
   }
}


template< typename GlF, typename IP>
void CurveGatherPackInfo<GlF,IP>::gatherFinished()
{
   sortReceivedData();
   dataProcessor_->process( receivedData );
   receivedData.clear();
}



template< typename GlF, typename IP>
void CurveGatherPackInfo<GlF,IP>::sortReceivedData()
{
   //sort according to "t" value, which is the first entry of the inner vector
   std::sort(receivedData.begin(), receivedData.end(), Compare() );

#ifndef NDEBUG
   // check that we have collected everything,
   // therefore all the t values have to have equal distance
   if (receivedData.size() > 1) {
      real_t dist = receivedData[1][0] - receivedData[0][0];

      for(uint_t i=1; i<receivedData.size(); ++i ) {
         WALBERLA_ASSERT_FLOAT_EQUAL( receivedData[i][0], receivedData[i-1][0] + dist );
      }
   }
#endif
}





} // namespace gather
} // namespace walberla



