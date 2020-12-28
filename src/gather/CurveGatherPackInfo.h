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
//! \file CurveGatherPackInfo.h
//! \ingroup gather
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Packs data along a Curve
//
//======================================================================================================================

#pragma once

#include "DataProcessor.h"
#include "GatherPackInfo.h"
#include "core/math/Vector3.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "field/interpolators/TrilinearInterpolatorFwd.h"

#include <functional>
#include <functional>
#include <map>
#include <vector>



namespace walberla {
namespace gather {


//**********************************************************************************************************************
/*! Packs data of a field along a curve
*
* \ingroup gather
*
*/
//**********************************************************************************************************************
template< typename GhostLayerField_T,
          typename Interpolator = field::TrilinearInterpolator< GhostLayerField_T, real_t >  >
class CurveGatherPackInfo : public GatherPackInfo
{
   public:
      typedef Vector3<real_t>  RealVec3;
      typedef std::function<RealVec3 (real_t t) > CurveCallback;


      //**Construction & Destruction*************************************************************************
      /*! \name Construction & Destruction */
      //@{


      /**
       * Construction using a callback function
       * @param bs             structured block storage, needed to map global to local coordinates
       * @param fieldID        BlockdataID pointing to a GhostLayerField of type GhostLayerField_T
       * @param curveCallback  curve defining function, taking one real_t parameter,
       *                       returning a point in "world" coordinates. Same coordinates as AABB of blockfield
       * @param tStart         Start value of free curve parameter t
       * @param tEnd           End value of free curve parameter t
       * @param numSamples     number of subdivisions of interval [tStart,tEnd]
       * @param dp             The data processor acts as output of gathered data.
       *                       When gathered data is available the process() method is called,
       *                       with an array containing data points. The first element of the data point is the
       *                       free curve parameter t, followed by the elements given by the data interpolator
       */
      CurveGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                           ConstBlockDataID fieldID,
                           CurveCallback curveCallback,
                           real_t tStart, real_t tEnd, uint_t numSamples,
                           const shared_ptr<DataProcessor> & dp );

      /**
       * Construction using curve definition as string
       *
       * see constructor above
       * @param curveX      string containing the function that defines the x value of the curve
       *                    the name of the free parameter is 't'
       *                    example: cos(t)+2*t+exp(t)
       */
      CurveGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                           ConstBlockDataID field,
                           const std::string & curveX, const std::string & curveY, const std::string & curveZ,
                           real_t tStart, real_t tEnd, uint_t numSamples,
                           const shared_ptr<DataProcessor> & dp);

      /**
       * Construction using vector of sample points
       *
       * @samplePoints  Curve definition using a vector of points in R^3.
       *                The points are expected to be in "world" coordinates
       */
      CurveGatherPackInfo( const shared_ptr<StructuredBlockStorage> & bs,
                           ConstBlockDataID field,
                           const std::vector<RealVec3 > & samplePoints,
                           const shared_ptr<DataProcessor> & dp);


      virtual ~CurveGatherPackInfo() {}

      //@}
      //****************************************************************************************************************



      //** Packing Interface  *******************************************************************************
      /*! \name Packing Interface  */
      //@{

      virtual void packData  ( const IBlock * sender,
                               mpi::SendBuffer & outBuffer );

      virtual void unpackData( mpi::RecvBuffer & buffer );


      virtual void gatherFinished();
      //@}
      //****************************************************************************************************************


   protected:

      //** Helper Functions ********************************************************************************************
      /*! \name Helper Functions */
      //@{

      /**
       * Helper function needed by all constructors
       * Computes conversion of world coordinates to local cell coordinates.
       * If sampled point @p lies in a locally allocated block, it is added to localSamplePoints_ .
       * Initializes also the totalNumberOfSamples_ member.
       */
      void addSample(real_t t, const RealVec3 & p);

      /**
       * Sorts received data according to the free curve parameter t
       * necessary because a curve may cross a block multiple times for multiple non connected
       * t intervals
       */
      void sortReceivedData();

      //@}
      //****************************************************************************************************************




      //** Members for Sending ******************************************************************************
      /*! \name  Members for Sending  */
      //@{

      shared_ptr<StructuredBlockStorage> blocks_;

      /// DataInterpolator acting as source, for the data that has to be packed
      ConstBlockDataID fieldID_;

      struct Samples {
            std::vector<real_t>    t;
            std::vector<RealVec3>  coord;
      };

      /// For every LOCAL block, where data has to be packed, a vector of
      /// sample points is stored. Points are stored in local cell coordinates
      std::map<const IBlock*, Samples > localSamplePoints_;
      //@}
      //****************************************************************************************************************



      //** Members for Receiving  ***************************************************************************
      /*! \name Members for Receiving  */
      //@{

      /// Two dimensional array of received data:
      /// the outer vector has one entry for every sample point
      /// the inner vector represents one sample point
      ///  - the first entry (receivedData[i][0] ) is the t-value
      ///  - subsequent entries are for the f values
      std::vector<std::vector<real_t> > receivedData;

      /// Channel for output of gathered data
      shared_ptr<DataProcessor> dataProcessor_;


      /// Helper class for sorting the receivedData array according to t-value
      struct Compare : public std::function<bool(std::vector<real_t> , std::vector<real_t>)>
      {
         inline bool operator()(const std::vector<real_t> & v1, const std::vector<real_t> & v2) const {
            return v1[0] < v2[0];
         }
      };
      //@}
      //****************************************************************************************************************
};


} // namespace gather
} // namespace walberla


#include "CurveGatherPackInfo.impl.h"
