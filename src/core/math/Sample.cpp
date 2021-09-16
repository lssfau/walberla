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
//! \file Sample.cpp
//! \ingroup core
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//! \brief Implementations for class Sample
//
//======================================================================================================================

#include "Sample.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/MPIManager.h"
#include "core/StringUtility.h"

#include <functional>
#include <iterator>
#include <ostream>


namespace walberla {
namespace math {

/*******************************************************************************************************************//**
 * \brief   Combines the samples from all processes and stores the result on each process.
 *
 * Note that this is a collective MPI operation. It has to be called by all processes!
 **********************************************************************************************************************/
void Sample::mpiAllGather()
{
   WALBERLA_MPI_SECTION()
   {
      std::vector< real_t > input( begin(), end() );

      std::vector< real_t > result = mpi::allGatherv( input );

      clear();
      insert( result.begin(), result.end() );
   }
}

/*******************************************************************************************************************//**
 * \brief   Combines the samples from all processes and stores the result on process 'rank'.
 *
 * Note that this is a collective MPI operation. It has to be called by all processes!
 *
 * \param   rank   The rank of the process the combined sample is stored on.
 **********************************************************************************************************************/
void Sample::mpiGather(int rank)
{
   WALBERLA_MPI_SECTION()
   {
      std::vector< real_t > input( begin(), end() );

      std::vector< real_t > result = mpi::gatherv( input, rank );

      clear();
      insert( result.begin(), result.end() );
   }
   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT( rank == 0 );
   }
}

/*******************************************************************************************************************//**
 * \brief   Combines the samples from all processes and stores the result on the root process.
 *
 * Note that this is a collective MPI operation. It has to be called by all processes!
 **********************************************************************************************************************/
void Sample::mpiGatherRoot()
{
   mpiGather( 0 );
}

/*******************************************************************************************************************//**
 * \brief   Calculates the median of the sample.
 *
 * In case of size() being an even number, the average of the two central elements is returned.
 *
 * \returns  The median.
 **********************************************************************************************************************/
real_t Sample::median() const
{
   WALBERLA_ASSERT(!empty());

   auto it = begin();
   std::advance(it, size() / 2u);

   if( size() % 2u != 0 )
      return *it;
   else
   {
      auto it2 = it;
      --it2;
      return real_c(0.5) * (*it + *it2);
   }
}

/*******************************************************************************************************************//**
 * \brief   Calculates the variance of the sample.
 *
 * The variance calculated here is the _uncorrected_ variance.
 * See: http://en.wikipedia.org/w/index.php?title=Bessel%27s_correction&oldid=526066331
 *
 * \param    theMean   the mean()
 *
 * \returns  The uncorrected variance.
 **********************************************************************************************************************/
real_t Sample::variance( real_t theMean ) const
{
   WALBERLA_ASSERT(!empty());

   KahanAccumulator< real_t > acc;
   for(auto it = begin(); it != end(); ++it)
   {
      real_t val = *it - theMean;
      acc += val*val;
   }
      
   return acc.get() / real_c(size());
}

/*******************************************************************************************************************//**
 * \brief   Calculates the relative standard deviation of the sample.
 *
 * Equals stdDeviation / mean()
 *
 * \returns  The standard deviation of the sample.
 **********************************************************************************************************************/
real_t Sample::relativeStdDeviation() const
{
   real_t theMean = mean();
   return std::sqrt( variance( theMean ) ) / theMean;
}

/*******************************************************************************************************************//**
 * \brief   Calculates the median absolute deviation (MAD) of the sample.
 *
 * MAD is a robust alternative to the standard deviation.
 *
 * See http://en.wikipedia.org/w/index.php?title=Median_absolute_deviation&oldid=608254065
 *
 * \returns  The MAD.
 **********************************************************************************************************************/
real_t Sample::mad() const
{
   WALBERLA_ASSERT( !empty() );

   real_t theMedian = median();

   std::vector<real_t> deviations( size() );

   auto valueIt     = begin();
   auto deviationIt = deviations.begin();

   while( valueIt != end() )
      *deviationIt++ = std::fabs( theMedian - *valueIt++ );

   WALBERLA_ASSERT_EQUAL( deviationIt, deviations.end() );

   if( size() % 2LU == 1LU )
   {
      auto medianPosition = deviations.begin() + size() / 2LU;
      std::nth_element( deviations.begin(), medianPosition, deviations.end() );
      return *medianPosition;
   }
   else
   {
      auto upperMedianPosition = deviations.begin() + size() / 2LU;
      std::nth_element( deviations.begin(), upperMedianPosition, deviations.end() );

      real_t upperMedian = *upperMedianPosition;
      real_t lowerMedian = *std::max_element( deviations.begin(), upperMedianPosition );
      
      return ( lowerMedian + upperMedian ) / real_t(2);
   }
}

/*******************************************************************************************************************//**
 * \brief   Calculates the Gini coefficient of the sample.
 *
 * http://en.wikipedia.org/w/index.php?title=Gini_coefficient&oldid=608263369
 *
 * \pre size() > 1
 * \pre mean() > 0
 *
 * \returns  Gini coefficient.
 **********************************************************************************************************************/
real_t Sample::giniCoefficient() const
{
   WALBERLA_ASSERT_GREATER( size(), 1LU );
   WALBERLA_ASSERT_GREATER( mean(), real_t(0) );

   real_t sum0 = 0;
   real_t sum1 = 0;
   uint_t i    = 1;

   for( auto it = begin(); it != end(); ++it )
   {
      sum0 += *it * real_t( i++ );
      sum1 += *it;
   }

   const real_t theSize = real_c( size() );
   return real_t(1) - ( real_t(2) / ( theSize - real_t(1) ) ) * ( theSize - sum0 / sum1 );
}

/*******************************************************************************************************************//**
 * \brief   Calculates a quantile of the sample.
 *
 * To understand how the quantiles are calculated please see http://tinyurl.com/d8vm37f. Quantiles
 * are rounded outwards.
 *
 * \returns  The quantile.
 **********************************************************************************************************************/
real_t Sample::quantile(const real_t p) const
{
   WALBERLA_ASSERT_GREATER_EQUAL(p, 0);
   WALBERLA_ASSERT_LESS_EQUAL(p, 1);
   WALBERLA_ASSERT(!empty());

   size_t idx;
   if(p > real_c(0.5))
      idx = numeric_cast<size_t>( std::ceil( p * real_c(size()-1) ) );
   else
      idx = numeric_cast<size_t>( std::floor( p * real_c(size()-1) ) );

   WALBERLA_ASSERT_LESS(idx, size());

   auto it = begin();
   std::advance(it, numeric_cast< std::multiset<real_t>::difference_type >(idx));
   WALBERLA_ASSERT(it != end());

   return *it;
}


/*******************************************************************************************************************//**
 * \brief   Generates a string with attributes of the sample.
 *
 * The following patters are replaced in the format string:
 *
 *    %min       by min()
 *    %max       by max()
 *    %sum       by sum()
 *    %mean      by mean()
 *    %med       by median()
 *    %var       by variance()
 *    %stddev    by stdDeviation()
 *    %relstddev by relativeStdDeviation()
 *    %mad       by mad()
 *    %size       by size()
 *
 * \returns  The formatted string.
 **********************************************************************************************************************/
std::string Sample::format(const std::string & formatString) const
{

   std::string result = formatString;

   if( !empty() )
   {
      string_replace_all( result, "%min", std::to_string( min() ) );
      string_replace_all( result, "%max", std::to_string( max() ) );
      string_replace_all( result, "%sum", std::to_string( sum() ) );
      string_replace_all( result, "%mean", std::to_string( mean() ) );
      string_replace_all( result, "%med", std::to_string( median() ) );
      string_replace_all( result, "%var", std::to_string( variance() ) );
      string_replace_all( result, "%stddev", std::to_string( stdDeviation() ) );
      string_replace_all( result, "%relstddev", std::to_string( relativeStdDeviation() ) );
      string_replace_all( result, "%mad", std::to_string( mad() ) );
      string_replace_all( result, "%size", std::to_string( size() ) );
   }
   else // empty()
   {
      string_replace_all( result, "%min",       "N/A" );
      string_replace_all( result, "%max",       "N/A" );
      string_replace_all( result, "%sum",       "N/A" );
      string_replace_all( result, "%mean",      "N/A" );
      string_replace_all( result, "%med",       "N/A" );
      string_replace_all( result, "%var",       "N/A" );
      string_replace_all( result, "%stddev",    "N/A" );
      string_replace_all( result, "%relstddev", "N/A" );
      string_replace_all( result, "%mad",       "N/A" );
      string_replace_all( result, "%size", "0" );
   }

   return result;
}

const std::string Sample::DEFAULT_FORMAT_STRING = "Sample has %size values in [%min, %max], "
                                                  "sum = %sum, mean = %mean, med = %med, "
                                                  "stddev = %stddev (relative: %relstddev), mad = %mad";

/*******************************************************************************************************************//**
 * \brief   Stream output operator for class Sample
 **********************************************************************************************************************/
std::ostream & operator<<( std::ostream & os, const Sample & statReal )
{
   os << "[" << statReal.size() << "]{";
   if(!statReal.empty())
   {
      auto it = statReal.begin();
      os << *it++;
      while(it != statReal.end())
         os << ", " << *it++;
   }
   os << "}";
   return os;
}

} // namespace math
} // namespace walberla
