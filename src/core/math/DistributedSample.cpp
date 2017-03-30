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
//! \file DistributedSample.cpp
//! \ingroup core
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//! \author Christian Godenschwager <christian.godenschwager@fau.de>
//
//======================================================================================================================

#include "DistributedSample.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Reduce.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/math/special_functions/pow.hpp>



namespace walberla {
namespace math {

/*******************************************************************************************************************//**
 * \brief   Combines the samples from all processes and stores the result on each process.
 *
 * Note that this is a collective MPI operation. It has to be called by all processes!
 **********************************************************************************************************************/
void DistributedSample::mpiAllGather()
{
   sum_ = real_t(0);
   min_ =  std::numeric_limits< real_t >::max();
   max_ = -std::numeric_limits< real_t >::max();
   size_ = uint_c( data_.size() );

   for( auto it = data_.begin(); it != data_.end(); ++it )
   {
      sum_ += *it;
      min_ = std::min( min_, *it );
      max_ = std::max( max_, *it );
   }

   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( sum_, mpi::SUM );
      mpi::allReduceInplace( min_, mpi::MIN );
      mpi::allReduceInplace( max_, mpi::MAX );
      mpi::allReduceInplace( size_, mpi::SUM );
   }

   mean_ = sum_ / real_c(size_);
   variance_ = real_t(0);

   for( auto it = data_.begin(); it != data_.end(); ++it )
      variance_ += boost::math::pow<2>( *it - mean_ );

   WALBERLA_MPI_SECTION()
   {
      mpi::allReduceInplace( variance_, mpi::SUM );
   }

   variance_ /= real_c(size_);
   
   if( size_ == uint_t(0) )
   {
      min_ = real_t(0);
      max_ = real_t(0);
      mean_ = real_t(0);
      variance_ = real_t(0);
   }
}

/*******************************************************************************************************************//**
 * \brief   Combines the samples from all processes and stores the result on process 'rank'.
 *
 * Note that this is a collective MPI operation. It has to be called by all processes!
 *
 * \param   rank   The rank of the process the combined sample is stored on.
 **********************************************************************************************************************/
void DistributedSample::mpiGather( int rank )
{
   WALBERLA_NON_MPI_SECTION()
   {
      WALBERLA_ASSERT( rank == 0 );
   }

   sum_ = real_t(0);
   min_ =  std::numeric_limits< real_t >::max();
   max_ = -std::numeric_limits< real_t >::max();
   size_ = uint_c( data_.size() );

   for( auto it = data_.begin(); it != data_.end(); ++it )
   {
      sum_ += *it;
      min_ = std::min( min_, *it );
      max_ = std::max( max_, *it );
   }

   WALBERLA_MPI_SECTION()
   {
      mpi::reduceInplace( sum_, mpi::SUM, rank );
      mpi::reduceInplace( min_, mpi::MIN, rank );
      mpi::reduceInplace( max_, mpi::MAX, rank );
      mpi::reduceInplace( size_, mpi::SUM, rank );
   }

   mean_ = sum_ / real_c(size_);
   variance_ = real_t(0);

   for( auto it = data_.begin(); it != data_.end(); ++it )
      variance_ += boost::math::pow<2>( *it - mean_ );

   WALBERLA_MPI_SECTION()
   {
      mpi::reduceInplace( variance_, mpi::SUM, rank );
   }

   variance_ /= real_c(size_);
   
   if( size_ == uint_t(0) )
   {
      min_ = real_t(0);
      max_ = real_t(0);
      mean_ = real_t(0);
      variance_ = real_t(0);
   }   
}

/*******************************************************************************************************************//**
 * \brief   Combines the samples from all processes and stores the result on the root process.
 *
 * Note that this is a collective MPI operation. It has to be called by all processes!
 **********************************************************************************************************************/
void DistributedSample::mpiGatherRoot()
{
   mpiGather( 0 );
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
 *    %var       by variance()
 *    %stddev    by stdDeviation()
 *    %relstddev by relativeStdDeviation()
 *    %size      by size()
 *
 * \returns  The formated string.
 **********************************************************************************************************************/
std::string DistributedSample::format( const std::string & formatString ) const
{
   using boost::algorithm::replace_all;
   using boost::lexical_cast;

   std::string result = formatString;

   if( size_ > uint_t(0) )
   {
      if( result.find( "%min" ) != std::string::npos )
         replace_all(result, "%min", lexical_cast<std::string>( min_ ) );
      if( result.find( "%max" ) != std::string::npos )
         replace_all(result, "%max", lexical_cast<std::string>( max_ ) );
      if( result.find( "%sum" ) != std::string::npos )
         replace_all(result, "%sum",lexical_cast<std::string>( sum_ ) );
      if( result.find( "%mean" ) != std::string::npos )
         replace_all(result, "%mean", lexical_cast<std::string>( mean_ ) );
      if( result.find( "%var" ) != std::string::npos )
         replace_all(result, "%var", lexical_cast<std::string>( variance_ ) );
      if( result.find( "%stddev" ) != std::string::npos )
         replace_all(result, "%stddev", lexical_cast<std::string>( stdDeviation() ) );
      if( result.find( "%relstddev" ) != std::string::npos )
         replace_all(result, "%relstddev", lexical_cast<std::string>( relativeStdDeviation() ) );
      if( result.find( "%size" ) != std::string::npos )
         replace_all(result, "%size", lexical_cast<std::string>( size_ ) );
   }
   else // empty()
   {
      replace_all(result, "%min",       "N/A" );
      replace_all(result, "%max",       "N/A" );
      replace_all(result, "%sum",       "N/A" );
      replace_all(result, "%mean",      "N/A" );
      replace_all(result, "%var",       "N/A" );
      replace_all(result, "%stddev",    "N/A" );
      replace_all(result, "%relstddev", "N/A" );
      replace_all(result, "%size",      "0"   );
   }

   return result;
}

const std::string DistributedSample::DEFAULT_FORMAT_STRING = "Sample has %size values in [%min, %max], "
                                                             "sum = %sum, mean = %mean, "
                                                             "stddev = %stddev (relative: %relstddev)";

} // namespace math
} // namespace walberla
