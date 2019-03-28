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
//! \file AccuracyEvaluationLinePlot.h
//! \ingroup field
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/DataTypes.h"
#include "core/Set.h"
#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIManager.h"
#include "core/mpi/Gatherv.h"
#include "core/mpi/Reduce.h"
#include "core/uid/SUID.h"

#include "domain_decomposition/StructuredBlockStorage.h"

#include "field/EvaluationFilter.h"
#include "field/iterators/IteratorMacros.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <string>
#include <type_traits>


namespace walberla {
namespace field {



namespace internal {

template< typename T >
struct AccuracyEvaluationPlotData
{
   // bug in clang 3.3, cannot use a member function to compare
   Vector3< real_t > center;
   T value;
   T solution;
};
// this struct only exists because of a bug in clang 3.3, see comment above
template< typename T >
struct AccuracyEvaluationPlotDataSorter
{
   inline bool operator()( const AccuracyEvaluationPlotData<T> & lhs,
                           const AccuracyEvaluationPlotData<T> & rhs ) const
   {
      return checkYAxis ? ( lhs.center[1] < rhs.center[1] ) : ( lhs.center[2] < rhs.center[2] );
   }
   bool checkYAxis;
};

template< typename T >
void accuracyEvaluationLinePlotIO( std::ofstream & file, const std::vector< AccuracyEvaluationPlotData<T> > & points )
{
   file << "# position [1] [2] [3], simulation [4], exact solution [5], "
           "| u - u_exact | / | u_exact | (rel. error) [10], | u - u_exact | (abs. error) [11]\n";
   for( auto point = points.begin(); point != points.end(); ++point )
   {
      T diff = std::abs( point->value - point->solution );
      file << point->center[0] << " "
           << point->center[1] << " "
           << point->center[2] << " "
           << point->value << " "
           << point->solution << " "
           << ( diff / std::abs( point->solution ) ) << " "
           << diff << "\n";
   }
}

template<>
inline void accuracyEvaluationLinePlotIO( std::ofstream & file, const std::vector< AccuracyEvaluationPlotData<Vector3<real_t>> > & points )
{
   file << "# position [1] [2] [3], simulation [4] [5] [6], exact solution [7] [8] [9], "
           "|| u - u_exact || / || u_exact || (rel. error) [10], || u - u_exact || (abs. error) [11]\n";
   for( auto point = points.begin(); point != points.end(); ++point )
   {
      Vector3< real_t > diff = point->value - point->solution;
      file << point->center[0] << " "
           << point->center[1] << " "
           << point->center[2] << " "
           << point->value[0] << " "
           << point->value[1] << " "
           << point->value[2] << " "
           << point->solution[0] << " "
           << point->solution[1] << " "
           << point->solution[2] << " "
           << ( diff.length() / point->solution.length() ) << " "
           << diff.length() << "\n";
   }
}

} // namespace internal



//**********************************************************************************************************************
/*!
*   \brief Class for plotting simulation (and solution) values along a line through the domain
*
*   \section docAccuracyEvaluationLinePlot Accuracy Evaluation Line Plot
*
*   Class for evaluating the accuracy of a simulation by comparing simulation values with values provided by a solution
*   function (can be the analytical solution, if one exists). Simulation and solution values are evaluated along a line
*   perpendicular to the x-axis. The data is saved to a file that can be used for plotting graphs.
*
*   Do not create objects of class AccuracyEvaluationLinePlot directly, better use one of the various
*   'makeAccuracyEvaluationLinePlot' functions below!
*
*   Template parameters:
*   - Field_T: the field storing the simulation values (also works if the field stores data of type Vector3)
*   - SolutionFunction_T: type of the solution function - must return Field_T::value_type and must take one parameter
*                         of type Vector3<real_t> that corresponds to a position inside the simulation domain
*                         (simulation domain = returned by calling getDomain() at the block storage)
*   - Filter_T: the type of the evaluation filter (see \ref docEvaluationFilter in 'EvaluationFilter.h')
*
*   Parameters for setting up and controlling accuracy evaluation:
*   - blocks: the block storage
*   - fieldId: block data ID of the field
*   - solution: the solution callback function - must return the solution when called with a position inside the domain
*   - filter: the evaluation filter that indicates which cells are processed
*   - y-axis: true = plot along y axis, false = plot along z axis
*   - line point: a point that specifies the position of the line - this point is always relative to the entire
*                 domain with [<0,0,0>,<1,1,1>], i.e., all coordinates of this point (x, y, and z) must be in [0;1]
*   - normalization factor: an optional factor the simulation values and the solution values are multiplied with
*   - domain normalization: By default, the data points in the file are relative to the simulation space as given by
*                           the domain bounding box (blocks->getDomain()). However, you can overwrite the size of the
*                           domain via the domain normalization parameter. If you do so, the data points in the file
*                           will correspond to this 'normalized' domain. This essentially scales the axis in the plot
*                           to a desired range.
*   - required and incompatible selectors
*
*   You do not have to specify an evaluation filter! If you do not specify any filter, _all_ cells on the line are
*   processed and no cell is excluded.
*
*   If you want to use a flag field as evaluation filter, fitting 'makeAccuracyEvaluationLinePlot' functions already
*   exist. These functions need an additional template parameter FlagField_T and you have to provide the block data ID
*   of the flag field together with a set of flag UIDs that specify which cells need to be processed.
*
*   There also exist 'makeAccuracyEvaluationLinePlot' functions that take configuration file data as an additional
*   parameter in order to parse the configuration file for setting up and controlling accuracy evaluation. The
*   configuration file block looks like as follows:
*
*   \code
*   AccuracyEvaluationLinePlot
*   {
*      y             [boolean]; // plot along y- or z-axis?
*      point         [Vector3: <real,real,real>]; // line point
*      normalization [floating point value]; // normalization factor
*      domain        [AABB: [<real,real,real>,<real,real,real>]]; // domain normalization
*   }
*   \endcode
*
*   Example:
*
*   \code
*   AccuracyEvaluationLinePlot
*   {
*      y             true;
*      point         <0.5,0.5,0.5>; // = line through the middle of the domain
*      normalization 1;
*      domain        [<-1,-1,-1>,<1,1,1>];
*   }
*   \endcode
*
*   Note that the shared pointer returned by all 'makeAccuracyEvaluationLinePlot' functions can be dereferenced and
*   called with a string (=filename) in order to generate a plot straight away or it can be passed to a
*   'makeAccuracyEvaluationLinePlotter' function (see \ref docAccuracyEvaluationLinePlotter) in order to generate a
*   series of plots.
*/
//**********************************************************************************************************************

template< typename Field_T,
          typename SolutionFunction_T = std::function< typename Field_T::value_type ( const Vector3< real_t > & ) >,
          typename Filter_T = DefaultEvaluationFilter >
class AccuracyEvaluationLinePlot
{
public:

   AccuracyEvaluationLinePlot( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                               const SolutionFunction_T & solution, const Filter_T & filter, const bool yAxis = true,
                               const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                               const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), fieldId_( fieldId ), solution_( solution ), filter_( filter ), yAxis_( yAxis ),
      relLinePoint_( Vector3<real_t>( real_c(0.5) ) ), normalizationFactor_( real_t(1) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      auto _blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'AccuracyEvaluationLinePlot' for a block storage object that doesn't exist anymore" );
      domainNormalization_ = _blocks->getDomain();
   }
   
   AccuracyEvaluationLinePlot( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                               const SolutionFunction_T & solution, const bool yAxis = true,
                               const Set<SUID> & requiredSelectors     = Set<SUID>::emptySet(),
                               const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() ) :
      blocks_( blocks ), fieldId_( fieldId ), solution_( solution ), filter_( Filter_T() ), yAxis_( yAxis ),
      relLinePoint_( Vector3<real_t>( real_c(0.5) ) ), normalizationFactor_( real_t(1) ),
      requiredSelectors_(requiredSelectors), incompatibleSelectors_( incompatibleSelectors )
   {
      static_assert( (std::is_same< Filter_T, DefaultEvaluationFilter >::value),
                     "This constructor is only available if DefaultEvaluationFilter is set as filter type!" );

      auto _blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'AccuracyEvaluationLinePlot' for a block storage object that doesn't exist anymore" );
      domainNormalization_ = _blocks->getDomain();
   }   

   void setLinePoint( const Vector3< real_t > & p ) { relLinePoint_ = p; }
   void setNormalizationFactor( const real_t f ) { normalizationFactor_ = f; }
   void setDomainNormalization( const math::AABB & d ) { domainNormalization_ = d; }
   
   void operator()( const std::string & filename );

private:
   
   weak_ptr< StructuredBlockStorage > blocks_;
   ConstBlockDataID fieldId_;
   
   SolutionFunction_T solution_;
   Filter_T filter_;
   
   bool yAxis_; // true = plot along y axis, false = plot along z axis
   Vector3< real_t > relLinePoint_; // all values are relative to the entire domain, i.e., all values must be in [0;1]
   
   real_t normalizationFactor_;
   math::AABB domainNormalization_;

   Set<SUID> requiredSelectors_;
   Set<SUID> incompatibleSelectors_;

}; // class AccuracyEvaluationLinePlot



template< typename Field_T, typename SolutionFunction_T, typename Filter_T >
void AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, Filter_T >::operator()( const std::string & filename )
{
   mpi::SendBuffer sendBuffer;

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'AccuracyEvaluationLinePlot' for a block storage object that doesn't exist anymore" );

   const auto & domainAABB = blocks->getDomain();
   Vector3<real_t> p = domainAABB.min();
   p[0] += relLinePoint_[0] * domainAABB.xSize();
   p[1] += relLinePoint_[1] * domainAABB.ySize();
   p[2] += relLinePoint_[2] * domainAABB.zSize();

   for( auto block = blocks->begin( requiredSelectors_, incompatibleSelectors_ ); block != blocks->end(); ++block )
   {
      const Field_T * field = block->template getData< const Field_T >( fieldId_ );

      filter_( *block );

#ifndef _OPENMP

      WALBERLA_FOR_ALL_CELLS_XYZ( field,

         if( filter_(x,y,z) )
         {
            const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
            const auto center = aabb.center();
            if( aabb.contains( p[0], yAxis_ ? center[1] : p[1], yAxis_ ? p[2] : center[2] ) )
            {
               sendBuffer << center << field->get(x,y,z) << solution_( center );
            }
         }
      )

#else

      // WALBERLA_FOR_ALL_CELLS macros cannot be used since they do not support additional omp pragmas inside the kernel.
      // The additional omp critical section, however, is required.

      const CellInterval size = field->xyzSize();

      if( size.zSize() >= size.ySize() )
      {
         const int izSize = int_c( size.zSize() );
         #pragma omp parallel for schedule(static)
         for( int iz = 0; iz < izSize; ++iz )
         {
            cell_idx_t z = cell_idx_c( iz );
            for( cell_idx_t y = size.yMin(); y <= size.yMax(); ++y ) {
               for( cell_idx_t x = size.xMin(); x <= size.xMax(); ++x )
               {
                  if( filter_(x,y,z) )
                  {
                     const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                     const auto center = aabb.center();
                     if( aabb.contains( p[0], yAxis_ ? center[1] : p[1], yAxis_ ? p[2] : center[2] ) )
                     {
                        #pragma omp critical (AccuracyEvaluationLinePlot)
                        {
                           sendBuffer << center << field->get(x,y,z) << solution_( center );
                        }
                     }
                  }
               }
            }
         }
      }
      else
      {
         const int iySize = int_c( size.ySize() );
         #pragma omp parallel for schedule(static)
         for( int iy = 0; iy < iySize; ++iy )
         {
            cell_idx_t y = cell_idx_c( iy );
            for( cell_idx_t z = size.zMin(); z <= size.zMax(); ++z ) {
               for( cell_idx_t x = size.xMin(); x <= size.xMax(); ++x )
               {
                  if( filter_(x,y,z) )
                  {
                     const auto aabb = blocks->getBlockLocalCellAABB( *block, Cell(x,y,z) );
                     const auto center = aabb.center();
                     if( aabb.contains( p[0], yAxis_ ? center[1] : p[1], yAxis_ ? p[2] : center[2] ) )
                     {
                        #pragma omp critical (AccuracyEvaluationLinePlot)
                        {
                           sendBuffer << center << field->get(x,y,z) << solution_( center );
                        }
                     }
                  }
               }
            }
         }
      }

#endif
   }

   mpi::RecvBuffer buffer;
   mpi::gathervBuffer( sendBuffer, buffer );

   WALBERLA_ROOT_SECTION()
   {
      typedef typename Field_T::value_type Value_T;
      
      std::vector< internal::AccuracyEvaluationPlotData<Value_T> > points;
      while( !buffer.isEmpty() )
      {
         internal::AccuracyEvaluationPlotData<Value_T> data;
         buffer >> data.center >> data.value >> data.solution;
         
         data.center[0] = domainNormalization_.xMin() + ( ( data.center[0] - domainAABB.xMin() ) / domainAABB.xSize() ) * domainNormalization_.xSize();
         data.center[1] = domainNormalization_.yMin() + ( ( data.center[1] - domainAABB.yMin() ) / domainAABB.ySize() ) * domainNormalization_.ySize();
         data.center[2] = domainNormalization_.zMin() + ( ( data.center[2] - domainAABB.zMin() ) / domainAABB.zSize() ) * domainNormalization_.zSize();
         
         data.value *= normalizationFactor_;
         data.solution *= normalizationFactor_;
 
         points.push_back( data );
      }

      internal::AccuracyEvaluationPlotDataSorter<Value_T> plotDataSorter;
      plotDataSorter.checkYAxis = yAxis_;
      std::sort( points.begin(), points.end(), plotDataSorter );

      std::ofstream file( filename.c_str() );
      file << "# accuracy evaluation of data '" << blocks->getBlockDataIdentifier( fieldId_ ) << "'\n";
      internal::accuracyEvaluationLinePlotIO<Value_T>( file, points );
      file.close();
   }
}



/////////////////////////////////////////////////////////////////////////
// makeAccuracyEvaluationLinePlot functions without configuration file //
/////////////////////////////////////////////////////////////////////////

template< typename Field_T, typename SolutionFunction_T >
shared_ptr< AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T > > makeAccuracyEvaluationLinePlot( const weak_ptr< StructuredBlockStorage > & blocks,
                                                                                                        const ConstBlockDataID & fieldId, const SolutionFunction_T & solution,
                                                                                                        const bool yAxis = true,
                                                                                                        const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                                                                                        const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   typedef AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T > AE_T;
   return shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, yAxis, requiredSelectors, incompatibleSelectors ) );
}

template< typename Field_T, typename FlagField_T, typename SolutionFunction_T >
shared_ptr< AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeAccuracyEvaluationLinePlot( const weak_ptr< StructuredBlockStorage > & blocks,
                                const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                                const SolutionFunction_T & solution,
                                const bool yAxis = true,
                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   typedef AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T> > AE_T;
   return shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                        yAxis, requiredSelectors, incompatibleSelectors ) );
}

template< typename Field_T, typename Filter_T, typename SolutionFunction_T >
shared_ptr< AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, Filter_T > >
makeAccuracyEvaluationLinePlot( const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                                const SolutionFunction_T & solution, const Filter_T & filter,
                                const bool yAxis = true,
                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   typedef AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, Filter_T > AE_T;
   return shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, filter, yAxis, requiredSelectors, incompatibleSelectors ) );
}



///////////////////////////////////////////////////////////////////
// makeAccuracyEvaluationLinePlot functions + configuration file //
///////////////////////////////////////////////////////////////////

namespace internal {

const std::string accuracyEvaluationLinePlotConfigBlock("AccuracyEvaluationLinePlot");

inline void accuracyEvaluationLinePlotConfigParser( const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                                                    bool & defaultYAxis, Vector3<real_t> & defaultRelLinePoint,
                                                    real_t & defaultNormalizationFactor, math::AABB & defaultDomainNormalization )
{
   if( parentBlockHandle )
   {
      Config::BlockHandle block = parentBlockHandle.getBlock( configBlockName );
      if( block )
      {
         defaultYAxis = block.getParameter< bool >( "y", defaultYAxis );
         defaultRelLinePoint = block.getParameter< Vector3<real_t> >( "point", defaultRelLinePoint );
         defaultNormalizationFactor = block.getParameter< real_t >( "normalization", defaultNormalizationFactor );
         defaultDomainNormalization = block.getParameter< math::AABB >( "domain", defaultDomainNormalization );
      }
   }
}

inline void accuracyEvaluationLinePlotConfigParser( const shared_ptr< Config > & config, const std::string & configBlockName,
                                                    bool & defaultYAxis, Vector3<real_t> & defaultRelLinePoint,
                                                    real_t & defaultNormalizationFactor, math::AABB & defaultDomainNormalization )
{
   if( !!config )
      accuracyEvaluationLinePlotConfigParser( config->getGlobalBlock(), configBlockName,
                                              defaultYAxis, defaultRelLinePoint, defaultNormalizationFactor, defaultDomainNormalization );
}

} // namespace internal

#define WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_CONFIG_PARSER( config ) \
   bool defaultYAxis( true ); \
   Vector3<real_t> defaultRelLinePoint( real_c(0.5) ); \
   real_t defaultNormalizationFactor( real_t(1) ); \
   auto _blocks = blocks.lock(); \
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to execute 'makeAccuracyEvaluationLinePlot' for a block storage object that doesn't exist anymore" ); \
   math::AABB defaultDomainNormalization( _blocks->getDomain() ); \
   internal::accuracyEvaluationLinePlotConfigParser( config, configBlockName, defaultYAxis, defaultRelLinePoint, defaultNormalizationFactor, defaultDomainNormalization );

#define WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_SET_AND_RETURN() \
   evaluation->setLinePoint( defaultRelLinePoint ); \
   evaluation->setNormalizationFactor( defaultNormalizationFactor ); \
   evaluation->setDomainNormalization( defaultDomainNormalization ); \
   return evaluation;

template< typename Field_T, typename SolutionFunction_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T > >
makeAccuracyEvaluationLinePlot( const Config_T & config,
                                const weak_ptr< StructuredBlockStorage > & blocks,
                                const ConstBlockDataID & fieldId, const SolutionFunction_T & solution,
                                const std::string & configBlockName = internal::accuracyEvaluationLinePlotConfigBlock,
                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_CONFIG_PARSER( config )
   typedef AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T > AE_T;
   auto evaluation = shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, defaultYAxis, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_SET_AND_RETURN()
}

template< typename Field_T, typename FlagField_T, typename SolutionFunction_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T> > >
makeAccuracyEvaluationLinePlot( const Config_T & config,
                                const weak_ptr< StructuredBlockStorage > & blocks,
                                const ConstBlockDataID & fieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                                const SolutionFunction_T & solution,
                                const std::string & configBlockName = internal::accuracyEvaluationLinePlotConfigBlock,
                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_CONFIG_PARSER( config )
   typedef AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, FlagFieldEvaluationFilter<FlagField_T> > AE_T;
   auto evaluation = shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                   defaultYAxis, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_SET_AND_RETURN()
}

template< typename Field_T, typename Filter_T, typename SolutionFunction_T, typename Config_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, Filter_T > >
makeAccuracyEvaluationLinePlot( const Config_T & config,
                                const weak_ptr< StructuredBlockStorage > & blocks, const ConstBlockDataID & fieldId,
                                const SolutionFunction_T & solution, const Filter_T & filter,
                                const std::string & configBlockName = internal::accuracyEvaluationLinePlotConfigBlock,
                                const Set<SUID> & requiredSelectors = Set<SUID>::emptySet(),
                                const Set<SUID> & incompatibleSelectors = Set<SUID>::emptySet() )
{
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_CONFIG_PARSER( config )
   typedef AccuracyEvaluationLinePlot< Field_T, SolutionFunction_T, Filter_T > AE_T;
   auto evaluation = shared_ptr< AE_T >( new AE_T( blocks, fieldId, solution, filter, defaultYAxis, requiredSelectors, incompatibleSelectors ) );
   WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_SET_AND_RETURN()
}



   


///////////////////////////////////
// AccuracyEvaluationLinePlotter //
///////////////////////////////////
   
namespace internal {

const std::string accuracyEvaluationLinePlotterFilename("linePlot");
const std::string accuracyEvaluationLinePlotterExtension("dat");
const std::string accuracyEvaluationLinePlotterConfigBlock("AccuracyEvaluationLinePlotter");

}

//**********************************************************************************************************************
/*!
*   \brief Class for plotting simulation (and solution) values along a line through the domain
*
*   \section docAccuracyEvaluationLinePlotter Accuracy Evaluation Line Plotter
*
*   Class for generating a series of plots. Build upon AccuracyEvaluationLinePlot,
*   see \ref docAccuracyEvaluationLinePlot !
*
*   Do not create objects of class AccuracyEvaluationLinePlotter directly, better use one of the
*   'makeAccuracyEvaluationLinePlotter' functions below!
*
*   Parameters for setting up and controlling the line plotter:
*   - plot: instance of AccuracyEvaluationLinePlot (see \ref docAccuracyEvaluationLinePlot)
*   - evaluation frequency: the plotting interval - if set to '0', no plotting data is created.
*   - filename: the filename base
*   - file extension: the file extension
*
*   There also exist 'makeAccuracyEvaluationLinePlotter' functions that take configuration file data as an additional
*   parameter in order to parse the configuration file for setting up and controlling the line plotter. The
*   configuration file block looks like as follows:
*
*   \code
*   AccuracyEvaluationLinePlotter
*   {
*      frequency [unsigned integer]; // the plot frequency
*      filename  [string]; // filename base
*      extension [string]; // file extension
*   }
*   \endcode
*
*   Example:
*
*   \code
*   AccuracyEvaluationLinePlotter
*   {
*      frequency 100;
*      filename  AccuracyPlot;
*      extension txt;
*   }
*   \endcode
*
*   Note that the shared pointer returned by all 'makeAccuracyEvaluationLinePlotter' functions can be captured by a
*   SharedFunctor for immediate registration at a time loop (see field::makeSharedFunctor).
*/
//**********************************************************************************************************************

template< typename AccuracyEvaluationLinePlot_T >
class AccuracyEvaluationLinePlotter
{
public:

   AccuracyEvaluationLinePlotter( const shared_ptr< AccuracyEvaluationLinePlot_T > & plot,
                                  const uint_t evaluationFrequency,
                                  const std::string & filename = internal::accuracyEvaluationLinePlotterFilename,
                                  const std::string & fileExtension = internal::accuracyEvaluationLinePlotterExtension ) :
      plot_( plot ), executionCounter_( uint_t(0) ), evaluationFrequency_( evaluationFrequency ),
      filename_( filename ), fileExtension_( fileExtension )
   {}

   void operator()()
   {
      ++executionCounter_;
      if( evaluationFrequency_ == uint_t(0) || ( executionCounter_ - uint_c(1) ) % evaluationFrequency_ != 0 )
         return;

      std::ostringstream oss;
      oss << filename_ << "_" << ( executionCounter_ - uint_c(1) ) << "." << fileExtension_;

      (*plot_)( oss.str() );
   }

private:

   shared_ptr< AccuracyEvaluationLinePlot_T > plot_;

   uint_t executionCounter_;
   uint_t evaluationFrequency_;

   std::string filename_;
   std::string fileExtension_;

}; // class AccuracyEvaluationLinePlotter



////////////////////////////////////////////////////////////////////////////
// makeAccuracyEvaluationLinePlotter functions without configuration file //
////////////////////////////////////////////////////////////////////////////

template< typename AccuracyEvaluationLinePlot_T >
shared_ptr< AccuracyEvaluationLinePlotter< AccuracyEvaluationLinePlot_T > >
makeAccuracyEvaluationLinePlotter( const shared_ptr< AccuracyEvaluationLinePlot_T > & plot,
                                   const uint_t evaluationFrequency,
                                   const std::string & filename = internal::accuracyEvaluationLinePlotterFilename,
                                   const std::string & fileExtension = internal::accuracyEvaluationLinePlotterExtension )
{
   typedef AccuracyEvaluationLinePlotter< AccuracyEvaluationLinePlot_T > AE_T;
   return shared_ptr< AE_T >( new AE_T( plot, evaluationFrequency, filename, fileExtension ) );
}

//////////////////////////////////////////////////////////////////////
// makeAccuracyEvaluationLinePlotter functions + configuration file //
//////////////////////////////////////////////////////////////////////

namespace internal {

inline void accuracyEvaluationLinePlotterConfigParser( const Config::BlockHandle & parentBlockHandle, const std::string & configBlockName,
                                                       uint_t & defaultEvaluationFrequency, std::string & defaultFilename, std::string & defaultFileExtension )
{
   if( parentBlockHandle )
   {
      Config::BlockHandle block = parentBlockHandle.getBlock( configBlockName );
      if( block )
      {
         defaultEvaluationFrequency = block.getParameter< uint_t >( "frequency", defaultEvaluationFrequency );
         defaultFilename = block.getParameter< std::string >( "filename", defaultFilename );
         defaultFileExtension = block.getParameter< std::string >( "extension", defaultFileExtension );
      }
   }
}

inline void accuracyEvaluationLinePlotterConfigParser( const shared_ptr< Config > & config, const std::string & configBlockName,
                                                       uint_t & defaultEvaluationFrequency, std::string & defaultFilename, std::string & defaultFileExtension )
{
   if( !!config )
      accuracyEvaluationLinePlotterConfigParser( config->getGlobalBlock(), configBlockName, defaultEvaluationFrequency, defaultFilename, defaultFileExtension );
}

} // namespace internal

template< typename Config_T, typename AccuracyEvaluationLinePlot_T > // Config_T may be 'shared_ptr< Config >' or 'Config::BlockHandle'
shared_ptr< AccuracyEvaluationLinePlotter< AccuracyEvaluationLinePlot_T > >
makeAccuracyEvaluationLinePlotter( const Config_T & config,
                                   const shared_ptr< AccuracyEvaluationLinePlot_T > & plot,
                                   const std::string & configBlockName = internal::accuracyEvaluationLinePlotterConfigBlock )
{
   uint_t defaultEvaluationFrequency( uint_t(0) );
   std::string defaultFilename( internal::accuracyEvaluationLinePlotterFilename );
   std::string defaultFileExtension( internal::accuracyEvaluationLinePlotterExtension );
   internal::accuracyEvaluationLinePlotterConfigParser( config, configBlockName, defaultEvaluationFrequency, defaultFilename, defaultFileExtension );
   typedef AccuracyEvaluationLinePlotter< AccuracyEvaluationLinePlot_T > AE_T;
   return shared_ptr< AE_T >( new AE_T( plot, defaultEvaluationFrequency, defaultFilename, defaultFileExtension ) );
}



#undef WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_CONFIG_PARSER
#undef WALBERLA_FIELD_MAKE_ACCURACY_EVALUATION_LINE_PLOT_SET_AND_RETURN

} // namespace field
} // namespace walberla
