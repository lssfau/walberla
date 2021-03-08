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
//! \file TimeStep.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once

#include "EqualLevelBorderStreamCorrection.h"
#include "LinearExplosion.h"
#include "PdfFieldPackInfo.h"
#include "TimeStepPdfPackInfo.h"

#include "blockforest/StructuredBlockForest.h"
#include "blockforest/communication/NonUniformBufferedScheme.h"

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/math/Uint.h"
#include "core/selectable/IsSetSelected.h"
#include "core/timing/TimingPool.h"

#include "lbm/lattice_model/NeighborsStencil.h"

#include "stencil/D2Q9.h"
#include "stencil/D3Q27.h"



namespace walberla {
namespace lbm {
namespace refinement {



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
class TimeStep
{
public:
   
   using CommunicationStencil_T = typename NeighborsStencil<LatticeModel_T>::type;

   using VoidFunction = std::function<void (const uint_t, const uint_t)>;  // parameters: level, execution count
   using BlockFunction = std::function<void (IBlock *, const uint_t, const uint_t)>; // parameters: block, level, execution count
   
private:

   class VoidFunctionWrappper
   {
   public:
      VoidFunctionWrappper( std::vector< std::pair< VoidFunction, std::string > > & functions, const uint_t index ) :
         functions_( functions ), index_( index ) {}
      void operator()( const uint_t level, const uint_t executionCount )
      {
         functions_[index_].first( level, executionCount );
      }
   private:
      std::vector< std::pair< VoidFunction, std::string > > & functions_;
      uint_t index_;
   };
   
   class BlockFunctionWrappper
   {
   public:
      BlockFunctionWrappper( std::vector< std::pair< BlockFunction, std::string > > & functions, const uint_t index ) :
         functions_( functions ), index_( index ) {}
      void operator()( IBlock * block, const uint_t level, const uint_t executionCount )
      {
         functions_[index_].first( block, level, executionCount );
      }
   private:
      std::vector< std::pair< BlockFunction, std::string > > & functions_;
      uint_t index_;
   };

   template< typename F >
   class SharedVoidFunctor
   {
   public:
      SharedVoidFunctor( const shared_ptr<F> & functorPtr ) : functorPtr_( functorPtr ) { }
      void operator()( const uint_t level, const uint_t executionCount ){ (*functorPtr_)( level, executionCount ); }
   private:
      shared_ptr<F> functorPtr_;
   };

   template< typename F >
   class SharedBlockFunctor
   {
   public:
      SharedBlockFunctor( const shared_ptr<F> & functorPtr ) : functorPtr_( functorPtr ) { }
      void operator()( IBlock * block, const uint_t level, const uint_t executionCount ){ (*functorPtr_)( block, level, executionCount ); }
   private:
      shared_ptr<F> functorPtr_;
   };

public:

   using Stencil = typename LatticeModel_T::Stencil;

   static const uint_t StreamIncludedGhostLayers = 2;

   // constructors

   TimeStep( weak_ptr< StructuredBlockForest > blocks, shared_ptr< Sweep_T > & sweep,
             const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
             const Set<SUID> & requiredBlockSelectors = Set<SUID>::emptySet(),
             const Set<SUID> & incompatibleBlockSelectors = Set<SUID>::emptySet() );

   TimeStep( weak_ptr< StructuredBlockForest > blocks, shared_ptr< Sweep_T > & sweep,
             const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
             const shared_ptr< TimeStepPdfPackInfo > & pdfPackInfo,
             const Set<SUID> & requiredBlockSelectors = Set<SUID>::emptySet(),
             const Set<SUID> & incompatibleBlockSelectors = Set<SUID>::emptySet() );
   
   // member functions
 
   bool asynchronousCommunicationIsUsed() const { return asynchronousCommunication_; }
   void asynchronousCommunication( const bool value = true ) { asynchronousCommunication_ = value; }

   bool optimizedCommunicationIsUsed() const { return optimizedCommunication_; }
   void optimizeCommunication( const bool value = true )
   {
      optimizedCommunication_ = value;
      pdfPackInfo_->optimizeEqualLevelCommunication( optimizedCommunication_ );
      pdfPackInfo_->optimizeForLinearExplosion( optimizedCommunication_ && performLinearExplosion_ );
   }

   bool equalLevelBorderStreamCorrectionIsPerformed() const { return performEqualLevelBorderStreamCorrection_; }
   void performEqualLevelBorderStreamCorrection( const bool value = true )
   {
      performEqualLevelBorderStreamCorrection_ = value;
   }
   
   bool linearExplosionIsPerformed() const { return performLinearExplosion_; }
   void performLinearExplosion( const bool value = true )
   {
      performLinearExplosion_ = value;
      pdfPackInfo_->optimizeForLinearExplosion( optimizedCommunication_ && performLinearExplosion_ );
   }

   void deactivateTiming() { timing_ = false; }
   void enableTiming( const shared_ptr<WcTimingPool> & timingPool, const shared_ptr<WcTimingPool> & levelwiseTimingPool )
   {
      timing_ = true;
      timingPool_ = timingPool;
      levelwiseTimingPool_ = levelwiseTimingPool;
      refresh( postCollideVoidFunctions_.size() );
   }

   void enableTiming()
   {
      timing_ = true;
      if( !timingPool_ )
         timingPool_ = make_shared<WcTimingPool>();
      if( !levelwiseTimingPool_ )
         levelwiseTimingPool_ = make_shared<WcTimingPool>();
      refresh( postCollideVoidFunctions_.size() );
   }

   const shared_ptr<WcTimingPool> & getTimingPool()          const { return timingPool_;          }
   const shared_ptr<WcTimingPool> & getLevelWiseTimingPool() const { return levelwiseTimingPool_; }

   void operator()()
   {
      auto blocks = blocks_.lock();
      WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );
      if( blocks->getNumberOfLevels() > postCollideVoidFunctions_.size() )
         refresh( blocks->getNumberOfLevels() );
      recursiveStep( uint_t(0), uint_t(0) );
   }

   void addPackInfo( const typename blockforest::communication::NonUniformBufferedScheme< CommunicationStencil_T >::PackInfo & packInfo )
   {
      communication_.addPackInfo( packInfo );
   }
   
   inline void addPostCollideVoidFunction ( const  VoidFunction & function, const std::string & identifier );
   inline void addPostCollideBlockFunction( const BlockFunction & function, const std::string & identifier );
   inline void addPostCollideVoidFunction ( const  VoidFunction & function, const std::string & identifier, const uint_t level );
   inline void addPostCollideBlockFunction( const BlockFunction & function, const std::string & identifier, const uint_t level );
   
   inline void addPostBoundaryHandlingVoidFunction ( const  VoidFunction & function, const std::string & identifier );
   inline void addPostBoundaryHandlingBlockFunction( const BlockFunction & function, const std::string & identifier );
   inline void addPostBoundaryHandlingVoidFunction ( const  VoidFunction & function, const std::string & identifier, const uint_t level );
   inline void addPostBoundaryHandlingBlockFunction( const BlockFunction & function, const std::string & identifier, const uint_t level );
   
   inline void addPostStreamVoidFunction ( const  VoidFunction & function, const std::string & identifier );
   inline void addPostStreamBlockFunction( const BlockFunction & function, const std::string & identifier );
   inline void addPostStreamVoidFunction ( const  VoidFunction & function, const std::string & identifier, const uint_t level );
   inline void addPostStreamBlockFunction( const BlockFunction & function, const std::string & identifier, const uint_t level );
   
   template< typename F >
   inline void addPostCollideVoidFunction ( const shared_ptr<F> & functorPtr, const std::string & identifier );
   template< typename F >
   inline void addPostCollideBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier );
   template< typename F >
   inline void addPostCollideVoidFunction ( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level );
   template< typename F >
   inline void addPostCollideBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level );

   template< typename F >
   inline void addPostBoundaryHandlingVoidFunction ( const shared_ptr<F> & functorPtr, const std::string & identifier );
   template< typename F >
   inline void addPostBoundaryHandlingBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier );
   template< typename F >
   inline void addPostBoundaryHandlingVoidFunction ( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level );
   template< typename F >
   inline void addPostBoundaryHandlingBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level );

   template< typename F >
   inline void addPostStreamVoidFunction ( const shared_ptr<F> & functorPtr, const std::string & identifier );
   template< typename F >
   inline void addPostStreamBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier );
   template< typename F >
   inline void addPostStreamVoidFunction ( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level );
   template< typename F >
   inline void addPostStreamBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level );

private:

   void init( const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId );
   void consistencyChecks( const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId ) const;
   void refresh( const uint_t levels );

   std::string getTimingPoolString( const std::string & name, const std::string & suffix = std::string("") ) const;
   std::string getLevelwiseTimingPoolString( const std::string & name, const uint_t level, const std::string & suffix = std::string("") ) const;

   void createTimers( const uint_t levels );
   
   template< typename Function >
   void addFunction( std::vector< std::vector< std::pair< Function, std::string > > > & functions,
                     const Function & function, const std::string & identifier );
   
   template< typename Function >
   void addFunction( std::vector< std::vector< std::pair< Function, std::string > > > & functions,
                     const Function & function, const std::string & identifier, const uint_t level );

   std::vector< Block * > selectedBlocks( const uint_t level ) const;

   inline void startTiming( const std::string & name, const uint_t level, const std::string & suffix = std::string("") );
   inline void  stopTiming( const std::string & name, const uint_t level, const std::string & suffix = std::string("") );

   void collide( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount );
   void stream( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount );
   void finishStream( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount );
   void streamCollide( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount );

   void startCommunicationEqualLevel( const uint_t level );
   void   endCommunicationEqualLevel( const uint_t level );
   void startCommunicationCoarseToFine( const uint_t level );
   void   endCommunicationCoarseToFine( const uint_t level );
   void startCommunicationFineToCoarse( const uint_t level );
   void   endCommunicationFineToCoarse( const uint_t level );

   void performLinearExplosion( std::vector< Block * > & blocks, const uint_t level );

   void recursiveStep( const uint_t level, const uint_t executionCount );
   
   
   
   weak_ptr< StructuredBlockForest > blocks_;

            shared_ptr< Sweep_T >          sweep_;                   // stream & collide
   typename BoundaryHandling_T::BlockSweep boundarySweep_;           // pre-stream boundary treatment
   typename BoundaryHandling_T::BlockSweep boundarySweepWithLayers_; // pre-stream boundary treatment (including ghost layers)

   bool asynchronousCommunication_;
   bool optimizedCommunication_;

   shared_ptr< TimeStepPdfPackInfo > pdfPackInfo_;

   blockforest::communication::NonUniformBufferedScheme< CommunicationStencil_T > communication_;

   bool performEqualLevelBorderStreamCorrection_;
   EqualLevelBorderStreamCorrection< LatticeModel_T > equalLevelBorderStreamCorrection_;

   bool performLinearExplosion_;
   LinearExplosion< LatticeModel_T, BoundaryHandling_T > linearExplosion_;

   bool timing_;
   shared_ptr<WcTimingPool> timingPool_;
   shared_ptr<WcTimingPool> levelwiseTimingPool_;
   
   Set<SUID> requiredBlockSelectors_;
   Set<SUID> incompatibleBlockSelectors_;
   
   std::vector< std::pair< VoidFunction,  std::string > >  globalPostCollideVoidFunctions_;
   std::vector< std::pair< BlockFunction, std::string > >  globalPostCollideBlockFunctions_;
   
   std::vector< std::pair< VoidFunction,  std::string > >  globalPostBoundaryHandlingVoidFunctions_;
   std::vector< std::pair< BlockFunction, std::string > >  globalPostBoundaryHandlingBlockFunctions_;
   
   std::vector< std::pair< VoidFunction,  std::string > >  globalPostStreamVoidFunctions_;
   std::vector< std::pair< BlockFunction, std::string > >  globalPostStreamBlockFunctions_;

   std::vector< std::vector< std::pair< VoidFunction,  std::string > > >  postCollideVoidFunctions_;
   std::vector< std::vector< std::pair< BlockFunction, std::string > > >  postCollideBlockFunctions_;

   std::vector< std::vector< std::pair< VoidFunction,  std::string > > >  postBoundaryHandlingVoidFunctions_;
   std::vector< std::vector< std::pair< BlockFunction, std::string > > >  postBoundaryHandlingBlockFunctions_;

   std::vector< std::vector< std::pair< VoidFunction,  std::string > > >  postStreamVoidFunctions_;
   std::vector< std::vector< std::pair< BlockFunction, std::string > > >  postStreamBlockFunctions_;

}; // class TimeStep

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
const uint_t TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::StreamIncludedGhostLayers;



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::TimeStep( weak_ptr< StructuredBlockForest > blocks, shared_ptr< Sweep_T > & sweep,
                                                                   const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
                                                                   const Set<SUID> & requiredBlockSelectors /*= Set<SUID>::emptySet()*/,
                                                                   const Set<SUID> & incompatibleBlockSelectors /*= Set<SUID>::emptySet()*/ ) :
   blocks_( blocks ), sweep_( sweep ),
   boundarySweep_( BoundaryHandling_T::getBlockSweep( boundaryHandlingId, uint_t(0) ) ),
   boundarySweepWithLayers_( BoundaryHandling_T::getBlockSweep( boundaryHandlingId, StreamIncludedGhostLayers ) ),
   asynchronousCommunication_( true ), optimizedCommunication_( true ),
#ifdef NDEBUG   
   pdfPackInfo_( make_shared< lbm::refinement::PdfFieldPackInfo< LatticeModel_T > >( pdfFieldId, true, true ) ),
#else
   pdfPackInfo_( make_shared< lbm::refinement::PdfFieldPackInfo< LatticeModel_T, BoundaryHandling_T > >( pdfFieldId, boundaryHandlingId, true, true ) ),
#endif
   communication_( blocks, requiredBlockSelectors, incompatibleBlockSelectors ),
   performEqualLevelBorderStreamCorrection_( true ), equalLevelBorderStreamCorrection_( pdfFieldId ),
   performLinearExplosion_( true ), linearExplosion_( pdfFieldId, boundaryHandlingId ),
   timing_( false ),
   requiredBlockSelectors_( requiredBlockSelectors ), incompatibleBlockSelectors_( incompatibleBlockSelectors )
{
   init( pdfFieldId, boundaryHandlingId );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::TimeStep( weak_ptr< StructuredBlockForest > blocks, shared_ptr< Sweep_T > & sweep,
                                                                   const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
                                                                   const shared_ptr< TimeStepPdfPackInfo > & pdfPackInfo,
                                                                   const Set<SUID> & requiredBlockSelectors /*= Set<SUID>::emptySet()*/,
                                                                   const Set<SUID> & incompatibleBlockSelectors /*= Set<SUID>::emptySet()*/ ) :
   blocks_( blocks ), sweep_( sweep ),
   boundarySweep_( BoundaryHandling_T::getBlockSweep( boundaryHandlingId, uint_t(0) ) ),
   boundarySweepWithLayers_( BoundaryHandling_T::getBlockSweep( boundaryHandlingId, StreamIncludedGhostLayers ) ),
   asynchronousCommunication_( true ), optimizedCommunication_( true ),
   pdfPackInfo_( pdfPackInfo ),
   communication_( blocks, requiredBlockSelectors, incompatibleBlockSelectors ),
   performEqualLevelBorderStreamCorrection_( true ), equalLevelBorderStreamCorrection_( pdfFieldId ),
   performLinearExplosion_( true ), linearExplosion_( pdfFieldId, boundaryHandlingId ),
   timing_( false ),
   requiredBlockSelectors_( requiredBlockSelectors ), incompatibleBlockSelectors_( incompatibleBlockSelectors )
{
   init( pdfFieldId, boundaryHandlingId );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::init( const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId )
{   
   communication_.addPackInfo( pdfPackInfo_ );
   communication_.setLocalMode( blockforest::WAIT );
   consistencyChecks( pdfFieldId, boundaryHandlingId );
   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );
   refresh( blocks->getNumberOfLevels() );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::consistencyChecks( const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId ) const
{
   // consistency checks / field validity checks

   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );
   for( auto block = blocks->begin(); block != blocks->end(); ++block )
   {
      if( !selectable::isSetSelected( block->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
         continue;

      auto *  pdfField = block->template getData< PdfField< LatticeModel_T > >( pdfFieldId );
      if( pdfField == nullptr )
      {
         WALBERLA_ABORT( "Could not get the PDF field from block " << block->getId() << ". Check if it is allocated on "
                         "the block and if the LatticeModel matches!" );
      }

      auto * boundaryHandling = block->template getData< BoundaryHandling_T >( boundaryHandlingId );
      if( boundaryHandling == nullptr )
      {
         WALBERLA_ABORT( "Could not get the boundary handling from block " << block->getId() << ". Check if it is "
                         "allocated on the block and if its type matches!" );
      }

      auto * flagField = boundaryHandling->getFlagField();
      WALBERLA_ASSERT_NOT_NULLPTR( flagField );

      if( LatticeModel_T::Stencil::D == uint_t(3) )
      {
         if( ( pdfField->xSize() & uint_t(1) ) == uint_t(1) || ( pdfField->ySize() & uint_t(1) ) == uint_t(1) ||
             ( pdfField->zSize() & uint_t(1) ) == uint_t(1) )
            WALBERLA_ABORT( "The x-, y-, and z-size of the PDF field on each block must be divisible by two!\n"
                            "(PDF field size on block " << block->getId() <<
                            " is " << pdfField->xSize() << " x " << pdfField->ySize() << " x " << pdfField->zSize() << ")" );

         if( pdfField->xSize() < uint_t(4) || pdfField->ySize() < uint_t(4) || pdfField->zSize() < uint_t(4) )
            WALBERLA_ABORT( "The size of the PDF field on each block must be at least 4x4x4 cells!\n"
                            "(PDF field size on block " << block->getId() <<
                            " is " << pdfField->xSize() << " x " << pdfField->ySize() << " x " << pdfField->zSize() << ")" );
      }
      else
      {
         WALBERLA_CHECK_EQUAL( LatticeModel_T::Stencil::D, uint_t(2) );

         if( ( pdfField->xSize() & uint_t(1) ) == uint_t(1) || ( pdfField->ySize() & uint_t(1) ) == uint_t(1) )
            WALBERLA_ABORT( "The x- and y-size of the PDF field on each block must be divisible by two!\n"
                            "(PDF field size on block " << block->getId() <<
                            " is " << pdfField->xSize() << " x " << pdfField->ySize() << " x " << pdfField->zSize() << ")" );

         if( pdfField->xSize() < uint_t(4) || pdfField->ySize() < uint_t(4) )
            WALBERLA_ABORT( "The size of the PDF field on each block must be at least 4x4x1 cells!\n"
                            "(PDF field size on block " << block->getId() <<
                            " is " << pdfField->xSize() << " x " << pdfField->ySize() << " x " << pdfField->zSize() << ")" );

         if( pdfField->zSize() != uint_t(1) )
            WALBERLA_ABORT( "The z-size of the PDF field on each block must be 1!\n"
                            "(PDF field size on block " << block->getId() <<
                            " is " << pdfField->xSize() << " x " << pdfField->ySize() << " x " << pdfField->zSize() << ")" );
      }

      if( pdfField->xSize() != flagField->xSize() || pdfField->ySize() != flagField->ySize() || pdfField->zSize() != flagField->zSize() )
         WALBERLA_ABORT( "The size of the PDF field must be identical to the size of the flag field!\n"
                         "(PDF field [" << pdfField->xSize() << " x " << pdfField->ySize() << " x " << pdfField->zSize() <<
                         "] vs. flag field [" << flagField->xSize() << " x " << flagField->ySize() << " x " << flagField->zSize() << "])" );

      if( pdfField->nrOfGhostLayers() < uint_t(4) )
         WALBERLA_ABORT( "The PDF field of each block needs at least 4 ghost layers!\n"
                         "(currently only possesses " << pdfField->nrOfGhostLayers() << " on block " << block->getId() << ")" );

      if( flagField->nrOfGhostLayers() < uint_t(4) )
         WALBERLA_ABORT( "The flag field of each block needs at least 4 ghost layers!\n"
                         "(currently only possesses " << flagField->nrOfGhostLayers() << " on block " << block->getId() << ")" );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::refresh( const uint_t levels )
{
   createTimers( levels );

   if( levels > postCollideVoidFunctions_.size() )
   {
      const uint_t previousLevels = postCollideVoidFunctions_.size();
      
      postCollideVoidFunctions_.resize( levels );
      postCollideBlockFunctions_.resize( levels );

      postBoundaryHandlingVoidFunctions_.resize( levels );
      postBoundaryHandlingBlockFunctions_.resize( levels );

      postStreamVoidFunctions_.resize( levels );
      postStreamBlockFunctions_.resize( levels );
      
      for( uint_t f = uint_t(0); f != globalPostCollideVoidFunctions_.size(); ++f )
      {
         VoidFunctionWrappper wrappedFunction( globalPostCollideVoidFunctions_, f );
         const std::string identifier = globalPostCollideVoidFunctions_[f].second;
         for( uint_t i = previousLevels; i != levels; ++i )
         {
            postCollideVoidFunctions_[i].push_back( std::make_pair( wrappedFunction, identifier ) );
            if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
               levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
         }
      }

      for( uint_t f = uint_t(0); f != globalPostCollideBlockFunctions_.size(); ++f )
      {
         BlockFunctionWrappper wrappedFunction( globalPostCollideBlockFunctions_, f );
         const std::string identifier = globalPostCollideBlockFunctions_[f].second;
         for( uint_t i = previousLevels; i != levels; ++i )
         {
            postCollideBlockFunctions_[i].push_back( std::make_pair( wrappedFunction, identifier ) );
            if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
               levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
         }
      }
      
      for( uint_t f = uint_t(0); f != globalPostBoundaryHandlingVoidFunctions_.size(); ++f )
      {
         VoidFunctionWrappper wrappedFunction( globalPostBoundaryHandlingVoidFunctions_, f );
         const std::string identifier = globalPostBoundaryHandlingVoidFunctions_[f].second;
         for( uint_t i = previousLevels; i != levels; ++i )
         {
            postBoundaryHandlingVoidFunctions_[i].push_back( std::make_pair( wrappedFunction, identifier ) );
            if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
               levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
         }
      }
      
      for( uint_t f = uint_t(0); f != globalPostBoundaryHandlingBlockFunctions_.size(); ++f )
      {
         BlockFunctionWrappper wrappedFunction( globalPostBoundaryHandlingBlockFunctions_, f );
         const std::string identifier = globalPostBoundaryHandlingBlockFunctions_[f].second;
         for( uint_t i = previousLevels; i != levels; ++i )
         {
            postBoundaryHandlingBlockFunctions_[i].push_back( std::make_pair( wrappedFunction, identifier ) );
            if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
               levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
         }
      }
      
      for( uint_t f = uint_t(0); f != globalPostStreamVoidFunctions_.size(); ++f )
      {
         VoidFunctionWrappper wrappedFunction( globalPostStreamVoidFunctions_, f );
         const std::string identifier = globalPostStreamVoidFunctions_[f].second;
         for( uint_t i = previousLevels; i != levels; ++i )
         {
            postStreamVoidFunctions_[i].push_back( std::make_pair( wrappedFunction, identifier ) );
            if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
               levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
         }
      }
      
      for( uint_t f = uint_t(0); f != globalPostStreamBlockFunctions_.size(); ++f )
      {
         BlockFunctionWrappper wrappedFunction( globalPostStreamBlockFunctions_, f );
         const std::string identifier = globalPostStreamBlockFunctions_[f].second;
         for( uint_t i = previousLevels; i != levels; ++i )
         {
            postStreamBlockFunctions_[i].push_back( std::make_pair( wrappedFunction, identifier ) );
            if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
               levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
         }
      }
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename Function >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addFunction( std::vector< std::vector< std::pair< Function, std::string > > > & functions,
                                                                           const Function & function, const std::string & identifier )
{
   for( uint_t i = uint_t(0); i < functions.size(); ++i )
   {
      functions[i].push_back( std::make_pair( function, identifier ) );
      if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, i ) ) )
         levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, i ) );
   }
   if( timing_ && ! timingPool_->timerExists( getTimingPoolString( identifier ) ) )
      timingPool_->registerTimer( getTimingPoolString( identifier ) );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename Function >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addFunction( std::vector< std::vector< std::pair< Function, std::string > > > & functions,
                                                                           const Function & function, const std::string & identifier,
                                                                           const uint_t level )
{
   if( level > ( postCollideVoidFunctions_.size() - uint_t(1) ) )
      refresh( level + uint_t(1) );
   functions[level].push_back( std::make_pair( function, identifier ) );
   if( timing_ && ! levelwiseTimingPool_->timerExists( getLevelwiseTimingPoolString( identifier, level ) ) )
      levelwiseTimingPool_->registerTimer( getLevelwiseTimingPoolString( identifier, level ) );
   if( timing_ && ! timingPool_->timerExists( getTimingPoolString( identifier ) ) )
      timingPool_->registerTimer( getTimingPoolString( identifier ) );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideVoidFunction( const VoidFunction & function, const std::string & identifier )
{
   VoidFunctionWrappper wrappedFunction( globalPostCollideVoidFunctions_, globalPostCollideVoidFunctions_.size() );
   globalPostCollideVoidFunctions_.emplace_back( function, identifier );
   addFunction< VoidFunction >( postCollideVoidFunctions_, wrappedFunction, identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideBlockFunction( const BlockFunction & function, const std::string & identifier )
{
   BlockFunctionWrappper wrappedFunction( globalPostCollideBlockFunctions_, globalPostCollideBlockFunctions_.size() );
   globalPostCollideBlockFunctions_.emplace_back( function, identifier );
   addFunction< BlockFunction >( postCollideBlockFunctions_, wrappedFunction, identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideVoidFunction( const VoidFunction & function, const std::string & identifier, const uint_t level )
{
   addFunction< VoidFunction >( postCollideVoidFunctions_, function, identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideBlockFunction( const BlockFunction & function, const std::string & identifier, const uint_t level )
{
   addFunction< BlockFunction >( postCollideBlockFunctions_, function, identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingVoidFunction( const VoidFunction & function, const std::string & identifier )
{
   VoidFunctionWrappper wrappedFunction( globalPostBoundaryHandlingVoidFunctions_, globalPostBoundaryHandlingVoidFunctions_.size() );
   globalPostBoundaryHandlingVoidFunctions_.emplace_back( function, identifier );
   addFunction< VoidFunction >( postBoundaryHandlingVoidFunctions_, wrappedFunction, identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingBlockFunction( const BlockFunction & function, const std::string & identifier )
{
   BlockFunctionWrappper wrappedFunction( globalPostBoundaryHandlingBlockFunctions_, globalPostBoundaryHandlingBlockFunctions_.size() );
   globalPostBoundaryHandlingBlockFunctions_.emplace_back( function, identifier );
   addFunction< BlockFunction >( postBoundaryHandlingBlockFunctions_, wrappedFunction, identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingVoidFunction( const VoidFunction & function, const std::string & identifier, const uint_t level )
{
   addFunction< VoidFunction >( postBoundaryHandlingVoidFunctions_, function, identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingBlockFunction( const BlockFunction & function, const std::string & identifier, const uint_t level )
{
   addFunction< BlockFunction >( postBoundaryHandlingBlockFunctions_, function, identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamVoidFunction( const VoidFunction & function, const std::string & identifier )
{
   VoidFunctionWrappper wrappedFunction( globalPostStreamVoidFunctions_, globalPostStreamVoidFunctions_.size() );
   globalPostStreamVoidFunctions_.emplace_back( function, identifier );
   addFunction< VoidFunction >( postStreamVoidFunctions_, wrappedFunction, identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamBlockFunction( const BlockFunction & function, const std::string & identifier )
{
   BlockFunctionWrappper wrappedFunction( globalPostStreamBlockFunctions_, globalPostStreamBlockFunctions_.size() );
   globalPostStreamBlockFunctions_.emplace_back( function, identifier );
   addFunction< BlockFunction >( postStreamBlockFunctions_, wrappedFunction, identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamVoidFunction( const VoidFunction & function, const std::string & identifier, const uint_t level )
{
   addFunction< VoidFunction >( postStreamVoidFunctions_, function, identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamBlockFunction( const BlockFunction & function, const std::string & identifier, const uint_t level )
{
   addFunction< BlockFunction >( postStreamBlockFunctions_, function, identifier, level );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideVoidFunction( const shared_ptr<F> & functorPtr, const std::string & identifier )
{
   addPostCollideVoidFunction( SharedVoidFunctor<F>(functorPtr), identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier )
{
   addPostCollideBlockFunction( SharedBlockFunctor<F>(functorPtr), identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideVoidFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level )
{
   addPostCollideVoidFunction( SharedVoidFunctor<F>(functorPtr), identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostCollideBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level )
{
   addPostCollideBlockFunction( SharedBlockFunctor<F>(functorPtr), identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingVoidFunction( const shared_ptr<F> & functorPtr, const std::string & identifier )
{
   addPostBoundaryHandlingVoidFunction( SharedVoidFunctor<F>(functorPtr), identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier )
{
   addPostBoundaryHandlingBlockFunction( SharedBlockFunctor<F>(functorPtr), identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingVoidFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level )
{
   addPostBoundaryHandlingVoidFunction( SharedVoidFunctor<F>(functorPtr), identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostBoundaryHandlingBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level )
{
   addPostBoundaryHandlingBlockFunction( SharedBlockFunctor<F>(functorPtr), identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamVoidFunction( const shared_ptr<F> & functorPtr, const std::string & identifier )
{
   addPostStreamVoidFunction( SharedVoidFunctor<F>(functorPtr), identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier )
{
   addPostStreamBlockFunction( SharedBlockFunctor<F>(functorPtr), identifier );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamVoidFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level )
{
   addPostStreamVoidFunction( SharedVoidFunctor<F>(functorPtr), identifier, level );
}

template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
template< typename F >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::addPostStreamBlockFunction( const shared_ptr<F> & functorPtr, const std::string & identifier, const uint_t level )
{
   addPostStreamBlockFunction( SharedBlockFunctor<F>(functorPtr), identifier, level );
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
std::string TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::getTimingPoolString( const std::string & name,
                                                                                          const std::string & suffix ) const
{
   std::ostringstream oss;
   oss  << name;
   if( !suffix.empty() )
      oss  << " " << suffix;
   return oss.str();
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
std::string TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::getLevelwiseTimingPoolString( const std::string & name, const uint_t level,
                                                                                                   const std::string & suffix ) const
{
   std::ostringstream oss;
   oss  << name << " (" << level << ")";
   if( !suffix.empty() )
      oss  << " " << suffix;
   return oss.str();
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::createTimers( const uint_t levels )
{
   // unify registered timers across all processes
   
   if( timing_ )
   {
      std::vector< std::string > timers;
      timers.push_back( getTimingPoolString( "boundary handling" ) );
      timers.push_back( getTimingPoolString( "collide" ) );
      timers.push_back( getTimingPoolString( "communication coarse to fine", "[pack & send]" ) );
      timers.push_back( getTimingPoolString( "communication coarse to fine", "[wait & unpack]" ) );
      timers.push_back( getTimingPoolString( "communication equal level", "[pack & send]" ) );
      timers.push_back( getTimingPoolString( "communication equal level", "[wait & unpack]" ) );
      timers.push_back( getTimingPoolString( "communication fine to coarse", "[pack & send]" ) );
      timers.push_back( getTimingPoolString( "communication fine to coarse", "[wait & unpack]" ) );
      timers.push_back( getTimingPoolString( "equal level border stream correction", "[prepare]" ) );
      timers.push_back( getTimingPoolString( "equal level border stream correction", "[apply]" ) );
      timers.push_back( getTimingPoolString( "linear explosion" ) );
      timers.push_back( getTimingPoolString( "stream" ) );
      timers.push_back( getTimingPoolString( "stream & collide" ) );
      
      for( auto timer = timers.begin(); timer != timers.end(); ++timer )
         if( ! timingPool_->timerExists( *timer ) )
            timingPool_->registerTimer( *timer );
            
      timers.clear();
      timers.push_back( getLevelwiseTimingPoolString( "stream & collide", levels - uint_t(1) ) );
      for( uint_t i = uint_t(0); i < levels; ++i )
      {
         timers.push_back( getLevelwiseTimingPoolString( "boundary handling", i ) );
         timers.push_back( getLevelwiseTimingPoolString( "collide", i ) );
         timers.push_back( getLevelwiseTimingPoolString( "communication equal level", i, "[pack & send]" ) );
         timers.push_back( getLevelwiseTimingPoolString( "communication equal level", i, "[wait & unpack]" ) );
         timers.push_back( getLevelwiseTimingPoolString( "linear explosion", i ) );
         timers.push_back( getLevelwiseTimingPoolString( "stream", i ) );
         if( i != uint_t(0) )
         {
            timers.push_back( getLevelwiseTimingPoolString( "communication coarse to fine", i, "[pack & send]" ) );
            timers.push_back( getLevelwiseTimingPoolString( "communication coarse to fine", i, "[wait & unpack]" ) );
            timers.push_back( getLevelwiseTimingPoolString( "communication fine to coarse", i, "[pack & send]" ) );
            timers.push_back( getLevelwiseTimingPoolString( "communication fine to coarse", i, "[wait & unpack]" ) );
         }
         if( i != ( levels - uint_t(1) ) )
         {
            timers.push_back( getLevelwiseTimingPoolString( "equal level border stream correction", i, "[prepare]" ) );
            timers.push_back( getLevelwiseTimingPoolString( "equal level border stream correction", i, "[apply]" ) );
         }
      }
      
      for( auto timer = timers.begin(); timer != timers.end(); ++timer )
         if( ! levelwiseTimingPool_->timerExists( *timer ) )
            levelwiseTimingPool_->registerTimer( *timer );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
std::vector< Block * > TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::selectedBlocks( const uint_t level ) const
{
   std::vector< Block * > levelBlocks;
   auto blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );
   blocks->getBlocks( levelBlocks, level );
   
   std::vector< Block * > sBlocks;
   
   for( auto block = levelBlocks.begin(); block != levelBlocks.end(); ++block )
      if( selectable::isSetSelected( (*block)->getState(), requiredBlockSelectors_, incompatibleBlockSelectors_ ) )
         sBlocks.push_back( *block );
   
   return sBlocks;
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::startTiming( const std::string & name, const uint_t level,
                                                                                  const std::string & suffix )
{
   (*timingPool_)[ getTimingPoolString( name, suffix ) ].start();
   (*levelwiseTimingPool_)[ getLevelwiseTimingPoolString( name, level, suffix ) ].start();
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
inline void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::stopTiming( const std::string & name, const uint_t level,
                                                                                 const std::string & suffix )
{
   (*timingPool_)[ getTimingPoolString( name, suffix ) ].end();
   (*levelwiseTimingPool_)[ getLevelwiseTimingPoolString( name, level, suffix ) ].end();
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::collide( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount )
{
   if( timing_ )
   {
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         startTiming( "collide", level );
         sweep_->collide( *block );
         stopTiming( "collide", level );

         for( auto func = postCollideBlockFunctions_[level].begin(); func != postCollideBlockFunctions_[level].end(); ++func )
         {
            startTiming( func->second, level );
            (func->first)( *block, level, executionCount );
            stopTiming( func->second, level );
         }
      }
      for( auto func = postCollideVoidFunctions_[level].begin(); func != postCollideVoidFunctions_[level].end(); ++func )
      {
         startTiming( func->second, level );
         (func->first)( level, executionCount );
         stopTiming( func->second, level );
      }
   }
   else
   {
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         sweep_->collide( *block );

         for( auto func = postCollideBlockFunctions_[level].begin(); func != postCollideBlockFunctions_[level].end(); ++func )
            (func->first)( *block, level, executionCount );
      }
      for( auto func = postCollideVoidFunctions_[level].begin(); func != postCollideVoidFunctions_[level].end(); ++func )
         (func->first)( level, executionCount );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::stream( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount )
{
   auto _blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );

   const uint_t finestLevel = _blocks->getNumberOfLevels() - uint_t(1);
   const bool doEqualLevelBorderStreamCorrection = ( performEqualLevelBorderStreamCorrection_ && level != finestLevel );

   if( timing_ )
   {
      if( postBoundaryHandlingVoidFunctions_[level].empty() )
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            bool coarseNeighbors = false;
            for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
               coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);

            if( coarseNeighbors )
            {
               startTiming( "boundary handling", level );
               boundarySweepWithLayers_( *block );
               stopTiming( "boundary handling", level );

               for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
               {
                  startTiming( func->second, level );
                  (func->first)( *block, level, executionCount );
                  stopTiming( func->second, level );
               }

               startTiming( "stream", level );
               sweep_->stream( *block, StreamIncludedGhostLayers );
               stopTiming( "stream", level );
            }
            else
            {
               startTiming( "boundary handling", level );
               boundarySweep_( *block );
               stopTiming( "boundary handling", level );

               for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
               {
                  startTiming( func->second, level );
                  (func->first)( *block, level, executionCount );
                  stopTiming( func->second, level );
               }

               startTiming( "stream", level );
               sweep_->stream( *block, uint_t(0) );
               stopTiming( "stream", level );
            }
            
            if( doEqualLevelBorderStreamCorrection )
            {
               startTiming( "equal level border stream correction", level, "[prepare]" );
               equalLevelBorderStreamCorrection_.prepare( *block );
               stopTiming( "equal level border stream correction", level, "[prepare]" );
            }
         }
      }
      else
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            bool coarseNeighbors = false;
            for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
               coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);

            if( coarseNeighbors )
            {
               startTiming( "boundary handling", level );
               boundarySweepWithLayers_( *block );
               stopTiming( "boundary handling", level );
            }
            else
            {
               startTiming( "boundary handling", level );
               boundarySweep_( *block );
               stopTiming( "boundary handling", level );
            }
            for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
            {
               startTiming( func->second, level );
               (func->first)( *block, level, executionCount );
               stopTiming( func->second, level );
            }
         }
         for( auto func = postBoundaryHandlingVoidFunctions_[level].begin(); func != postBoundaryHandlingVoidFunctions_[level].end(); ++func )
         {
            startTiming( func->second, level );
            (func->first)( level, executionCount );
            stopTiming( func->second, level );
         }
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            bool coarseNeighbors = false;
            for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
               coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);
               
            if( coarseNeighbors )
            {
               startTiming( "stream", level );
               sweep_->stream( *block, StreamIncludedGhostLayers );
               stopTiming( "stream", level );
            }
            else
            {
               startTiming( "stream", level );
               sweep_->stream( *block, uint_t(0) );
               stopTiming( "stream", level );
            }
            
            if( doEqualLevelBorderStreamCorrection )
            {
               startTiming( "equal level border stream correction", level, "[prepare]" );
               equalLevelBorderStreamCorrection_.prepare( *block );
               stopTiming( "equal level border stream correction", level, "[prepare]" );
            }
         }
      }
   }
   else
   {
      if( postBoundaryHandlingVoidFunctions_[level].empty() )
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            bool coarseNeighbors = false;
            for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
               coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);

            if( coarseNeighbors )
            {
               boundarySweepWithLayers_( *block );
               for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                  (func->first)( *block, level, executionCount );
               sweep_->stream( *block, StreamIncludedGhostLayers );
            }
            else
            {
               boundarySweep_( *block );
               for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                  (func->first)( *block, level, executionCount );
               sweep_->stream( *block, uint_t(0) );
            }
            
            if( doEqualLevelBorderStreamCorrection )
               equalLevelBorderStreamCorrection_.prepare( *block );
         }
      }
      else
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            bool coarseNeighbors = false;
            for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
               coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);

            if( coarseNeighbors )
            {
               boundarySweepWithLayers_( *block );
            }
            else
            {
               boundarySweep_( *block );
            }
            for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
               (func->first)( *block, level, executionCount );
         }
         for( auto func = postBoundaryHandlingVoidFunctions_[level].begin(); func != postBoundaryHandlingVoidFunctions_[level].end(); ++func )
            (func->first)( level, executionCount );
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            bool coarseNeighbors = false;
            for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
               coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);
                     
            if( coarseNeighbors )
            {
               sweep_->stream( *block, StreamIncludedGhostLayers );
            }
            else
            {
               sweep_->stream( *block, uint_t(0) );
            }
            
            if( doEqualLevelBorderStreamCorrection )
               equalLevelBorderStreamCorrection_.prepare( *block );            
         }
      }
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::finishStream( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount )
{
   auto _blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );

   const uint_t finestLevel = _blocks->getNumberOfLevels() - uint_t(1);
   const bool doEqualLevelBorderStreamCorrection = ( performEqualLevelBorderStreamCorrection_ && level != finestLevel );

   if( timing_ )
   {
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         if( doEqualLevelBorderStreamCorrection )
         {
            startTiming( "equal level border stream correction", level, "[apply]" );
            equalLevelBorderStreamCorrection_.apply( *block );
            stopTiming( "equal level border stream correction", level, "[apply]" );
         }
         for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
         {
            startTiming( func->second, level );
            (func->first)( *block, level, executionCount );
            stopTiming( func->second, level );
         }
      }
      for( auto func = postStreamVoidFunctions_[level].begin(); func != postStreamVoidFunctions_[level].end(); ++func )
      {
         startTiming( func->second, level );
         (func->first)( level, executionCount );
         stopTiming( func->second, level );
      }
   }
   else
   {
      for( auto block = blocks.begin(); block != blocks.end(); ++block )
      {
         if( doEqualLevelBorderStreamCorrection )
            equalLevelBorderStreamCorrection_.apply( *block );
            
         for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
            (func->first)( *block, level, executionCount );
      }
      for( auto func = postStreamVoidFunctions_[level].begin(); func != postStreamVoidFunctions_[level].end(); ++func )
         (func->first)( level, executionCount );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::streamCollide( std::vector< Block * > & blocks, const uint_t level, const uint_t executionCount )
{
#ifndef NDEBUG
   auto _blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );
   WALBERLA_ASSERT_EQUAL( level, _blocks->getNumberOfLevels() - uint_t(1) );
#endif

   if( postStreamVoidFunctions_[level].empty() )
   {
      if( timing_ )
      {
         if( postBoundaryHandlingVoidFunctions_[level].empty() )
         {
            for( auto block = blocks.begin(); block != blocks.end(); ++block )
            {
               bool coarseNeighborsOrPostStreamFunctions = !(postStreamBlockFunctions_[level].empty());
               for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighborsOrPostStreamFunctions; ++i )
                  coarseNeighborsOrPostStreamFunctions = (*block)->neighborhoodSectionHasLargerBlock(i);

               if( coarseNeighborsOrPostStreamFunctions )
               {
                  startTiming( "boundary handling", level );
                  boundarySweepWithLayers_( *block );
                  stopTiming( "boundary handling", level );

                  for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                  {
                     startTiming( func->second, level );
                     (func->first)( *block, level, executionCount );
                     stopTiming( func->second, level );
                  }

                  startTiming( "stream", level );
                  sweep_->stream( *block, StreamIncludedGhostLayers );
                  stopTiming( "stream", level );

                  for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
                  {
                     startTiming( func->second, level );
                     (func->first)( *block, level, executionCount );
                     stopTiming( func->second, level );
                  }

                  startTiming( "collide", level );
                  sweep_->collide( *block );
                  stopTiming( "collide", level );
               }
               else
               {
                  startTiming( "boundary handling", level );
                  boundarySweep_( *block );
                  stopTiming( "boundary handling", level );

                  for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                  {
                     startTiming( func->second, level );
                     (func->first)( *block, level, executionCount );
                     stopTiming( func->second, level );
                  }

                  startTiming( "stream & collide", level );
                  (*sweep_)( *block );
                  stopTiming( "stream & collide", level );
               }
               for( auto func = postCollideBlockFunctions_[level].begin(); func != postCollideBlockFunctions_[level].end(); ++func )
               {
                  startTiming( func->second, level );
                  (func->first)( *block, level, executionCount + uint_t(1) );
                  stopTiming( func->second, level );
               }
            }
            for( auto func = postCollideVoidFunctions_[level].begin(); func != postCollideVoidFunctions_[level].end(); ++func )
            {
               startTiming( func->second, level );
               (func->first)( level, executionCount + uint_t(1) );
               stopTiming( func->second, level );
            }
         }
         else
         {
            for( auto block = blocks.begin(); block != blocks.end(); ++block )
            {
               bool coarseNeighbors = false;
               for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
                  coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);

               if( coarseNeighbors )
               {
                  startTiming( "boundary handling", level );
                  boundarySweepWithLayers_( *block );
                  stopTiming( "boundary handling", level );
               }
               else
               {
                  startTiming( "boundary handling", level );
                  boundarySweep_( *block );
                  stopTiming( "boundary handling", level );
               }
               for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
               {
                  startTiming( func->second, level );
                  (func->first)( *block, level, executionCount );
                  stopTiming( func->second, level );
               }
            }
            for( auto func = postBoundaryHandlingVoidFunctions_[level].begin(); func != postBoundaryHandlingVoidFunctions_[level].end(); ++func )
            {
               startTiming( func->second, level );
               (func->first)( level, executionCount );
               stopTiming( func->second, level );
            }
            for( auto block = blocks.begin(); block != blocks.end(); ++block )
            {
               bool coarseNeighborsOrPostStreamFunctions = !(postStreamBlockFunctions_[level].empty());
               for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighborsOrPostStreamFunctions; ++i )
                  coarseNeighborsOrPostStreamFunctions = (*block)->neighborhoodSectionHasLargerBlock(i);

               if( coarseNeighborsOrPostStreamFunctions )
               {
                  startTiming( "stream", level );
                  sweep_->stream( *block, StreamIncludedGhostLayers );
                  stopTiming( "stream", level );

                  for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
                  {
                     startTiming( func->second, level );
                     (func->first)( *block, level, executionCount );
                     stopTiming( func->second, level );
                  }

                  startTiming( "collide", level );
                  sweep_->collide( *block );
                  stopTiming( "collide", level );
               }
               else
               {
                  startTiming( "stream & collide", level );
                  (*sweep_)( *block );
                  stopTiming( "stream & collide", level );
               }
               for( auto func = postCollideBlockFunctions_[level].begin(); func != postCollideBlockFunctions_[level].end(); ++func )
               {
                  startTiming( func->second, level );
                  (func->first)( *block, level, executionCount + uint_t(1) );
                  stopTiming( func->second, level );
               }
            }
            for( auto func = postCollideVoidFunctions_[level].begin(); func != postCollideVoidFunctions_[level].end(); ++func )
            {
               startTiming( func->second, level );
               (func->first)( level, executionCount + uint_t(1) );
               stopTiming( func->second, level );
            }
         }
      }
      else
      {
         if( postBoundaryHandlingVoidFunctions_[level].empty() )
         {
            for( auto block = blocks.begin(); block != blocks.end(); ++block )
            {
               bool coarseNeighborsOrPostStreamFunctions = !(postStreamBlockFunctions_[level].empty());
               for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighborsOrPostStreamFunctions; ++i )
                  coarseNeighborsOrPostStreamFunctions = (*block)->neighborhoodSectionHasLargerBlock(i);

               if( coarseNeighborsOrPostStreamFunctions )
               {
                  boundarySweepWithLayers_( *block );
                  for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                     (func->first)( *block, level, executionCount );
                  sweep_->stream( *block, StreamIncludedGhostLayers );
                  for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
                     (func->first)( *block, level, executionCount );
                  sweep_->collide( *block );
               }
               else
               {
                  boundarySweep_( *block );
                  for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                     (func->first)( *block, level, executionCount );
                  (*sweep_)( *block );
               }
               for( auto func = postCollideBlockFunctions_[level].begin(); func != postCollideBlockFunctions_[level].end(); ++func )
                  (func->first)( *block, level, executionCount + uint_t(1) );
            }
            for( auto func = postCollideVoidFunctions_[level].begin(); func != postCollideVoidFunctions_[level].end(); ++func )
               (func->first)( level, executionCount + uint_t(1) );
         }
         else
         {
            for( auto block = blocks.begin(); block != blocks.end(); ++block )
            {
               bool coarseNeighbors = false;
               for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighbors; ++i )
                  coarseNeighbors = (*block)->neighborhoodSectionHasLargerBlock(i);

               if( coarseNeighbors )
               {
                  boundarySweepWithLayers_( *block );
               }
               else
               {
                  boundarySweep_( *block );
               }
               for( auto func = postBoundaryHandlingBlockFunctions_[level].begin(); func != postBoundaryHandlingBlockFunctions_[level].end(); ++func )
                  (func->first)( *block, level, executionCount );
            }
            for( auto func = postBoundaryHandlingVoidFunctions_[level].begin(); func != postBoundaryHandlingVoidFunctions_[level].end(); ++func )
               (func->first)( level, executionCount );
            for( auto block = blocks.begin(); block != blocks.end(); ++block )
            {
               bool coarseNeighborsOrPostStreamFunctions = !(postStreamBlockFunctions_[level].empty());
               for( uint_t i = uint_t(0); i < uint_t(26) && !coarseNeighborsOrPostStreamFunctions; ++i )
                  coarseNeighborsOrPostStreamFunctions = (*block)->neighborhoodSectionHasLargerBlock(i);

               if( coarseNeighborsOrPostStreamFunctions )
               {
                  sweep_->stream( *block, StreamIncludedGhostLayers );
                  for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
                     (func->first)( *block, level, executionCount );
                  sweep_->collide( *block );
               }
               else
               {
                  (*sweep_)( *block );
               }
               for( auto func = postCollideBlockFunctions_[level].begin(); func != postCollideBlockFunctions_[level].end(); ++func )
                  (func->first)( *block, level, executionCount + uint_t(1) );
            }
            for( auto func = postCollideVoidFunctions_[level].begin(); func != postCollideVoidFunctions_[level].end(); ++func )
               (func->first)( level, executionCount + uint_t(1) );
         }
      }
   }
   else
   {
      stream( blocks, level, executionCount );

      if( timing_ )
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
            {
               startTiming( func->second, level );
               (func->first)( *block, level, executionCount );
               stopTiming( func->second, level );
            }
         }
         for( auto func = postStreamVoidFunctions_[level].begin(); func != postStreamVoidFunctions_[level].end(); ++func )
         {
            startTiming( func->second, level );
            (func->first)( level, executionCount );
            stopTiming( func->second, level );
         }
      }
      else
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            for( auto func = postStreamBlockFunctions_[level].begin(); func != postStreamBlockFunctions_[level].end(); ++func )
               (func->first)( *block, level, executionCount );
         }
         for( auto func = postStreamVoidFunctions_[level].begin(); func != postStreamVoidFunctions_[level].end(); ++func )
            (func->first)( level, executionCount );
      }

      collide( blocks, level, executionCount + uint_t(1) );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::startCommunicationEqualLevel( const uint_t level )
{
   if( timing_ )
   {
      startTiming( "communication equal level", level, "[pack & send]" );
      communication_.startCommunicateEqualLevel( level );
      stopTiming( "communication equal level", level, "[pack & send]" );
   }
   else
   {
      communication_.startCommunicateEqualLevel( level );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::endCommunicationEqualLevel( const uint_t level )
{
   if( timing_ )
   {
      startTiming( "communication equal level", level, "[wait & unpack]" );
      communication_.waitCommunicateEqualLevel( level );
      stopTiming( "communication equal level", level, "[wait & unpack]" );
   }
   else
   {
      communication_.waitCommunicateEqualLevel( level );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::startCommunicationCoarseToFine( const uint_t level )
{
   if( timing_ )
   {
      startTiming( "communication coarse to fine", level, "[pack & send]" );
      communication_.startCommunicateCoarseToFine( level );
      stopTiming( "communication coarse to fine", level, "[pack & send]" );
   }
   else
   {
      communication_.startCommunicateCoarseToFine( level );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::endCommunicationCoarseToFine( const uint_t level )
{
   if( timing_ )
   {
      startTiming( "communication coarse to fine", level, "[wait & unpack]" );
      communication_.waitCommunicateCoarseToFine( level );
      stopTiming( "communication coarse to fine", level, "[wait & unpack]" );
   }
   else
   {
      communication_.waitCommunicateCoarseToFine( level );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::startCommunicationFineToCoarse( const uint_t level )
{
   if( timing_ )
   {
      startTiming( "communication fine to coarse", level, "[pack & send]" );
      communication_.startCommunicateFineToCoarse( level );
      stopTiming( "communication fine to coarse", level, "[pack & send]" );
   }
   else
   {
      communication_.startCommunicateFineToCoarse( level );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::endCommunicationFineToCoarse( const uint_t level )
{
   if( timing_ )
   {
      startTiming( "communication fine to coarse", level, "[wait & unpack]" );
      communication_.waitCommunicateFineToCoarse( level );
      stopTiming( "communication fine to coarse", level, "[wait & unpack]" );
   }
   else
   {
      communication_.waitCommunicateFineToCoarse( level );
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::performLinearExplosion( std::vector< Block * > & blocks, const uint_t level )
{
   if( performLinearExplosion_ )
   {
      if( timing_ )
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
         {
            startTiming( "linear explosion", level );
            linearExplosion_( *block );
            stopTiming( "linear explosion", level );
         }
      }
      else
      {
         for( auto block = blocks.begin(); block != blocks.end(); ++block )
            linearExplosion_( *block );
      }
   }
}



template< typename LatticeModel_T, typename Sweep_T, typename BoundaryHandling_T >
void TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T >::recursiveStep( const uint_t level, const uint_t executionCount )
{
   // -> "A generic, mass conservative local grid refinement technique for lattice-Boltzmann schemes", Rhode et al., 2005

   auto _blocks = blocks_.lock();
   WALBERLA_CHECK_NOT_NULLPTR( _blocks, "Trying to access 'TimeStep' (refinement) for a block storage object that doesn't exist anymore" );

   const uint_t coarsestLevel = uint_t(0);
   const uint_t   finestLevel = _blocks->getNumberOfLevels() - uint_t(1);

   const uint_t executionCount1st = (executionCount + uint_t(1)) * uint_t(2) - uint_t(2);
   const uint_t executionCount2nd = (executionCount + uint_t(1)) * uint_t(2) - uint_t(1);

   std::vector< Block * > blocks = selectedBlocks( level );

   WALBERLA_LOG_DETAIL("Starting recursive step with level " << level << " and execution count " << executionCount);

   if( asynchronousCommunication_ && level != coarsestLevel )
   {
      WALBERLA_LOG_DETAIL("Start communication coarse to fine, initiated by fine level " << level );
      startCommunicationCoarseToFine( level ); // [start] explosion (initiated by fine level, involves "level" and "level-1")
   }

   WALBERLA_LOG_DETAIL("Colliding on level " << level);
   collide( blocks, level, executionCount1st );

   if( asynchronousCommunication_ )
   {
      WALBERLA_LOG_DETAIL("Start communication equal level, initiated by level " << level );
      startCommunicationEqualLevel( level ); // [start] equal level communication
   }

   if( level != finestLevel )
   {
      WALBERLA_LOG_DETAIL("Calling recursive step with level " << level + uint_t(1) << " and execution count " << executionCount1st );
      recursiveStep( level + uint_t(1), executionCount1st );

      if( asynchronousCommunication_ ) {
         WALBERLA_LOG_DETAIL("Start communication fine to coarse, initiated by coarse level " << level );
         startCommunicationFineToCoarse(level + uint_t(1)); // [start] coalescence (initiated by coarse level)
      }
   }

   if( level != coarsestLevel )
   {
      if( !asynchronousCommunication_ ) {
         WALBERLA_LOG_DETAIL("Start communication coarse to fine, initiated by fine level " << level );
         startCommunicationCoarseToFine(level); // [start] explosion (initiated by fine level, involves "level" and "level-1")
      }
      WALBERLA_LOG_DETAIL("End communication coarse to fine, initiated by fine level " << level );
      endCommunicationCoarseToFine( level ); // [end] explosion (initiated by fine level, involves "level" and "level-1")
      WALBERLA_LOG_DETAIL("Perform linear explosion on level " << level );
      performLinearExplosion( blocks, level );
   }

   if( !asynchronousCommunication_ ) {
      WALBERLA_LOG_DETAIL("Start communication equal level, initiated by level " << level );
      startCommunicationEqualLevel(level); // [start] equal level communication
   }
   WALBERLA_LOG_DETAIL("End communication equal level, initiated by level " << level );
   endCommunicationEqualLevel( level ); // [end] equal level communication

   // performLinearExplosion( blocks, level ); // if equal level neighbor values are needed, linear explosion should be performed here

   if( level == finestLevel && level != coarsestLevel )
   {
      WALBERLA_LOG_DETAIL("Stream + collide on level " << level );
      streamCollide( blocks, level, executionCount1st );
   }
   else
   {
      WALBERLA_LOG_DETAIL("Stream on level " << level );
      stream( blocks, level, executionCount1st );

      if( level != finestLevel )
      {
         if( !asynchronousCommunication_ ) {
            WALBERLA_LOG_DETAIL("Start communication fine to coarse, initiated by coarse level " << level );
            startCommunicationFineToCoarse(level + uint_t(1)); // [start] coalescence (initiated by coarse level)
         }
         WALBERLA_LOG_DETAIL("End communication fine to coarse, initiated by coarse level " << level );
         endCommunicationFineToCoarse( level + uint_t(1) ); // [end] coalescence (initiated by coarse level)
      }

      WALBERLA_LOG_DETAIL("Finish stream on level " << level );
      finishStream( blocks, level, executionCount1st );

      if( level == coarsestLevel )
      {
         WALBERLA_LOG_DETAIL("End recursive step on level " << level );
         return;
      }

      WALBERLA_LOG_DETAIL("Colliding on level " << level);
      collide( blocks, level, executionCount2nd );
   }

   if( asynchronousCommunication_ ) {
      WALBERLA_LOG_DETAIL("Start communication equal level, initiated by level " << level );
      startCommunicationEqualLevel(level); // [start] equal level communication
   }

   if( level != finestLevel )
   {
      WALBERLA_LOG_DETAIL("Calling recursive step with level " << level + uint_t(1) << " and execution count " << executionCount2nd );
      recursiveStep( level + uint_t(1), executionCount2nd );
   
      if( asynchronousCommunication_ ) {
         WALBERLA_LOG_DETAIL("Start communication fine to coarse, initiated by coarse level " << level );
         startCommunicationFineToCoarse(level + uint_t(1)); // [start] coalescence (initiated by coarse level)
      }
   }

   if( !asynchronousCommunication_ ) {
      WALBERLA_LOG_DETAIL("Start communication equal level, initiated by level " << level );
      startCommunicationEqualLevel(level); // [start] equal level communication
   }
   WALBERLA_LOG_DETAIL("End communication equal level, initiated by level " << level );
   endCommunicationEqualLevel( level ); // [end] equal level communication

   WALBERLA_LOG_DETAIL("Stream on level " << level );
   stream( blocks, level, executionCount2nd );

   if( level != finestLevel )
   {
      if( !asynchronousCommunication_ )
      {
         WALBERLA_LOG_DETAIL("Start communication fine to coarse, initiated by coarse level " << level );
         startCommunicationFineToCoarse( level + uint_t(1) ); // [start] coalescence (initiated by coarse level)
      }
      WALBERLA_LOG_DETAIL("End communication fine to coarse, initiated by coarse level " << level );
      endCommunicationFineToCoarse( level + uint_t(1) ); // [end] coalescence (initiated by coarse level)
   }

   WALBERLA_LOG_DETAIL("Finish stream on level " << level );
   finishStream( blocks, level, executionCount2nd );

   WALBERLA_LOG_DETAIL("End recursive step on level " << level );
}




template< typename LatticeModel_T, typename BoundaryHandling_T, typename Sweep_T >
shared_ptr< TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T > >
makeTimeStep( weak_ptr< StructuredBlockForest > blocks, shared_ptr< Sweep_T > & sweep,
              const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
              const Set<SUID> & requiredBlockSelectors = Set<SUID>::emptySet(),
              const Set<SUID> & incompatibleBlockSelectors = Set<SUID>::emptySet() )
{
   using TS_T = TimeStep<LatticeModel_T, Sweep_T, BoundaryHandling_T>;
   return shared_ptr< TS_T >( new TS_T( blocks, sweep, pdfFieldId, boundaryHandlingId, requiredBlockSelectors, incompatibleBlockSelectors ) );
}



template< typename LatticeModel_T, typename BoundaryHandling_T, typename Sweep_T >
shared_ptr< TimeStep< LatticeModel_T, Sweep_T, BoundaryHandling_T > >
makeTimeStep( weak_ptr< StructuredBlockForest > blocks, shared_ptr< Sweep_T > & sweep,
              const BlockDataID & pdfFieldId, const BlockDataID & boundaryHandlingId,
              const shared_ptr< TimeStepPdfPackInfo > & pdfPackInfo,
              const Set<SUID> & requiredBlockSelectors = Set<SUID>::emptySet(),
              const Set<SUID> & incompatibleBlockSelectors = Set<SUID>::emptySet() )
{
   using TS_T = TimeStep<LatticeModel_T, Sweep_T, BoundaryHandling_T>;
   return shared_ptr< TS_T >( new TS_T( blocks, sweep, pdfFieldId, boundaryHandlingId, pdfPackInfo, requiredBlockSelectors, incompatibleBlockSelectors ) );
}



} // namespace refinement
} // namespace lbm
} // namespace walberla
