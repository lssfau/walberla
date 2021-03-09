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
//! \file CellwiseSweep.h
//! \ingroup lbm
//! \author Florian Schornbaum <florian.schornbaum@fau.de>
//
//======================================================================================================================

#pragma once



namespace walberla {
namespace lbm {



//**********************************************************************************************************************
/*!
*   \brief Class for performing the streaming and collision step of the LBM
*
*   \section docCellwiseSweep Cellwise Sweep for the LBM
*
*   Do not create objects of class CellwiseSweep directly, better use one of the various 'makeCellwiseSweep'
*   functions below! Advantage of using these make functions: typically, only the type of the lattice model must be
*   specified, all other templates can be deduced automatically. If you need a sweep for the LBM that performs
*   advection diffusion, look at the 'makeCellwiseAdvectionDiffusionSweep' functions.
*
*   Template parameters:
*   - LatticeModel_T: type of the LB lattice model
*   - Filter_T: the type of the evaluation filter (see \ref docEvaluationFilter in 'EvaluationFilter.h')
*   - DensityVelocityIn_T: type of the density/velocity calculation functor that specifies how the density and
*                          equilibrium velocity are calculated after the streaming step/prior to the collision step
*                          (see \ref docDensityVelocityCalculation in 'DensityVelocityCallback.h')
*   - DensityVelocityOut_T: type of the density/velocity callback functor - can be used to store the density and/or
*                           velocity (both of which are calculated prior to the collision step) to a separate field.
*                           (see \ref docDensityVelocityCallback in 'DensityVelocityCallback.h')
*
*   Parameters passed to the 'makeCellwiseSweep' functions:
*   - pdfFieldId: the block data ID of the PDF field
*                 (the temporary field needed during the streaming step is managed automatically)
*   - src + dst: block data IDs to a source and a destination PDF field
*                (for manually controlling the temporary field that is used during the streaming step)
*   - filter: the evaluation filter that indicates which cells are processed
*   - densityVelocityIn: the density/velocity calculation functor that specifies how the density and
*                        equilibrium velocity are calculated after the streaming step/prior to the collision step
*   - densityVelocityOut: the density/velocity callback functor that can be used to store the density and/or
*                         velocity (both of which are calculated prior to the collision step) to a separate field.
*
*   The 'makeCellwiseAdvectionDiffusionSweep' functions need a block data ID for and the type of a velocity field
*   that is used to fetch the velocity that is required for the equilibrium distribution in the collision step.
*
*   You do not have to specify an evaluation filter or a density/velocity callback function! If you do not specify
*   any filter and if you do not provide a callback function, _all_ cells are processed (no cell is excluded) and
*   the density and velocity are not stored in a separate field.
*
*   If you want to use a flag field as evaluation filter, fitting 'makeCellwiseSweep' functions already exist.
*   These functions need an additional template parameter FlagField_T and you have to provide the block data ID of the
*   flag field together with a set of flag UIDs that specify which cells need to be processed.
*   Similarly, if you want to store the density and/or velocity in a separate field, fitting 'makeCellwiseSweep'
*   functions exist. They need the type of the velocity field and/or the type of the density field as additional
*   template parameters, and you must provide block data IDs that correspond to these fields.
*
*   Note that the shared pointer returned by all 'makeCellwiseSweep' functions can be captured by a SharedSweep
*   for immediate registration at a time loop (see domain_decomposition::makeSharedSweep).
*/
//**********************************************************************************************************************



template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T, class Enable = void >
class CellwiseSweep
{
  static_assert( never_true<LatticeModel_T>::value, "Instantiating 'lbm::CellwiseSweep' failed, possible reasons:\n"
                                                     " - For your current LB lattice model, there is yet no implementation for class 'lbm::CellwiseSweep'.\n"
                                                     " - Your custom (?) collision model and/or force model does not conform to the concept of a collision/force model.\n"
                                                     " - You are providing a wrong number of template parameters to class 'lbm::CellwiseSweep'.\n"
                                                     "   'lbm::CellwiseSweep' needs 4 template parameters:\n"
                                                     "   1. The LB lattice model\n"
                                                     "   2. The type of the cell evaluation filter\n"
                                                     "   3. The type of the density/velocity calculation (in) function\n"
                                                     "   4. The type of the density/velocity callback (out) function\n" );
   
};



///////////////////////////////////
// HELPER MACROS FOR COMMON CODE //
///////////////////////////////////

#define WALBERLA_LBM_CELLWISE_SWEEP_CLASS_HEAD_AND_STREAM( specialization ) \
   template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T > \
   class CellwiseSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, typename std::enable_if< specialization >::type > \
      : public SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T > \
   { \
   public: \
   \
      typedef typename SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >::PdfField_T  PdfField_T; \
      typedef typename LatticeModel_T::Stencil  Stencil_T; \
      \
      CellwiseSweep( const BlockDataID & pdfFieldId, \
                     const Filter_T & _filter = walberla::field::DefaultEvaluationFilter(), \
                     const DensityVelocityIn_T & _densityVelocityIn = DefaultDensityEquilibriumVelocityCalculation(), \
                     const DensityVelocityOut_T & _densityVelocityOut = DefaultDensityVelocityCallback() ) : \
         SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >( pdfFieldId, _filter, _densityVelocityIn, _densityVelocityOut ) {} \
      \
      CellwiseSweep( const BlockDataID & src, const BlockDataID & dst, \
                     const Filter_T & _filter = walberla::field::DefaultEvaluationFilter(), \
                     const DensityVelocityIn_T & _densityVelocityIn = DefaultDensityEquilibriumVelocityCalculation(), \
                     const DensityVelocityOut_T & _densityVelocityOut = DefaultDensityVelocityCallback() ) : \
         SweepBase< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T >( src, dst, _filter, _densityVelocityIn, _densityVelocityOut ) {} \
      \
      void operator()( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) ) \
      { \
         streamCollide( block, numberOfGhostLayersToInclude ); \
      } \
      \
      void streamCollide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) ); \
      \
      void stream ( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) ); \
      void collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude = uint_t(0) ); \
   }; \
   \
   template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T > \
   void CellwiseSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, typename std::enable_if< specialization >::type \
      >::stream( IBlock * const block, const uint_t numberOfGhostLayersToInclude ) \
   { \
      PdfField_T * src( NULL ); \
      PdfField_T * dst( NULL ); \
      this->getFields( block, src, dst ); \
      StreamPull< LatticeModel_T >::execute( src, dst, block, this->filter_, numberOfGhostLayersToInclude ); \
   }



#define WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_HEAD( specialization) \
   template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T > \
   void CellwiseSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, typename std::enable_if< specialization >::type \
      >::streamCollide( IBlock * const block, const uint_t numberOfGhostLayersToInclude ) \
   { \
      PdfField_T * src( NULL ); \
      PdfField_T * dst( NULL ); \
      \
      this->getFields( block, src, dst ); \
      \
      WALBERLA_ASSERT_GREATER( src->nrOfGhostLayers(), numberOfGhostLayersToInclude ); \
      WALBERLA_ASSERT_GREATER_EQUAL( dst->nrOfGhostLayers(), numberOfGhostLayersToInclude ); \
      \
      const auto & lm = src->latticeModel(); \
      dst->resetLatticeModel( lm ); /* required so that member functions for getting density and equilibrium velocity can be called for dst! */ \
      \
      this->filter( *block ); \
      this->densityVelocityIn( *block ); \
      this->densityVelocityOut( *block );

#define WALBERLA_LBM_CELLWISE_SWEEP_STREAM_COLLIDE_FOOT() src->swapDataPointers( dst ); }



#define WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_HEAD( specialization ) \
   template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T > \
   void CellwiseSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T, typename std::enable_if< specialization >::type \
      >::collide( IBlock * const block, const uint_t numberOfGhostLayersToInclude ) \
   { \
      PdfField_T * src = this->getSrcField( block ); \
      WALBERLA_ASSERT_GREATER_EQUAL( src->nrOfGhostLayers(), numberOfGhostLayersToInclude ); \
      \
      const auto & lm = src->latticeModel(); \
      \
      this->filter( *block ); \
      this->densityVelocityIn( *block ); \
      this->densityVelocityOut( *block );

#define WALBERLA_LBM_CELLWISE_SWEEP_COLLIDE_FOOT() }



#define WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_STREAM_COLLIDE_PULL() \
         const real_t vC  = src->get( x  , y  , z, Stencil_T::idx[C]  ); \
         const real_t vN  = src->get( x  , y-1, z, Stencil_T::idx[N]  ); \
         const real_t vS  = src->get( x  , y+1, z, Stencil_T::idx[S]  ); \
         const real_t vW  = src->get( x+1, y  , z, Stencil_T::idx[W]  ); \
         const real_t vE  = src->get( x-1, y  , z, Stencil_T::idx[E]  ); \
         const real_t vNW = src->get( x+1, y-1, z, Stencil_T::idx[NW] ); \
         const real_t vNE = src->get( x-1, y-1, z, Stencil_T::idx[NE] ); \
         const real_t vSW = src->get( x+1, y+1, z, Stencil_T::idx[SW] ); \
         const real_t vSE = src->get( x-1, y+1, z, Stencil_T::idx[SE] );

#define WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_COLLIDE_GET() \
         const real_t vC  = src->get( x, y, z, Stencil_T::idx[C]  ); \
         const real_t vN  = src->get( x, y, z, Stencil_T::idx[N]  ); \
         const real_t vS  = src->get( x, y, z, Stencil_T::idx[S]  ); \
         const real_t vW  = src->get( x, y, z, Stencil_T::idx[W]  ); \
         const real_t vE  = src->get( x, y, z, Stencil_T::idx[E]  ); \
         const real_t vNW = src->get( x, y, z, Stencil_T::idx[NW] ); \
         const real_t vNE = src->get( x, y, z, Stencil_T::idx[NE] ); \
         const real_t vSW = src->get( x, y, z, Stencil_T::idx[SW] ); \
         const real_t vSE = src->get( x, y, z, Stencil_T::idx[SE] );

#define WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_DENSITY_VELOCITY_INCOMP() \
         const real_t velXTerm = vE + vNE + vSE; \
         const real_t velYTerm = vN + vNW; \
         \
         const real_t rho = vC + vS + vW + vSW + velXTerm + velYTerm; \
         \
         const real_t velX = velXTerm - vW  - vNW - vSW; \
         const real_t velY = velYTerm + vNE - vS  - vSW - vSE;

#define WALBERLA_LBM_CELLWISE_SWEEP_D2Q9_DENSITY_VELOCITY_COMP() \
         const real_t velXTerm = vE + vNE + vSE; \
         const real_t velYTerm = vN + vNW; \
         \
         const real_t rho = vC + vS + vW + vB + vSW + velXTerm + velYTerm; \
         const real_t invRho = real_t(1.0) / rho; \
         \
         const real_t velX = invRho * ( velXTerm - vW  - vNW - vSW ); \
         const real_t velY = invRho * ( velYTerm + vNE - vS  - vSW - vSE );



#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_STREAM_COLLIDE_PULL() \
         const real_t vC  = src->get( x  , y  , z  , Stencil_T::idx[C]  ); \
         const real_t vN  = src->get( x  , y-1, z  , Stencil_T::idx[N]  ); \
         const real_t vS  = src->get( x  , y+1, z  , Stencil_T::idx[S]  ); \
         const real_t vW  = src->get( x+1, y  , z  , Stencil_T::idx[W]  ); \
         const real_t vE  = src->get( x-1, y  , z  , Stencil_T::idx[E]  ); \
         const real_t vT  = src->get( x  , y  , z-1, Stencil_T::idx[T]  ); \
         const real_t vB  = src->get( x  , y  , z+1, Stencil_T::idx[B]  ); \
         const real_t vNW = src->get( x+1, y-1, z  , Stencil_T::idx[NW] ); \
         const real_t vNE = src->get( x-1, y-1, z  , Stencil_T::idx[NE] ); \
         const real_t vSW = src->get( x+1, y+1, z  , Stencil_T::idx[SW] ); \
         const real_t vSE = src->get( x-1, y+1, z  , Stencil_T::idx[SE] ); \
         const real_t vTN = src->get( x  , y-1, z-1, Stencil_T::idx[TN] ); \
         const real_t vTS = src->get( x  , y+1, z-1, Stencil_T::idx[TS] ); \
         const real_t vTW = src->get( x+1, y  , z-1, Stencil_T::idx[TW] ); \
         const real_t vTE = src->get( x-1, y  , z-1, Stencil_T::idx[TE] ); \
         const real_t vBN = src->get( x  , y-1, z+1, Stencil_T::idx[BN] ); \
         const real_t vBS = src->get( x  , y+1, z+1, Stencil_T::idx[BS] ); \
         const real_t vBW = src->get( x+1, y  , z+1, Stencil_T::idx[BW] ); \
         const real_t vBE = src->get( x-1, y  , z+1, Stencil_T::idx[BE] );

#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_COLLIDE_GET() \
         const real_t vC  = src->get( x, y, z, Stencil_T::idx[C]  ); \
         const real_t vN  = src->get( x, y, z, Stencil_T::idx[N]  ); \
         const real_t vS  = src->get( x, y, z, Stencil_T::idx[S]  ); \
         const real_t vW  = src->get( x, y, z, Stencil_T::idx[W]  ); \
         const real_t vE  = src->get( x, y, z, Stencil_T::idx[E]  ); \
         const real_t vT  = src->get( x, y, z, Stencil_T::idx[T]  ); \
         const real_t vB  = src->get( x, y, z, Stencil_T::idx[B]  ); \
         const real_t vNW = src->get( x, y, z, Stencil_T::idx[NW] ); \
         const real_t vNE = src->get( x, y, z, Stencil_T::idx[NE] ); \
         const real_t vSW = src->get( x, y, z, Stencil_T::idx[SW] ); \
         const real_t vSE = src->get( x, y, z, Stencil_T::idx[SE] ); \
         const real_t vTN = src->get( x, y, z, Stencil_T::idx[TN] ); \
         const real_t vTS = src->get( x, y, z, Stencil_T::idx[TS] ); \
         const real_t vTW = src->get( x, y, z, Stencil_T::idx[TW] ); \
         const real_t vTE = src->get( x, y, z, Stencil_T::idx[TE] ); \
         const real_t vBN = src->get( x, y, z, Stencil_T::idx[BN] ); \
         const real_t vBS = src->get( x, y, z, Stencil_T::idx[BS] ); \
         const real_t vBW = src->get( x, y, z, Stencil_T::idx[BW] ); \
         const real_t vBE = src->get( x, y, z, Stencil_T::idx[BE] );

#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_INCOMP() \
         const real_t velXTerm = vE + vNE + vSE + vTE + vBE; \
         const real_t velYTerm = vN + vNW + vTN + vBN; \
         const real_t velZTerm = vT + vTS + vTW; \
         \
         const real_t rho = vC + vS + vW + vB + vSW + vBS + vBW + velXTerm + velYTerm + velZTerm; \
         \
         const real_t velX = velXTerm - vW  - vNW - vSW - vTW - vBW; \
         const real_t velY = velYTerm + vNE - vS  - vSW - vSE - vTS - vBS; \
         const real_t velZ = velZTerm + vTN + vTE - vB  - vBN - vBS - vBW - vBE;

#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q19_DENSITY_VELOCITY_COMP() \
         const real_t velXTerm = vE + vNE + vSE + vTE + vBE; \
         const real_t velYTerm = vN + vNW + vTN + vBN; \
         const real_t velZTerm = vT + vTS + vTW; \
         \
         const real_t rho = vC + vS + vW + vB + vSW + vBS + vBW + velXTerm + velYTerm + velZTerm; \
         const real_t invRho = real_t(1.0) / rho; \
         \
         const real_t velX = invRho * ( velXTerm - vW  - vNW - vSW - vTW - vBW ); \
         const real_t velY = invRho * ( velYTerm + vNE - vS  - vSW - vSE - vTS - vBS ); \
         const real_t velZ = invRho * ( velZTerm + vTN + vTE - vB  - vBN - vBS - vBW - vBE );



#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_STREAM_COLLIDE_PULL() \
         const real_t vC   = src->get( x  , y  , z  , Stencil_T::idx[C]   ); \
         const real_t vN   = src->get( x  , y-1, z  , Stencil_T::idx[N]   ); \
         const real_t vS   = src->get( x  , y+1, z  , Stencil_T::idx[S]   ); \
         const real_t vW   = src->get( x+1, y  , z  , Stencil_T::idx[W]   ); \
         const real_t vE   = src->get( x-1, y  , z  , Stencil_T::idx[E]   ); \
         const real_t vT   = src->get( x  , y  , z-1, Stencil_T::idx[T]   ); \
         const real_t vB   = src->get( x  , y  , z+1, Stencil_T::idx[B]   ); \
         const real_t vNW  = src->get( x+1, y-1, z  , Stencil_T::idx[NW]  ); \
         const real_t vNE  = src->get( x-1, y-1, z  , Stencil_T::idx[NE]  ); \
         const real_t vSW  = src->get( x+1, y+1, z  , Stencil_T::idx[SW]  ); \
         const real_t vSE  = src->get( x-1, y+1, z  , Stencil_T::idx[SE]  ); \
         const real_t vTN  = src->get( x  , y-1, z-1, Stencil_T::idx[TN]  ); \
         const real_t vTS  = src->get( x  , y+1, z-1, Stencil_T::idx[TS]  ); \
         const real_t vTW  = src->get( x+1, y  , z-1, Stencil_T::idx[TW]  ); \
         const real_t vTE  = src->get( x-1, y  , z-1, Stencil_T::idx[TE]  ); \
         const real_t vBN  = src->get( x  , y-1, z+1, Stencil_T::idx[BN]  ); \
         const real_t vBS  = src->get( x  , y+1, z+1, Stencil_T::idx[BS]  ); \
         const real_t vBW  = src->get( x+1, y  , z+1, Stencil_T::idx[BW]  ); \
         const real_t vBE  = src->get( x-1, y  , z+1, Stencil_T::idx[BE]  ); \
         const real_t vTNE = src->get( x-1, y-1, z-1, Stencil_T::idx[TNE] ); \
         const real_t vTNW = src->get( x+1, y-1, z-1, Stencil_T::idx[TNW] ); \
         const real_t vTSE = src->get( x-1, y+1, z-1, Stencil_T::idx[TSE] ); \
         const real_t vTSW = src->get( x+1, y+1, z-1, Stencil_T::idx[TSW] ); \
         const real_t vBNE = src->get( x-1, y-1, z+1, Stencil_T::idx[BNE] ); \
         const real_t vBNW = src->get( x+1, y-1, z+1, Stencil_T::idx[BNW] ); \
         const real_t vBSE = src->get( x-1, y+1, z+1, Stencil_T::idx[BSE] ); \
         const real_t vBSW = src->get( x+1, y+1, z+1, Stencil_T::idx[BSW] );

#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_COLLIDE_GET() \
         const real_t vC   = src->get( x, y, z, Stencil_T::idx[C]   ); \
         const real_t vN   = src->get( x, y, z, Stencil_T::idx[N]   ); \
         const real_t vS   = src->get( x, y, z, Stencil_T::idx[S]   ); \
         const real_t vW   = src->get( x, y, z, Stencil_T::idx[W]   ); \
         const real_t vE   = src->get( x, y, z, Stencil_T::idx[E]   ); \
         const real_t vT   = src->get( x, y, z, Stencil_T::idx[T]   ); \
         const real_t vB   = src->get( x, y, z, Stencil_T::idx[B]   ); \
         const real_t vNW  = src->get( x, y, z, Stencil_T::idx[NW]  ); \
         const real_t vNE  = src->get( x, y, z, Stencil_T::idx[NE]  ); \
         const real_t vSW  = src->get( x, y, z, Stencil_T::idx[SW]  ); \
         const real_t vSE  = src->get( x, y, z, Stencil_T::idx[SE]  ); \
         const real_t vTN  = src->get( x, y, z, Stencil_T::idx[TN]  ); \
         const real_t vTS  = src->get( x, y, z, Stencil_T::idx[TS]  ); \
         const real_t vTW  = src->get( x, y, z, Stencil_T::idx[TW]  ); \
         const real_t vTE  = src->get( x, y, z, Stencil_T::idx[TE]  ); \
         const real_t vBN  = src->get( x, y, z, Stencil_T::idx[BN]  ); \
         const real_t vBS  = src->get( x, y, z, Stencil_T::idx[BS]  ); \
         const real_t vBW  = src->get( x, y, z, Stencil_T::idx[BW]  ); \
         const real_t vBE  = src->get( x, y, z, Stencil_T::idx[BE]  ); \
         const real_t vTNE = src->get( x, y, z, Stencil_T::idx[TNE] ); \
         const real_t vTNW = src->get( x, y, z, Stencil_T::idx[TNW] ); \
         const real_t vTSE = src->get( x, y, z, Stencil_T::idx[TSE] ); \
         const real_t vTSW = src->get( x, y, z, Stencil_T::idx[TSW] ); \
         const real_t vBNE = src->get( x, y, z, Stencil_T::idx[BNE] ); \
         const real_t vBNW = src->get( x, y, z, Stencil_T::idx[BNW] ); \
         const real_t vBSE = src->get( x, y, z, Stencil_T::idx[BSE] ); \
         const real_t vBSW = src->get( x, y, z, Stencil_T::idx[BSW] );

#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_INCOMP() \
         const real_t velXTerm = vE + vNE + vSE + vTE + vBE + vTNE + vTSE + vBNE + vBSE; \
         const real_t velYTerm = vN + vNW + vTN + vBN + vTNW + vBNW; \
         const real_t velZTerm = vT + vTS + vTW + vTSW; \
         \
         const real_t rho = vC + vS + vW + vB + vSW + vBS + vBW + vBSW + velXTerm + velYTerm + velZTerm; \
         \
         const real_t velX = velXTerm - vW  - vNW - vSW - vTW - vBW - vTNW - vTSW - vBNW - vBSW ; \
         const real_t velY = velYTerm + vNE + vTNE + vBNE - vS  - vSW - vSE - vTS - vBS  - vTSE - vTSW - vBSE - vBSW; \
         const real_t velZ = velZTerm + vTN + vTE + vTNE + vTNW + vTSE - vB  - vBN - vBS - vBW  - vBE  - vBNE - vBNW - vBSE - vBSW;

#define WALBERLA_LBM_CELLWISE_SWEEP_D3Q27_DENSITY_VELOCITY_COMP() \
         const real_t velXTerm = vE + vNE + vSE + vTE + vBE + vTNE + vTSE + vBNE + vBSE; \
         const real_t velYTerm = vN + vNW + vTN + vBN + vTNW + vBNW; \
         const real_t velZTerm = vT + vTS + vTW + vTSW; \
         \
         const real_t rho = vC + vS + vW + vB + vSW + vBS + vBW + vBSW + velXTerm + velYTerm + velZTerm; \
         const real_t invRho = real_t(1.0) / rho; \
         \
         const real_t velX = invRho * ( velXTerm - vW  - vNW - vSW - vTW - vBW - vTNW - vTSW - vBNW - vBSW ); \
         const real_t velY = invRho * ( velYTerm + vNE + vTNE + vBNE - vS  - vSW - vSE - vTS - vBS  - vTSE - vTSW - vBSE - vBSW ); \
         const real_t velZ = invRho * ( velZTerm + vTN + vTE + vTNE + vTNW + vTSE - vB  - vBN - vBS - vBW  - vBE  - vBNE - vBNW - vBSE - vBSW );

} // namespace lbm
} // namespace walberla



///////////////////////////////////
// LATTICE MODEL SPECIALIZATIONS //
///////////////////////////////////

#include "lbm/srt/CellwiseSweep.impl.h"
#include "lbm/trt/CellwiseSweep.impl.h"
#include "lbm/mrt/CellwiseSweep.impl.h"
#include "lbm/cumulant/CellwiseSweep.impl.h"



/////////////////////////////////
// makeCellwiseSweep FUNCTIONS //
/////////////////////////////////

namespace walberla {
namespace lbm {

// block data IDs of PDF data + arbitrary filter + arbitrary density/velocity calculation (in) and callback (out)

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T > >
makeCellwiseSweep( const BlockDataID & pdfFieldId, const Filter_T & filter,
                   const DensityVelocityIn_T & densityVelocityIn, const DensityVelocityOut_T & densityVelocityOut )
{
   using CS_T = CellwiseSweep<LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T>;
   return shared_ptr< CS_T >( new CS_T( pdfFieldId, filter, densityVelocityIn, densityVelocityOut ) );
}

template< typename LatticeModel_T, typename Filter_T, typename DensityVelocityIn_T, typename DensityVelocityOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T > >
makeCellwiseSweep( const BlockDataID & src, const BlockDataID & dst, const Filter_T & filter,
                   const DensityVelocityIn_T & densityVelocityIn, const DensityVelocityOut_T & densityVelocityOut )
{
   using CS_T = CellwiseSweep<LatticeModel_T, Filter_T, DensityVelocityIn_T, DensityVelocityOut_T>;
   return shared_ptr< CS_T >( new CS_T( src, dst, filter, densityVelocityIn, densityVelocityOut ) );
}

// only block data IDs of PDF data

template< typename LatticeModel_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter, DefaultDensityEquilibriumVelocityCalculation, DefaultDensityVelocityCallback > >
makeCellwiseSweep( const BlockDataID & pdfFieldId )
{
   return makeCellwiseSweep< LatticeModel_T >( pdfFieldId, walberla::field::DefaultEvaluationFilter(),
                                               DefaultDensityEquilibriumVelocityCalculation(), DefaultDensityVelocityCallback() );
}

template< typename LatticeModel_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter, DefaultDensityEquilibriumVelocityCalculation, DefaultDensityVelocityCallback > >
makeCellwiseSweep( const BlockDataID & src, const BlockDataID & dst )
{
   return makeCellwiseSweep< LatticeModel_T >( src, dst, walberla::field::DefaultEvaluationFilter(),
                                               DefaultDensityEquilibriumVelocityCalculation(), DefaultDensityVelocityCallback() );
}

// block data IDs of PDF data + flag field as filter

template< typename LatticeModel_T, typename FlagField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            DefaultDensityEquilibriumVelocityCalculation, DefaultDensityVelocityCallback > >
makeCellwiseSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate )
{
   return makeCellwiseSweep< LatticeModel_T >( pdfFieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                               DefaultDensityEquilibriumVelocityCalculation(), DefaultDensityVelocityCallback() );
}

template< typename LatticeModel_T, typename FlagField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            DefaultDensityEquilibriumVelocityCalculation, DefaultDensityVelocityCallback > >
makeCellwiseSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate )
{
   return makeCellwiseSweep< LatticeModel_T >( src, dst, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                               DefaultDensityEquilibriumVelocityCalculation(), DefaultDensityVelocityCallback() );
}

// block data IDs of PDF data + flag field as filter + block data ID of velocity field (out)

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            DefaultDensityEquilibriumVelocityCalculation, VelocityCallback<VelocityField_T> > >
makeCellwiseSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                   const BlockDataID & velocityFieldId )
{
   return makeCellwiseSweep< LatticeModel_T >( pdfFieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                               DefaultDensityEquilibriumVelocityCalculation(), VelocityCallback<VelocityField_T>( velocityFieldId ) );
}

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            DefaultDensityEquilibriumVelocityCalculation, VelocityCallback<VelocityField_T> > >
makeCellwiseSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                   const BlockDataID & velocityFieldId )
{
   return makeCellwiseSweep< LatticeModel_T >( src, dst, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                               DefaultDensityEquilibriumVelocityCalculation(), VelocityCallback<VelocityField_T>( velocityFieldId ) );
}

// block data IDs of PDF data + flag field as filter + block data IDs of velocity and density field (out)

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, typename DensityField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            DefaultDensityEquilibriumVelocityCalculation, DensityVelocityCallback<VelocityField_T,DensityField_T> > >
makeCellwiseSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                   const BlockDataID & velocityFieldId, const BlockDataID & densityFieldId  )
{
   return makeCellwiseSweep< LatticeModel_T >( pdfFieldId, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                               DefaultDensityEquilibriumVelocityCalculation(),
                                               DensityVelocityCallback<VelocityField_T,DensityField_T>( velocityFieldId, densityFieldId ) );
}

template< typename LatticeModel_T, typename FlagField_T, typename VelocityField_T, typename DensityField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            DefaultDensityEquilibriumVelocityCalculation, DensityVelocityCallback<VelocityField_T,DensityField_T> > >
makeCellwiseSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                   const BlockDataID & velocityFieldId, const BlockDataID & densityFieldId )
{
   return makeCellwiseSweep< LatticeModel_T >( src, dst, walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                               DefaultDensityEquilibriumVelocityCalculation(),
                                               DensityVelocityCallback<VelocityField_T,DensityField_T>( velocityFieldId, densityFieldId ) );
}



/////////////////////////
// ADVECTION DIFFUSION //
/////////////////////////

// block data IDs of PDF data + block data ID of velocity field + arbitrary filter + arbitrary density/velocity callback (out)

template< typename LatticeModel_T, typename VelocityField_T, typename Filter_T, typename DensityVelocityOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, Filter_T, AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DensityVelocityOut_T > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & velocityFieldId, const Filter_T & filter,
                                     const DensityVelocityOut_T & densityVelocityOut )
{
   using CS_T = CellwiseSweep<LatticeModel_T, Filter_T, AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DensityVelocityOut_T>;
   return shared_ptr< CS_T >( new CS_T( pdfFieldId, filter,
                                        AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>( velocityFieldId ), densityVelocityOut ) );
}

template< typename LatticeModel_T, typename VelocityField_T, typename Filter_T, typename DensityVelocityOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, Filter_T, AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DensityVelocityOut_T > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & velocityFieldId,
                                     const Filter_T & filter, const DensityVelocityOut_T & densityVelocityOut )
{
   using CS_T = CellwiseSweep<LatticeModel_T, Filter_T, AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DensityVelocityOut_T>;
   return shared_ptr< CS_T >( new CS_T( src, dst, filter,
                                        AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>( velocityFieldId ), densityVelocityOut ) );
}

// only block data IDs of PDF data + block data ID of velocity field

template< typename LatticeModel_T, typename VelocityField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DefaultDensityVelocityCallback > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & velocityFieldId )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( pdfFieldId, velocityFieldId, walberla::field::DefaultEvaluationFilter(),
                                                                                  DefaultDensityVelocityCallback() );
}

template< typename LatticeModel_T, typename VelocityField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::DefaultEvaluationFilter,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DefaultDensityVelocityCallback > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & velocityFieldId )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( src, dst, velocityFieldId, walberla::field::DefaultEvaluationFilter(),
                                                                                  DefaultDensityVelocityCallback() );
}

// block data IDs of PDF data + block data ID of velocity field + flag field as filter

template< typename LatticeModel_T, typename VelocityField_T, typename FlagField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DefaultDensityVelocityCallback > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & velocityFieldId,
                                     const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( pdfFieldId, velocityFieldId,
                                                                                  walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                                                  DefaultDensityVelocityCallback() );
}

template< typename LatticeModel_T, typename VelocityField_T, typename FlagField_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DefaultDensityVelocityCallback > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & velocityFieldId,
                                     const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( src, dst, velocityFieldId,
                                                                                  walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                                                  DefaultDensityVelocityCallback() );
}

// block data IDs of PDF data + block data ID of velocity field + flag field as filter + block data ID of velocity field (out)

template< typename LatticeModel_T, typename VelocityField_T, typename FlagField_T, typename VelocityFieldOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, VelocityCallback<VelocityFieldOut_T> > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & velocityFieldId, const ConstBlockDataID & flagFieldId,
                                     const Set< FlagUID > & cellsToEvaluate, const BlockDataID & velocityFieldOutId )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( pdfFieldId, velocityFieldId,
                                                                                  walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                                                  VelocityCallback<VelocityFieldOut_T>( velocityFieldOutId ) );
}

template< typename LatticeModel_T, typename VelocityField_T, typename FlagField_T, typename VelocityFieldOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, VelocityCallback<VelocityFieldOut_T> > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & velocityFieldId,
                                     const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate, const BlockDataID & velocityFieldOutId )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( src, dst, velocityFieldId,
                                                                                  walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                                                  VelocityCallback<VelocityFieldOut_T>( velocityFieldOutId ) );
}

// block data IDs of PDF data + block data ID of velocity field + flag field as filter + block data IDs of velocity and density field (out)

template< typename LatticeModel_T, typename VelocityField_T, typename FlagField_T, typename VelocityFieldOut_T, typename DensityFieldOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DensityVelocityCallback<VelocityFieldOut_T,DensityFieldOut_T> > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & pdfFieldId, const ConstBlockDataID & velocityFieldId, const ConstBlockDataID & flagFieldId,
                                     const Set< FlagUID > & cellsToEvaluate, const BlockDataID & velocityFieldOutId, const BlockDataID & densityFieldOutId )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( pdfFieldId, velocityFieldId,
                                                                                  walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                                                  DensityVelocityCallback<VelocityFieldOut_T,DensityFieldOut_T>( velocityFieldOutId, densityFieldOutId ) );
}

template< typename LatticeModel_T, typename VelocityField_T, typename FlagField_T, typename VelocityFieldOut_T, typename DensityFieldOut_T >
shared_ptr< CellwiseSweep< LatticeModel_T, walberla::field::FlagFieldEvaluationFilter<FlagField_T>,
            AdvectionDiffusionDensityEquilibriumVelocityCalculation<VelocityField_T>, DensityVelocityCallback<VelocityFieldOut_T,DensityFieldOut_T> > >
makeCellwiseAdvectionDiffusionSweep( const BlockDataID & src, const BlockDataID & dst, const ConstBlockDataID & velocityFieldId,
                                     const ConstBlockDataID & flagFieldId, const Set< FlagUID > & cellsToEvaluate,
                                     const BlockDataID & velocityFieldOutId, const BlockDataID & densityFieldOutId )
{
   return makeCellwiseAdvectionDiffusionSweep< LatticeModel_T, VelocityField_T >( src, dst, velocityFieldId,
                                                                                  walberla::field::FlagFieldEvaluationFilter<FlagField_T>( flagFieldId, cellsToEvaluate ),
                                                                                  DensityVelocityCallback<VelocityFieldOut_T,DensityFieldOut_T>( velocityFieldOutId, densityFieldOutId ) );
}



} // namespace lbm
} // namespace walberla
