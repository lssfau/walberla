//====================================================================================================================
//  Caution: This file has been generated automatically. All manual changes are lost when file is regenerated!
//           Changes should be done in Stencil.in.h,and then all stencils classes can be generated again.
//====================================================================================================================
#ifndef DOXY_SKIP_INTERNAL
#pragma once

#include "Directions.h"
#include "Iterator.h"


namespace walberla {
namespace stencil {


   namespace internal
   {
      //***************************************************************************************************************
      /*! \brief Defines a D-dimensional stencil info with Q directions.
      *
      * \ingroup stencil
      *
      * A stencil is defined by picking a subset of all possible directions. These directions are listed in the
      * dir array. So when iterating a stencil, one iterates the dir array which contains the directions.
      * To get the properties for these directions, the global arrays of Directions.h can be used:
      *
      \code
        using namespace stencil;
        for( uint_t i = 0; i < D3EdgeCornerStencil::Size; ++i )
        {
           int cx = cx[ D3EdgeCornerStencil::dir[i] ];
           Direction inverse = inverseDir[ D3EdgeCornerStencil::dir[i] ];
        }
      \endcode
      *
      * Since all these arrays are constant, the lookup can be resolved at compile time and no performance overhead
      * should be generated.
      *
      * For more convenient iteration, use the provided iterator.
      * The following code is equivalent to the code example above:
      \code
        using namespace stencil;
        for( auto dir = D3EdgeCornerStencil::begin(); dir != D3EdgeCornerStencil::end(); ++dir )
        {
          int cx = dir.cx();
          Direction inverse = dir.inverseDir();
        }
      \endcode
      *
      */
      //***************************************************************************************************************
      template<typename Dummy = int>
      struct D3EdgeCornerStencil
      {
         //**Member Variables******************************************************************************************
         /*! \name Member Variables*/
         //@{

         static const char * NAME;

         static const uint_t D = 3;
         static const uint_t Q = 20;

         static const uint_t POS_Q = 20 / 2;

         static const uint_t Dimension = 3;
         static const uint_t Size      = 20;

         static const bool   containsCenter   = false;
         static const uint_t noCenterFirstIdx = 0;

         static const Direction dir           [20];
         static const Direction dir_pos       [POS_Q];
         static const uint_t    idx           [NR_OF_DIRECTIONS];
         static const Direction d_per_d       [NR_OF_DIRECTIONS][20/2];
         static const uint_t    d_per_d_length[NR_OF_DIRECTIONS];

         static const Direction dir_neighbors        [NR_OF_DIRECTIONS][NR_OF_DIRECTIONS];
         static const uint_t    dir_neighbors_length [NR_OF_DIRECTIONS];

         static bool   containsDir(Direction d) { return idx[d] < NR_OF_DIRECTIONS; }
         static uint_t invDirIdx  (Direction d) { return idx[stencil::inverseDir[d]]; }
         //}
         //************************************************************************************************************


         //** Iteration   *********************************************************************************************
         /*! \name Iteration*/
         //@{

         using iterator = stencil::Iterator<D3EdgeCornerStencil>;

         static iterator begin()           { return iterator(0); }
         static iterator beginNoCenter()   { return iterator(noCenterFirstIdx); }
         static iterator end()             { return iterator(20); }

         //}
         //************************************************************************************************************
      };


      template<typename Dummy> const uint_t D3EdgeCornerStencil<Dummy>::D;
      template<typename Dummy> const uint_t D3EdgeCornerStencil<Dummy>::Q;
      template<typename Dummy> const uint_t D3EdgeCornerStencil<Dummy>::POS_Q;
      template<typename Dummy> const uint_t D3EdgeCornerStencil<Dummy>::Dimension;
      template<typename Dummy> const uint_t D3EdgeCornerStencil<Dummy>::Size;
      template<typename Dummy> const bool   D3EdgeCornerStencil<Dummy>::containsCenter;
      template<typename Dummy> const uint_t D3EdgeCornerStencil<Dummy>::noCenterFirstIdx;

      template<typename Dummy> const char * D3EdgeCornerStencil<Dummy>::NAME = "D3EdgeCornerStencil";

      /// Subset of directions. Defines the stencil
      template<typename Dummy>
      const Direction D3EdgeCornerStencil<Dummy>::dir[20] = { NW,NE,SW,SE,TN,TS,TW,TE,BN,BS,BW,BE,TNE,TNW,TSE,TSW,BNE,BNW,BSE,BSW };


      /**
       * \brief Contains only half of the directions ( the positive ones )
       *
       * Use this f.e. in situations where direction have to be handled pairwise,
       * the direction together with its inverse direction.
       */
      template<typename Dummy>
      const Direction D3EdgeCornerStencil<Dummy>::dir_pos[POS_Q] = { NE,SE,TN,TE,BN,BE,TNE,TSE,BNE,BSE };


      /**
       * \brief Maps direction enums, to an index
       *
       *  Use this when working with fields: The direction enum
       *  cannot be used as field index directly. Consider having a D2Q4
       *  stencil, then the fourth field index f ranges from 0 to 3, but the
       *  direction enums are S=1,N=2,E=3,W=4 start at 1.
       *  So a stencil class is needed to map back from direction to field index
       */
      template<typename Dummy>
      const uint_t D3EdgeCornerStencil<Dummy>::idx[NR_OF_DIRECTIONS] = { INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19 };


      /**
       * \brief Maps a direction to a set of sub-directions
       *
       * A sub-direction is defined as a direction that contains the letters of the original direction
       * i.e. SW is a sub-direction of S and W.
       * the following code shows how to iterate over all sub-directions, that are contained in the stencil of W
       *
       \code
       template<typename Stencil>
       //...
       stencil::Direction dir = stencil::W;
       for (int i=0; i < Stencil::d_per_d_length[dir]; ++i)
           stencil::Direction subDirection = Stencil::d_per_d[dir][i];
       \endcode
       *
       */
      template<typename Dummy>
      const Direction D3EdgeCornerStencil<Dummy>::d_per_d[NR_OF_DIRECTIONS][20/2] = { {},
								{NW,NE,TN,BN,TNE,TNW,BNE,BNW},
								{SW,SE,TS,BS,TSE,TSW,BSE,BSW},
								{NW,SW,TW,BW,TNW,TSW,BNW,BSW},
								{NE,SE,TE,BE,TNE,TSE,BNE,BSE},
								{TN,TS,TW,TE,TNE,TNW,TSE,TSW},
								{BN,BS,BW,BE,BNE,BNW,BSE,BSW},
								{NW,TNW,BNW},
								{NE,TNE,BNE},
								{SW,TSW,BSW},
								{SE,TSE,BSE},
								{TN,TNE,TNW},
								{TS,TSE,TSW},
								{TW,TNW,TSW},
								{TE,TNE,TSE},
								{BN,BNE,BNW},
								{BS,BSE,BSW},
								{BW,BNW,BSW},
								{BE,BNE,BSE},
								{TNE},
								{TNW},
								{TSE},
								{TSW},
								{BNE},
								{BNW},
								{BSE},
								{BSW} };


      /**
       * \brief Length of the d_per_d array
       * For usage see documentation of d_per_d
       */
      template<typename Dummy>
      const uint_t D3EdgeCornerStencil<Dummy>::d_per_d_length [NR_OF_DIRECTIONS] = { 0,8,8,8,8,8,8,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1 };


      /**
       * \brief Views directions as cells in a 3x3x3 grid. Describes neighborhood between cells/directions.
       *
       * Basis a 3x3x3 grid where every cell corresponds to a direction.
       * The following array provides answer to the question: What are the neighboring cells of a given
       * cell using the current stencil. If a neighbor is not in this 3x3x3 grid it can not be described
       * by a direction and is not returned. Therefore N , for example, has more neighbors that NW.
       * By definition, the neighbors of C are identical to the dir[] array without C itself.
       */
      template<typename Dummy>
      const Direction D3EdgeCornerStencil<Dummy>::dir_neighbors[NR_OF_DIRECTIONS][NR_OF_DIRECTIONS] = { {NW,NE,SW,SE,TN,TS,TW,TE,BN,BS,BW,BE,TNE,TNW,TSE,TSW,BNE,BNW,BSE,BSW},
								{W,E,T,B,TW,TE,BW,BE,TNE,TNW,BNE,BNW},
								{W,E,T,B,TW,TE,BW,BE,TSE,TSW,BSE,BSW},
								{N,S,T,B,TN,TS,BN,BS,TNW,TSW,BNW,BSW},
								{N,S,T,B,TN,TS,BN,BS,TNE,TSE,BNE,BSE},
								{N,S,W,E,NW,NE,SW,SE,TNE,TNW,TSE,TSW},
								{N,S,W,E,NW,NE,SW,SE,BNE,BNW,BSE,BSW},
								{C,T,B,TN,TW,BN,BW},
								{C,T,B,TN,TE,BN,BE},
								{C,T,B,TS,TW,BS,BW},
								{C,T,B,TS,TE,BS,BE},
								{C,W,E,NW,NE,TW,TE},
								{C,W,E,SW,SE,TW,TE},
								{C,N,S,NW,SW,TN,TS},
								{C,N,S,NE,SE,TN,TS},
								{C,W,E,NW,NE,BW,BE},
								{C,W,E,SW,SE,BW,BE},
								{C,N,S,NW,SW,BN,BS},
								{C,N,S,NE,SE,BN,BS},
								{C,N,E,T},
								{C,N,W,T},
								{C,S,E,T},
								{C,S,W,T},
								{C,N,E,B},
								{C,N,W,B},
								{C,S,E,B},
								{C,S,W,B} };


      /**
       * \brief Length of the dir_neighbors array
       * For usage see documentation of dir_neighbors
       */
      template<typename Dummy>
      const uint_t D3EdgeCornerStencil<Dummy>::dir_neighbors_length[NR_OF_DIRECTIONS] = { 20,12,12,12,12,12,12,7,7,7,7,7,7,7,7,7,7,7,7,4,4,4,4,4,4,4,4 };

   } // namespace internal

   using D3EdgeCornerStencil = internal::D3EdgeCornerStencil<>;

} // namespace stencil
} // namespace walberla


#endif // DOXY_SKIP_INTERNAL
