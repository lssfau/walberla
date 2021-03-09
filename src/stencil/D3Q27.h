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
        for( uint_t i = 0; i < D3Q27::Size; ++i )
        {
           int cx = cx[ D3Q27::dir[i] ];
           Direction inverse = inverseDir[ D3Q27::dir[i] ];
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
        for( auto dir = D3Q27::begin(); dir != D3Q27::end(); ++dir )
        {
          int cx = dir.cx();
          Direction inverse = dir.inverseDir();
        }
      \endcode
      *
      */
      //***************************************************************************************************************
      template<typename Dummy = int>
      struct D3Q27
      {
         //**Member Variables******************************************************************************************
         /*! \name Member Variables*/
         //@{

         static const char * NAME;

         static const uint_t D = 3;
         static const uint_t Q = 27;

         static const uint_t POS_Q = 27 / 2;

         static const uint_t Dimension = 3;
         static const uint_t Size      = 27;

         static const bool   containsCenter   = true;
         static const uint_t noCenterFirstIdx = 1;

         static const Direction dir           [27];
         static const Direction dir_pos       [POS_Q];
         static const uint_t    idx           [NR_OF_DIRECTIONS];
         static const Direction d_per_d       [NR_OF_DIRECTIONS][27/2];
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

         using iterator = stencil::Iterator<D3Q27>;

         static iterator begin()           { return iterator(0); }
         static iterator beginNoCenter()   { return iterator(noCenterFirstIdx); }
         static iterator end()             { return iterator(27); }

         //}
         //************************************************************************************************************
      };


      template<typename Dummy> const uint_t D3Q27<Dummy>::D;
      template<typename Dummy> const uint_t D3Q27<Dummy>::Q;
      template<typename Dummy> const uint_t D3Q27<Dummy>::POS_Q;
      template<typename Dummy> const uint_t D3Q27<Dummy>::Dimension;
      template<typename Dummy> const uint_t D3Q27<Dummy>::Size;
      template<typename Dummy> const bool   D3Q27<Dummy>::containsCenter;
      template<typename Dummy> const uint_t D3Q27<Dummy>::noCenterFirstIdx;

      template<typename Dummy> const char * D3Q27<Dummy>::NAME = "D3Q27";

      /// Subset of directions. Defines the stencil
      template<typename Dummy>
      const Direction D3Q27<Dummy>::dir[27] = { C,N,S,W,E,T,B,NW,NE,SW,SE,TN,TS,TW,TE,BN,BS,BW,BE,TNE,TNW,TSE,TSW,BNE,BNW,BSE,BSW };


      /**
       * \brief Contains only half of the directions ( the positive ones )
       *
       * Use this f.e. in situations where direction have to be handled pairwise,
       * the direction together with its inverse direction.
       */
      template<typename Dummy>
      const Direction D3Q27<Dummy>::dir_pos[POS_Q] = { N,E,T,NE,SE,TN,TE,BN,BE,TNE,TSE,BNE,BSE };


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
      const uint_t D3Q27<Dummy>::idx[NR_OF_DIRECTIONS] = { 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };


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
      const Direction D3Q27<Dummy>::d_per_d[NR_OF_DIRECTIONS][27/2] = { {C},
								{N,NW,NE,TN,BN,TNE,TNW,BNE,BNW},
								{S,SW,SE,TS,BS,TSE,TSW,BSE,BSW},
								{W,NW,SW,TW,BW,TNW,TSW,BNW,BSW},
								{E,NE,SE,TE,BE,TNE,TSE,BNE,BSE},
								{T,TN,TS,TW,TE,TNE,TNW,TSE,TSW},
								{B,BN,BS,BW,BE,BNE,BNW,BSE,BSW},
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
      const uint_t D3Q27<Dummy>::d_per_d_length [NR_OF_DIRECTIONS] = { 1,9,9,9,9,9,9,3,3,3,3,3,3,3,3,3,3,3,3,1,1,1,1,1,1,1,1 };


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
      const Direction D3Q27<Dummy>::dir_neighbors[NR_OF_DIRECTIONS][NR_OF_DIRECTIONS] = { {N,S,W,E,T,B,NW,NE,SW,SE,TN,TS,TW,TE,BN,BS,BW,BE,TNE,TNW,TSE,TSW,BNE,BNW,BSE,BSW},
								{C,W,E,T,B,NW,NE,TN,TW,TE,BN,BW,BE,TNE,TNW,BNE,BNW},
								{C,W,E,T,B,SW,SE,TS,TW,TE,BS,BW,BE,TSE,TSW,BSE,BSW},
								{C,N,S,T,B,NW,SW,TN,TS,TW,BN,BS,BW,TNW,TSW,BNW,BSW},
								{C,N,S,T,B,NE,SE,TN,TS,TE,BN,BS,BE,TNE,TSE,BNE,BSE},
								{C,N,S,W,E,NW,NE,SW,SE,TN,TS,TW,TE,TNE,TNW,TSE,TSW},
								{C,N,S,W,E,NW,NE,SW,SE,BN,BS,BW,BE,BNE,BNW,BSE,BSW},
								{C,N,W,T,B,TN,TW,BN,BW,TNW,BNW},
								{C,N,E,T,B,TN,TE,BN,BE,TNE,BNE},
								{C,S,W,T,B,TS,TW,BS,BW,TSW,BSW},
								{C,S,E,T,B,TS,TE,BS,BE,TSE,BSE},
								{C,N,W,E,T,NW,NE,TW,TE,TNE,TNW},
								{C,S,W,E,T,SW,SE,TW,TE,TSE,TSW},
								{C,N,S,W,T,NW,SW,TN,TS,TNW,TSW},
								{C,N,S,E,T,NE,SE,TN,TS,TNE,TSE},
								{C,N,W,E,B,NW,NE,BW,BE,BNE,BNW},
								{C,S,W,E,B,SW,SE,BW,BE,BSE,BSW},
								{C,N,S,W,B,NW,SW,BN,BS,BNW,BSW},
								{C,N,S,E,B,NE,SE,BN,BS,BNE,BSE},
								{C,N,E,T,NE,TN,TE},
								{C,N,W,T,NW,TN,TW},
								{C,S,E,T,SE,TS,TE},
								{C,S,W,T,SW,TS,TW},
								{C,N,E,B,NE,BN,BE},
								{C,N,W,B,NW,BN,BW},
								{C,S,E,B,SE,BS,BE},
								{C,S,W,B,SW,BS,BW} };


      /**
       * \brief Length of the dir_neighbors array
       * For usage see documentation of dir_neighbors
       */
      template<typename Dummy>
      const uint_t D3Q27<Dummy>::dir_neighbors_length[NR_OF_DIRECTIONS] = { 26,17,17,17,17,17,17,11,11,11,11,11,11,11,11,11,11,11,11,7,7,7,7,7,7,7,7 };

   } // namespace internal

   using D3Q27 = internal::D3Q27<>;

} // namespace stencil
} // namespace walberla


#endif // DOXY_SKIP_INTERNAL
