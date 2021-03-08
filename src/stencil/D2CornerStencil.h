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
        for( uint_t i = 0; i < D2CornerStencil::Size; ++i )
        {
           int cx = cx[ D2CornerStencil::dir[i] ];
           Direction inverse = inverseDir[ D2CornerStencil::dir[i] ];
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
        for( auto dir = D2CornerStencil::begin(); dir != D2CornerStencil::end(); ++dir )
        {
          int cx = dir.cx();
          Direction inverse = dir.inverseDir();
        }
      \endcode
      *
      */
      //***************************************************************************************************************
      template<typename Dummy = int>
      struct D2CornerStencil
      {
         //**Member Variables******************************************************************************************
         /*! \name Member Variables*/
         //@{

         static const char * NAME;

         static const uint_t D = 2;
         static const uint_t Q = 4;

         static const uint_t POS_Q = 4 / 2;

         static const uint_t Dimension = 2;
         static const uint_t Size      = 4;

         static const bool   containsCenter   = false;
         static const uint_t noCenterFirstIdx = 0;

         static const Direction dir           [4];
         static const Direction dir_pos       [POS_Q];
         static const uint_t    idx           [NR_OF_DIRECTIONS];
         static const Direction d_per_d       [NR_OF_DIRECTIONS][4/2];
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

         using iterator = stencil::Iterator<D2CornerStencil>;

         static iterator begin()           { return iterator(0); }
         static iterator beginNoCenter()   { return iterator(noCenterFirstIdx); }
         static iterator end()             { return iterator(4); }

         //}
         //************************************************************************************************************
      };


      template<typename Dummy> const uint_t D2CornerStencil<Dummy>::D;
      template<typename Dummy> const uint_t D2CornerStencil<Dummy>::Q;
      template<typename Dummy> const uint_t D2CornerStencil<Dummy>::POS_Q;
      template<typename Dummy> const uint_t D2CornerStencil<Dummy>::Dimension;
      template<typename Dummy> const uint_t D2CornerStencil<Dummy>::Size;
      template<typename Dummy> const bool   D2CornerStencil<Dummy>::containsCenter;
      template<typename Dummy> const uint_t D2CornerStencil<Dummy>::noCenterFirstIdx;

      template<typename Dummy> const char * D2CornerStencil<Dummy>::NAME = "D2CornerStencil";

      /// Subset of directions. Defines the stencil
      template<typename Dummy>
      const Direction D2CornerStencil<Dummy>::dir[4] = { NW,NE,SW,SE };


      /**
       * \brief Contains only half of the directions ( the positive ones )
       *
       * Use this f.e. in situations where direction have to be handled pairwise,
       * the direction together with its inverse direction.
       */
      template<typename Dummy>
      const Direction D2CornerStencil<Dummy>::dir_pos[POS_Q] = { NE,SE };


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
      const uint_t D2CornerStencil<Dummy>::idx[NR_OF_DIRECTIONS] = { INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,0,1,2,3,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR,INVALID_DIR };


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
      const Direction D2CornerStencil<Dummy>::d_per_d[NR_OF_DIRECTIONS][4/2] = { {},
								{NW,NE},
								{SW,SE},
								{NW,SW},
								{NE,SE},
								{},
								{},
								{NW},
								{NE},
								{SW},
								{SE},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{},
								{} };


      /**
       * \brief Length of the d_per_d array
       * For usage see documentation of d_per_d
       */
      template<typename Dummy>
      const uint_t D2CornerStencil<Dummy>::d_per_d_length [NR_OF_DIRECTIONS] = { 0,2,2,2,2,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };


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
      const Direction D2CornerStencil<Dummy>::dir_neighbors[NR_OF_DIRECTIONS][NR_OF_DIRECTIONS] = { {NW,NE,SW,SE},
								{W,E},
								{W,E},
								{N,S},
								{N,S},
								{TNE,TNW,TSE,TSW},
								{BNE,BNW,BSE,BSW},
								{C},
								{C},
								{C},
								{C},
								{TW,TE},
								{TW,TE},
								{TN,TS},
								{TN,TS},
								{BW,BE},
								{BW,BE},
								{BN,BS},
								{BN,BS},
								{T},
								{T},
								{T},
								{T},
								{B},
								{B},
								{B},
								{B} };


      /**
       * \brief Length of the dir_neighbors array
       * For usage see documentation of dir_neighbors
       */
      template<typename Dummy>
      const uint_t D2CornerStencil<Dummy>::dir_neighbors_length[NR_OF_DIRECTIONS] = { 4,2,2,2,2,4,4,1,1,1,1,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1 };

   } // namespace internal

   using D2CornerStencil = internal::D2CornerStencil<>;

} // namespace stencil
} // namespace walberla


#endif // DOXY_SKIP_INTERNAL
