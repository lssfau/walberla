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
        for( uint_t i = 0; i < $name::Size; ++i )
        {
           int cx = cx[ $name::dir[i] ];
           Direction inverse = inverseDir[ $name::dir[i] ];
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
        for( auto dir = $name::begin(); dir != $name::end(); ++dir )
        {
          int cx = dir.cx();
          Direction inverse = dir.inverseDir();
        }
      \endcode
      *
      */
      //***************************************************************************************************************
      template<typename Dummy = int>
      struct $name
      {
         //**Member Variables******************************************************************************************
         /*! \name Member Variables*/
         //@{

         static const char * NAME;

         static const uint_t D = $D;
         static const uint_t Q = $Q;

         static const uint_t POS_Q = $Q / 2;

         static const uint_t Dimension = $D;
         static const uint_t Size      = $Q;

         static const bool   containsCenter   = $containsCenter;
         static const uint_t noCenterFirstIdx = $noCenterFirstIndex;

         static const std::array< Direction, $Q > dir;
         static const std::array< Direction, POS_Q > dir_pos;
         static const std::array< uint_t, NR_OF_DIRECTIONS > idx;
         static const std::array< std::array<Direction, $Q/2>, NR_OF_DIRECTIONS> d_per_d;
         static const std::array< uint_t, NR_OF_DIRECTIONS > d_per_d_length;

         static const std::array< std::array< Direction, NR_OF_DIRECTIONS >, NR_OF_DIRECTIONS > dir_neighbors;
         static const std::array< uint_t, NR_OF_DIRECTIONS > dir_neighbors_length;

         static bool   containsDir(Direction d) { return idx[d] < NR_OF_DIRECTIONS; }
         static uint_t invDirIdx  (Direction d) { return idx[stencil::inverseDir[d]]; }
         //}
         //************************************************************************************************************


         //** Iteration   *********************************************************************************************
         /*! \name Iteration*/
         //@{

         using iterator = stencil::Iterator<$name>;

         static iterator begin()           { return iterator(0); }
         static iterator beginNoCenter()   { return iterator(noCenterFirstIdx); }
         static iterator end()             { return iterator($Q); }

         //}
         //************************************************************************************************************
      };


      template<typename Dummy> const uint_t $name<Dummy>::D;
      template<typename Dummy> const uint_t $name<Dummy>::Q;
      template<typename Dummy> const uint_t $name<Dummy>::POS_Q;
      template<typename Dummy> const uint_t $name<Dummy>::Dimension;
      template<typename Dummy> const uint_t $name<Dummy>::Size;
      template<typename Dummy> const bool   $name<Dummy>::containsCenter;
      template<typename Dummy> const uint_t $name<Dummy>::noCenterFirstIdx;

      template<typename Dummy> const char * $name<Dummy>::NAME = "$name";

      /// Subset of directions. Defines the stencil
      template<typename Dummy>
      const std::array< Direction, $Q > $name<Dummy>::dir{ $dirs };


      /**
       * \brief Contains only half of the directions ( the positive ones )
       *
       * Use this f.e. in situations where direction have to be handled pairwise,
       * the direction together with its inverse direction.
       */
      template<typename Dummy>
      const std::array< Direction, $name<Dummy>::POS_Q > $name<Dummy>::dir_pos{ $dir_pos };


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
      const std::array< uint_t, NR_OF_DIRECTIONS > $name<Dummy>::idx{ $indexFromDir };


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
      const std::array< std::array< Direction, $Q/2 >, NR_OF_DIRECTIONS > $name<Dummy>::d_per_d{ { $d_per_d } };


      /**
       * \brief Length of the d_per_d array
       * For usage see documentation of d_per_d
       */
      template<typename Dummy>
      const std::array< uint_t, NR_OF_DIRECTIONS > $name<Dummy>::d_per_d_length{ $d_per_d_length };


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
      const std::array< std::array< Direction, NR_OF_DIRECTIONS >, NR_OF_DIRECTIONS > $name<Dummy>::dir_neighbors{ { $dir_neighbors } };


      /**
       * \brief Length of the dir_neighbors array
       * For usage see documentation of dir_neighbors
       */
      template<typename Dummy>
      const std::array<uint_t, NR_OF_DIRECTIONS> $name<Dummy>::dir_neighbors_length{ $dir_neighbors_length };

   } // namespace internal

   using $name = internal::$name<>;

} // namespace stencil
} // namespace walberla

