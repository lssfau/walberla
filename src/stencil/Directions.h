//=====================================================================================================================
/*!
 *  \file   Directions.h
 *  \brief  Defines all stencil directions and their properties, and a general DxQy stencil class
 *  \author Martin Bauer <martin.bauer@fau.de>
 */
//=====================================================================================================================

#pragma once

// core includes
#include "core/DataTypes.h"
#include "core/cell/Cell.h"
#include "core/debug/Debug.h"

// STL includes
#include <string>
#include <cmath>

namespace walberla {

namespace stencil {


   const uint_t NR_OF_DIRECTIONS = 27;

   /*******************************************************************************************************************
    * Enumeration of all possible directions
    * \ingroup stencil
    *
    * Stencils hold a subset of these directions. Changes in this enum will affect the complete stencil module.
    * All arrays that define direction properties (dx,dy,dz, mirrorX, ...) have to be adapted. Also the direction list
    * in the stencil generation code (generate.py) has to be adapted.
    ******************************************************************************************************************/
   enum Direction
   {
      C   =  0,  //!< Center
      N   =  1,  //!< North
      S   =  2,  //!< South
      W   =  3,  //!< West
      E   =  4,  //!< East
      T   =  5,  //!< Top
      B   =  6,  //!< Bottom
      NW  =  7,  //!< North-West
      NE  =  8,  //!< North-East
      SW  =  9,  //!< South-West
      SE  = 10,  //!< South-East
      TN  = 11,  //!< Top-North
      TS  = 12,  //!< Top-South
      TW  = 13,  //!< Top-West
      TE  = 14,  //!< Top-East
      BN  = 15,  //!< Bottom-North
      BS  = 16,  //!< Bottom-South
      BW  = 17,  //!< Bottom-West
      BE  = 18,  //!< Bottom-East
      TNE = 19,  //!< Top-North-East
      TNW = 20,  //!< Top-North-West
      TSE = 21,  //!< Top-South-East
      TSW = 22,  //!< Top-South-West
      BNE = 23,  //!< Bottom-North-East
      BNW = 24,  //!< Bottom-North-West
      BSE = 25,  //!< Bottom-South-East
      BSW = 26,  //!< Bottom-South-West
      INVALID_DIR = 27  //!< Invalid direction
   };


   /*******************************************************************************************************************
    * Binary direction enumeration where each bit stands for a direction
    * \ingroup stencil
    *
    * Every direction is represented by a single bit. It is NOT the case that for example NW has to bits set,
    * and BNW three bits!
    ******************************************************************************************************************/
   enum BinaryDirection
   {
      Bin_C   =  1<<0,  //!< Center
      Bin_N   =  1<<1,  //!< North
      Bin_S   =  1<<2,  //!< South
      Bin_W   =  1<<3,  //!< West
      Bin_E   =  1<<4,  //!< East
      Bin_T   =  1<<5,  //!< Top
      Bin_B   =  1<<6,  //!< Bottom
      Bin_NW  =  1<<7,  //!< North-West
      Bin_NE  =  1<<8,  //!< North-East
      Bin_SW  =  1<<9,  //!< South-West
      Bin_SE  = 1<<10,  //!< South-East
      Bin_TN  = 1<<11,  //!< Top-North
      Bin_TS  = 1<<12,  //!< Top-South
      Bin_TW  = 1<<13,  //!< Top-West
      Bin_TE  = 1<<14,  //!< Top-East
      Bin_BN  = 1<<15,  //!< Bottom-North
      Bin_BS  = 1<<16,  //!< Bottom-South
      Bin_BW  = 1<<17,  //!< Bottom-West
      Bin_BE  = 1<<18,  //!< Bottom-East
      Bin_TNE = 1<<19,  //!< Top-North-East
      Bin_TNW = 1<<20,  //!< Top-North-West
      Bin_TSE = 1<<21,  //!< Top-South-East
      Bin_TSW = 1<<22,  //!< Top-South-West
      Bin_BNE = 1<<23,  //!< Bottom-North-East
      Bin_BNW = 1<<24,  //!< Bottom-North-West
      Bin_BSE = 1<<25,  //!< Bottom-South-East
      Bin_BSW = 1<<26   //!< Bottom-South-West
   };

   /// The x component for each direction  \ingroup stencil
   const int cx[NR_OF_DIRECTIONS] =  {
   // C   N   S   W   E   T   B  NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE TNE TNW TSE TSW BNE BNW BSE BSW
      0,  0,  0, -1,  1,  0,  0, -1,  1, -1,  1,  0,  0, -1,  1,  0,  0, -1,  1,  1, -1,  1, -1,  1, -1,  1, -1
   };

   /// The y component for each direction \ingroup stencil
   const int cy[NR_OF_DIRECTIONS] =  {
   // C   N   S   W   E   T   B  NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE TNE TNW TSE TSW BNE BNW BSE BSW
      0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1, -1,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1
   };

   /// The z component for each direction \ingroup stencil
   const int cz[NR_OF_DIRECTIONS] =  {
   // C   N   S   W   E   T   B  NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE TNE TNW TSE TSW BNE BNW BSE BSW
      0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1
   };

   /// The x,y,z component for each direction \ingroup stencil
   const int c[3][NR_OF_DIRECTIONS] = {
      {
   // C   N   S   W   E   T   B  NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE TNE TNW TSE TSW BNE BNW BSE BSW
      0,  0,  0, -1,  1,  0,  0, -1,  1, -1,  1,  0,  0, -1,  1,  0,  0, -1,  1,  1, -1,  1, -1,  1, -1,  1, -1
      }, {
   // C   N   S   W   E   T   B  NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE TNE TNW TSE TSW BNE BNW BSE BSW
      0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1, -1,  0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1
      }, {
   // C   N   S   W   E   T   B  NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE TNE TNW TSE TSW BNE BNW BSE BSW
      0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1,  1,  1, -1, -1, -1, -1,  1,  1,  1,  1, -1, -1, -1, -1
      }
   };

   /// The x,y,z component for each normalized direction \ingroup stencil
   const real_t cNorm[3][NR_OF_DIRECTIONS] = {
      {
         real_t(0), real_t(0), real_t(0), real_t(-1), real_t(1), real_t(0), real_t(0), real_t(-1) / std::sqrt( real_t(2) ),
         real_t(1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ), real_t(1) / std::sqrt( real_t(2) ), real_t(0), real_t(0),
         real_t(-1) / std::sqrt( real_t(2) ), real_t(1) / std::sqrt( real_t(2) ), real_t(0), real_t(0), real_t(-1) / std::sqrt( real_t(2) ),
         real_t(1) / std::sqrt( real_t(2) ), real_t(1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) ),
         real_t(1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) ), real_t(1) / std::sqrt( real_t(3) ),
         real_t(-1) / std::sqrt( real_t(3) ), real_t(1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) )
      }, {
         real_t(0), real_t(1), real_t(-1), real_t(0), real_t(0), real_t(0), real_t(0), real_t(1) / std::sqrt( real_t(2) ),
         real_t(1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ),
         real_t(1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ), real_t(0), real_t(0), real_t(1) / std::sqrt( real_t(2) ),
         real_t(-1) / std::sqrt( real_t(2) ), real_t(0), real_t(0), real_t(1) / std::sqrt( real_t(3) ), real_t(1) / std::sqrt( real_t(3) ),
         real_t(-1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) ), real_t(1) / std::sqrt( real_t(3) ),
         real_t(1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) )
      }, {
         real_t(0), real_t(0), real_t(0), real_t(0), real_t(0), real_t(1), real_t(-1), real_t(0), real_t(0), real_t(0), real_t(0),
         real_t(1) / std::sqrt( real_t(2) ), real_t(1) / std::sqrt( real_t(2) ), real_t(1) / std::sqrt( real_t(2) ),
         real_t(1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ),
         real_t(-1) / std::sqrt( real_t(2) ), real_t(-1) / std::sqrt( real_t(2) ), real_t(1) / std::sqrt( real_t(3) ),
         real_t(1) / std::sqrt( real_t(3) ), real_t(1) / std::sqrt( real_t(3) ), real_t(1) / std::sqrt( real_t(3) ),
         real_t(-1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) ), real_t(-1) / std::sqrt( real_t(3) ),
         real_t(-1) / std::sqrt( real_t(3) )
      }
   };

   /// String representation for each direction \ingroup stencil
   const std::string dirToString[NR_OF_DIRECTIONS] =  {
      "C", "N", "S", "W", "E", "T", "B",
      "NW", "NE", "SW", "SE", "TN", "TS", "TW", "TE", "BN", "BS", "BW","BE",
      "TNE", "TNW", "TSE", "TSW", "BNE", "BNW", "BSE", "BSW",
   };

   /// Binary encoded direction for each direction \ingroup stencil
   const BinaryDirection dirToBinary[27] = {
      Bin_C, Bin_N, Bin_S, Bin_W, Bin_E, Bin_T, Bin_B,
      Bin_NW, Bin_NE, Bin_SW, Bin_SE, Bin_TN, Bin_TS, Bin_TW, Bin_TE, Bin_BN, Bin_BS, Bin_BW, Bin_BE,
      Bin_TNE, Bin_TNW, Bin_TSE, Bin_TSW, Bin_BNE, Bin_BNW, Bin_BSE, Bin_BSW,
   };

   /// Inverse directions  \ingroup stencil
   const Direction inverseDir[NR_OF_DIRECTIONS] = {
      C, S, N, E, W, B, T,
      SE, SW, NE, NW, BS, BN, BE, BW, TS, TN, TE, TW,
      BSW, BSE, BNW, BNE, TSW, TSE, TNW, TNE
   };

   /// Length for each direction \ingroup stencil
   const real_t dirLength [NR_OF_DIRECTIONS] = {
        real_t(0), real_t(1), real_t(1), real_t(1), real_t(1), real_t(1), real_t(1),
        std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), 
        std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), std::sqrt( real_t(2) ),
        std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), std::sqrt( real_t(2) ), std::sqrt( real_t(2) ),
        std::sqrt( real_t(3) ), std::sqrt( real_t(3) ), std::sqrt( real_t(3) ), std::sqrt( real_t(3) ),
        std::sqrt( real_t(3) ), std::sqrt( real_t(3) ), std::sqrt( real_t(3) ), std::sqrt( real_t(3) )
   };


   const real_t gaussianWeights [NR_OF_DIRECTIONS] =
   {
      //C  N   S   W   E   T   B   NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE  TNE TNW TSE TSW BNE BNW BSE BSW
      real_t(8) / real_t(64),
      real_t(4) / real_t(64), real_t(4) / real_t(64), real_t(4) / real_t(64), real_t(4) / real_t(64),
      real_t(4) / real_t(64), real_t(4) / real_t(64),
      real_t(2) / real_t(64), real_t(2) / real_t(64), real_t(2) / real_t(64), real_t(2) / real_t(64),
      real_t(2) / real_t(64), real_t(2) / real_t(64), real_t(2) / real_t(64), real_t(2) / real_t(64), 
      real_t(2) / real_t(64), real_t(2) / real_t(64), real_t(2) / real_t(64), real_t(2) / real_t(64),
      real_t(1) / real_t(64), real_t(1) / real_t(64), real_t(1) / real_t(64), real_t(1) / real_t(64), 
      real_t(1) / real_t(64), real_t(1) / real_t(64), real_t(1) / real_t(64), real_t(1) / real_t(64)
   };


   const uint_t gaussianMultipliers [NR_OF_DIRECTIONS] =
   {
      //C  N   S   W   E   T   B   NW  NE  SW  SE  TN  TS  TW  TE  BN  BS  BW  BE  TNE TNW TSE TSW BNE BNW BSE BSW
      uint_t(8u),
      uint_t(4u), uint_t(4u), uint_t(4u), uint_t(4u) ,
      uint_t(4u), uint_t(4u),
      uint_t(2u), uint_t(2u), uint_t(2u), uint_t(2u) ,
      uint_t(2u), uint_t(2u), uint_t(2u), uint_t(2u) ,
      uint_t(2u), uint_t(2u), uint_t(2u), uint_t(2u) ,
      uint_t(1u), uint_t(1u), uint_t(1u), uint_t(1u) ,
      uint_t(1u), uint_t(1u), uint_t(1u), uint_t(1u)
   };


   /// The mirrored directions (flip W-E)  \ingroup stencil
   const Direction mirrorX[NR_OF_DIRECTIONS] = {
      C, N, S, E, W, T, B,
      NE, NW, SE, SW, TN, TS, TE, TW, BN, BS, BE, BW,
      TNW, TNE, TSW, TSE, BNW, BNE, BSW, BSE
   };

   /// The mirrored directions (flip N-S) \ingroup stencil
   const Direction mirrorY[NR_OF_DIRECTIONS] = {
      C, S, N, W, E, T, B,
      SW, SE, NW, NE, TS, TN, TW, TE, BS, BN, BW, BE,
      TSE, TSW, TNE, TNW, BSE, BSW, BNE, BNW
   };

   /// The mirrored directions (flip T-B) \ingroup stencil
   const Direction mirrorZ[NR_OF_DIRECTIONS] = {
      C, N, S, W, E, B, T,
      NW, NE, SW, SE, BN, BS, BW, BE, TN, TS, TW, TE,
      BNE, BNW, BSE, BSW, TNE, TNW, TSE, TSW
   };



   /// Maps from 2D directions (C, N, S, W, E, NW, NE, SW, SE) to 3D directions, by slicing through x,y or z coordinate \ingroup stencil
   /// The first array index represents the slice dimension ( 0 for x, 1 for y, 2 for z)
   /// Example: printing a slice through x coordinate (keeping x fixed) of a D3Q19 field:
   /** \code
    *   GhostLayerField<real_t,19,1> pdfField;
        for(auto i = D2Q9::begin(); i != D2Q9::end(); ++i)
            cout << pdfField.get(x,y,z, D3Q19::idx[ map2Dto3D[0][*i] ] ) << endl;
       \endcode
   */
   const Direction map2Dto3D[3][NR_OF_DIRECTIONS] =
   {
            { C, T, B, S, N, INVALID_DIR, INVALID_DIR, TS, TN, BS, BN,
              INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR,
              INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR },

            { C, T, B, W, E, INVALID_DIR, INVALID_DIR, TW, TE, BW, BE,
              INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR,
              INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR },

            { C, N, S, W, E, INVALID_DIR, INVALID_DIR, NW, NE, SW, SE,
              INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR,
              INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR, INVALID_DIR }
   };

   /// Maps (direction,axis) pair to direction
   /// \param axis     0,1 or 2 standing for x,y,z
   /// \param minOrMax if true, the direction pointing in the negative axis direction is returned,
   ///                 if false, the positive axis direction
   inline Direction directionFromAxis( int axis, bool minOrMax )
   {
      WALBERLA_ASSERT_LESS( axis, 3 );
      WALBERLA_ASSERT_GREATER_EQUAL( axis, 0 );
           if ( axis==0 &&  minOrMax ) return W;
      else if ( axis==0 && !minOrMax ) return E;
      else if ( axis==1 &&  minOrMax ) return S;
      else if ( axis==1 && !minOrMax ) return N;
      else if ( axis==2 &&  minOrMax ) return B;
      else if ( axis==2 && !minOrMax ) return T;
      else                             return INVALID_DIR;
   }

   /// Maps (direction,axis) pair to direction
   /// \param axis     0,1 or 2 standing for x,y,z
   /// \param minOrMax if true, the direction pointing in the negative axis direction is returned,
   ///                 if false, the positive axis direction
   inline Direction directionFromAxis( uint_t axis, bool minOrMax )
   {
      WALBERLA_ASSERT_LESS( axis, 3 );
           if ( axis==0 &&  minOrMax ) return W;
      else if ( axis==0 && !minOrMax ) return E;
      else if ( axis==1 &&  minOrMax ) return S;
      else if ( axis==1 && !minOrMax ) return N;
      else if ( axis==2 &&  minOrMax ) return B;
      else if ( axis==2 && !minOrMax ) return T;
      else                             return INVALID_DIR;
   }



   /// Computes neighbor from cell in direction d
   /// \param cell   origin cell
   /// \param d      direction pointing towards the computed neighbor
   inline Cell operator+( const Cell & cell, const Direction d )
   {
      return Cell( cell.x() + cx[d], cell.y() + cy[d], cell.z() + cz[d] );
   }



   /// Computes neighbor from cell in direction inverseDir[d]
   /// \param cell   origin cell
   /// \param d      direction pointing away from the computed neighbor
   inline Cell operator-( const Cell & cell, const Direction d )
   {
      return Cell( cell.x() - cx[d], cell.y() - cy[d], cell.z() - cz[d] );
   }



   /// Shifts cell to its neighbor in direction d
   /// \param cell   shifted cell
   /// \param d      direction in which to shift the cell
   inline Cell & operator+=( Cell & cell, const Direction d )
   {
      cell.x() += cx[d];
      cell.y() += cy[d];
      cell.z() += cz[d];

      return cell;
   }



   /// Shifts cell to its neighbor in direction inverseDir[d]
   /// \param cell   shifted cell
   /// \param d      direction opposite to which the cell gets shifted
   inline Cell & operator-=( Cell & cell, const Direction d )
   {
      cell.x() -= cx[d];
      cell.y() -= cy[d];
      cell.z() -= cz[d];

      return cell;
   }


} // namespace stencil
} // namespace walberla


