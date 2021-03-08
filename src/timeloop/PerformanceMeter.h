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
//! \file PerformanceMeter.h
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/timing/Timer.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>


namespace walberla {
namespace timeloop {


   /********************************************************************************************************************
    * \brief Class for measuring performance characteristics like for example MFLUPS
    *
    * \ingroup timeloop
    *
    * This class measures performance numbers related to a cell count. For example the number of fluid cell updates
    * per second, or the number of interface cell updates per second.
    * This cell number can either be the total number of cells in the domain or can be counted using a custom
    * function ( every n'th timestep or just once )
    *
      There are two usage modes:
      1) The class can be linked together with a timeloop using getBeforeFunction() and getAfterFunction()
         Everything that is executed between these functions is measured.
          \code
          PerformanceMeter perfMeter ( blockStorage );
          timeloop.addFuncBeforeTimeStep( perfMeter.getBeforeFunction() );
          timeloop.addFuncAfterTimeStep ( perfMeter.getAfterFunction() );
          perfMeter.addMeasurement ( "MLUPS",  1e-6 ); // uses total number of cells ( no counting )
          // counts the number of cells using custom count function, counting is done every 20th step
          perfMeter.addMeasurement ( "MFLUPS", myCountFunction, 1e-6, 20);
          //(...)
          perfMeter.logResultOnRoot();
          \endcode

      2) Measure everything between two calls of operator()
         In this case the PerformanceMeter should be added as first entry to the timeloop:
          \code
          PerformanceMeter perfMeter ( blockStorage );
          timeloop.addFuncBeforeTimeStep ( perfMeter );
          // ... (same as above)
          \endcode
         If the operator() is called N times, (N-1) timing measurements are done.
    *******************************************************************************************************************/
   class PerformanceMeter
   {
   public:
      using CountFunction = std::function<uint_t (const IBlock &)>;


      PerformanceMeter( StructuredBlockStorage & blockStorage );


      //** Timeloop Wiring  ********************************************************************************************
      /*! \name Timeloop Wiring */
      //@{

      void operator() ();
      void timingStart();
      void timingEnd();

      std::function<void () > getBeforeFunction();
      std::function<void () > getAfterFunction();

      void clear();
      //@}
      //****************************************************************************************************************


      //** Measurement Definitions *************************************************************************************
      /*! \name Measurement Definitions */
      //@{
      void addMeasurement( const std::string & name, const CountFunction& countFunction,
                           uint_t countFreq = 0, real_t scaling = 1 );

      void addMeasurement( const std::string & name, real_t scaling = 1 );

      template<typename FField>
      void addMeasurementBasedOnFlagField( const std::string & name, ConstBlockDataID flagFieldID,
                                           typename FField::flag_t activeMask,
                                           uint_t countFreq = 0, real_t scaling = 1 );
      //@}
      //****************************************************************************************************************


      //** Output ******************************************************************************************************
      /*! \name Output */
      //@{
      void logResultOnRoot( );
      void print( std::ostream & os, int targetRank=0 );

      shared_ptr< std::map<std::string, real_t> >  getReduced ( int targetRank = 0 );
      //@}
      //****************************************************************************************************************

   private:

      /// This struct holds a measurement definition
      struct Measurement
      {
         Measurement ( CountFunction cf, const std::string & n, real_t sc, uint_t countFreq )
            : countFunction ( cf), name(n), scaling(sc), countingFreq(countFreq), counts(0), avgCellsPerTimeStep(0)
         {}

         CountFunction countFunction;
         std::string   name;                 //!< name of the measured performance number
         real_t        scaling;              //!< scaling factor, usually power of ten, for (M)FLUPS
         uint_t        countingFreq;         //!< every n'th timestep cells are counted, if 0 cells are counted only once
         uint_t        counts;               //!< how often was the number of cells counted until now
         real_t        avgCellsPerTimeStep;  //!< the average number of cells per timestep
      };

      void   timingStartStop();
      void   reduce ( std::vector<real_t> & totalCells, int targetRank );

      void   updateCellCounts();

      WcTimer timer_;

      StructuredBlockStorage & blockStorage_;

      std::vector<Measurement> measurements_;

      /// Used for operator() to indicate if called the first time
      bool firstTimingStartStopCall_;
   };










   template<typename FlagField_T>
   uint_t flagFieldCountFunction( const IBlock & block, ConstBlockDataID flagFieldID,
                                  typename FlagField_T::flag_t mask )
   {
      uint_t cellCounter = 0;

      const FlagField_T * ffield = block.template getData<FlagField_T>( flagFieldID );
      for( auto cell = ffield->begin(); cell != ffield->end(); ++cell )
         if ( isPartOfMaskSet( cell,mask ) )
            ++cellCounter;

      return cellCounter;
   }

   template<typename FField>
   void PerformanceMeter::addMeasurementBasedOnFlagField( const std::string & name, ConstBlockDataID flagFieldID,
                                                          typename FField::flag_t activeMask,
                                                          uint_t countFreq , real_t scaling )
   {
      this->addMeasurement( name, std::bind( flagFieldCountFunction<FField>,  std::placeholders::_1, flagFieldID, activeMask ),
                           countFreq, scaling );
   }



} // namespace timeloop
} // namespace walberla


