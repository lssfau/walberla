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
//! \file Timeloop.h
//! \ingroup timeloop
//! \author Martin Bauer <martin.bauer@fau.de>
//! \brief Header file for Timeloop
//
//======================================================================================================================

#pragma once

#include "ITimeloop.h"

#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/selectable/SetSelectableObject.h"
#include "core/timing/TimingPool.h"
#include "core/uid/SUID.h"

#include <functional>


namespace walberla {
namespace timeloop {

typedef std::function<void ()> VoidFctNoArguments;


//*******************************************************************************************************************
/*! Abstract base class for time loops.
*
* Supports registration of selectable functions that are run before/after a timestep.
* The function doTimeStep() runs the timestep itself and must be implemented by derived classes.
*
* \ingroup timeloop
*/
//*******************************************************************************************************************
class Timeloop : public ITimeloop
{
private:

   class LoggingStamp;
   friend class LoggingStamp;
   class LoggingStamp : public logging::Logging::CustomStamp
   {
   public:
      LoggingStamp( const Timeloop & timeloop ) : timeloop_( timeloop ) {}
      std::string stamp() override
      {
         std::ostringstream oss;
         int indention;

         if( timeloop_.nrOfTimeSteps_ > 0 )
            indention = int_c( std::ceil( std::log10( real_c( timeloop_.nrOfTimeSteps_ ) ) ) );
         else if( timeloop_.curTimeStep_ > 0 )
            indention = int_c( std::ceil( std::log10( real_c( timeloop_.curTimeStep_ ) ) ) );
         else
            indention = 0;

         oss << std::setw( indention )
             << std::setfill(' ') << std::right << timeloop_.curTimeStep_;
         return std::string("[") + oss.str() + std::string("]");
      }
      uint_t maxStampWidth() override
      {
         if( timeloop_.nrOfTimeSteps_ > 0 )
            return uint_c( std::ceil( std::log10( real_c( timeloop_.nrOfTimeSteps_ ) ) ) ) + uint_c(2);
         else if( timeloop_.curTimeStep_ > 0 )
            return uint_c( std::ceil( std::log10( real_c( timeloop_.curTimeStep_ ) ) ) ) + uint_c(2);
         else
            return uint_c(2);
      }
   private:
      const Timeloop & timeloop_;
   };

   class LoggingStampManager
   {
   public:
      LoggingStampManager( const shared_ptr< LoggingStamp > & stamp, const bool useCustomStamp )
         : useCustomStamp_( useCustomStamp )
      {
         if( useCustomStamp_ )
            logging::Logging::instance()->addCustomStamp( stamp );
      }
      ~LoggingStampManager()
      {
         if( useCustomStamp_ )
            logging::Logging::instance()->clearCustomStamp();
      }
   private:
      const bool useCustomStamp_;
   };

public:
   //**Construction & Destruction************************************************************************************
   /*! \name Construction & Destruction */
   //@{
   Timeloop( uint_t nrOfTimeSteps );

   ~Timeloop() override = default;
   //@}
   //****************************************************************************************************************


   //** Execution Control *******************************************************************************************
   /*! \name Execution Control*/
   //@{
   void run() override                  { run(true); }
   void run( const bool logTimeStep );
   void run( WcTimingPool & timing, const bool logTimeStep = true );

   void singleStep() override { singleStep(true); }
   void singleStep( const bool logTimeStep );
   void singleStep( WcTimingPool & timing, const bool logTimeStep = true );

   void stop() override;
   void synchronizedStop( bool stop ) override;

    void setCurrentTimeStepToZero()     { curTimeStep_ = 0;  }
    void setCurrentTimeStep( uint_t ts) override { curTimeStep_ = ts; }

    //@}
   //****************************************************************************************************************


   //** Registration Functions **************************************************************************************
   /*! \name Registration Functions */
   //@{
   typedef size_t FctHandle;


    FctHandle addFuncBeforeTimeStep(const VoidFctNoArguments & f,
                                    const std::string & identifier = std::string(),
                                    const Set<SUID> & require      = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatible = Set<SUID>::emptySet() );

    void      addFuncBeforeTimeStep(const FctHandle & fctToBindTo,
                                    const VoidFctNoArguments & f,
                                    const std::string & identifier = std::string(),
                                    const Set<SUID> & require      = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatible = Set<SUID>::emptySet() );


    FctHandle addFuncAfterTimeStep (const VoidFctNoArguments & f,
                                    const std::string & identifier     = std::string(),
                                    const Set<SUID> & require          = Set<SUID>::emptySet(),
                                    const Set<SUID> & exludingSelector = Set<SUID>::emptySet() );

    void      addFuncAfterTimeStep (const FctHandle & fctToBindTo,
                                    const VoidFctNoArguments & f,
                                    const std::string & identifier = std::string(),
                                    const Set<SUID> & require      = Set<SUID>::emptySet(),
                                    const Set<SUID> & incompatible = Set<SUID>::emptySet() );

   //@}
   //****************************************************************************************************************


   //** Timestep ****************************************************************************************************
   /*! \name Timestep */
   //@{
   uint_t getCurrentTimeStep() const override   { return curTimeStep_;   }
   uint_t getNrOfTimeSteps()   const override   { return nrOfTimeSteps_; }

   //@}
   //****************************************************************************************************************


protected:

   virtual void doTimeStep(const Set<SUID> &selectors) = 0;
   virtual void doTimeStep(const Set<SUID> &selectors, WcTimingPool &timing) = 0;


   void executeSelectable(const selectable::SetSelectableObject<VoidFctNoArguments,SUID> & selectable,
                          const Set<SUID> & selector,
                          const std::string & what );
   void executeSelectable(const selectable::SetSelectableObject<VoidFctNoArguments,SUID> & selectable,
                          const Set<SUID> & selector,
                          const std::string & what,
                          WcTimingPool & tp);


   uint_t curTimeStep_;   ///< current time step
   uint_t nrOfTimeSteps_; ///< total number of time steps

   typedef selectable::SetSelectableObject<VoidFctNoArguments, SUID> SelectableFunc;
   std::vector<SelectableFunc> beforeFunctions_;
   std::vector<SelectableFunc> afterFunctions_;

   bool stop_;
};



} // namespace timeloop
} // namespace walberla



//======================================================================================================================
//
//  EXPORT
//
//======================================================================================================================

namespace walberla {
   using timeloop::Timeloop;
}

