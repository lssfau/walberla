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
//! \file PoiseuilleInitializer.h
//! \ingroup geometry
//! \author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once

#include "core/config/Config.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include "geometry/initializer/Initializer.h"

namespace walberla {
namespace lbm {
namespace initializer {


   //*******************************************************************************************************************
   /*! Initializes the complete domain with a Poiseuille (Channel) Flow
   *
   * \ingroup geometry
   *
   * Configuration file syntax:
   *    <blockName>
   *     {
   *       scenario          rect2D|pipe;
   *       boundary          pressureDriven|velocityDriven|forceDriven
   *
   *       pressureDiff 0.1; // either pressureDiff or velocity is required
   *       velocity 1.0      // either one can be calculated from the other
   *
   *       flowAxis 0;       // specifying the direction of the channel, i.e. the axis along which the fluid is flowing
   *       parabolaAxis 1;   // only required for rect2D setup. Axis where no-slip boundaries are set -
   *                         // the remaining third axis has to be periodic. By default chosen as a non-periodic axis
   *                         // that is not the flowAxis
   *     }
   *
   *  scenario:
   *     - rect2D: rectangular channel where boundary in direction of periodicAxis are set periodic
   *               and the remaining boundary is set to no-slip
   *     - pipe:   channel with circular cross section
   *  boundary:
   *     - pressureDriven:  at the end of flowAxis pressure boundaries are used
   *     - velocityDriven:  at the inflow a UBB prescribing a parabolic velocity profile is set
   *                        at the outflow a Outlet boundary is set
   *     - forceDriven:     flowAxis has to be periodic. No boundaries are set, but an external force term is set
   *                        in the Lattice model
   *
   *  flowAxis:
   *     specifying the direction of the channel, i.e. the axis along which the fluid is flowing
   *     0 for x, 1 for y, 2 = z
   */
   //*******************************************************************************************************************
   template<typename BoundaryHandling_T, typename LM,
            typename SimplePressure = lbm::SimplePressure <LM, typename BoundaryHandling_T::flag_t>,
            typename UBB            = lbm::UBB            <LM, typename BoundaryHandling_T::flag_t>  >
   class Poiseuille : public geometry::initializer::Initializer
   {
   public:
      enum Scenario       {  RECT_2D, PIPE  };
      enum BoundaryType   {  PRESSURE_DRIVEN, VELOCITY_DRIVEN, FORCE_DRIVEN };

      /*
      enum Axis {
         X_AXIS=0,
         Y_AXIS=1,
         Z_AXIS=2,
         INVALID_AXIS=3
      };*/
      typedef uint_t Axis;
      static const uint_t X_AXIS;
      static const uint_t Y_AXIS;
      static const uint_t Z_AXIS;
      static const uint_t INVALID_AXIS;

      Poiseuille( StructuredBlockStorage & blocks, BlockDataID & handlerID, BlockDataID  & pdfFieldID,
                  field::FlagUID noSlipFlag, field::FlagUID ubbFlag,
                  field::FlagUID pressureFlag1, field::FlagUID pressureFlag2 );


      void init( BlockStorage & , const Config::BlockHandle & blockHandle ) override    {  init( blockHandle );   }

      void init( const Config::BlockHandle & blockHandle );
      void init( Scenario scenario, BoundaryType boundaryType, real_t pressureDiff, Axis flowAxis, Axis parabolaAxis = INVALID_AXIS );

   protected:
      void initVelocityBoundary( Scenario scenario, Axis flowAxis, real_t maxVelocity, Axis parabolaAxis );
      void initPressureBoundary( Axis flowAxis, real_t pressureDiff );
      void initNoSlipBoundaries( Scenario scenario, Axis flowAxis, Axis parabolaAxis );

      void setInitVelocity( Scenario scenario, Axis flowAxis, real_t maxVelocity, Axis parabolaAxis );

      real_t getVelocity  ( const Cell & globalCell, Scenario scenario, Axis flowAxis, real_t maxVelocity, Axis parabolaAxis );
      real_t getPipeRadius( Scenario scenario, Axis flowAxis, Axis parabolaAxis ) const;

      real_t velocityFromPressureDiff( Scenario scenario, real_t pressureDiff, Axis flowAxis, Axis parabolaAxis );
      real_t pressureDiffFromVelocity( Scenario scenario, real_t velocity, Axis flowAxis, Axis parabolaAxis );

      Axis checkParabolaAxis( Axis parabolaAxis, Axis flowAxis );

      StructuredBlockStorage & storage_;
      BlockDataID handlerID_;
      BlockDataID pdfFieldID_;

      field::FlagUID noSlipFlag_;
      field::FlagUID ubbFlag_;
      field::FlagUID pressureFlag1_;
      field::FlagUID pressureFlag2_;

      real_t latticeViscosity_;

      Vector3<real_t> midPoint_;    // point in top-north-west corner
      Vector3<real_t> maxPoint_;    // midpoint of domain
   };


} // namespace initializer
} // namespace lbm
} // namespace walberla


#include "PoiseuilleInitializer.impl.h"


