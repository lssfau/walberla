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
//! \file Evaluation.cpp
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================
# pragma once

#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/Filesystem.h"
#include "core/math/Sample.h"

#include "lbm_generated/field/PdfField.h"

#include "Types.h"
#include "Setup.h"
#include "FlowAroundSphereInfoHeader.h"

#include <deque>
#include <fstream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace walberla;

using FlagField_T          = FlagField< uint8_t >;
using BoundaryCollection_T = lbm::FlowAroundSphereBoundaryCollection< FlagField_T >;

namespace walberla
{

struct DragCoefficient {
   uint_t timestep;
   double Fx;
   double Fy;
   double Fz;
   double cDRealArea;
   double cLRealArea;
   double cDDiscreteArea;
   double cLDiscreteArea;
   DragCoefficient(uint_t t, Vector3<double> f, double cdR, double clR, double cdD, double clD) : timestep(t), Fx(f[0]), Fy(f[1]), Fz(f[2]), cDRealArea(cdR), cLRealArea(clR), cDDiscreteArea(cdD), cLDiscreteArea(clD) {}
};

class Evaluation
{
 public:
   Evaluation(std::shared_ptr< StructuredBlockForest >& blocks, const uint_t checkFrequency, const uint_t rampUpTime,
              BoundaryCollection_T & boundaryCollection,
              const IDs& ids, const Setup& setup,
              const bool logToStream = true, const bool logToFile = true)
      : blocks_(blocks), checkFrequency_(checkFrequency), rampUpTime_(rampUpTime),
        boundaryCollection_(boundaryCollection), ids_(ids), setup_(setup),
        logToStream_(logToStream), logToFile_(logToFile)
   {
      WALBERLA_ASSERT_NOT_NULLPTR(blocks)
      maxLevel_ = blocks->getDepth();
      blocks->getBlocks(finestBlocks_, maxLevel_);

      coefficients_.resize(uint_c(4));
      coefficientExtremas_.resize(uint_c(4));

      const double factor = setup_.dxC / setup_.dxF;
      diameterSphere      = double_c(2.0) * (double_c(setup_.sphereRadius) * factor);
      surfaceAreaSphere   = math::pi * diameterSphere * diameterSphere;
      meanVelocity        = setup_.inflowVelocity; // (double_c(4.0) * setup_.inflowVelocity) / double_c(9.0);

      dragFilename_ = std::string("dragSphereRe_") + std::to_string(uint_t(setup_.reynoldsNumber)) + std::string("_meshLevels_") +  std::to_string(uint_t(setup_.refinementLevels + 1)) + std::string(".csv");

      WALBERLA_ROOT_SECTION(){
         if (logToFile_){
            filesystem::path dataFile( dragFilename_.c_str() );
            if( filesystem::exists( dataFile ) )
               std::remove( dataFile.string().c_str() );

            std::ofstream outfile( dragFilename_.c_str() );
            outfile << "SEP=," << "\n";

            setup_.writeParameterHeader(outfile);

            outfile << "timestep," << "Fx," << "Fy," << "Fz," << "cDRealArea," << "cLRealArea," << "cDDiscreteArea," << "cLDiscreteArea" << "\n";
            outfile.close();
         }
      }

      refresh();

      WALBERLA_LOG_INFO_ON_ROOT(
         "Evaluation initialised:"
            "\n   + Sphere real area:                      " << surfaceAreaSphere
         << "\n   + Sphere real diameter:                  " << diameterSphere
         << "\n   + Sphere discrete area drag coefficient: " << AD_
         << "\n   + Sphere discrete area lift coefficient: " << AL_
      )

   }

   void operator()();
   void forceCalculation(const uint_t level); // for calculating the force
   void resetForce();

   std::function<void (const uint_t)> forceCalculationFunctor()
   {
      return [this](uint_t level) { forceCalculation(level); };
   }

   std::function<void()> resetForceFunctor()
   {
      return [this]() { resetForce(); };
   }

   void refresh();

 protected:
   bool initialized_{false};

   std::shared_ptr< StructuredBlockForest > blocks_;
   uint_t maxLevel_;
   std::vector<Block *> finestBlocks_;

   uint_t coarseExecutionCounter_{ uint_c(0) };
   uint_t fineExecutionCounter_{ uint_c(0) };
   uint_t checkFrequency_;
   uint_t rampUpTime_;

   BoundaryCollection_T & boundaryCollection_;

   IDs ids_;
   Setup setup_;

   double diameterSphere;
   double surfaceAreaSphere;
   double meanVelocity;

   Vector3< double > force_;
   double AD_;
   double AL_;
   std::vector<DragCoefficient> dragResults;

   std::vector< std::deque< double > > coefficients_;
   std::vector< std::pair< double, double > > coefficientExtremas_;

   bool logToStream_;
   bool logToFile_;
   std::string dragFilename_;

}; // class Evaluation

}