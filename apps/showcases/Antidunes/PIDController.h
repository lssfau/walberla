//======================================================================================================================
/*!
 *  \file   PIDController.h
 */
//======================================================================================================================

#pragma once

#include <core/DataTypes.h>
#include <core/logging/Logging.h>

using namespace walberla;
using walberla::real_t;

class PIDController
{
 public:
   PIDController();

   PIDController(const real_t commandVariable, const real_t initialActuatingVariable, const real_t proportionalGain,
                 const real_t derivateGain, const real_t integralGain, const real_t maxRamp,
                 const real_t minActuatingVariable, const real_t maxActuatingVariable);

   PIDController(const real_t commandVariable, const real_t initialActuatingVariable, const real_t proportionalGain,
                 const real_t derivateGain, const real_t integralGain);

   real_t update(const real_t controlledVariable);

   real_t getProportionalGain() const { return proportionalGain_; }
   real_t getDerivateGain() const { return derivateGain_; }
   real_t getIntegralGain() const { return integralGain_; }

   void writeStateToFile(std::string filename) const;

   void readStateFromFile(std::string filename);

 private:
   real_t commandVariable_;
   real_t actuatingVariable_;

   real_t proportionalGain_;
   real_t derivateGain_;
   real_t integralGain_;

   real_t maxRamp_;
   real_t minActuatingVariable_;
   real_t maxActuatingVariable_;

   real_t errorHistory_[3];
   real_t errorIntegral_;
};
