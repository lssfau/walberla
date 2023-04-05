//======================================================================================================================
/*!
 *  \file   PIDController.cpp
 */
//======================================================================================================================

#include "PIDController.h"

#include <algorithm>
#include <fstream>

using namespace walberla;
using walberla::real_t;

PIDController::PIDController()
   : commandVariable_(0), actuatingVariable_(0), proportionalGain_(0), derivateGain_(0), integralGain_(0), maxRamp_(0),
     minActuatingVariable_(0), maxActuatingVariable_(0), errorIntegral_(0)
{
   std::fill(errorHistory_, errorHistory_ + sizeof(errorHistory_) / sizeof(real_t), real_t(0));
}

PIDController::PIDController(const real_t commandVariable, const real_t initialActuatingVariable,
                             const real_t proportionalGain, const real_t derivateGain, const real_t integralGain,
                             const real_t maxRamp, const real_t minActuatingVariable, const real_t maxActuatingVariable)
   : commandVariable_(commandVariable), actuatingVariable_(initialActuatingVariable),
     proportionalGain_(proportionalGain), derivateGain_(derivateGain), integralGain_(integralGain), maxRamp_(maxRamp),
     minActuatingVariable_(minActuatingVariable), maxActuatingVariable_(maxActuatingVariable), errorIntegral_(0)
{
   std::fill(errorHistory_, errorHistory_ + sizeof(errorHistory_) / sizeof(real_t), real_t(0));

   if (integralGain_ > real_t(0))
      errorIntegral_ = initialActuatingVariable / integralGain_;
   else
      errorIntegral_ = real_t(0);
}

PIDController::PIDController(const real_t commandVariable, const real_t initialActuatingVariable,
                             const real_t proportionalGain, const real_t derivateGain, const real_t integralGain)
   : commandVariable_(commandVariable), actuatingVariable_(initialActuatingVariable),
     proportionalGain_(proportionalGain), derivateGain_(derivateGain), integralGain_(integralGain),
     maxRamp_(std::numeric_limits< real_t >::max()), minActuatingVariable_(-std::numeric_limits< real_t >::max()),
     maxActuatingVariable_(std::numeric_limits< real_t >::max()), errorIntegral_(0)
{
   std::fill(errorHistory_, errorHistory_ + sizeof(errorHistory_) / sizeof(real_t), real_t(0));

   if (integralGain_ > real_t(0))
      errorIntegral_ = initialActuatingVariable / integralGain_;
   else
      errorIntegral_ = real_t(0);
}

real_t PIDController::update(const real_t controlledVariable)
{
   static const real_t ONE_OVER_SIX = real_t(1) / real_t(6);
   const real_t error               = commandVariable_ - controlledVariable;

   const real_t d =
      (error + real_t(3) * errorHistory_[0] - real_t(3) * errorHistory_[1] - errorHistory_[2]) * ONE_OVER_SIX;
   std::rotate(errorHistory_, errorHistory_ + 1, errorHistory_ + sizeof(errorHistory_) / sizeof(real_t));
   errorHistory_[sizeof(errorHistory_) / sizeof(real_t) - size_t(1)] = error;

   real_t newActuationVariable = proportionalGain_ * error + derivateGain_ * d + integralGain_ * errorIntegral_;

   if (std::fabs(actuatingVariable_ - newActuationVariable) < maxRamp_)
   {
      errorIntegral_ += error;
      newActuationVariable = proportionalGain_ * error + derivateGain_ * d + integralGain_ * errorIntegral_;
   }

   const real_t maxValue = std::min(actuatingVariable_ + maxRamp_, maxActuatingVariable_);
   const real_t minValue = std::max(actuatingVariable_ - maxRamp_, minActuatingVariable_);

   actuatingVariable_ = std::min(std::max(minValue, newActuationVariable), maxValue);

   return actuatingVariable_;
}

void PIDController::writeStateToFile(std::string filename) const
{
   std::ofstream file;
   file.open(filename);
   file << std::setprecision(16);
   file << commandVariable_ << std::endl;
   file << actuatingVariable_ << std::endl;
   file << proportionalGain_ << std::endl;
   file << derivateGain_ << std::endl;
   file << integralGain_ << std::endl;
   file << maxRamp_ << std::endl;
   file << minActuatingVariable_ << std::endl;
   file << maxActuatingVariable_ << std::endl;
   file << errorHistory_[0] << std::endl;
   file << errorHistory_[1] << std::endl;
   file << errorHistory_[2] << std::endl;
   file << errorIntegral_ << std::endl;
   file.close();
}

void PIDController::readStateFromFile(std::string filename)
{
   std::ifstream file;
   file.open(filename);
   file >> commandVariable_;
   file >> actuatingVariable_;
   file >> proportionalGain_;
   file >> derivateGain_;
   file >> integralGain_;
   file >> maxRamp_;
   file >> minActuatingVariable_;
   file >> maxActuatingVariable_;
   file >> errorHistory_[0];
   file >> errorHistory_[1];
   file >> errorHistory_[2];
   file >> errorIntegral_;
   file.close();
}
