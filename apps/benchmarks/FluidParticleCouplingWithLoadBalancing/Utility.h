#pragma once

#include "lbm_mesapd_coupling/amr/BlockInfo.h"

namespace walberla {
namespace lbm_mesapd_coupling  {
namespace amr {

/*
 * Result from the workload evaluation as described in
 *  Rettinger, Ruede - "Dynamic Load Balancing Techniques for Particulate Flow Simulations", 2019, Computation
 */
real_t fittedLBMWeightEvaluationFunction(const BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   real_t weight = real_t(7.597476065046571e-06) * real_c(Ce) + real_t(8.95723566283202e-05) * real_c(F) + real_t(-0.1526111388616016);
   return std::max(weight,real_t(0));
}
real_t fittedBHWeightEvaluationFunction(const BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t NB = blockInfo.numberOfNearBoundaryCells;
   real_t weight = real_t(1.3067711379655123e-07) * real_c(Ce) + real_t(0.0007289549127205142) * real_c(NB) + real_t(-0.1575698071795788);
   return std::max(weight,real_t(0));
}
real_t fittedRPDWeightEvaluationFunction(const BlockInfo& blockInfo)
{
   uint_t Pl = blockInfo.numberOfLocalParticles;
   uint_t Pg = blockInfo.numberOfGhostParticles;
   uint_t Sc = blockInfo.numberOfRPDSubCycles;
   real_t cPlPg2 = real_t(2.402288635599054e-05);
   real_t cPl    = real_t(0.00040932622363097144);
   real_t cPg    = real_t(0.0007268941363125683);
   real_t c      = real_t(2.01883028312316e-19);
   real_t weight = real_c(Sc) * ( cPlPg2 * real_c(Pl+Pg) * real_c(Pl+Pg) + cPl * real_c(Pl) + cPg * real_c(Pg) + c );
   return std::max(weight,real_t(0));
}
real_t fittedCoup1WeightEvaluationFunction(const BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   uint_t Pl = blockInfo.numberOfLocalParticles;
   uint_t Pg = blockInfo.numberOfGhostParticles;
   real_t weight = real_t(5.610203409278647e-06) * real_c(Ce) + real_t(-7.257635845636656e-07) * real_c(F) + real_t(0.02049703546054693) * real_c(Pl) + real_t(0.04248208493809902) * real_c(Pg) + real_t(-0.26609470510074784);
   return std::max(weight,real_t(0));
}
real_t fittedCoup2WeightEvaluationFunction(const BlockInfo& blockInfo)
{
   uint_t Ce = blockInfo.numberOfCells;
   uint_t F  = blockInfo.numberOfFluidCells;
   uint_t Pl = blockInfo.numberOfLocalParticles;
   uint_t Pg = blockInfo.numberOfGhostParticles;
   real_t weight = real_t(7.198479654682179e-06) * real_c(Ce) + real_t(1.178247475854302e-06) * real_c(F) + real_t(-0.0026401549115124628) * real_c(Pl) + real_t(0.008459646786179298) * real_c(Pg) + real_t(-0.001077320113275954);
   return std::max(weight,real_t(0));
}
real_t fittedTotalWeightEvaluationFunction(const BlockInfo& blockInfo)
{
   return fittedLBMWeightEvaluationFunction(blockInfo) + fittedBHWeightEvaluationFunction(blockInfo) +
          fittedRPDWeightEvaluationFunction(blockInfo) + fittedCoup1WeightEvaluationFunction(blockInfo) +
          fittedCoup2WeightEvaluationFunction(blockInfo);
}

} //namespace amr
} //namespace lbm_mesapd_coupling
} //namespace walberla

