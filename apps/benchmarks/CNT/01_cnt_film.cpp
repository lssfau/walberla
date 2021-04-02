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
//! \file
//! \author Igor Ostanin <i.ostanin@skoltech.ru>
//! \author Grigorii Drozdov <drozd013@umn.edu>
//! \author Sebastian Eibl <sebastian.eibl@fau.de>
//
//======================================================================================================================

#include "Config.h"
#include "InitializeCNTs.h"
#include "SphericalSegmentAccessor.h"
#include "SQLProperties.h"
#include "Statistics.h"
#include "TerminalColors.h"

#include "mesa_pd/domain/BlockForestDomain.h"
#include "mesa_pd/data/Flags.h"
#include "mesa_pd/data/LinkedCells.h"
#include "mesa_pd/data/ParticleStorage.h"
#include "mesa_pd/data/SparseLinkedCells.h"
#include <mesa_pd/kernel/AssocToBlock.h>
#include "mesa_pd/kernel/InsertParticleIntoLinkedCells.h"
#include "mesa_pd/kernel/InsertParticleIntoSparseLinkedCells.h"
#include "mesa_pd/kernel/ParticleSelector.h"
#include "mesa_pd/kernel/cnt/AnisotropicVDWContact.h"
#include "mesa_pd/kernel/cnt/VBondContact.h"
#include "mesa_pd/kernel/cnt/ViscousDamping.h"
#include "mesa_pd/kernel/cnt/Parameters.h"
#include "mesa_pd/kernel/PFCDamping.h"
#include "mesa_pd/kernel/VelocityVerlet.h"
#include "mesa_pd/mpi/ContactFilter.h"
#include "mesa_pd/mpi/notifications/ForceTorqueNotification.h"
#include "mesa_pd/mpi/ReduceProperty.h"
#include "mesa_pd/mpi/SyncNextNeighbors.h"
#include "mesa_pd/sorting/LinearizedCompareFunctor.h"
#include "mesa_pd/vtk/ParticleVtkOutput.h"

#include "blockforest/Initialization.h"
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/timing/TimingPool.h"
#include "core/timing/TimingTree.h"
#include "sqlite/SQLite.h"
#include "vtk/VTKOutput.h"

namespace walberla {
namespace mesa_pd {

int main(int argc, char **argv)
{
   Environment env(argc, argv);
   walberla::mpi::MPIManager::instance()->useWorldComm();

   logging::Logging::instance()->setStreamLogLevel(logging::Logging::INFO);
   logging::Logging::instance()->setFileLogLevel(logging::Logging::INFO);

   WcTimingPool globalTP;
   WcTimingPool loopTP;
   WcTimingTree tt;

   WALBERLA_LOG_INFO_ON_ROOT("loading configuration parameters");
   globalTP["load config"].start();
   //---------------Generation of CNT configuration-----------------------------
   FilmSpecimen filmSpecimen; // Data structure for CNT film
   auto cfg = env.config();
   if (cfg == nullptr)
   WALBERLA_ABORT("No config specified!");
   const Config::BlockHandle cntConf = cfg->getBlock("CNT");
   loadFromConfig(filmSpecimen, cntConf);
   //===========================================================================
   globalTP["load config"].end();

   WALBERLA_LOG_INFO_ON_ROOT("creating domain setup");
   globalTP["domain setup"].start();
   Vector3<bool> periodicity = Vector3<bool>(true, true, false);
   if (filmSpecimen.oopp) periodicity[2] = true;
   real_t maxZ = filmSpecimen.sizeZ;
   real_t minZ = 0;
   if (!filmSpecimen.oopp)
   {
      maxZ += real_c(filmSpecimen.numSegs) * filmSpecimen.spacing;
      minZ -= real_c(filmSpecimen.numSegs) * filmSpecimen.spacing;
   }

   auto forest = blockforest::createBlockForest(
         AABB(0, 0, minZ, filmSpecimen.sizeX, filmSpecimen.sizeY, maxZ), // simulation domain
         Vector3<uint_t>(filmSpecimen.numBlocksX, filmSpecimen.numBlocksY, filmSpecimen.numBlocksZ), // blocks in each direction
         periodicity // periodicity
   );
   if (!forest) WALBERLA_ABORT("No BlockForest created ... exiting!");

   auto domain = std::make_shared<domain::BlockForestDomain> (forest);

   auto simulationDomain = forest->getDomain();
   WALBERLA_UNUSED(simulationDomain);
   auto localDomain = forest->begin()->getAABB();
   for (auto& blk : *forest)
   {
      localDomain.merge(blk.getAABB());
   }
   globalTP["domain setup"].end();

   WALBERLA_LOG_INFO_ON_ROOT("creating initial particle setup");
   globalTP["generate CNT"].start();
   auto ps = std::make_shared<data::ParticleStorage>(10);
   auto ac = SphericalSegmentAccessor(ps);
   auto lc = data::SparseLinkedCells(localDomain.getExtended(4_r * kernel::cnt::outer_radius), 4_r * kernel::cnt::outer_radius );

   int64_t numParticles = 0;
   if (filmSpecimen.initialConfigurationFile.empty())
   {
      numParticles = generateCNTs(filmSpecimen, ps, *domain);
   } else
   {
      numParticles = loadCNTs(filmSpecimen.initialConfigurationFile, ps);
   }
   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   WALBERLA_LOG_INFO_ON_ROOT("#particles created: " << numParticles);
   globalTP["generate CNT"].end();

   WALBERLA_LOG_INFO_ON_ROOT("setting up VTK output");
   auto vtkOutput = make_shared<mesa_pd::vtk::ParticleVtkOutput>(ps);
   vtkOutput->addOutput<data::SelectParticlePosition>("position");
   auto vtkWriter = walberla::vtk::createVTKOutput_PointData(vtkOutput,
                                                             "cnt",
                                                             1,
                                                             filmSpecimen.vtkFolder,
                                                             "simulation_step",
                                                             true,
                                                             true,
                                                             true,
                                                             filmSpecimen.useMPIIO);

   WALBERLA_LOG_INFO_ON_ROOT("setting up interaction models");
   kernel::InsertParticleIntoSparseLinkedCells ipilc;
   kernel::cnt::AnisotropicVDWContact anisotropicVdwContact;
   kernel::cnt::VBondContact vBondContact;
   kernel::cnt::ViscousDamping viscous_damping(filmSpecimen.viscousDamping * 1052.0_r,
                                               filmSpecimen.viscousDamping * 1052.0_r);
   kernel::PFCDamping pfcDamping(filmSpecimen.localDamping);
   kernel::VelocityVerletPreForceUpdate  vv_pre(kernel::cnt::dT);
   kernel::VelocityVerletPostForceUpdate vv_post(kernel::cnt::dT);
   mpi::ContactFilter                    contactFilter;
   mpi::ReduceProperty                   RP;
   mpi::SyncNextNeighbors                SNN;

   WALBERLA_LOG_INFO_ON_ROOT("warmup");
   SNN(*ps, *domain);
   sorting::LinearizedCompareFunctor linear(lc.domain_, Vector3<uint_t>(lc.numCellsPerDim_[0]));
   ps->sort(linear);
   vtkWriter->write(true);

   WALBERLA_LOG_INFO_ON_ROOT("running simulation");
   int64_t contactsChecked  = 0;
   int64_t contactsDetected = 0;
   int64_t contactsTreated  = 0;

   std::vector<Statistics> statistics;

   globalTP["sim loop"].start();
   for (auto i = 0; i < filmSpecimen.simulationSteps; ++i)
   {
      Statistics stats;

      loopTP["VV pre"].start();
      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_pre,
                          ac);
      loopTP["VV pre"].end();

      loopTP["SNN"].start();
      SNN(*ps, *domain);
      loopTP["SNN"].end();

      loopTP["GenerateSparseLinkedCells"].start();
      lc.clear();
      ps->forEachParticle(false, kernel::SelectAll(), ac, ipilc, ac, lc);
      loopTP["GenerateSparseLinkedCells"].end();

      loopTP["interaction"].start();
      contactsChecked  = 0;
      contactsDetected = 0;
      contactsTreated  = 0;
      lc.forEachParticlePairHalf(false,
                                 kernel::SelectAll(),
                                 ac,
                                 [&](size_t p_idx1, size_t p_idx2)
                                 {
                                    //ensure ordering
                                    if (ac.getUid(p_idx1) > ac.getUid(p_idx2))
                                       std::swap(ac.getUidRef(p_idx1), ac.getUidRef(p_idx1));
                                    
                                    ++contactsChecked;
                                    auto dist = ac.getPosition(p_idx1) - ac.getPosition(p_idx2);
                                    auto dist2 = math::dot(dist, dist);
                                    auto cutoff = ac.getInteractionRadius(p_idx1) + ac.getInteractionRadius(p_idx2);
                                    if (dist2 < cutoff * cutoff)
                                    {
                                       auto contactPoint = ac.getPosition(p_idx2) + 0.5_r * dist;
                                       ++contactsDetected;
                                       if (contactFilter(p_idx1, p_idx2, ac, contactPoint, *domain))
                                       {
                                          ++contactsTreated;
                                          constexpr bool period_X = false;
                                          constexpr bool period_Y = false;
                                          constexpr int max_segsX = 0;
                                          constexpr int max_segsY = 0;
                                          auto sameCluster = ac.getClusterID(p_idx2) == ac.getClusterID(p_idx1);
                                          auto segmentDistance = std::abs(ac.getSegmentID(p_idx2) - ac.getSegmentID(p_idx1));
                                          if ((sameCluster) && ((segmentDistance == 1) ||
                                                                ((period_X) && (segmentDistance == max_segsX - 1)) ||
                                                                ((period_Y) && (segmentDistance == max_segsY - 1))
                                          )
                                                )
                                          {
                                             vBondContact(p_idx1, p_idx2, ac);
                                             stats.bendingEnergy += vBondContact.getLastBendingEnergy();
                                             stats.shearEnergy += vBondContact.getLastShearEnergy();
                                             stats.tensileEnergy += vBondContact.getLastTensileEnergy();
                                             stats.twistingEnergy += vBondContact.getLastTwistingEnergy();
                                             stats.stretchEnergy += vBondContact.getLastBendingEnergy() +
                                                   vBondContact.getLastShearEnergy() +
                                                   vBondContact.getLastTensileEnergy() +
                                                   vBondContact.getLastTwistingEnergy();
                                          }
                                          else if ((sameCluster) && (segmentDistance < 20))
                                          {
                                             //do nothing
                                          } else
                                          {
                                             anisotropicVdwContact(p_idx1, p_idx2, ac);
                                             stats.vdwEnergy += anisotropicVdwContact.getLastEnergy();
                                             stats.numAlignedPairs += anisotropicVdwContact.isParallel() ? 1 : 0;
                                          }
                                          viscous_damping(p_idx1, p_idx2, ac);
                                       }
                                    }
                                 });
      loopTP["interaction"].end();

      loopTP["ReduceForce"].start();
      RP.operator()<ForceTorqueNotification>(*ps);
      loopTP["ReduceForce"].end();

      loopTP["PFDDamping"].start();
      ps->forEachParticle(false, kernel::SelectAll(), ac, pfcDamping, ac);
      loopTP["PFDDamping"].end();

      loopTP["VV post"].start();
      ps->forEachParticle(false,
                          kernel::SelectAll(),
                          ac,
                          vv_post,
                          ac);
      loopTP["VV post"].end();

      if ((i + 1) % filmSpecimen.saveVTKEveryNthStep == 0)
      {
         loopTP["save vtk"].start();
         WALBERLA_LOG_INFO_ON_ROOT(CYAN << "            SAVING VTK                  " << RESET);
         vtkWriter->write(true);
         loopTP["save vtk"].end();
      }
      if ((i + 1) % filmSpecimen.saveEnergyEveryNthStep == 0)
      {
         loopTP["save energy"].start();
         WALBERLA_LOG_INFO_ON_ROOT(CYAN << "            SAVING ENERGY                  " << RESET);
         ps->forEachParticle(false,
                             kernel::SelectLocal(),
                             ac,
                             [&](const size_t p_idx)
                             {
                                stats.kineticEnergy += 0.5_r * ac.getMass(p_idx) * ac.getLinearVelocity(p_idx) * ac.getLinearVelocity(p_idx);
                                stats.kineticEnergy += 0.5_r * ac.getAngularVelocity(p_idx) * ac.getInertiaBF(p_idx) * ac.getAngularVelocity(p_idx);
                             });
         walberla::mpi::reduceInplace(stats, walberla::mpi::SUM);
         WALBERLA_ROOT_SECTION()
         {
            statistics.push_back(stats);
            saveStatistics(filmSpecimen.energyFolder + "/statistics.txt", statistics);
         }
         loopTP["save energy"].end();
      }
      if ((i + 1) % filmSpecimen.saveConfEveryNthStep == 0)
      {
         loopTP["save config"].start();
         WALBERLA_LOG_INFO_ON_ROOT(CYAN << "            SAVING CONFIG                  " << RESET);
         saveConfig(ps,
                    filmSpecimen.confFolder + "/conf_" + std::to_string((i + 1) * 0.00002) + "ns",
                    IO_MODE::USE_MPI_IO);
         loopTP["save config"].end();
      }
   }
   globalTP["sim loop"].end();

   auto loopTPReduced = loopTP.getReduced();
   WALBERLA_LOG_DEVEL_ON_ROOT(*loopTPReduced);
   auto globalTPReduced = globalTP.getReduced();
   WALBERLA_LOG_DEVEL_ON_ROOT(*globalTPReduced);
   auto ttReduced = tt.getReduced();
   WALBERLA_LOG_DEVEL_ON_ROOT(ttReduced);

   numParticles = 0;
   int64_t numGhostParticles = 0;
   ps->forEachParticle(false,
                       kernel::SelectAll(),
                       ac,
                       [&](const size_t idx)
                       {
                          if (data::particle_flags::isSet( ac.getFlagsRef(idx), data::particle_flags::GHOST))
                          {
                             ++numGhostParticles;
                          } else
                          {
                             ++numParticles;
                          }
                       });

   walberla::mpi::reduceInplace(numParticles, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(numGhostParticles, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(contactsDetected, walberla::mpi::SUM);
   walberla::mpi::reduceInplace(contactsTreated, walberla::mpi::SUM);


   WALBERLA_ROOT_SECTION()
   {
      std::map<std::string, walberla::int64_t> integerProperties;
      std::map<std::string, double> realProperties;
      std::map<std::string, std::string> stringProperties;

      addBuildInfoToSQL(integerProperties, realProperties, stringProperties);
      addDomainPropertiesToSQL(*forest, integerProperties, realProperties, stringProperties);
      addSlurmPropertiesToSQL(integerProperties, realProperties, stringProperties);

      integerProperties["num_particles"] = numParticles;
      integerProperties["num_ghost_particles"] = numGhostParticles;

      integerProperties["contacts_detected"] = contactsDetected;
      integerProperties["contacts_treated"] = contactsTreated;

      auto runId = sqlite::storeRunInSqliteDB(filmSpecimen.sqlFile,
                                              integerProperties,
                                              stringProperties,
                                              realProperties);
      sqlite::storeTimingPoolInSqliteDB(filmSpecimen.sqlFile, runId, *loopTPReduced, "loop");
      sqlite::storeTimingPoolInSqliteDB(filmSpecimen.sqlFile, runId, *globalTPReduced, "global");
      sqlite::storeTimingTreeInSqliteDB(filmSpecimen.sqlFile, runId, tt, "tt");
   }

   return EXIT_SUCCESS;
}

} // namespace mesa_pd
} // namespace walberla

int main(int argc, char *argv[])
{
   return walberla::mesa_pd::main(argc, argv);
}

