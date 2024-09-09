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
//! \file ShiftedPeriodicity.h
//! \ingroup boundary
//! \author Helen Schottenhamml <helen.schottenhamml@fau.de>
//
//======================================================================================================================

#include "blockforest/Block.h"
#include "blockforest/BlockID.h"
#include "blockforest/StructuredBlockForest.h"

#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "core/math/AABBFwd.h"
#include "core/math/Vector3.h"
#include "core/mpi/MPIWrapper.h"
#include "core/mpi/MPIManager.h"

#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/IBlockID.h"
#include "domain_decomposition/MapPointToPeriodicDomain.h"

#include <array>
#include <cstdlib>
#include <iterator>
#include <map>
#include <memory>
#include <tuple>
#include <vector>


namespace walberla {
namespace boundary {

//*******************************************************************************************************************
/*!
 * A periodicity boundary condition that adds a user-defined spatial shift to the field when applied.
 * This shift can prevent the locking of large-scale turbulent features in the flow direction, see e.g.,
 * Munters et al. (https://doi.org/10.1063/1.4941912).
 *
 * Periodicity defined in the blockforest must be turned off in the normal-direction. Only uniform grids are supported
 * for the moment. Normal direction and shift direction must not coincide.
 * Moreover, this boundary condition is only applied at the minimum and maximum of the domain in a given direction, but
 * cannot be applied in the interior of the domain for now.
 *
 * This base class handles the logistics of splitting the field and communication.
 *
 * @tparam Field_T Type of the ghost-layer field that is shifted periodically
 */
//*******************************************************************************************************************
template<typename Derived_T, typename Field_T>
class ShiftedPeriodicityBase {

 public:

   using FieldType = Field_T;
   using ValueType = typename FieldType::value_type;
   using ShiftType = int;

   ShiftedPeriodicityBase( const std::weak_ptr<StructuredBlockForest> & blockForest,
                           const BlockDataID & fieldID, const uint_t fieldGhostLayers,
                           const uint_t normalDir, const uint_t shiftDir, const ShiftType shiftValue )
      : blockForest_(blockForest), normalDir_(normalDir), shiftDir_(shiftDir),
        fieldID_( fieldID ), fieldGhostLayers_(fieldGhostLayers)
   {
      auto sbf = blockForest_.lock();

      WALBERLA_ASSERT_NOT_NULLPTR( sbf )

      WALBERLA_CHECK(sbf->storesUniformBlockGrid(), "Periodic shift is currently only implemented for uniform grids.")
      WALBERLA_CHECK(sbf->containsGlobalBlockInformation(), "For the periodic shift, the blockforest must be constructed to retain global information.")

      shift_ = Vector3<ShiftType>{};
      auto adaptedShiftValue = shiftValue % ShiftType(sbf->getNumberOfCells(shiftDir_));
      shift_[shiftDir_] = ShiftType(adaptedShiftValue >= 0 ? adaptedShiftValue : adaptedShiftValue + ShiftType(sbf->getNumberOfCells(shiftDir_)));

      // sanity checks
      WALBERLA_CHECK_UNEQUAL( shiftDir_, normalDir_, "Direction of periodic shift and boundary normal must not coincide." )

      WALBERLA_CHECK( sbf->isPeriodic(shiftDir_), "Blockforest must be periodic in direction " << shiftDir_ << "!" )
      WALBERLA_CHECK( !sbf->isPeriodic(normalDir_), "Blockforest must NOT be periodic in direction " << normalDir_ << "!" )

   }

   Vector3<ShiftType> shift() const { return shift_; }
   uint_t shiftDirection() const { return shiftDir_; }
   uint_t normalDirection() const { return normalDir_; }

   void operator()() {

      auto mpiInstance = mpi::MPIManager::instance();
      const auto currentRank = numeric_cast<mpi::MPIRank>(mpiInstance->rank());

      const auto sbf = blockForest_.lock();
      WALBERLA_ASSERT_NOT_NULLPTR( sbf )

      // only setup send and receive information once at the beginning
      if(!setupPeriodicity_){
         setupPeriodicity();
         setupPeriodicity_ = true;
      }

      // set up local to avoid temporary fields; key is unique tag to ensure correct communication
      // if several messages have the same source and destination
      // thanks to the unique tag, the same buffer can be used for both local and MPI communication
      std::map<int, std::vector<ValueType>> buffer;
      std::map<BlockID, std::map<int, MPI_Request>> sendRequests;
      std::map<BlockID, std::map<int, MPI_Request>> recvRequests;

      std::vector<IBlock*> localBlocks{};
      sbf->getBlocks(localBlocks);

      /// packing and schedule sends

      for( auto block : localBlocks ) {

         const auto blockID = dynamic_cast<Block*>(block)->getId();

         WALBERLA_ASSERT_EQUAL(sendAABBs_[blockID].size(), recvAABBs_[blockID].size() )

         for(uint_t i = 0; i < sendAABBs_[blockID].size(); ++i) {

            auto & sendInfo = sendAABBs_[blockID][i];

            const auto sendRank = std::get<0>(sendInfo);
            const auto sendTag = std::get<3>(sendInfo);

            const auto sendCI = std::get<2>(sendInfo);

            const auto sendSize = sendCI.numCells() * uint_c(fSize_);

            buffer[sendTag].resize(sendSize);
            static_cast<Derived_T*>(this)->packBuffer(block, sendCI, buffer[sendTag]);

            WALBERLA_MPI_SECTION() {
               // schedule sends if MPI
               if (sendRank != currentRank)
               {
                  MPI_Isend(buffer[sendTag].data(), mpi::MPISize(buffer[sendTag].size() * sizeof(ValueType)), MPI_BYTE,
                            sendRank, sendTag, mpiInstance->comm(), &sendRequests[blockID][sendTag]);
               }
            }

         }

      }

      /// schedule receives

      for( auto block : localBlocks ) {

         const auto blockID = dynamic_cast<Block*>(block)->getId();

         WALBERLA_ASSERT_EQUAL(sendAABBs_[blockID].size(), recvAABBs_[blockID].size() )

         for(uint_t i = 0; i < recvAABBs_[blockID].size(); ++i) {

            auto & recvInfo = recvAABBs_[blockID][i];

            const auto recvRank = std::get<0>(recvInfo);
            const auto recvTag   = std::get<3>(recvInfo);
            const auto recvCI = std::get<2>(recvInfo);

            WALBERLA_MPI_SECTION() {
               // do not schedule receives for local communication
               if (recvRank != currentRank) {
                  const auto recvSize = recvCI.numCells() * uint_c(fSize_);
                  buffer[recvTag].resize(recvSize);

                  // Schedule receives
                  MPI_Irecv(buffer[recvTag].data(), mpi::MPISize(buffer[recvTag].size() * sizeof(ValueType)), MPI_BYTE,
                            recvRank, recvTag, mpiInstance->comm(), &recvRequests[blockID][recvTag]);
               }
            }

         }

      }

      /// unpacking

      for( auto block : localBlocks ) {

         const auto blockID = dynamic_cast<Block*>(block)->getId();

         for(uint_t i = 0; i < recvAABBs_[blockID].size(); ++i) {

            auto & recvInfo = recvAABBs_[blockID][i];

            const auto recvRank = std::get<0>(recvInfo);
            const auto recvTag = std::get<3>(recvInfo);

            const auto recvCI = std::get<2>(recvInfo);

            if(recvRank == currentRank) {
               WALBERLA_ASSERT_GREATER(buffer.count(recvTag), 0)
               static_cast<Derived_T*>(this)->unpackBuffer(block, recvCI, buffer[recvTag]);
            } else {
               WALBERLA_MPI_SECTION() {
                  MPI_Status status;
                  MPI_Wait(&recvRequests[blockID][recvTag], &status);

                  WALBERLA_ASSERT_GREATER(buffer.count(recvTag), 0)
                  static_cast< Derived_T* >(this)->unpackBuffer(block, recvCI, buffer[recvTag]);
               }
            }

         }

      }

      WALBERLA_MPI_SECTION() {
         // wait for all communication to be finished
         for (auto& block : localBlocks)
         {
            const auto blockID = dynamic_cast< Block* >(block)->getId();

            if (sendRequests[blockID].empty()) return;

            std::vector< MPI_Request > v;
            std::transform(sendRequests[blockID].begin(), sendRequests[blockID].end(), std::back_inserter(v),
                           second(sendRequests[blockID]));
            MPI_Waitall(int_c(sendRequests[blockID].size()), v.data(), MPI_STATUSES_IGNORE);
         }
      }

   }

 private:

   /*!
    * Creates send and receive information (source, target rank, source / target blocks, source and destination CI
    * (here in form of an AABB), and unique communication tag) per block if the shift does not coincide with the
    * block size.
    *
    * If the shift does not coincide with the block size, each block is split into two send sub-blocks and two receive
    * sub-blocks whose source / target block will differ (maybe also the rank).
    * Afterwards, all the information necessary for the communication is determined and saved in the corresponding
    * data structure.
    */
   void processDoubleAABB( const AABB & originalAABB, const std::shared_ptr<StructuredBlockForest> & sbf,
                          const real_t nGL, const BlockID & blockID, const int normalShift ) {

      WALBERLA_ASSERT(normalShift == -1 || normalShift == 1)

      const auto localShift = normalShift == 1 ? shift_ : -shift_;

      const auto offset = ShiftType(int_c(shift_[shiftDir_]) % int_c(sbf->getNumberOfCellsPerBlock(shiftDir_)));
      Vector3<ShiftType> normalShiftVector{};
      normalShiftVector[normalDir_] = ShiftType(normalShift * int_c(sbf->getNumberOfCellsPerBlock(normalDir_)));

      const auto remDir = uint_t(3) - normalDir_ - shiftDir_;

      AABB aabb1 = originalAABB;
      AABB aabb2 = originalAABB;

      if( normalShift == 1 ) {
         aabb1.setAxisBounds(shiftDir_, aabb1.min(shiftDir_), aabb1.max(shiftDir_) - real_c(offset));
         aabb2.setAxisBounds(shiftDir_, aabb2.max(shiftDir_) - real_c(offset), aabb2.max(shiftDir_));
      } else {
         aabb1.setAxisBounds(shiftDir_, aabb1.min(shiftDir_), aabb1.min(shiftDir_) + real_c(offset));
         aabb2.setAxisBounds(shiftDir_, aabb2.min(shiftDir_) + real_c(offset), aabb2.max(shiftDir_));
      }

      auto center1 = aabb1.center();
      auto center2 = aabb2.center();

      // account for ghost layers
      aabb1.setAxisBounds(remDir, aabb1.min(remDir) - nGL, aabb1.max(remDir) + nGL);
      aabb2.setAxisBounds(remDir, aabb2.min(remDir) - nGL, aabb2.max(remDir) + nGL);

      auto localSendAABB1 = aabb1;
      auto localRecvAABB1 = aabb1;
      auto localSendAABB2 = aabb2;
      auto localRecvAABB2 = aabb2;

      localSendAABB1.setAxisBounds(shiftDir_, localSendAABB1.min(shiftDir_), localSendAABB1.max(shiftDir_) + nGL);
      localSendAABB2.setAxisBounds(shiftDir_, localSendAABB2.min(shiftDir_) - nGL, localSendAABB2.max(shiftDir_));
      localRecvAABB1.setAxisBounds(shiftDir_, localRecvAABB1.min(shiftDir_) - nGL, localRecvAABB1.max(shiftDir_));
      localRecvAABB2.setAxisBounds(shiftDir_, localRecvAABB2.min(shiftDir_), localRecvAABB2.max(shiftDir_) + nGL);

      if(normalShift == 1) { // at maximum of domain -> send inner layers, receive ghost layers
         localSendAABB1.setAxisBounds(normalDir_, localSendAABB1.max(normalDir_) - nGL, localSendAABB1.max(normalDir_));
         localSendAABB2.setAxisBounds(normalDir_, localSendAABB2.max(normalDir_) - nGL, localSendAABB2.max(normalDir_));
         localRecvAABB1.setAxisBounds(normalDir_, localRecvAABB1.max(normalDir_), localRecvAABB1.max(normalDir_) + nGL);
         localRecvAABB2.setAxisBounds(normalDir_, localRecvAABB2.max(normalDir_), localRecvAABB2.max(normalDir_) + nGL);
      } else { // at minimum of domain -> send inner layers, receive ghost layers
         localSendAABB1.setAxisBounds(normalDir_, localSendAABB1.min(normalDir_), localSendAABB1.min(normalDir_) + nGL);
         localSendAABB2.setAxisBounds(normalDir_, localSendAABB2.min(normalDir_), localSendAABB2.min(normalDir_) + nGL);
         localRecvAABB1.setAxisBounds(normalDir_, localRecvAABB1.min(normalDir_) - nGL, localRecvAABB1.min(normalDir_));
         localRecvAABB2.setAxisBounds(normalDir_, localRecvAABB2.min(normalDir_) - nGL, localRecvAABB2.min(normalDir_));
      }

      WALBERLA_ASSERT_FLOAT_EQUAL(localSendAABB1.volume(), localRecvAABB1.volume())
      WALBERLA_ASSERT_FLOAT_EQUAL(localSendAABB2.volume(), localRecvAABB2.volume())

      center1 += (localShift + normalShiftVector);
      center2 += (localShift + normalShiftVector);

      std::array<bool, 3> virtualPeriodicity{false};
      virtualPeriodicity[normalDir_] = true;
      virtualPeriodicity[shiftDir_] = true;

      domain_decomposition::mapPointToPeriodicDomain( virtualPeriodicity, sbf->getDomain(), center1 );
      domain_decomposition::mapPointToPeriodicDomain( virtualPeriodicity, sbf->getDomain(), center2 );

      uint_t rank1{};
      uint_t rank2{};

      BlockID id1{};
      BlockID id2{};

      const auto blockInformation = sbf->getBlockInformation();
      WALBERLA_ASSERT(blockInformation.active())

      blockInformation.getProcess(rank1, center1[0], center1[1], center1[2]);
      blockInformation.getProcess(rank2, center2[0], center2[1], center2[2]);

      blockInformation.getId(id1, center1[0], center1[1], center1[2]);
      blockInformation.getId(id2, center2[0], center2[1], center2[2]);

      // define unique send / receive tags for communication

      const int atMaxTagSend = normalShift < 0 ? 0 : 1;
      const int atMaxTagRecv = normalShift < 0 ? 1 : 0;

      const int sendTag1 = ((int_c(blockID.getID()) + int_c(id1.getID() * numGlobalBlocks_)) * 2 + atMaxTagSend) * 2 + 0;
      const int sendTag2 = ((int_c(blockID.getID()) + int_c(id2.getID() * numGlobalBlocks_)) * 2 + atMaxTagSend) * 2 + 1;
      const int recvTag2 = ((int_c(id2.getID()) + int_c(blockID.getID() * numGlobalBlocks_)) * 2 + atMaxTagRecv) * 2 + 0;
      const int recvTag1 = ((int_c(id1.getID()) + int_c(blockID.getID() * numGlobalBlocks_)) * 2 + atMaxTagRecv) * 2 + 1;

      auto block = sbf->getBlock(blockID);

      sendAABBs_[blockID].emplace_back(mpi::MPIRank(rank1), id1, globalAABBToLocalCI(localSendAABB1, sbf, block), sendTag1);
      sendAABBs_[blockID].emplace_back(mpi::MPIRank(rank2), id2, globalAABBToLocalCI(localSendAABB2, sbf, block), sendTag2);
      recvAABBs_[blockID].emplace_back(mpi::MPIRank(rank2), id2, globalAABBToLocalCI(localRecvAABB2, sbf, block), recvTag2);
      recvAABBs_[blockID].emplace_back(mpi::MPIRank(rank1), id1, globalAABBToLocalCI(localRecvAABB1, sbf, block), recvTag1);

      WALBERLA_LOG_DETAIL("blockID = " << blockID.getID() << ", normalShift = " << normalShift
                                       << "\n\tsendRank1 = " << rank1 << "\tsendID1 = " << id1.getID() << "\tsendTag1 = " << sendTag1 << "\taabb = " << localSendAABB1
                                       << "\n\tsendRank2 = " << rank2 << "\tsendID2 = " << id2.getID() << "\tsendTag2 = " << sendTag2 << "\taabb = " << localSendAABB2
                                       << "\n\trecvRank1 = " << rank1 << "\trecvID1 = " << id1.getID() << "\trecvTag1 = " << recvTag1 << "\taabb = " << localRecvAABB1
                                       << "\n\trecvRank2 = " << rank2 << "\trecvID2 = " << id2.getID() << "\trecvTag2 = " << recvTag2 << "\taabb = " << localRecvAABB2
      )
   }

   /*!
    * Does the same things as `processDoubleAABB` but assuming that the shift is a multiple of the block size, i.e.,
    * every field slice is shifted one or several entire blocks. In this case, every block only has ONE send and ONE
    * receive block whose information must be stored.
    */
   void processSingleAABB( const AABB & originalAABB, const std::shared_ptr<StructuredBlockForest> & sbf,
                          const real_t nGL, const BlockID & blockID, const int normalShift ) {

      WALBERLA_ASSERT(normalShift == -1 || normalShift == 1)

      const auto localShift = normalShift == 1 ? shift_ : -shift_;

      Vector3<ShiftType> normalShiftVector{};
      normalShiftVector[normalDir_] = ShiftType(normalShift * int_c(sbf->getNumberOfCellsPerBlock(normalDir_)));

      // aabb
      auto aabb = originalAABB.getExtended(nGL);
      auto center = aabb.center();

      // receive from the interior of domain in ghost layers
      auto localSendAABB = aabb;
      auto localRecvAABB = aabb;

      if(normalShift == 1) { // at maximum of domain -> send inner layers, receive ghost layers
         localSendAABB.setAxisBounds(normalDir_, localSendAABB.max(normalDir_) - 2 * nGL, localSendAABB.max(normalDir_) - nGL);
         localRecvAABB.setAxisBounds(normalDir_, localRecvAABB.max(normalDir_) - nGL, localRecvAABB.max(normalDir_));
      } else { // at minimum of domain -> send inner layers, receive ghost layers
         localSendAABB.setAxisBounds(normalDir_, localSendAABB.min(normalDir_) + nGL, localSendAABB.min(normalDir_) + 2 * nGL);
         localRecvAABB.setAxisBounds(normalDir_, localRecvAABB.min(normalDir_), localRecvAABB.min(normalDir_) + nGL);
      }

      WALBERLA_ASSERT_FLOAT_EQUAL(localSendAABB.volume(), localRecvAABB.volume())

      center += (localShift + normalShiftVector);

      std::array<bool, 3> virtualPeriodicity{false};
      virtualPeriodicity[normalDir_] = true;
      virtualPeriodicity[shiftDir_] = true;

      domain_decomposition::mapPointToPeriodicDomain( virtualPeriodicity, sbf->getDomain(), center );

      uint_t rank{};

      BlockID id{};

      const auto blockInformation = sbf->getBlockInformation();
      WALBERLA_ASSERT(blockInformation.active())

      blockInformation.getProcess(rank, center[0], center[1], center[2]);

      blockInformation.getId(id, center[0], center[1], center[2]);

      // define unique send / receive tags for communication

      const int atMaxTagSend = normalShift < 0 ? 0 : 1;
      const int atMaxTagRecv = normalShift < 0 ? 1 : 0;

      const int sendTag = ((int_c(blockID.getID()) + int_c(id.getID() * numGlobalBlocks_)) * 2 + atMaxTagSend) * 2 + 0;
      const int recvTag = ((int_c(id.getID()) + int_c(blockID.getID() * numGlobalBlocks_)) * 2 + atMaxTagRecv) * 2 + 0;

      auto block = sbf->getBlock(blockID);

      sendAABBs_[blockID].emplace_back(mpi::MPIRank(rank), id, globalAABBToLocalCI(localSendAABB, sbf, block), sendTag);
      recvAABBs_[blockID].emplace_back(mpi::MPIRank(rank), id, globalAABBToLocalCI(localRecvAABB, sbf, block), recvTag);

      WALBERLA_LOG_DETAIL("blockID = " << blockID.getID() << ", normalShift = " << normalShift
                                       << "\n\tsendRank = " << rank << "\tsendID = " << id.getID() << "\tsendTag = " << sendTag << "\taabb = " << localSendAABB
                                       << "\n\trecvRank = " << rank << "\trecvID = " << id.getID() << "\trecvTag = " << recvTag << "\taabb = " << localRecvAABB
      )
   }

   /*!
    * Determines which blocks are participating in the boundary condition / communication, and calls the appropriate
    * functions to extract and save the communication information.
    */
   void setupPeriodicity() {

      const auto sbf = blockForest_.lock();
      WALBERLA_ASSERT_NOT_NULLPTR( sbf )

      auto & blockInformation = sbf->getBlockInformation();
      WALBERLA_ASSERT(blockInformation.active())
      std::vector<std::shared_ptr<IBlockID>> globalBlocks;
      blockInformation.getAllBlocks(globalBlocks);
      numGlobalBlocks_ = globalBlocks.size();

      // get blocks of current processor
      std::vector<IBlock*> localBlocks;
      sbf->getBlocks(localBlocks);

      const auto nGL = real_c(fieldGhostLayers_);

      const bool shiftWholeBlock = (shift_[shiftDir_] % ShiftType(sbf->getNumberOfCells(*localBlocks[0], shiftDir_))) == 0;

      // get f-size
      {
         auto* field = localBlocks[0]->getData<FieldType>(fieldID_);
         WALBERLA_ASSERT_NOT_NULLPTR(field)
         fSize_ = cell_idx_c(field->fSize());
      }

      for( auto block : localBlocks ) {

         // get minimal ghost layer region (in normal direction)
         const auto blockAABB = block->getAABB();
         const auto blockID = dynamic_cast<Block*>(block)->getId();

         const bool atMin = sbf->atDomainMinBorder(normalDir_, *block);
         const bool atMax = sbf->atDomainMaxBorder(normalDir_, *block);

         if(atMin) {
            if(shiftWholeBlock) {
               processSingleAABB(blockAABB, sbf, nGL, blockID, -1);
            } else {
               processDoubleAABB(blockAABB, sbf, nGL, blockID, -1);
            }
         }
         if(atMax)  {
            if(shiftWholeBlock) {
               processSingleAABB(blockAABB, sbf, nGL, blockID, +1);
            } else {
               processDoubleAABB(blockAABB, sbf, nGL, blockID, +1);
            }
         }

      }

   }

   Vector3<cell_idx_t> toCellVector( const Vector3<real_t> & point ) {
      return Vector3<cell_idx_t >{ cell_idx_c(point[0]), cell_idx_c(point[1]), cell_idx_c(point[2]) };
   }

   CellInterval globalAABBToLocalCI( const AABB & aabb, const std::shared_ptr<StructuredBlockForest> & sbf, IBlock * const block ) {
      auto globalCI = CellInterval{toCellVector(aabb.min()), toCellVector(aabb.max()) - Vector3<cell_idx_t>(1, 1, 1)};
      CellInterval localCI;
      sbf->transformGlobalToBlockLocalCellInterval(localCI, *block, globalCI);

      return localCI;
   }

   template< typename tPair >
   struct second_t {
      typename tPair::second_type operator()( const tPair& p ) const { return p.second; }
   };
   template< typename tMap >
   second_t< typename tMap::value_type > second( const tMap& ) { return second_t< typename tMap::value_type >(); }


   std::weak_ptr<StructuredBlockForest> blockForest_;

   uint_t normalDir_;
   uint_t shiftDir_;
   Vector3<ShiftType> shift_;

   // for each local block, stores the ranks where to send / receive, the corresponding block IDs,
   // the local AABBs that need to be packed / unpacked, and a unique tag for communication
   std::map<BlockID, std::vector<std::tuple<mpi::MPIRank, BlockID, CellInterval, int>>> sendAABBs_{};
   std::map<BlockID, std::vector<std::tuple<mpi::MPIRank, BlockID, CellInterval, int>>> recvAABBs_{};

   bool setupPeriodicity_{false};
   uint_t numGlobalBlocks_{};

 protected:

   const BlockDataID fieldID_;

   const uint_t fieldGhostLayers_;
   cell_idx_t fSize_;

}; // class ShiftedPeriodicityBase


//*******************************************************************************************************************
/*!
 * A periodicity boundary condition that adds a user-defined spatial shift to the field when applied.
 * This shift can prevent the locking of large-scale turbulent features in the flow direction, see e.g.,
 * Munters et al. (https://doi.org/10.1063/1.4941912).
 *
 * Periodicity defined in the blockforest must be turned off in the normal-direction.
 *
 * This class handles the CPU-specific packing and unpacking of the communication buffers.
 *
 * @tparam GhostLayerField_T Type of the ghost-layer field that is shifted periodically
 */
//*******************************************************************************************************************
template<typename GhostLayerField_T>
class ShiftedPeriodicity : public ShiftedPeriodicityBase<ShiftedPeriodicity<GhostLayerField_T>, GhostLayerField_T> {

   using Base = ShiftedPeriodicityBase<ShiftedPeriodicity<GhostLayerField_T>, GhostLayerField_T>;
   friend Base;

 public:

   using ValueType = typename GhostLayerField_T::value_type;
   using ShiftType = typename Base::ShiftType;

   ShiftedPeriodicity( const std::weak_ptr<StructuredBlockForest> & blockForest,
                       const BlockDataID & fieldID, const uint_t fieldGhostLayers,
                      const uint_t normalDir, const uint_t shiftDir, const ShiftType shiftValue )
      : Base( blockForest, fieldID, fieldGhostLayers, normalDir, shiftDir, shiftValue )
   {}

 private:

   void packBuffer(IBlock * const block, const CellInterval & ci,
                   std::vector<ValueType> & buffer) {

      // get field
      auto field = block->getData<GhostLayerField_T>(this->fieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(field)

      auto bufferIt = buffer.begin();

      // forward iterate over ci and add values to value vector
      for(auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt) {
         for(cell_idx_t f = 0; f < this->fSize_; ++f, ++bufferIt) {
            WALBERLA_ASSERT(field->coordinatesValid(cellIt->x(), cellIt->y(), cellIt->z(), f))
            *bufferIt = field->get(*cellIt, f);
         }
      }

   }

   void unpackBuffer(IBlock* const block, const CellInterval & ci,
                     const std::vector<ValueType> & buffer) {

      // get field
      auto field = block->getData<GhostLayerField_T>(this->fieldID_);
      WALBERLA_ASSERT_NOT_NULLPTR(field)

      auto bufferIt = buffer.begin();

      // forward iterate over ci and add values to value vector
      for(auto cellIt = ci.begin(); cellIt != ci.end(); ++cellIt) {
         for(cell_idx_t f = 0; f < this->fSize_; ++f, ++bufferIt) {
            WALBERLA_ASSERT(field->coordinatesValid(cellIt->x(), cellIt->y(), cellIt->z(), f))
            field->get(*cellIt, f) = *bufferIt;
         }
      }

   }

}; // class ShiftedPeriodicity

} // namespace lbm
} // namespace walberla
