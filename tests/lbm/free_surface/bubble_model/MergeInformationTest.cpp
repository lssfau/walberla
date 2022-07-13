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
//! \file MergeInformationTest.cpp
//! \ingroup lbm/free_surface/bubble_model
//! \author Martin Bauer
//! \author Christoph Schwarzmeier <christoph.schwarzmeier@fau.de>
//! \brief Test bubble merge and reordering, merge registering, and merge communication.
//
//======================================================================================================================

#include "lbm/free_surface/bubble_model/MergeInformation.h"

#include "core/Environment.h"
#include "core/debug/TestSubsystem.h"
#include "core/logging/Logging.h"
#include "core/mpi/MPIManager.h"

#include <iostream>
#include <vector>

namespace std
{
// overload std::cout for printing vector content
template< typename T >
std::ostream& operator<<(std::ostream& os, const std::vector< T >& vec)
{
   os << "Vector (" << vec.size() << "):  ";
   for (auto i = vec.begin(); i != vec.end(); ++i)
   {
      os << *i << " ";
   }
   os << std::endl;

   return os;
}
} // namespace std

namespace walberla
{
namespace free_surface
{
namespace bubble_model
{
// test if renameVec_ and vecCompare are equal
void checkRenameVec(const bubble_model::MergeInformation& mi, const std::vector< bubble_model::BubbleID >& vecCompare)
{
   // renameVec_ and vecCompare must have the same size
   WALBERLA_ASSERT_EQUAL(mi.renameVec_.size(), vecCompare.size());

   bool equal = true;
   for (size_t i = 0; i < mi.renameVec_.size(); ++i)
      if (mi.renameVec_[i] != vecCompare[i])
      {
         equal = false;
         break;
      }

   WALBERLA_CRITICAL_SECTION_START
   if (!equal)
   {
      WALBERLA_LOG_WARNING("Process " << MPIManager::instance()->rank()
                                      << ": rename vector different than expected result.\n"
                                      << "\tExpected:\n"
                                      << "\t\t" << vecCompare << "\tResult:\n"
                                      << "\t\t" << mi.renameVec_);
   }
   WALBERLA_CRITICAL_SECTION_END

   WALBERLA_CHECK(equal);
}
} // namespace bubble_model

namespace MergeInformationTest
{
void mergeAndReorderTest1()
{
   std::vector< bubble_model::Bubble > bubbles;
   bubbles.emplace_back(real_c(10));                       // BubbleID 0
   bubbles.emplace_back(bubble_model::Bubble(real_c(11))); // BubbleID 1
   bubbles.emplace_back(bubble_model::Bubble(real_c(12))); // BubbleID 2
   bubbles.emplace_back(bubble_model::Bubble(real_c(13))); // BubbleID 3

   bubble_model::MergeInformation mi(bubbles.size());
   mi.registerMerge(0, 1);

   // merge bubbles with IDs 0 and 1 to ID 0:
   // - ID 4 gets assigned ID 1 (highest ID, i.e., last position is copied to free space position 1)
   // - highest ID is now 3 and bubbles vector must have reduced to size 2
   mi.mergeAndReorderBubbleVector(bubbles);

   WALBERLA_CHECK_EQUAL(bubbles.size(), 3);
   WALBERLA_CHECK_FLOAT_EQUAL(bubbles[0].getInitVolume(), real_c(10 + 11));
   WALBERLA_CHECK_FLOAT_EQUAL(bubbles[1].getInitVolume(), real_c(13));
   WALBERLA_CHECK_FLOAT_EQUAL(bubbles[2].getInitVolume(), real_c(12));
}

void mergeAndReorderTest2()
{
   std::vector< bubble_model::Bubble > bubbles;
   bubbles.emplace_back(bubble_model::Bubble(real_c(10))); // BubbleID 0
   bubbles.emplace_back(bubble_model::Bubble(real_c(11))); // BubbleID 1
   bubbles.emplace_back(bubble_model::Bubble(real_c(12))); // BubbleID 2
   bubbles.emplace_back(bubble_model::Bubble(real_c(13))); // BubbleID 3
   bubbles.emplace_back(bubble_model::Bubble(real_c(14))); // BubbleID 4
   bubbles.emplace_back(real_c(15));                       // BubbleID 5

   bubble_model::MergeInformation mi(bubbles.size());
   mi.registerMerge(4, 3);
   mi.registerMerge(3, 2);
   mi.registerMerge(0, 1);

   // merge bubbles with IDs 4, 3, 2 to ID 2, merge bubbles with IDs 0 and 1 to ID 0:
   // - ID 5 gets assigned ID 1 (highest ID, i.e., last position is copied to free space position 1)
   // - highest ID is now 3 and bubbles vector must have reduced to size 2
   mi.mergeAndReorderBubbleVector(bubbles);

   WALBERLA_CHECK_EQUAL(bubbles.size(), 3);
   WALBERLA_CHECK_FLOAT_EQUAL(bubbles[0].getInitVolume(), real_c(10 + 11));
   WALBERLA_CHECK_FLOAT_EQUAL(bubbles[1].getInitVolume(), real_c(15));
   WALBERLA_CHECK_FLOAT_EQUAL(bubbles[2].getInitVolume(), real_c(12 + 13 + 14));
}

void mergeRegisterTest()
{
   // create new MergeInformation with 6 bubbles
   bubble_model::MergeInformation mi(6);

   // initial ID distribution
   std::vector< bubble_model::BubbleID > correctRenameVec = { 0, 1, 2, 3, 4, 5 };
   checkRenameVec(mi, correctRenameVec);

   // Function registerMerge() only registers the merge by (temporarily) renaming bubble IDs appropriately. Therefore,
   // in below comments, it is more meaningful to write "position x (in renameVec_) is renamed to the ID at position
   // y" than "ID x is renamed to ID y".

   // position 5 is renamed to the ID at position 3
   mi.registerMerge(5, 3);
   correctRenameVec = { 0, 1, 2, 3, 4, 3 };
   checkRenameVec(mi, correctRenameVec);

   // position 3 is renamed to the ID at position 1
   mi.registerMerge(3, 1);
   correctRenameVec = { 0, 1, 2, 1, 4, 3 };
   checkRenameVec(mi, correctRenameVec);

   // since position 3 was already renamed to the ID at position 1, registerMerge(1, 0) is called internally and
   // position 1 is renamed to the ID at position 0; position 3 is renamed to the ID at position 0
   mi.registerMerge(3, 0);
   correctRenameVec = { 0, 0, 2, 0, 4, 3 };
   checkRenameVec(mi, correctRenameVec);

   // since position 5 was already renamed to the ID of position 3, registerMerge(3, 2) is called internally; since
   // position 3 was already renamed to the ID of position 0, registerMerge(2, 0) is called internally; position 2 is
   // renamed to the ID at position 0, position 5 is renamed to the ID at position 0
   mi.registerMerge(5, 2);
   correctRenameVec = { 0, 0, 0, 0, 4, 0 };
   checkRenameVec(mi, correctRenameVec);
}

void mergeCommunicationTest()
{
   auto mpiManager = MPIManager::instance();

   // this test is only meaningful with multiple processes
   if (!mpiManager->isMPIInitialized() || mpiManager->numProcesses() < 3) { return; }

   // create new MergeInformation with 5 bubbles and renameVec_={ 0, 1, 2, 3, 4 }
   bubble_model::MergeInformation mi(5);

   if (mpiManager->rank() == 0) { mi.registerMerge(1, 0); }
   else
   {
      if (mpiManager->rank() == 1)
      {
         mi.registerMerge(3, 2);
         mi.registerMerge(2, 1);
      }
      else
      {
         if (mpiManager->rank() == 2) { mi.registerMerge(4, 1); }
      }
   }

   // before communication:
   // process 0: renameVec_={ 0, 0, 2, 3, 4 }
   // process 1: renameVec_={ 0, 1, 1, 1, 4 }
   // process 2: renameVec_={ 0, 1, 2, 3, 1 }

   mi.communicateMerges();

   // after communication:
   // process 0 to 2: renameVec_={ 0, 0, 0, 0, 0 }
   std::vector< bubble_model::BubbleID > res = { 0, 0, 0, 0, 0 };
   checkRenameVec(mi, res);
}

int main(int argc, char** argv)
{
   debug::enterTestMode();
   Environment env(argc, argv);

   // the MPI communicator is normally created with the block forest; we will not use a block forest in this test and
   // thus choose the MPI communicator manually
   WALBERLA_MPI_SECTION() { MPIManager::instance()->useWorldComm(); }

   mergeAndReorderTest1();
   mergeAndReorderTest2();

   mergeRegisterTest();
   mergeCommunicationTest();

   return EXIT_SUCCESS;
}
} // namespace MergeInformationTest
} // namespace free_surface
} // namespace walberla

int main(int argc, char** argv) { return walberla::free_surface::MergeInformationTest::main(argc, argv); }
