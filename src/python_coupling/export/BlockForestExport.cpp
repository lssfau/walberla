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
//! \file Exports.cpp
//! \ingroup blockforest
//! \author Martin Bauer <martin.bauer@fau.de>
//! \author Markus Holzer <markus.holzer@fau.de>
//
//======================================================================================================================

// Do not reorder includes - the include order is important
#include "python_coupling/PythonWrapper.h"

#ifdef WALBERLA_BUILD_WITH_PYTHON

#   include "blockforest/Initialization.h"
#   include "blockforest/SetupBlock.h"
#   include "blockforest/SetupBlockForest.h"
#   include "blockforest/StructuredBlockForest.h"

#   include "core/StringUtility.h"
#   include "core/mpi/MPIIO.h"

#   include "stencil/D3Q19.h"

#   include <memory>
#   include <pybind11/stl.h>
#   include <sstream>

#   include "BlockForestExport.h"
#   include "python_coupling/helper/OwningIterator.h"

namespace walberla
{
namespace blockforest
{
std::string printSetupBlock(const SetupBlock& b)
{
   std::stringstream out;
   out << "SetupBlock at " << b.getAABB();
   return out.str();
}

namespace py = pybind11;

//======================================================================================================================
//
//  StructuredBlockForest
//
//======================================================================================================================

#ifdef WALBERLA_BUILD_WITH_PYTHON

void NoSuchBlockData::translate(  const NoSuchBlockData & e ) {
   throw py::cast_error(e.what());
}

void BlockDataNotConvertible::translate(  const BlockDataNotConvertible & e ) {
   throw py::cast_error(e.what());
}
#else

void NoSuchBlockData::translate(  const NoSuchBlockData &  ) {}

void BlockDataNotConvertible::translate(  const BlockDataNotConvertible &  ) {}

#endif

BlockDataID blockDataIDFromString(BlockStorage& bs, const std::string& stringID)
{
   auto ids = bs.getBlockDataIdentifiers();

   for (uint_t i = 0; i < ids.size(); ++i)
      if (ids[i] == stringID) return BlockDataID(i);

   throw NoSuchBlockData();
}

BlockDataID blockDataIDFromString(IBlock& block, const std::string& stringID)
{
   return blockDataIDFromString(block.getBlockStorage(), stringID);
}

BlockDataID blockDataIDFromString(StructuredBlockForest& bs, const std::string& stringID)
{
   return blockDataIDFromString(bs.getBlockStorage(), stringID);
}

py::iterator StructuredBlockForest_iter(const shared_ptr< StructuredBlockForest >& bf) // NOLINT
{
   // shared_ptr<StructuredBlockForest> s = py::cast< shared_ptr<StructuredBlockForest> > ( StructuredBlockForest );

   std::vector< const IBlock* > blocks;
   bf->getBlocks(blocks);
   std::vector<py::object> resultList;
   resultList.reserve(blocks.size());

   for (auto it = blocks.begin(); it != blocks.end(); ++it)
   {
      py::object theObject = py::cast(*it);
      resultList.push_back(theObject);
   }

   return python_coupling::make_owning_iterator(resultList);
}

py::object StructuredBlockForest_getItem(const shared_ptr< StructuredBlockForest >& bf, uint_t i) // NOLINT
{
   if (i >= bf->size()) { throw py::value_error("Index out of bounds"); }

   std::vector< const IBlock* > blocks;
   bf->getBlocks(blocks);

   py::object theObject = py::cast(blocks[i]);
   return theObject;
}

std::vector<py::object> StructuredBlockForest_blocksOverlappedByAABB(StructuredBlockForest& s, const AABB& aabb)
{
   std::vector< IBlock* > blocks;
   s.getBlocksOverlappedByAABB(blocks, aabb);

   std::vector<py::object> resultList;
   for (auto it = blocks.begin(); it != blocks.end(); ++it)
      resultList.push_back(py::cast(*it));
   return resultList;
}

std::vector<py::object> StructuredBlockForest_blocksContainedWithinAABB(StructuredBlockForest& s, const AABB& aabb)
{
   std::vector< IBlock* > blocks;
   s.getBlocksContainedWithinAABB(blocks, aabb);

   std::vector<py::object> resultList;
   for (auto it = blocks.begin(); it != blocks.end(); ++it)
      resultList.push_back(py::cast(*it));
   return resultList;
}

py::object SbF_transformGlobalToLocal(StructuredBlockForest& s, IBlock& block, const py::object& global)
{
   if (py::isinstance< CellInterval >(global))
   {
      CellInterval ret;
      s.transformGlobalToBlockLocalCellInterval(ret, block, py::cast< CellInterval >(global));
      return py::cast(ret);
   }
   else if (py::isinstance< Cell >(global))
   {
      Cell ret;
      s.transformGlobalToBlockLocalCell(ret, block, py::cast< Cell >(global));
      return py::cast(ret);
   }

   throw py::value_error("Only CellIntervals and cells can be transformed");
}

py::object SbF_transformLocalToGlobal(StructuredBlockForest& s, IBlock& block, const py::object& local)
{
   if (py::isinstance< CellInterval >(local))
   {
      CellInterval ret;
      s.transformBlockLocalToGlobalCellInterval(ret, block, py::cast< CellInterval >(local));
      return py::cast(ret);
   }
   else if (py::isinstance< Cell >(local))
   {
      Cell ret;
      s.transformBlockLocalToGlobalCell(ret, block, py::cast< Cell >(local));
      return py::cast(ret);
   }
   throw py::value_error("Only CellIntervals and cells can be transformed");
}

void SbF_writeBlockData(StructuredBlockForest& s, const std::string& blockDataId, const std::string& file)
{
   mpi::SendBuffer buffer;
   s.serializeBlockData(blockDataIDFromString(s, blockDataId), buffer);
   mpi::writeMPIIO(file, buffer);
}

void SbF_readBlockData(StructuredBlockForest& s, const std::string& blockDataId, const std::string& file)
{
   mpi::RecvBuffer buffer;
   mpi::readMPIIO(file, buffer);

   s.deserializeBlockData(blockDataIDFromString(s, blockDataId), buffer);
   if (!buffer.isEmpty())
   { throw py::cast_error("Reading failed - file does not contain matching data for this type."); }
}

CellInterval SbF_getBlockCellBB(StructuredBlockForest& s, const IBlock* block) { return s.getBlockCellBB(*block); }

std::array<real_t , 3> SbF_mapToPeriodicDomain1(StructuredBlockForest& s, real_t x, real_t y, real_t z)
{
   Vector3< real_t > res(x, y, z);
   s.mapToPeriodicDomain(res);
   return std::array< real_t, 3 >{ res[0], res[1], res[2] };
}

std::array<real_t , 3> SbF_mapToPeriodicDomain2(StructuredBlockForest& s, const std::array<real_t, 3>& in)
{
   Vector3< real_t > tmp(in[0], in[1], in[2]);
   s.mapToPeriodicDomain(tmp);
   return std::array< real_t, 3 >{ tmp[0], tmp[1], tmp[2] };
}

Cell SbF_mapToPeriodicDomain3(StructuredBlockForest& s, Cell in, uint_t level = 0)
{
   s.mapToPeriodicDomain(in, level);
   return in;
}

py::object SbF_getBlock1(StructuredBlockForest& s, const real_t x, const real_t y, const real_t z)
{
   return py::cast(s.getBlock(x, y, z));
}

py::object SbF_getBlock2(StructuredBlockForest& s, const std::array<real_t, 3>& v)
{
   return py::cast(s.getBlock(Vector3<real_t>(v[0], v[1], v[2])));

}

py::tuple SbF_periodic(StructuredBlockForest& s)
{
   return py::make_tuple(s.isXPeriodic(), s.isYPeriodic(), s.isZPeriodic());
}

bool p_blockExists1(StructuredBlockForest& s, const std::array<real_t, 3>& v)
{
   return s.blockExists(Vector3<real_t>(v[0], v[1], v[2]));

}

bool p_blockExistsLocally1(StructuredBlockForest& s, const std::array<real_t, 3>& v)
{
   return s.blockExistsLocally(Vector3<real_t>(v[0], v[1], v[2]));

}

bool p_blockExistsRemotely1(StructuredBlockForest& s, const std::array<real_t, 3>& v)
{
   return s.blockExistsRemotely(Vector3<real_t>(v[0], v[1], v[2]));

}

bool SbF_atDomainXMinBorder(StructuredBlockForest& s, const IBlock* b) { return s.atDomainXMinBorder(*b); }
bool SbF_atDomainXMaxBorder(StructuredBlockForest& s, const IBlock* b) { return s.atDomainXMaxBorder(*b); }
bool SbF_atDomainYMinBorder(StructuredBlockForest& s, const IBlock* b) { return s.atDomainYMinBorder(*b); }
bool SbF_atDomainYMaxBorder(StructuredBlockForest& s, const IBlock* b) { return s.atDomainYMaxBorder(*b); }
bool SbF_atDomainZMinBorder(StructuredBlockForest& s, const IBlock* b) { return s.atDomainZMinBorder(*b); }
bool SbF_atDomainZMaxBorder(StructuredBlockForest& s, const IBlock* b) { return s.atDomainZMaxBorder(*b); }

void exportBlockForest(py::module_& m)
{
   using namespace pybind11::literals;

   bool (StructuredBlockForest::*p_blockExists2)(const real_t, const real_t, const real_t) const =
      &StructuredBlockForest::blockExists;
   bool (StructuredBlockForest::*p_blockExistsLocally2)(const real_t, const real_t, const real_t) const =
      &StructuredBlockForest::blockExistsLocally;
   bool (StructuredBlockForest::*p_blockExistsRemotely2)(const real_t, const real_t, const real_t) const =
      &StructuredBlockForest::blockExistsRemotely;

   py::class_< StructuredBlockForest, std::shared_ptr< StructuredBlockForest > >(m, "StructuredBlockForest")
      .def("getNumberOfLevels", &StructuredBlockForest::getNumberOfLevels)
      .def_property_readonly("getDomain", &StructuredBlockForest::getDomain)
      .def("mapToPeriodicDomain", &SbF_mapToPeriodicDomain1)
      .def("mapToPeriodicDomain", &SbF_mapToPeriodicDomain2)
      .def("mapToPeriodicDomain", &SbF_mapToPeriodicDomain3)
      .def("__getitem__", &StructuredBlockForest_getItem, py::keep_alive< 0, 1 >())
      .def("__len__", &StructuredBlockForest::size)
      .def("getBlock", SbF_getBlock1, py::keep_alive< 0, 1 >())
      .def("getBlock", SbF_getBlock2, py::keep_alive< 0, 1 >())
      .def("containsGlobalBlockInformation", &StructuredBlockForest::containsGlobalBlockInformation)
      .def("blocksOverlappedByAABB", &StructuredBlockForest_blocksOverlappedByAABB, py::keep_alive< 0, 1 >())
      .def("blocksContainedWithinAABB", &StructuredBlockForest_blocksContainedWithinAABB, py::keep_alive< 0, 1 >())
      .def("blockExists", &p_blockExists1)
      .def("blockExists", p_blockExists2)
      .def("blockExistsLocally", &p_blockExistsLocally1)
      .def("blockExistsLocally", p_blockExistsLocally2)
      .def("blockExistsRemotely", &p_blockExistsRemotely1)
      .def("blockExistsRemotely", p_blockExistsRemotely2)
      .def("atDomainXMinBorder", &SbF_atDomainXMinBorder)
      .def("atDomainXMaxBorder", &SbF_atDomainXMaxBorder)
      .def("atDomainYMinBorder", &SbF_atDomainYMinBorder)
      .def("atDomainYMaxBorder", &SbF_atDomainYMaxBorder)
      .def("atDomainZMinBorder", &SbF_atDomainZMinBorder)
      .def("atDomainZMaxBorder", &SbF_atDomainZMaxBorder)
      .def("dx", &StructuredBlockForest::dx)
      .def("dy", &StructuredBlockForest::dy)
      .def("dz", &StructuredBlockForest::dz)
      .def("getDomainCellBB", &StructuredBlockForest::getDomainCellBB, "level"_a=0)
      .def("getBlockCellBB", &SbF_getBlockCellBB)
      .def("transformGlobalToLocal", &SbF_transformGlobalToLocal)
      .def("transformLocalToGlobal", &SbF_transformLocalToGlobal)
      .def("writeBlockData", &SbF_writeBlockData)
      .def("readBlockData", &SbF_readBlockData)
      .def("__iter__", &StructuredBlockForest_iter, py::keep_alive< 0, 1 >())
      .def_property_readonly("containsGlobalBlockInformation", &StructuredBlockForest::containsGlobalBlockInformation)
      .def_property_readonly("periodic", &SbF_periodic);

   py::class_< SetupBlock, shared_ptr< SetupBlock > >(m, "SetupBlock")
      .def("get_level", &SetupBlock::getLevel)
      .def("set_workload", &SetupBlock::setWorkload)
      .def("get_workload", &SetupBlock::getWorkload)
      .def("set_memory", &SetupBlock::setMemory)
      .def("get_memory", &SetupBlock::getMemory)
      .def("get_aabb", &SetupBlock::getAABB)
      .def("__repr__", &printSetupBlock);

   m.def(
      "createUniformBlockGrid",
      [](std::array< uint_t, 3 > blocks, std::array< uint_t, 3 > cellsPerBlock, real_t dx,
         bool oneBlockPerProcess, std::array< bool, 3 > periodic, bool keepGlobalBlockInformation) {
         return blockforest::createUniformBlockGrid(blocks[0], blocks[1], blocks[2], cellsPerBlock[0], cellsPerBlock[1],
                                                    cellsPerBlock[2], dx, oneBlockPerProcess, periodic[0], periodic[1], periodic[2],
                                                    keepGlobalBlockInformation);
      },
      "blocks"_a, "cellsPerBlock"_a, "dx"_a = real_t(1), "oneBlockPerProcess"_a = true,
      "periodic"_a = std::array< bool, 3 >{ false, false, false }, "keepGlobalBlockInformation"_a = false);

    m.def(
            "createUniformBlockGrid",
            [](std::array< uint_t, 3 > cells, real_t dx,
               bool oneBlockPerProcess, std::array< bool, 3 > periodic, bool keepGlobalBlockInformation) {
                Vector3<uint_t> cellsVec(cells[0], cells[1], cells[2]);
                Vector3<uint_t> cellsPerBlock;
                Vector3<uint_t> blocks;
                uint_t nrOfProcesses = uint_c( MPIManager::instance()->numProcesses() );

                calculateCellDistribution( cellsVec, nrOfProcesses, blocks, cellsPerBlock );


                return blockforest::createUniformBlockGrid(blocks[0], blocks[1], blocks[2], cellsPerBlock[0], cellsPerBlock[1],
                                                           cellsPerBlock[2], dx, oneBlockPerProcess, periodic[0], periodic[1], periodic[2],
                                                           keepGlobalBlockInformation);
            },
            "cells"_a,"dx"_a = real_t(1), "oneBlockPerProcess"_a = true,
            "periodic"_a = std::array< bool, 3 >{ false, false, false }, "keepGlobalBlockInformation"_a = false);
}

} // namespace blockforest
} // namespace walberla

#endif // WALBERLA_BUILD_WITH_PYTHON