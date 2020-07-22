#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from mesa_pd import Module
import mesa_pd.data as data
import mesa_pd.kernel as kernel
import mesa_pd.mpi as mpi

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate all necessary files for the waLBerla mesa_pd module.')
    parser.add_argument('path', help='Where should the files be created?')
    args = parser.parse_args()

    mpd = Module(args.path)
    ps = mpd.add(data.ParticleStorage())

    ps.set_shapes("Sphere", "HalfSpace")

    ps.add_property("position",         "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("linearVelocity",   "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("invMass",          "walberla::real_t",        defValue="real_t(1)", syncMode="ON_GHOST_CREATION")
    ps.add_property("force",            "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")

    ps.add_property("shapeID",          "size_t",                  defValue="",          syncMode="ON_GHOST_CREATION")
    ps.add_property("rotation",         "walberla::mesa_pd::Rot3", defValue="",          syncMode="ALWAYS")
    ps.add_property("angularVelocity",  "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("torque",           "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")

    ps.add_property("type",             "uint_t",                  defValue="0",         syncMode="ON_GHOST_CREATION")

    ps.add_property("flags",            "walberla::mesa_pd::data::particle_flags::FlagT", defValue="",
                    syncMode="ON_GHOST_CREATION")
    ps.add_property("nextParticle",     "int",                     defValue="-1",        syncMode="NEVER")

    ps.add_include("blockforest/BlockForest.h")
    ps.add_property("currentBlock",     "blockforest::BlockID",    defValue="",          syncMode="NEVER")

    ps.add_property("oldContactHistory", "std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>",
                    defValue="", syncMode="ALWAYS")
    ps.add_property("newContactHistory", "std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>",
                    defValue="", syncMode="NEVER")

    ch = mpd.add(data.ContactHistory())
    ch.add_property("tangentialSpringDisplacement", "walberla::mesa_pd::Vec3", defValue="real_t(0)")

    mpd.add(data.LinkedCells())
    mpd.add(data.SparseLinkedCells())
    mpd.add(data.ShapeStorage(ps))

    mpd.add(kernel.DoubleCast(ps))
    mpd.add(kernel.ExplicitEuler())
    mpd.add(kernel.ForceLJ())
    mpd.add(kernel.HeatConduction())
    mpd.add(kernel.InsertParticleIntoLinkedCells())
    mpd.add(kernel.LinearSpringDashpot())
    mpd.add(kernel.NonLinearSpringDashpot())
    mpd.add(kernel.SingleCast(ps))
    mpd.add(kernel.SpringDashpot())
    mpd.add(kernel.TemperatureIntegration())
    mpd.add(kernel.VelocityVerlet())
    mpd.add(kernel.VelocityVerlet())

    mpd.add(mpi.BroadcastProperty())
    mpd.add(mpi.ClearGhostOwnerSync())
    mpd.add(mpi.ClearNextNeighborSync())
    mpd.add(mpi.Notifications(ps))
    ftn = mpd.add(mpi.PropertyNotification('ForceTorqueNotification'))
    ftn.add_property('force', 'mesa_pd::Vec3', 'Vec3(real_t(0))')
    ftn.add_property('torque', 'mesa_pd::Vec3', 'Vec3(real_t(0))')
    mpd.add(mpi.ReduceContactHistory())
    mpd.add(mpi.ReduceProperty())
    mpd.add(mpi.ShapePackUnpack(ps))
    mpd.add(mpi.SyncGhostOwners(ps))
    mpd.add(mpi.SyncNextNeighbors(ps))
    mpd.add(mpi.SyncNextNeighborsNoGhosts(ps))

    mpd.generate()
