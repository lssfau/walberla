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
    parser.add_argument('-y', action='store_true', help='Silent mode. Accept all questions with yes.')
    args = parser.parse_args()

    mpd = Module(args.path)
    mpd.enable_openmp(False)
    ps = mpd.add(data.ParticleStorage())
    ps.set_shapes('Sphere', 'HalfSpace', 'CylindricalBoundary', 'Box', 'Ellipsoid', 'ConvexPolyhedron')
    ps.add_property("position", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("linearVelocity", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("invMass", "walberla::real_t", defValue="real_t(1)", syncMode="ON_GHOST_CREATION")
    ps.add_property("force", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")
    ps.add_property("oldForce", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")

    ps.add_property("shapeID", "size_t", defValue="", syncMode="ON_GHOST_CREATION")
    ps.add_property("rotation", "walberla::mesa_pd::Rot3", defValue="", syncMode="ALWAYS")
    ps.add_property("angularVelocity", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("torque", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")
    ps.add_property("oldTorque", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ON_OWNERSHIP_CHANGE")

    ps.add_property("radiusAtTemperature", "walberla::real_t", defValue="real_t(0)", syncMode="ALWAYS")

    ps.add_include("blockforest/BlockForest.h")
    ps.add_property("currentBlock", "blockforest::BlockID", defValue="", syncMode="NEVER")

    ps.add_property("type", "uint_t", defValue="0", syncMode="ON_GHOST_CREATION")

    ps.add_property("flags", "walberla::mesa_pd::data::particle_flags::FlagT", defValue="",
                    syncMode="ON_GHOST_CREATION")
    ps.add_property("nextParticle", "int", defValue="-1", syncMode="NEVER")

    ps.add_property("oldContactHistory", "std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>",
                    defValue="", syncMode="ALWAYS")
    ps.add_property("newContactHistory", "std::map<walberla::id_t, walberla::mesa_pd::data::ContactHistory>",
                    defValue="", syncMode="NEVER")

    ps.add_property("temperature", "walberla::real_t", defValue="real_t(0)", syncMode="ALWAYS")
    ps.add_property("heatFlux", "walberla::real_t", defValue="real_t(0)", syncMode="NEVER")

    # Properties for HCSITS
    ps.add_property("dv", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")
    ps.add_property("dw", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")

    # Properties for lbm_mesapd_coupling:
    ps.add_property("hydrodynamicForce", "walberla::mesa_pd::Vec3", defValue="real_t(0)",
                    syncMode="ON_OWNERSHIP_CHANGE")
    ps.add_property("hydrodynamicTorque", "walberla::mesa_pd::Vec3", defValue="real_t(0)",
                    syncMode="ON_OWNERSHIP_CHANGE")
    ps.add_property("oldHydrodynamicForce", "walberla::mesa_pd::Vec3", defValue="real_t(0)",
                    syncMode="ON_OWNERSHIP_CHANGE")
    ps.add_property("oldHydrodynamicTorque", "walberla::mesa_pd::Vec3", defValue="real_t(0)",
                    syncMode="ON_OWNERSHIP_CHANGE")

    ch = mpd.add(data.ContactHistory())
    ch.add_property("tangentialSpringDisplacement", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    ch.add_property("isSticking", "bool", defValue="false")
    ch.add_property("impactVelocityMagnitude", "real_t", defValue="real_t(0)")

    cs = mpd.add(data.ContactStorage())
    cs.add_property("id1", "walberla::id_t", defValue="walberla::id_t(-1)")
    cs.add_property("id2", "walberla::id_t", defValue="walberla::id_t(-1)")
    cs.add_property("distance", "real_t", defValue="real_t(1)")
    cs.add_property("normal", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("position", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("t", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("o", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("r1", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("r2", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("mu", "real_t", defValue="real_t(0)")
    cs.add_property("p", "walberla::mesa_pd::Vec3", defValue="real_t(0)")
    cs.add_property("diag_nto", "walberla::mesa_pd::Mat3", defValue="real_t(0)")
    cs.add_property("diag_nto_inv", "walberla::mesa_pd::Mat3", defValue="real_t(0)")
    cs.add_property("diag_to_inv", "walberla::mesa_pd::Mat2", defValue="real_t(0)")
    cs.add_property("diag_n_inv", "real_t", defValue="real_t(0)")
    cs.add_property("p", "walberla::mesa_pd::Vec3", defValue="real_t(0)")

    mpd.add(data.LinkedCells())
    mpd.add(data.SparseLinkedCells())
    mpd.add(data.ShapeStorage(ps))

    mpd.add(kernel.DetectAndStoreContacts())
    mpd.add(kernel.DoubleCast(ps))
    mpd.add(kernel.ExplicitEuler())
    mpd.add(kernel.ForceLJ())
    mpd.add(kernel.HCSITSRelaxationStep())
    mpd.add(kernel.HeatConduction())
    mpd.add(kernel.InitParticlesForHCSITS())
    mpd.add(kernel.InitContactsForHCSITS())
    mpd.add(kernel.IntegrateParticlesHCSITS())
    mpd.add(kernel.InsertParticleIntoLinkedCells())
    mpd.add(kernel.InsertParticleIntoSparseLinkedCells())
    mpd.add(kernel.LinearSpringDashpot())
    mpd.add(kernel.NonLinearSpringDashpot())
    mpd.add(kernel.PFCDamping())
    mpd.add(kernel.SemiImplicitEuler())
    mpd.add(kernel.SingleCast(ps))
    mpd.add(kernel.SpringDashpot())
    mpd.add(kernel.SpringDashpotSpring())
    mpd.add(kernel.TemperatureIntegration())
    mpd.add(kernel.VelocityVerlet())

    mpd.add(mpi.BroadcastProperty())
    mpd.add(mpi.ClearGhostOwnerSync())
    mpd.add(mpi.ClearNextNeighborSync())
    mpd.add(mpi.Notifications(ps))
    ftn = mpd.add(mpi.PropertyNotification('ForceTorqueNotification'))
    ftn.add_property('force', 'mesa_pd::Vec3', 'Vec3(real_t(0))')
    ftn.add_property('torque', 'mesa_pd::Vec3', 'Vec3(real_t(0))')
    hftn = mpd.add(mpi.PropertyNotification('HydrodynamicForceTorqueNotification'))
    hftn.add_property('hydrodynamicForce', 'mesa_pd::Vec3', 'Vec3(real_t(0))')
    hftn.add_property('hydrodynamicTorque', 'mesa_pd::Vec3', 'Vec3(real_t(0))')
    hfn = mpd.add(mpi.PropertyNotification('HeatFluxNotification'))
    hfn.add_property('heatFlux', 'real_t', 'real_t(0)')
    mpd.add(mpi.ReduceContactHistory())
    mpd.add(mpi.ReduceProperty())
    mpd.add(mpi.ShapePackUnpack(ps))
    mpd.add(mpi.SyncGhostOwners(ps))
    mpd.add(mpi.SyncNextNeighbors(ps))
    mpd.add(mpi.SyncNextNeighborsNoGhosts(ps))

    mpd.generate(not args.y)
