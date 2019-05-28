#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from mesa_pd.accessor import Accessor
import mesa_pd.data as data
import mesa_pd.kernel as kernel
import mesa_pd.mpi as mpi

import argparse
import numpy as np
import os

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Generate all necessary files for the waLBerla mesa_pd module.')
   parser.add_argument('path', help='Where should the files be created?')
   parser.add_argument("-f", "--force", help="Generate the files even if not inside a waLBerla directory.",
                       action="store_true")
   args = parser.parse_args()

   if ((not os.path.isfile(args.path + "/src/walberla.h")) and (not args.force)):
      raise RuntimeError(args.path + " is not the path to a waLBerla root directory! Specify -f to generate the files anyway.")

   os.makedirs(args.path + "/src/mesa_pd/common", exist_ok = True)
   os.makedirs(args.path + "/src/mesa_pd/data", exist_ok = True)
   os.makedirs(args.path + "/src/mesa_pd/domain", exist_ok = True)
   os.makedirs(args.path + "/src/mesa_pd/kernel", exist_ok = True)
   os.makedirs(args.path + "/src/mesa_pd/mpi/notifications", exist_ok = True)
   os.makedirs(args.path + "/src/mesa_pd/vtk", exist_ok = True)

   shapes = ["Sphere", "HalfSpace"]

   ps    = data.ParticleStorage()
   ch    = data.ContactHistory()
   lc    = data.LinkedCells()
   ss    = data.ShapeStorage(ps, shapes)

   ps.addProperty("position",         "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
   ps.addProperty("linearVelocity",   "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
   ps.addProperty("invMass",          "walberla::real_t",        defValue="real_t(1)", syncMode="COPY")
   ps.addProperty("force",            "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")

   ps.addProperty("shapeID",          "size_t",                  defValue="",          syncMode="COPY")
   ps.addProperty("rotation",         "walberla::mesa_pd::Rot3", defValue="",          syncMode="ALWAYS")
   ps.addProperty("angularVelocity",  "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
   ps.addProperty("torque",           "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="NEVER")

   ps.addProperty("type",             "uint_t",                  defValue="0",         syncMode="COPY")

   ps.addProperty("flags",            "walberla::mesa_pd::data::particle_flags::FlagT", defValue="", syncMode="COPY")
   ps.addProperty("nextParticle",     "int",                     defValue="-1",        syncMode="NEVER")

   kernels = []
   kernels.append( kernel.DoubleCast(shapes) )
   kernels.append( kernel.ExplicitEuler() )
   kernels.append( kernel.ExplicitEulerWithShape() )
   kernels.append( kernel.ForceLJ() )
   kernels.append( kernel.HeatConduction() )
   kernels.append( kernel.InsertParticleIntoLinkedCells() )
   kernels.append( kernel.LinearSpringDashpot() )
   kernels.append( kernel.NonLinearSpringDashpot() )
   kernels.append( kernel.SingleCast(shapes) )
   kernels.append( kernel.SpringDashpot() )
   kernels.append( kernel.TemperatureIntegration() )
   kernels.append( kernel.VelocityVerlet() )
   kernels.append( kernel.VelocityVerletWithShape() )

   ac = Accessor()
   for k in kernels:
      ac.mergeRequirements(k.getRequirements())
   ac.printSummary()

   comm = []
   comm.append(mpi.BroadcastProperty())
   comm.append(mpi.ClearNextNeighborSync())
   comm.append(mpi.ReduceContactHistory())
   comm.append(mpi.ReduceProperty())
   comm.append(mpi.SyncNextNeighbors(ps))


   ps.generate(args.path + "/src/mesa_pd/")
   ch.generate(args.path + "/src/mesa_pd/")
   lc.generate(args.path + "/src/mesa_pd/")
   ss.generate(args.path + "/src/mesa_pd/")

   for k in kernels:
      k.generate(args.path + "/src/mesa_pd/")

   for c in comm:
      c.generate(args.path + "/src/mesa_pd/")
