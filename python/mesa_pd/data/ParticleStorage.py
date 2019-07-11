# -*- coding: utf-8 -*-

import numpy as np
from ..Container import Container
from ..utility import generateFile

class ParticleStorage(Container):
   def __init__(self):
      super().__init__()

      self.addInclude("mesa_pd/data/Flags.h")

      self.addProperty("uid",               "walberla::id_t",       defValue = "UniqueID<data::Particle>::invalidID()", syncMode="ALWAYS")
      self.addProperty("position",          "walberla::mesa_pd::Vec3", defValue = "real_t(0)", syncMode="ALWAYS")
      self.addProperty("interactionRadius", "walberla::real_t",     defValue = "real_t(0)", syncMode="COPY")
      self.addProperty("flags",             "walberla::mesa_pd::data::particle_flags::FlagT", defValue = "", syncMode="COPY")
      self.addProperty("owner",             "int",                  defValue = "-1", syncMode="COPY")
      self.addProperty("ghostOwners",       "std::unordered_set<walberla::mpi::MPIRank>", defValue = "", syncMode="MIGRATION")

   def generate(self, path):
      self.unrollDimension()

      print("="*90)
      print("Creating Particle Datastructure:")
      print("")
      print("{0: <20}{1: <30}{2: <20}{3: <10}".format("Type", "Name", "Def. Value", "SyncMode"))
      print("="*90)
      for prop in self.properties:
         print("{0: <20.19}{1: <30.29}{2: <20.19}{3: <10.9}".format(prop.type, prop.name, prop.defValue, prop.syncMode))
      print("="*90)

      context = dict()
      context["includes"]    = self.includes
      context["properties"]  = self.properties

      generateFile(path, 'data/ParticleStorage.templ.h', context, filename='data/ParticleStorage.h')
      generateFile(path, 'data/ParticleAccessor.templ.h', context, filename='data/ParticleAccessor.h')

      generateFile(path, 'mpi/notifications/ForceTorqueNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/HeatFluxNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/ParseMessage.templ.h', context)
      generateFile(path, 'mpi/notifications/ParticleCopyNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/NewGhostParticleNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/ParticleMigrationNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/ParticleRemoteMigrationNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/ParticleRemovalInformationNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/ParticleRemovalNotification.templ.h', context)
      generateFile(path, 'mpi/notifications/ParticleUpdateNotification.templ.h', context)
