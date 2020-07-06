# -*- coding: utf-8 -*-

from ..utility import generate_file


class Notifications:
    def __init__(self, particle_storage):
        self.ps = particle_storage

    def generate(self, module):
        ctx = {'module': module, 'particle': self.ps.get_context()}

        generate_file(module['module_path'], 'mpi/notifications/ParseMessage.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleCopyNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleGhostCopyNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/NewGhostParticleNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleMigrationNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleRemoteMigrationNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleRemovalInformationNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleRemovalNotification.templ.h', ctx)
        generate_file(module['module_path'], 'mpi/notifications/ParticleUpdateNotification.templ.h', ctx)
