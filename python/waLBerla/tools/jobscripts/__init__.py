"""This module can creates job scripts for various supercomputers
    See :func:createJobscript
"""

from __future__ import print_function, absolute_import, division, unicode_literals
from datetime import timedelta
from waLBerla.tools.jobscripts.hornet import createJobscript as _cr_hornet
from waLBerla.tools.jobscripts.supermuc import createJobscript as _cr_supermuc
from waLBerla.tools.jobscripts.supermuc_phase2 import createJobscript as _cr_supermuc2
from waLBerla.tools.jobscripts.supermucng import createJobscript as _cr_supermuc_ng
from waLBerla.tools.jobscripts.pizdaint_hybrid import createJobscript as _cr_pizdainth


def createJobscript(*args, **kwargs):
    """
        :param machine:     Currently supported target machines are  ``supermuc``, ``supermuc_phase2``,
                            ``juqueen`` and ``hornet``
        :param nodes:       Number of nodes to run on. You can either specify nodes or cores.
        :param cores:       specify eiter nodes or cores. If using more than one node the nodes
                            have to be filled completely
        :param job_class:   optional, the jobclass is usually computed depending on number of nodes and wall_time,
                            this parameter overrides this
        :param initial_dir: initial working directory of the job, optional, defaults to home directory
        :param job_name:    name of the job in the queuing system, defaults to 'waLBerla'
        :param output_file: file where stdout will be redirected to by the queueing system
        :param input_file:  file where stderr will be redirected to by the queueing system
        :param energy_tag:  energy tag for SuperMUC[1,2]

        Use one of the following options:

        Run single program with different parameter files ( mpirun is prepended with correct number of processes )

        :param exe_name:         executable name, if not specified only the jobscript header is generated
        :param parameter_files:  list of parameter files to simulate

        Run multiple programs:

        :param commands:  can be either a list of two-tuples with (executableName, configFile), which are then run
                          in this order with mpirun or a list of string which are just appended to the jobscript file
    """
    funcs = {
        'supermuc': _cr_supermuc,
        'supermuc_phase2': _cr_supermuc2,
        'supermuc_ng': _cr_supermuc_ng,
        'hornet': _cr_hornet,
        'pizdaint_hybrid': _cr_pizdainth,
    }
    if 'machine' not in kwargs or kwargs['machine'] not in funcs.keys():
        raise ValueError("Specify which machine to use with 'machine={}'".format(list(funcs.keys())))

    if 'wall_time' in kwargs and isinstance(kwargs['wall_time'], int):
        kwargs['wall_time'] = timedelta(seconds=kwargs['wall_time'])

    return funcs[kwargs['machine']](*args, **kwargs)
