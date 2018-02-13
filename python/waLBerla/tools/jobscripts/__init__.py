"""This module can creates job scripts for various supercomputers
    See :func:createJobscript
"""

from __future__ import print_function, absolute_import, division, unicode_literals
from datetime import timedelta
from waLBerla.tools.jobscripts.hornet          import createJobscript as _cr_hornet
from waLBerla.tools.jobscripts.supermuc        import createJobscript as _cr_supermuc
from waLBerla.tools.jobscripts.supermuc_phase2 import createJobscript as _cr_supermuc2
from waLBerla.tools.jobscripts.juqueen         import createJobscript as _cr_juqueen
from waLBerla.tools.jobscripts.pizdaint_hybrid import createJobscript as _cr_pizdainth


def createJobscript(*args, **kwargs):
    """
        :param machine:     Currently supported target machines are  ``supermuc``, ``supermuc_phase2``, ``juqueen`` and ``hornet``
        :param nodes:       Number of nodes to run on. You can either specify nodes or cores. 
        :param cores:       specify eiter nodes or cores. If using more than one node the nodes have to be filled completely
        :param job_class:   optional, the jobclass is usually computed depending on number of nodes and wall_time, this parameter overrides this
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
        
        :param commands:  can be either a list of two-tuples with (executableName, configFile), which are then run in this order with mpirun
                          or a list of string which are just appended to the jobscript file
    """
    if 'machine' not in kwargs:
        raise ValueError("Specify which machine to use with 'machine=<supermuc,juqueen,hornet>'")

    if 'wall_time' in kwargs and isinstance(kwargs['wall_time'], int):
        kwargs['wall_time'] = timedelta(seconds=kwargs['wall_time'])

    if kwargs['machine'].lower() == 'supermuc':        return _cr_supermuc  ( *args, **kwargs )
    if kwargs['machine'].lower() == 'supermuc_phase2': return _cr_supermuc2 ( *args, **kwargs )
    if kwargs['machine'].lower() == 'juqueen' :        return _cr_juqueen   ( *args, **kwargs )
    if kwargs['machine'].lower() == 'hornet'  :        return _cr_hornet    ( *args, **kwargs )
    if kwargs['machine'].lower() == 'pizdaint_hybrid': return _cr_pizdainth ( *args, **kwargs )
    raise ValueError( "Unknown Machine: supported machines <supermuc,supermuc_phase2,juqueen,hornet,pizdaint_hybrid>" )
    