from __future__ import print_function, absolute_import, division, unicode_literals

import os
import math


def createJobscript(wall_time=None, nodes=None, cores=None, job_class=None, island_count=None,
                    initial_dir='~', job_name="waLBerla", energy_tag="",
                    exe_name=None, parameter_files=[], commands=[], hyperthreading=1,
                    output_file=None, error_file=None, **kwargs):
    if type(hyperthreading) == bool:
        hyperthreading = 2 if hyperthreading else 1

    CORES_PER_NODE = 28 * hyperthreading
    NODES_PER_ISLAND = 512

    if wall_time and wall_time.total_seconds() > 48 * 3600:
        raise ValueError("No jobs longer that 48h allowed")

    if hyperthreading > 2:
        raise ValueError("SuperMUC supports only two way hyperthreading (requested %d)" % (hyperthreading,))

    if nodes is not None and cores is not None:
        raise ValueError("You can either specify nodes or cores - not both.")

    if nodes is None and cores is None:
        raise ValueError('Specify either cores or nodes.')

    if nodes is None:
        nodes = math.ceil(cores / CORES_PER_NODE)
    if cores is None:
        cores = nodes * CORES_PER_NODE

    if cores > CORES_PER_NODE and cores % CORES_PER_NODE != 0:
        raise ValueError("When using more than one node, the number of cores has to be a multiple of %d",
                         (CORES_PER_NODE,))

    if island_count is None:
        island_count = math.ceil(nodes / NODES_PER_ISLAND)

    if not output_file:
        output_file = job_name
    if not error_file:
        error_file = job_name

    if not job_class:
        if nodes <= 20:
            if wall_time.total_seconds() < 30 * 60:
                job_class = 'test'
            else:
                job_class = 'micro'
        elif nodes <= 512:
            job_class = 'general'
        elif nodes <= 2048:
            job_class = 'big'
        else:
            job_class = 'special'

    tasks_per_node = min(CORES_PER_NODE, cores)

    task_affinity = "core" if hyperthreading == 1 else "cpu"

    energy_tag_statements = ""
    if len(energy_tag) > 0:
        energy_tag_statements = "\n#@ energy_policy_tag = " + energy_tag + "\n#@ minimize_time_to_solution = yes"

    template_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "supermuc.job")
    result = open(template_file).read().format(cores=cores,
                                               nodes=nodes,
                                               initial_dir=initial_dir,
                                               tasks_per_node=tasks_per_node,
                                               job_class=job_class,
                                               job_name=job_name,
                                               island_count=island_count,
                                               energy_tag_statements=energy_tag_statements,
                                               wall_time=wall_time,
                                               output_file=output_file,
                                               error_file=error_file,
                                               task_affinity=task_affinity)

    exec_line = "mpiexec -n %d %s %s \n"

    if exe_name is not None:
        for param_file in parameter_files:
            result += exec_line % (cores, exe_name, param_file)

    for exe_paramfile_pair in commands:
        if type(exe_paramfile_pair) is not tuple:
            result += exe_paramfile_pair + "\n"
        else:
            result += exec_line % (cores, exe_paramfile_pair[0], exe_paramfile_pair[1])

    return result
