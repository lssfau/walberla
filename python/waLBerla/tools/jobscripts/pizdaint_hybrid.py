from __future__ import print_function, absolute_import, division, unicode_literals

import os
import math


def createJobscript(wall_time=None, nodes=None, cores=None, initial_dir=None, job_name="waLBerla",
                    exe_name=None, parameter_files=[], commands=[], hyperthreading=1,
                    output_file=None, error_file=None, account=None, **kwargs):
    if type(hyperthreading) is bool:
        hyperthreading = 2 if hyperthreading else 1

    CORES_PER_NODE = 1  # 12 * hyperthreading

    if wall_time and wall_time.total_seconds() > 24 * 3600:
        raise ValueError("No jobs longer that 24h allowed")

    if hyperthreading > 2:
        raise ValueError("PizDaint supports only two way hyperthreading (requested %d)" % (hyperthreading,))

    if nodes is not None and cores is not None:
        raise ValueError("You can either specify nodes or cores - not both.")

    if nodes is None and cores is None:
        raise ValueError('Specify either cores or nodes.')

    if nodes is None:
        nodes = math.ceil(cores / CORES_PER_NODE)
    if cores is None:
        cores = nodes * CORES_PER_NODE

    if cores > CORES_PER_NODE and cores % CORES_PER_NODE != 0:
        raise ValueError("When using more than one node, the number of cores has to be a multiple of 12")

    if not output_file:
        output_file = job_name
    if not error_file:
        error_file = job_name

    partition = 'normal'
    if nodes <= 4 and wall_time.total_seconds() < 30 * 60:
        partition = 'debug'

    tasks_per_node = min(CORES_PER_NODE, cores)
    additional_lines = ""
    if account:
        additional_lines += '#SBATCH --account=%s\n' % (account,)

    template_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "pizdaint_hybrid.job")
    result = open(template_file).read().format(cores=cores,
                                               nodes=nodes,
                                               tasks_per_core=hyperthreading,
                                               tasks_per_node=tasks_per_node,
                                               cpus_per_task=1,  # OpenMP num threads would go here
                                               initial_dir=initial_dir,
                                               output_file=output_file,
                                               additional_lines=additional_lines,
                                               error_file=error_file,
                                               partition=partition,
                                               job_name=job_name,
                                               wall_time=wall_time)

    exec_line = "srun -n %d %s %s \n"

    if exe_name is not None:
        for param_file in parameter_files:
            result += exec_line % (cores, exe_name, param_file)

    for exe_paramfile_pair in commands:
        if type(exe_paramfile_pair) is not tuple:
            result += exe_paramfile_pair + "\n"
        else:
            result += exec_line % ((cores,) + exe_paramfile_pair)

    return result
