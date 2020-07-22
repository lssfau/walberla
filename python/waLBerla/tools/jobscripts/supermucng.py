from __future__ import print_function, absolute_import, division, unicode_literals

import os
import math


def createJobscript_supermucng(wall_time=None, nodes=None, cores=None, job_class=None,
                               initial_dir='~', job_name="waLBerla", exe_name=None,
                               parameter_files=[], commands=[], hyperthreading=1,
                               output_file=None, error_file=None, account=None,
                               fixed_freq=True, omp_num_threads=1, **_):
    if type(hyperthreading) == bool:
        hyperthreading = 2 if hyperthreading else 1

    cores_per_node = 48 * hyperthreading

    if wall_time and wall_time.total_seconds() > 48 * 3600:
        raise ValueError("No jobs longer that 48h allowed")

    if hyperthreading > 2:
        raise ValueError("SuperMUC supports only two way hyperthreading (requested %d)" % (hyperthreading,))

    if nodes is not None and cores is not None:
        raise ValueError("You can either specify nodes or cores - not both.")

    if nodes is None and cores is None:
        raise ValueError('Specify either cores or nodes.')

    if nodes is None:
        nodes = math.ceil(cores / cores_per_node)
    if cores is None:
        cores = nodes * cores_per_node

    if cores > cores_per_node and cores % cores_per_node != 0:
        raise ValueError("When using more than one node, the number of cores has to be a multiple of %d",
                         (cores_per_node,))

    if not output_file:
        output_file = job_name
    if not error_file:
        error_file = job_name

    if not job_class:
        if nodes <= 16:
            if wall_time.total_seconds() < 30 * 60:
                job_class = 'test'
            else:
                job_class = 'micro'
        elif nodes <= 792:
            job_class = 'general'
        elif nodes <= 3168:
            job_class = 'big'
        else:
            job_class = 'special'

    if cores_per_node % omp_num_threads != 0:
        raise ValueError("Could not divede cores_per_node %d to omp_num_threads %d", (cores_per_node, omp_num_threads))
    tasks_per_node = min(cores_per_node // omp_num_threads, cores)

    omp_places = "cores" if hyperthreading == 1 else "threads"

    template_file = os.path.join(os.path.dirname(os.path.realpath(__file__)), "supermucng.job")
    additional_lines = ""
    if account:
        additional_lines += '#SBATCH --account=%s' % (account,)
    if fixed_freq:
        additional_lines += '#SBATCH --ear=off\n'

    result = open(template_file).read().format(cores=cores,
                                               nodes=nodes,
                                               initial_dir=initial_dir,
                                               tasks_per_node=tasks_per_node,
                                               tasks_per_core=hyperthreading,
                                               cpus_per_task=omp_num_threads,
                                               partition=job_class,
                                               job_name=job_name,
                                               wall_time=wall_time,
                                               omp_places=omp_places,
                                               output_file=output_file,
                                               additional_lines=additional_lines,
                                               error_file=error_file)

    exec_line = "srun %s %s \n"

    if exe_name is not None:
        for param_file in parameter_files:
            result += exec_line % (exe_name, param_file)

    for exe_paramfile_pair in commands:
        if type(exe_paramfile_pair) is not tuple:
            result += exe_paramfile_pair + "\n"
        else:
            result += exec_line % (exe_paramfile_pair[0], exe_paramfile_pair[1])

    return result


if __name__ == '__main__':
    from waLBerla.tools.jobscripts import createJobscript

    print(createJobscript(wall_time=60 * 60, nodes=4, exe_name='grandchem', parameter_files=['a.cfg', 'b.cfg'],
                          machine='supermuc_ng'))
