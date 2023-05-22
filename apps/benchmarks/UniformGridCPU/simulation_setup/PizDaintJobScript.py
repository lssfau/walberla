#!/usr/bin/env python3
import os
from waLBerla.tools.config import block_decomposition


job_script_header = """
#!/bin/bash -l
#SBATCH --job-name=scaling
#SBATCH --time=01:00:00
#SBATCH --nodes={nodes}
#SBATCH -o out_scaling_{nodes}_%j.txt
#SBATCH -e err_scaling_{nodes}_%j.txt
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account=s1042

source ~/env.sh

export MPICH_RDMA_ENABLED_CUDA=1  # allow GPU-GPU data transfer
export CRAY_CUDA_MPS=1            # allow GPU sharing
export MPICH_G2G_PIPELINE=256     # adapt maximum number of concurrent in-flight messages

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

export MPICH_RANK_REORDER_METHOD=3
export PMI_MMAP_SYNC_WAIT_TIME=300

cd {folder}
# grid_order -R -H -c 1,1,8 -g 16,16,8

ulimit -c 0
"""

job_script_exe_part = """

export WALBERLA_SCENARIO_IDX=0
while srun -n {nodes} ./{app} {config}
do
 ((WALBERLA_SCENARIO_IDX++))
done
"""

streaming_patterns = ['pull', 'push', 'aa', 'esotwist']
stencils = ['d3q27', 'd3q19']
methods = ['srt', 'mrt', 'cumulant', 'entropic']

all_executables = []

for stencil in stencils:
    for streaming_pattern in streaming_patterns:
        for method in methods:
            all_executables.append(f"UniformGridGPU_{stencil}_{streaming_pattern}_{method}")

all_executables = tuple(all_executables)


def generate_jobscripts(exe_names=all_executables):
    for node_count in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 2400]:
        folder_name = "scaling_{:04d}".format(node_count)
        os.makedirs(folder_name, exist_ok=True)

        # run grid_order
        import subprocess
        decomposition = block_decomposition(node_count)
        decomposition_str = ",".join(str(e) for e in decomposition)
        subprocess.check_call(['grid_order', '-R', '-H', '-g', decomposition_str])

        job_script = job_script_header.format(nodes=node_count, folder=os.path.join(os.getcwd(), folder_name))
        for exe in exe_names:
            job_script += job_script_exe_part.format(app="../" + exe, nodes=node_count,
                                                     config='../communication_compare.py')

        with open(os.path.join(folder_name, 'job.sh'), 'w') as f:
            f.write(job_script)


if __name__ == '__main__':
    print("Called without waLBerla - generating job scripts for PizDaint")
    generate_jobscripts()
