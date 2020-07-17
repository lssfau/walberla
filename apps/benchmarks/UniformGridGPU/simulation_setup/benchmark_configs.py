#!/usr/bin/env python3
"""
This is a waLBerla parameter file that tests (almost) all parameter combinations for GPU communication.
Build waLBerla with -DWALBERLA_BUILD_WITH_PYTHON=1  then run e.g.
 ./UniformGridBenchmarkGPU_AA_trt simulation_setup/benchmark_configs.py

Look at the end of the file to select the benchmark to run
"""

import os
import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars, checkAndUpdateSchema, storeSingle
from copy import deepcopy
import sys
import sqlite3

# Number of time steps run for a workload of 128^3 per GPU
# if double as many cells are on the GPU, half as many time steps are run etc.
# increase this to get more reliable measurements
TIME_STEPS_FOR_128_BLOCK = 200
DB_FILE = "gpu_benchmark.sqlite3"

BASE_CONFIG = {
    'DomainSetup': {
        'cellsPerBlock': (256, 128, 128),
        'periodic': (1, 1, 1),
    },
    'Parameters': {
        'omega': 1.8,
        'cudaEnabledMPI': False,
        'warmupSteps': 5,
        'outerIterations': 3,
    }
}


def num_time_steps(block_size):
    cells = block_size[0] * block_size[1] * block_size[2]
    time_steps = (128 ** 3 / cells) * TIME_STEPS_FOR_128_BLOCK
    return int(time_steps)


class Scenario:
    def __init__(self, cells_per_block=(256, 128, 128), **kwargs):
        self.config_dict = deepcopy(BASE_CONFIG)
        self.config_dict['Parameters'].update(kwargs)
        self.config_dict['DomainSetup']['blocks'] = block_decomposition(wlb.mpi.numProcesses())
        self.config_dict['DomainSetup']['cellsPerBlock'] = cells_per_block

    @wlb.member_callback
    def config(self, **kwargs):
        from pprint import pformat
        wlb.log_info_on_root("Scenario:\n" + pformat(self.config_dict))
        # Write out the configuration as text-based prm:
        # print(toPrm(self.config_dict))
        return self.config_dict

    @wlb.member_callback
    def results_callback(self, **kwargs):
        data = {}
        data.update(self.config_dict['Parameters'])
        data.update(self.config_dict['DomainSetup'])
        data.update(kwargs)
        data['executable'] = sys.argv[0]
        data['compile_flags'] = wlb.build_info.compiler_flags
        data['walberla_version'] = wlb.build_info.version
        data['build_machine'] = wlb.build_info.build_machine
        sequenceValuesToScalars(data)

        result = data
        sequenceValuesToScalars(result)
        num_tries = 4
        # check multiple times e.g. may fail when multiple benchmark processes are running
        for num_try in range(num_tries):
            try:
                checkAndUpdateSchema(result, "runs", DB_FILE)
                storeSingle(result, "runs", DB_FILE)
                break
            except sqlite3.OperationalError as e:
                wlb.log_warning("Sqlite DB writing failed: try {}/{}  {}".format(num_try + 1, num_tries, str(e)))


# -------------------------------------- Functions trying different parameter sets -----------------------------------


def overlap_benchmark():
    """Tests different communication overlapping strategies"""
    wlb.log_info_on_root("Running different communication overlap strategies")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    inner_outer_splits = [(1, 1, 1), (4, 1, 1), (8, 1, 1), (16, 1, 1), (32, 1, 1),
                          (4, 4, 1), (8, 8, 1), (16, 16, 1), (32, 32, 1),
                          (4, 4, 4), (8, 8, 8), (16, 16, 16), (32, 32, 32)]

    # 'GPUPackInfo_Baseline', 'GPUPackInfo_Streams'
    for comm_strategy in ['UniformGPUScheme_Baseline', 'UniformGPUScheme_Memcpy']:
        # no overlap
        scenarios.add(Scenario(timeStepStrategy='noOverlap', communicationScheme=comm_strategy,
                               innerOuterSplit=(1, 1, 1)))

        # overlap
        for overlap_strategy in ['simpleOverlap', 'complexOverlap']:
            for inner_outer_split in inner_outer_splits:
                scenario = Scenario(timeStepStrategy=overlap_strategy,
                                    communicationScheme=comm_strategy,
                                    innerOuterSplit=inner_outer_split,
                                    timesteps=num_time_steps((256, 128, 128)))
                scenarios.add(scenario)


def communication_compare():
    """Tests different communication strategies"""
    wlb.log_info_on_root("Running benchmarks to compare communication strategies")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    for block_size in [(128, 128, 128), (32, 32, 32), (64, 64, 64), (256, 256, 256)]:
        for comm_strategy in ['UniformGPUScheme_Baseline', 'UniformGPUScheme_Memcpy']:

            sc = Scenario(cells_per_block=block_size,
                          gpuBlockSize=(128, 1, 1),
                          timeStepStrategy='noOverlap',
                          communicationScheme=comm_strategy,
                          timesteps=num_time_steps(block_size))
            scenarios.add(sc)
            for inner_outer_split in [(4, 1, 1), (8, 1, 1), (16, 1, 1), (32, 1, 1)]:
                # ensure that the inner part of the domain is still large enough
                if 3 * inner_outer_split[0] > block_size[0]:
                    continue
                sc = Scenario(cells_per_block=block_size,
                              gpuBlockSize=(128, 1, 1),
                              timeStepStrategy='simpleOverlap',
                              innerOuterSplit=inner_outer_split,
                              communicationScheme=comm_strategy,
                              timesteps=num_time_steps(block_size))
                scenarios.add(sc)


def single_gpu_benchmark():
    """Benchmarks only the LBM compute kernel"""
    wlb.log_info_on_root("Running single GPU benchmarks")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    block_sizes = [(i, i, i) for i in (64, 128, 256, 384)] + [(512, 512, 128)]
    cuda_blocks = [(32, 1, 1), (64, 1, 1), (128, 1, 1), (256, 1, 1), (512, 1, 1),
                   (32, 2, 1), (64, 2, 1), (128, 2, 1), (256, 2, 1),
                   (32, 4, 1), (64, 4, 1), (128, 4, 1),
                   (32, 8, 1), (64, 8, 1),
                   (32, 16, 1)]
    for block_size in block_sizes:
        for cuda_block_size in cuda_blocks:
            scenario = Scenario(cells_per_block=block_size,
                                gpuBlockSize=cuda_block_size,
                                timeStepStrategy='kernelOnly',
                                timesteps=num_time_steps(block_size))
            scenarios.add(scenario)


# -------------------------------------- Optional job script generation for PizDaint ---------------------------------


job_script_header = """
#!/bin/bash -l
#SBATCH --job-name=scaling
#SBATCH --time=0:30:00
#SBATCH --nodes={nodes}
#SBATCH -o out_scaling_{nodes}_%j.txt
#SBATCH -e err_scaling_{nodes}_%j.txt
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account=d105

cd {folder}

source ~/env.sh

module load daint-gpu
module load craype-accel-nvidia60
export MPICH_RDMA_ENABLED_CUDA=1  # allow GPU-GPU data transfer
export CRAY_CUDA_MPS=1            # allow GPU sharing
export MPICH_G2G_PIPELINE=256     # adapt maximum number of concurrent in-flight messages

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export CRAY_CUDA_MPS=1

export MPICH_RANK_REORDER_METHOD=3
export PMI_MMAP_SYNC_WAIT_TIME=300


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


all_executables = ('UniformGridBenchmarkGPU_mrt_d3q27',
                   'UniformGridBenchmarkGPU_smagorinsky_d3q27',
                   'UniformGridBenchmarkGPU_cumulant'
                   'UniformGridBenchmarkGPU_cumulant_d3q27')


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
else:
    wlb.log_info_on_root("Batch run of benchmark scenarios, saving result to {}".format(DB_FILE))
    # Select the benchmark you want to run
    single_gpu_benchmark()
    # benchmarks different CUDA block sizes and domain sizes and measures single
    # GPU performance of compute kernel (no communication)
    # communication_compare(): benchmarks different communication routines, with and without overlap
    # overlap_benchmark(): benchmarks different communication overlap options
