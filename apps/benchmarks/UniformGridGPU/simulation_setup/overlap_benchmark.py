#!/usr/bin/env python3

import os
import pandas as pd
import waLBerla as wlb
from waLBerla.tools.config import block_decomposition
from waLBerla.tools.sqlitedb import sequenceValuesToScalars
from os import getcwd
from waLBerla.tools.jobscripts import createJobscript
from datetime import timedelta
from copy import deepcopy
import sys

CSV_FILE = "overlap_benchmark.csv"

BASE_CONFIG = {
    'DomainSetup': {
        'cellsPerBlock': (256, 128, 128),
        'periodic': (1, 1, 1),
    },
    'Parameters': {
        'omega': 1.8,
        'timesteps': 400,
        'cudaEnabledMPI': False,
        'warmupSteps': 5,
        'outerIterations': 3,
        'initShearFlow': True,
    }
}


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

        df = pd.DataFrame.from_records([data])
        if not os.path.isfile(CSV_FILE):
            df.to_csv(CSV_FILE, index=False)
        else:
            df.to_csv(CSV_FILE, index=False, mode='a', header=False)


def overlap_benchmark():
    scenarios = wlb.ScenarioManager()
    inner_outer_splits = [(1, 1, 1), (4, 1, 1), (8, 1, 1), (16, 1, 1), (32, 1, 1),
                          (4, 4, 1), (8, 8, 1), (16, 16, 1), (32, 32, 1),
                          (4, 4, 4), (8, 8, 8), (16, 16, 16), (32, 32, 32)]

    for comm_strategy in ['UniformGPUScheme_Baseline', 'UniformGPUScheme_Memcpy']:  # 'GPUPackInfo_Baseline', 'GPUPackInfo_Streams'
        # no overlap
        scenarios.add(Scenario(timeStepStrategy='noOverlap', communicationScheme=comm_strategy, innerOuterSplit=(1, 1, 1)))

        # overlap
        for overlap_strategy in ['simpleOverlap', 'complexOverlap']:
            for inner_outer_split in inner_outer_splits:
                scenario = Scenario(timeStepStrategy=overlap_strategy,
                                    communicationScheme=comm_strategy,
                                    innerOuterSplit=inner_outer_split)
                scenarios.add(scenario)


def single_gpu_benchmark():
    scenarios = wlb.ScenarioManager()
    block_sizes = [(i, i, i) for i in (64, 128, 256, 384)] + [(512, 512, 128)]
    cuda_blocks = [(32, 1, 1), (64, 1, 1), (128, 1, 1), (256, 1, 1), (512, 1, 1),
                   (32, 2, 1), (64, 2, 1), (128, 2, 1), (256, 2, 1),
                   (32, 4, 1), (64, 4, 1), (128, 4, 1),
                   (32, 8, 1), (64, 8, 1),
                   (32, 16, 1)]
    for block_size in block_sizes:
        for cuda_block_size in cuda_blocks:
            cells = block_size[0] * block_size[1] * block_size[2]
            time_steps_for_128_cubed = 1000
            time_steps = (128 ** 3 / cells) * time_steps_for_128_cubed
            scenario = Scenario(cells_per_block=block_size,
                                gpuBlockSize=cuda_block_size,
                                timeStepStrategy='kernelOnly',
                                timesteps=int(time_steps))
            scenarios.add(scenario)


all_executables = ('UniformGridBenchmarkGPU_AA_entropic',
                   'UniformGridBenchmarkGPU_AA_mrt',
                   'UniformGridBenchmarkGPU_AA_smagorinsky',
                   'UniformGridBenchmarkGPU_AA_srt',
                   'UniformGridBenchmarkGPU_AA_trt',
                   'UniformGridBenchmarkGPU_entropic',
                   'UniformGridBenchmarkGPU_mrt',
                   'UniformGridBenchmarkGPU_smagorinsky',
                   'UniformGridBenchmarkGPU_srt',
                   'UniformGridBenchmarkGPU_trt')


def generate_jobscripts(machine='pizdaint_hybrid',
                        exe_names=all_executables):
    for node_count in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 2400]:
        with open("job_overlap_benchmark_{:04d}.sh".format(node_count), 'w') as f:

            js = createJobscript(nodes=node_count,
                                 output_file='overlap_bench_{:04d}_%j.txt'.format(node_count),
                                 error_file='overlap_bench_{:04d}_%j.txt'.format(node_count),
                                 initial_dir=getcwd(),
                                 commands=list(("./" + exe, 'overlap_benchmark.py') for exe in exe_names),
                                 wall_time=timedelta(minutes=25),
                                 machine=machine,
                                 account='d105',
                                 )
            f.write(js)


if __name__ == '__main__':
    print("Called without waLBerla - generating job scripts for PizDaint")
    generate_jobscripts()
else:
    single_gpu_benchmark()
