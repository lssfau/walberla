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
        'outerIterations': 1,
    }
}


class Scenario:
    def __init__(self, **kwargs):
        self.config_dict = deepcopy(BASE_CONFIG)
        self.config_dict['Parameters'].update(kwargs)
        self.config_dict['DomainSetup']['blocks'] = block_decomposition(wlb.mpi.numProcesses())

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


def generate_jobscripts(machine='pizdaint_hybrid'):
    for node_count in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 2400]:
        with open("job_overlap_benchmark_{:04d}.sh".format(node_count), 'w') as f:
            js = createJobscript(nodes=node_count,
                                 output_file='overlap_bench_{:04d}_%j.txt'.format(node_count),
                                 error_file='overlap_bench_{:04d}_%j.txt'.format(node_count),
                                 initial_dir=getcwd(),
                                 exe_name='UniformGridBenchmarkGPU',
                                 parameter_files=['overlap_benchmark.py'],
                                 wall_time=timedelta(minutes=25),
                                 machine=machine,
                                 account='d105',
                                 )
            f.write(js)


if __name__ == '__main__':
    print("Called without waLBerla - generating job scripts for PizDaint")
    generate_jobscripts()
else:
    overlap_benchmark()
