# encoding: utf-8

import os
import pandas as pd
import waLBerla as wlb
import copy
from datetime import datetime


CommunicationSchemeType = {
    'GPUPackInfo_Baseline': 0,
    'GPUPackInfo_Streams': 1,
    'UniformGPUScheme_Baseline': 2,
    'UniformGPUScheme_Memcpy': 3,
}

CommunicationSchemeName = {
    0: 'GPUPackInfo_Baseline',
    1: 'GPUPackInfo_Streams',
    2: 'UniformGPUScheme_Baseline',
    3: 'UniformGPUScheme_Memcpy',
}

# Base configuration for the benchmark
BASE_CONFIG = {
    'DomainSetup' : {
        'cellsPerBlock': (64, 64, 64),
        'blocks': (1, 1, 1),
        'nrOfProcesses': (1, 1, 1),
        'periodic': (0, 0, 1),
        'dx': 1.0
    },
    'Parameters': {
        'omega': 1.8,
        'timesteps': 1001,
        'remainingTimeLoggerFrequency': 250,
        'vtkWriteFrequency': 0,
        'overlapCommunication': False,
        'cudaEnabledMPI': False,
        'initialVelocity': (0, 0, 0),
        'performanceReportFrequency': 250,
        'communicationScheme': CommunicationSchemeType['UniformGPUScheme_Baseline'],
    },
    'Boundaries': {
        'Border': [
            {'direction': 'W', 'walldistance': -1, 'flag': 'NoSlip'},
            {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
            {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
            {'direction': 'N', 'walldistance': -1, 'flag': 'UBB'},
        ]
    }
}


class BenchmarkScenario:
    def __init__(self, testcase, decomposition_axes=None):
        self.testcase = testcase
        self.scenario_config = copy.deepcopy(BASE_CONFIG)
        self.decomposition_axes = decomposition_axes

        now = datetime.now().replace(second=0, microsecond=0)
        self.output_filename = f'{self.testcase}_{now.strftime("%Y-%m-%d_%H-%M")}.csv'

    @wlb.member_callback
    def config(self, **kwargs):
        return self.scenario_config

    @wlb.member_callback
    def results_callback(self, **kwargs):
        block_setup = self.scenario_config.get('DomainSetup')
        params = self.scenario_config.get('Parameters')

        data = [{
            'processesX': block_setup.get('nrOfProcesses')[0],
            'processesY': block_setup.get('nrOfProcesses')[1],
            'processesZ': block_setup.get('nrOfProcesses')[2],
            'blocksX': block_setup.get('blocks')[0],
            'blocksY': block_setup.get('blocks')[1],
            'blocksZ': block_setup.get('blocks')[2],
            'cellsPerBlockX': block_setup.get('cellsPerBlock')[0],
            'cellsPerBlockY': block_setup.get('cellsPerBlock')[1],
            'cellsPerBlockZ': block_setup.get('cellsPerBlock')[2],
            'cudaEnabledMPI': params.get('cudaEnabledMPI'),
            'overlapCommunication': params.get('overlapCommunication'),
            'domainDecomposition': self.decomposition_axes,
            'communicationScheme': CommunicationSchemeName[params.get('communicationScheme')],
            'mlupsTotal': kwargs.get('mlups_total'),
            'mlupsProcess': kwargs.get('mlups_process'),
            'mflupsTotal': kwargs.get('mflups_total'),
            'mflupsProcess': kwargs.get('mflups_process'),
        }]

        self.save_data(data)

    def save_data(self, data):
        df = pd.DataFrame(data)
        if not os.path.isfile(self.output_filename):
            df.to_csv(self.output_filename, index=False)
        else:
            df.to_csv(self.output_filename, index=False, mode='a', header=False)

