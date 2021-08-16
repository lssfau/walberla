import os

import waLBerla as wlb
from waLBerla.tools.sqlitedb import sequenceValuesToScalars
from waLBerla.core_extension import makeSlice

import numpy as np
import pandas as pd
from lbmpy.phasefield_allen_cahn.parameter_calculation import calculate_parameters_rti


class Scenario:
    def __init__(self, cuda_enabled_mpi=False):
        # output frequencies
        self.vtkWriteFrequency = 1000

        # simulation parameters
        self.time = 2  # physical time in seconds

        self.cells = (128, 512, 128)
        self.blocks = (1, 1, 1)
        self.periodic = (1, 0, 1)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        # physical parameters
        self.density_heavy = 1.0
        self.reference_time = 4000
        self.dbWriteFrequency = self.reference_time // 20
        self.timesteps = int(self.reference_time * self.time) + 1

        self.capillary_number = 8.7
        self.reynolds_number = 3000
        self.atwood_number = 1
        self.peclet_number = 744
        self.density_ratio = 1000
        self.viscosity_ratio = 100

        self.parameters = calculate_parameters_rti(reference_length=self.cells[0],
                                                   reference_time=self.reference_time,
                                                   density_heavy=self.density_heavy,
                                                   capillary_number=self.capillary_number,
                                                   reynolds_number=self.reynolds_number,
                                                   atwood_number=self.atwood_number,
                                                   peclet_number=self.peclet_number,
                                                   density_ratio=self.density_ratio,
                                                   viscosity_ratio=self.viscosity_ratio)

        self.interface_thickness = 5
        self.tube = False

        # everything else
        self.scenario = 2  # 1 rising bubble or droplet, 2 RTI, 3 bubble field, 4 taylor bubble set up

        self.counter = 0
        self.yPositions = []

        self.cudaEnabledMpi = cuda_enabled_mpi
        self.cuda_blocks = (64, 2, 2)

        self.config_dict = self.config()

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
                'tube': self.tube
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'dbWriteFrequency': self.dbWriteFrequency,
                'remainingTimeLoggerFrequency': 10.0,
                'scenario': self.scenario,
                'cudaEnabledMpi': self.cudaEnabledMpi,
                'gpuBlockSize': self.cuda_blocks
            },
            'PhysicalParameters': {
                'density_liquid': self.density_heavy,
                'density_gas': self.parameters.get("density_light"),
                'surface_tension': self.parameters.get("surface_tension"),
                'mobility': self.parameters.get("mobility"),
                'gravitational_acceleration': self.parameters.get("gravitational_acceleration"),
                'relaxation_time_liquid': self.parameters.get("relaxation_time_heavy"),
                'relaxation_time_gas': self.parameters.get("relaxation_time_light"),
                'interface_thickness': self.interface_thickness
            },
            'Boundaries': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                ]
            },
        }

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, **kwargs):
        t = kwargs['timeStep']
        if t % self.dbWriteFrequency == 0:
            phase = wlb.field.gather(blocks, 'phase', makeSlice[:, :, :])
            if phase:
                data = {'timestep': t}
                data.update(self.config_dict['PhysicalParameters'])
                data.update({'total_timesteps': self.timesteps})
                data.update({'normalized_time': t / self.reference_time})
                data.update({'tube_setup': self.tube})
                data.update({'interface_thickness': self.interface_thickness})
                data.update({'capillary_number': self.capillary_number})
                data.update({'reynolds_number': self.reynolds_number})
                data.update({'atwood_number': self.atwood_number})
                data.update({'peclet_number': self.peclet_number})
                data.update({'density_ratio': self.density_ratio})
                data.update({'viscosity_ratio': self.viscosity_ratio})
                data.update({'reference_time': self.reference_time})
                data.update(kwargs)

                phase_field = np.asarray(phase).squeeze()
                stable = np.isfinite(np.sum(phase_field))
                mass = np.sum(phase_field)
                rel_max = np.max(phase_field) - 1
                rel_min = np.min(phase_field)
                data.update({'mass': mass})
                data.update({'rel_max': rel_max})
                data.update({'rel_min': rel_min})
                data.update({'stable': stable})

                if self.tube:
                    location_of_spike = self.get_interface_location(
                        phase_field[self.size[0] // 2, :, self.size[2] // 2])
                    a = np.where(phase_field < 0.5)
                    value = np.argmax(a[1])
                    location_of_bubble = self.get_interface_location(
                        phase_field[a[0][value], a[1][value] - 10:a[1][value] + 10, a[2][value]], a[1][value] - 10)

                    data.update({'location_of_spike': location_of_spike})
                    data.update({'location_of_bubble': location_of_bubble})
                else:
                    location_of_spike = self.get_interface_location(
                        phase_field[self.size[0] // 2, :, self.size[2] // 2])
                    location_of_bubble = self.get_interface_location(phase_field[0, :, 0])
                    location_of_saddle = self.get_interface_location(phase_field[0, :, self.size[2] // 2])

                    data.update({'location_of_spike': location_of_spike})
                    data.update({'location_of_bubble': location_of_bubble})
                    data.update({'location_of_saddle': location_of_saddle})

                sequenceValuesToScalars(data)

                csv_file = f"RTI_{data['stencil_phase']}_{data['stencil_hydro']}_Re_{self.reynolds_number}"
                csv_file += "_tube.csv" if self.tube else ".csv"

                df = pd.DataFrame.from_records([data])
                if not os.path.isfile(csv_file):
                    df.to_csv(csv_file, index=False)
                else:
                    df.to_csv(csv_file, index=False, mode='a', header=False)

    def get_interface_location(self, one_dimensional_array, shift=None):
        ny = self.size[1]
        l0 = self.size[0]

        index = np.argmax(one_dimensional_array > 0.5)

        if index > 0:
            zw1 = one_dimensional_array[index]
            zw2 = one_dimensional_array[index - 1]
            absolute_location = (index - 1) + (zw2 - 0.5) / (zw2 - zw1)
            if shift:
                absolute_location += shift
            return (absolute_location - ny // 2) / l0
        else:
            return -100


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
