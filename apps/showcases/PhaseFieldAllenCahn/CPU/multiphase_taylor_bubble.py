import os
import waLBerla as wlb
import numpy as np
import pandas as pd

from waLBerla.core_extension import makeSlice
from waLBerla.tools.sqlitedb import sequenceValuesToScalars


class Scenario:
    def __init__(self):
        # physical parameters
        self.tube_diameter = 64

        self.bond_number = 100
        self.morton_number = 0.015

        self.relaxation_rate_heavy = 1.76
        self.density_heavy = 1
        self.density_ratio = 744
        self.dynamic_viscosity_ratio = 4236

        self.interface_width = 3
        self.mobility = 0.08

        self.relaxation_time_heavy = 1 / self.relaxation_rate_heavy
        self.kinematic_viscosity_heavy = 1 / 3 * (self.relaxation_time_heavy - 0.5)
        self.dynamic_viscosity_heavy = self.density_heavy * self.kinematic_viscosity_heavy

        self.density_light = self.density_heavy / self.density_ratio
        self.dynamic_viscosity_light = self.dynamic_viscosity_heavy / self.dynamic_viscosity_ratio
        self.kinematic_viscosity_light = self.dynamic_viscosity_light / self.density_light
        self.relaxation_time_light = 3 * self.kinematic_viscosity_light + 0.5

        self.surface_tension = (self.bond_number * self.dynamic_viscosity_heavy**4 / self.morton_number /
                                self.density_heavy**2 / self.tube_diameter**2)**0.5

        self.gravitational_acceleration = - (self.morton_number * self.density_heavy * self.surface_tension**3 /
                                             self.dynamic_viscosity_heavy**4)

        self.reference_time = (self.tube_diameter / abs(self.gravitational_acceleration))**0.5
        self.timesteps = 20 * self.reference_time

        self.vtkWriteFrequency = self.reference_time
        self.dbWriteFrequency = self.reference_time
        self.meshWriteFrequency = self.reference_time

        # simulation parameters
        self.size = (self.tube_diameter, self.tube_diameter * 10, self.tube_diameter)
        self.blocks = (1, 1, 1)
        self.cells = (self.size[0] // self.blocks[0], self.size[1] // self.blocks[1], self.size[2] // self.blocks[2])
        self.periodic = (0, 0, 0)
        self.inner_radius = self.tube_diameter

        self.center_x = self.size[0] / 2
        self.center_y = self.size[1] / 2
        self.center_z = self.size[2] / 2

        self.scenario = 4   # 1 rising bubble, 2 RTI, 3 drop, 4 taylor bubble set up

        self.counter = 0
        self.yPositions = []

        self.eccentricity_or_pipe_ratio = False  # if True eccentricity is conducted otherwise pipe ratio
        self.ratio = 0

        self.start_transition = (self.size[1] // 2) - 2 * self.tube_diameter
        self.length_transition = 4 * self.tube_diameter

        setup = "eccentricity" if self.eccentricity_or_pipe_ratio else "ratio"

        self.csv_file = f"Taylor_bubble_D_{self.tube_diameter}_C_{setup}_{self.ratio}_W_" \
                        f"{self.interface_width}_M_{self.mobility}.csv"

        d = self.tube_diameter / 2
        dh = self.tube_diameter - d

        resolution = self.tube_diameter / 128

        self.Donut_D = 0.1 * self.tube_diameter / resolution
        self.Donut_h = dh / 6
        self.DonutTime = 0.5 * (self.tube_diameter + d) / 2

        self.config_dict = self.config()

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'domainSize': self.size,
                'cellsPerBlock': self.cells,
                'diameter': self.tube_diameter,
                'periodic': self.periodic,
                'inner_radius': self.inner_radius,
                'ratio': self.ratio,
                'start_transition': self.start_transition,
                'length_transition': self.length_transition,
                'eccentricity_or_pipe_ration': self.eccentricity_or_pipe_ratio,
                'tube': True
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'dbWriteFrequency': self.dbWriteFrequency,
                'meshWriteFrequency': self.meshWriteFrequency,
                'remainingTimeLoggerFrequency': 60.0,
                'scenario': self.scenario,
            },
            'PhysicalParameters': {
                'density_liquid': self.density_heavy,
                'density_gas': self.density_light,
                'surface_tension': self.surface_tension,
                'mobility': self.mobility,
                'gravitational_acceleration': self.gravitational_acceleration,
                'relaxation_time_liquid': self.relaxation_time_heavy - 0.5,
                'relaxation_time_gas': self.relaxation_time_light - 0.5,
                'interface_thickness': self.interface_width
            },
            'Boundaries': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'W', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'E', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
            },
            'TaylorBubble': {
                'bubble_radius': 0.75 * 0.5 * self.tube_diameter,
                'initial_height': 1 * self.tube_diameter,
                'length': 3 * self.tube_diameter,
            }
        }

    @wlb.member_callback
    def at_end_of_time_step(self, blocks, **kwargs):
        t = kwargs["timeStep"]
        if t % self.dbWriteFrequency == 0:
            wlb_field = wlb.field.gather(blocks, 'phase', makeSlice[:, :, self.size[2] // 2])
            if wlb_field:
                data = {'timestep': t}
                data.update(self.config_dict['Parameters'])
                data.update(self.config_dict['DomainSetup'])
                data.update(self.config_dict['PhysicalParameters'])
                data.update(self.config_dict['TaylorBubble'])
                data.update(kwargs)

                phase_field = np.asarray(wlb_field).squeeze()
                location_of_gas = np.where(phase_field < 0.5)
                center_of_mass = np.mean(location_of_gas, axis=1)

                assert np.isfinite(np.sum(phase_field)), "NaN detected in the phasefield"

                self.yPositions.append(center_of_mass[1])
                if len(self.yPositions) > 1:
                    speed = self.yPositions[-1] - self.yPositions[-2]
                else:
                    speed = 0

                data['center_of_mass_x'] = center_of_mass[0]
                data['center_of_mass_y'] = center_of_mass[1]
                # data['center_of_mass_z'] = center_of_mass[2]

                data['xCells'] = self.size[0]
                data['yCells'] = self.size[1]
                data['zCells'] = self.size[2]

                data['rising_velocity'] = speed

                data['StencilHydroLB'] = kwargs["stencil_hydro"]
                data['StencilPhaseLB'] = kwargs["stencil_phase"]
                sequenceValuesToScalars(data)

                df = pd.DataFrame.from_records([data])
                if not os.path.isfile(self.csv_file):
                    df.to_csv(self.csv_file, index=False)
                else:
                    df.to_csv(self.csv_file, index=False, mode='a', header=False)


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
