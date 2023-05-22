import waLBerla as wlb
from lbmpy.relaxationrates import relaxation_rate_from_lattice_viscosity


class Scenario:
    def __init__(self):
        self.timesteps = 10
        self.vtkWriteFrequency = 100

        self.cells = (64, 32, 32)
        self.blocks = (1, 1, 1)
        self.periodic = (0, 0, 0)

        self.constant_inflow = True

        self.diameter_sphere = min(self.cells) // 2
        self.u_max = 0.1
        self.reynolds_number = 1000000

        kinematic_viscosity = (self.diameter_sphere * self.u_max) / self.reynolds_number

        self.omega = relaxation_rate_from_lattice_viscosity(kinematic_viscosity)

        self.total_cells = (self.cells[0] * self.blocks[0],
                            self.cells[1] * self.blocks[1],
                            self.cells[2] * self.blocks[2])

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
                'oneBlockPerProcess': True
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'omega': self.omega,
                'u_max': self.u_max,
                'reynolds_number': self.reynolds_number,
                'diameter_sphere': self.diameter_sphere,
                'constant_inflow': self.constant_inflow
            },
            'Boundaries': {
                'Border': [
                    {'direction': 'N', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'S', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'W', 'walldistance': -1, 'flag': 'UBB'},
                    {'direction': 'E', 'walldistance': -1, 'flag': 'Outflow'},
                    {'direction': 'T', 'walldistance': -1, 'flag': 'NoSlip'},
                    {'direction': 'B', 'walldistance': -1, 'flag': 'NoSlip'},
                ],
                'Body': [
                    {'shape': "sphere",
                     'midpoint': (int(0.40 * self.total_cells[0]), self.total_cells[1] // 2, self.total_cells[2] // 2),
                     'radius': self.diameter_sphere // 2,
                     'flag': 'NoSlip'}
                ],
            },
        }


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
