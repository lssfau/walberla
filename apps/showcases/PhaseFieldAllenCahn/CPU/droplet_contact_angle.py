import waLBerla as wlb


class Scenario:
    def __init__(self):
        # output frequencies
        self.vtkWriteFrequency = 1000

        # simulation parameters
        self.timesteps = 10001
        self.cells = (32, 64, 32)
        self.blocks = (4, 1, 4)
        self.periodic = (0, 0, 0)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        self.overlappingWidth = (8, 1, 1)
        self.timeStepStrategy = 'normal'

        # bubble parameters
        self.dropletRadius = 24.0
        self.dropletMidPoint = (64, 24, 64)

        # everything else
        self.scenario = 1  # 1 rising bubble or droplet, 2 RTI, 3 bubble field, 4 taylor bubble set up

        self.counter = 0
        self.yPositions = []

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
                'tube': False
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'timeStepStrategy': self.timeStepStrategy,
                'overlappingWidth': self.overlappingWidth,
                'remainingTimeLoggerFrequency': 10.0,
                'scenario': self.scenario,
            },
            'PhysicalParameters': {
                'density_liquid': 1.0,
                'density_gas': 0.001,
                'surface_tension': 5e-5,
                'mobility': 0.05,
                'gravitational_acceleration': 0.0,
                'relaxation_time_liquid': 3 * 0.166,
                'relaxation_time_gas': 3 * 0.0166,
                'interface_thickness': 5
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
            'Bubble': {
                'bubbleMidPoint': self.dropletMidPoint,
                'bubbleRadius': self.dropletRadius,
                'bubble': False  # this means we are simulating a droplet rather than a bubble
            },
        }


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
