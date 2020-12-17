import waLBerla as wlb


class Scenario:
    def __init__(self):
        # output frequencies
        self.vtkWriteFrequency = 100

        # simulation parameters
        self.timesteps = 3
        edge_size = 200
        self.cells = (edge_size, edge_size, edge_size)
        self.blocks = (1, 1, 1)
        self.periodic = (1, 1, 1)
        self.size = (self.cells[0] * self.blocks[0],
                     self.cells[1] * self.blocks[1],
                     self.cells[2] * self.blocks[2])

        self.timeStepStrategy = 'phase_only'
        self.overlappingWidth = (1, 1, 1)
        self.cuda_block_size = (128, 1, 1)

        # bubble parameters
        self.bubbleRadius = min(self.size) // 4
        self.bubbleMidPoint = (self.size[0] / 2, self.size[1] / 2, self.size[2] / 2)

        self.scenario = 1  # 1 rising bubble, 2 RTI
        self.config_dict = self.config()

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'blocks': self.blocks,
                'cellsPerBlock': self.cells,
                'periodic': self.periodic,
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'useGui': 0,
                'remainingTimeLoggerFrequency': -1,
                'timeStepStrategy': self.timeStepStrategy,
                'overlappingWidth': self.overlappingWidth,
                'gpuBlockSize': self.cuda_block_size,
                'warmupSteps': 0,
                'scenario': self.scenario,
            },
            'Boundaries_GPU': {
                'Border': []
            },
            'Boundaries_CPU': {
                'Border': []
            },
            'Bubble': {
                'bubbleMidPoint': self.bubbleMidPoint,
                'bubbleRadius': self.bubbleRadius,
            },
        }


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario())
