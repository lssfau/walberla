import waLBerla as wlb


class Scenario:
    def __init__(self, cuda_enabled_mpi=False):
        # output frequencies
        self.vtkWriteFrequency = 1000

        # simulation parameters
        self.timesteps = 10000

        # domain decomposition can be specified manually by specifying the number of cells per block and the
        # number of blocks. The number of blocks must be equal to the MPI processes used. If only the total domain size
        # is specified with 'cells' waLBerla will take care of the decomposition depending on the number of MPI
        # processes at runtime

        # self.cell_per_block = (32, 64, 32)
        # self.blocks = (4, 1, 4)
        self.cells = (128, 64, 128)

        self.periodic = (0, 0, 0)

        # bubble parameters
        self.dropletRadius = 24.0
        self.dropletMidPoint = (64, 24, 64)

        # everything else
        self.scenario = 1  # 1 rising bubble or droplet, 2 RTI, 3 bubble field, 4 taylor bubble set up

        self.cudaEnabledMpi = cuda_enabled_mpi
        self.cuda_blocks = (128, 1, 1)

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                # 'blocks': self.blocks,
                # 'cellsPerBlock': self.cell_per_block,
                'cells': self.cells,
                'periodic': self.periodic,
                'tube': False
            },
            'Parameters': {
                'timesteps': self.timesteps,
                'vtkWriteFrequency': self.vtkWriteFrequency,
                'remainingTimeLoggerFrequency': 10.0,
                'scenario': self.scenario,
                'cudaEnabledMpi': self.cudaEnabledMpi,
                'gpuBlockSize': self.cuda_blocks
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
scenarios.add(Scenario(cuda_enabled_mpi=False))
