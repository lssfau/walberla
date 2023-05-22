import waLBerla as wlb


class Scenario:
    def __init__(self, domain_size=(64, 64, 64), root_blocks=(2, 2, 2),
                 cells_per_block=(32, 32, 32), refinement_depth=0):

        self.domain_size = domain_size
        self.root_blocks = root_blocks
        self.cells_per_block = cells_per_block
        self.refinement_depth = refinement_depth

        self.periodic = (0, 0, 0)

        self.config_dict = self.config(print_dict=False)

    @wlb.member_callback
    def config(self, print_dict=True):
        from pprint import pformat
        config_dict = {
            'DomainSetup': {
                'domainSize': self.domain_size,
                'rootBlocks': self.root_blocks,
                'cellsPerBlock': self.cells_per_block,
                'periodic': self.periodic
            },
            'Parameters': {
                'omega': 1.95,
                'timesteps': 1501,

                'refinementDepth': self.refinement_depth,
                'writeSetupForestAndReturn': False,
                'numProcesses': 1,

                'cudaEnabledMPI': False,
                'benchmarkKernelOnly': False,

                'remainingTimeLoggerFrequency': 3,

                'vtkWriteFrequency': 500,
            }
        }

        if print_dict and config_dict["Parameters"]["writeSetupForestAndReturn"] is False:
            wlb.log_info_on_root("Scenario:\n" + pformat(config_dict))
        return config_dict


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")

    domain_size = (64, 64, 64)
    cells_per_block = (32, 32, 32)

    root_blocks = tuple([d // c for d, c in zip(domain_size, cells_per_block)])

    scenarios = wlb.ScenarioManager()
    scenario = Scenario(domain_size=domain_size,
                        root_blocks=root_blocks,
                        cells_per_block=cells_per_block,
                        refinement_depth=1)
    scenarios.add(scenario)


validation_run()
