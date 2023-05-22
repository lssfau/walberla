import waLBerla as wlb


class Scenario:
    def __init__(self, domain_size=(32, 32, 32), root_blocks=(2, 2, 2),
                 cells_per_block=(16, 16, 16)):

        self.domain_size = domain_size
        self.root_blocks = root_blocks
        self.cells_per_block = cells_per_block

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
                'timesteps': 101,

                'refinementDepth': 1,
                'writeSetupForestAndReturn': False,
                'numProcesses': 1,

                'benchmarkKernelOnly': False,

                'remainingTimeLoggerFrequency': 3,

                'vtkWriteFrequency': 50,
            }
        }

        if print_dict:
            wlb.log_info_on_root("Scenario:\n" + pformat(config_dict))
        return config_dict


def validation_run():
    """Run with full periodic shear flow or boundary scenario (ldc) to check if the code works"""
    wlb.log_info_on_root("Validation run")
    wlb.log_info_on_root("")

    scenarios = wlb.ScenarioManager()
    scenario = Scenario()
    scenarios.add(scenario)


validation_run()
