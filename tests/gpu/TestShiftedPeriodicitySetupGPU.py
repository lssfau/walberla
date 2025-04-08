import waLBerla as wlb


class Scenario:
    def __init__(self, normal_dir, shift_dir, shift_value, periodicity, zero_centered_domain):
        self.normal_dir = normal_dir
        self.shift_dir = shift_dir
        self.shift_value = shift_value
        self.periodicity = tuple(periodicity)
        self.zero_centered_domain = zero_centered_domain

    @wlb.member_callback
    def config(self):

        return {
            'DomainSetup': {
                'blocks': (3, 3, 3),
                'cellsPerBlock': (4, 4, 4),
                'cartesianSetup': 0,
                'periodic': self.periodicity,
                'zeroCenteredDomain': self.zero_centered_domain
            },
            'Boundaries': {
                'ShiftedPeriodicity': {
                    'shiftDir': self.shift_dir,
                    'shiftValue': self.shift_value,
                    'normalDir': self.normal_dir
                }
            },
            'Logging': {
                'logLevel': "Info"
            }
        }


scenarios = wlb.ScenarioManager()

for normal_dir in (0, 1, 2):
    for shift_dir in (0, 1, 2):
        if normal_dir == shift_dir:
            continue
        periodicity = 3 * [0]
        periodicity[shift_dir] = 1
        for shift_value in (-3, 0, 2, 5, 8, 15):
            for zero_centered_domain in (True, False):
                scenarios.add(Scenario(normal_dir, shift_dir, shift_value, periodicity, zero_centered_domain))
