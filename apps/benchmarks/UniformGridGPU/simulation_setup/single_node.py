# encoding: utf-8

import itertools
import waLBerla as wlb
from base import cells_per_block, num_processes
from benchmark import BenchmarkScenario


scenarios = wlb.ScenarioManager()

testcase_name = "single-node"

assert num_processes == 1

for num_cells_per_block in cells_per_block:
    # Create a benchmark scenario
    scenario = BenchmarkScenario(testcase=testcase_name)
    scenario.scenario_config['DomainSetup']['cellsPerBlock'] = 3 * (num_cells_per_block,)
    # Add scenario for execution
    scenarios.add(scenario)
