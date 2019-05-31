import waLBerla
import subprocess
import re
from collections import defaultdict
import sys
from datetime import timedelta

base = (32, 16, 2, 64)
BlOCK_SIZES_SQ = [(i, i, i) for i in base]
BLOCK_SIZES_RECT = [(i, i, i // 2) for i in base] + [(i, i // 2, i // 2) for i in base]

time_for_benchmark = 0.25
outer_iterations = 2
db_file = 'FieldCommunication.sqlite'


def supermuc_network_spread():
    try:
        node_list = subprocess.check_output("scontrol show hostname $SLURM_JOB_NODELIST", shell=True, encoding='utf8')
    except subprocess.CalledProcessError:
        return defaultdict(lambda: 0)

    spread = defaultdict(set)
    for s in node_list.split("\n"):
        m = re.search("i(\d\d)r(\d\d)c(\d\d)s(\d\d)", s)
        if m:
            for name, idx in zip(['i', 'r', 'c', 's'], range(1, 5)):
                spread[name].add(m.group(idx))
    return {k: len(v) for k, v in spread.items()}


sng_network = supermuc_network_spread()


class AlreadySimulated:

    def __init__(self, db_file, properties=('processes0*processes1*processes2', 'layout', 'ghostLayers', 'cartesianCommunicator', 'stencil',
                                            'cellsPerBlock0', 'cellsPerBlock1', 'cellsPerBlock2',
                                            'blocksPerProcess', 'localCommunicationMode', 'singleMessage',
                                            'fieldsPdf', 'fieldsPdfOpt', 'fieldsVector', 'fieldsScalar',
                                            'buffered')):
        self.properties = properties
        import sqlite3
        conn = sqlite3.connect(db_file)
        self.data = set()

        try:
            for row in conn.execute("SELECT {} FROM runs;".format(",".join(self.properties))):
                self.data.add(row)
        except sqlite3.OperationalError:
            pass
        waLBerla.log_info_on_root("Loaded {} scenarios".format(len(self.data)))

    def in_db(self, args):
        return args in self.data


@waLBerla.callback("config")
def config(processes=None):
    simulated_db = AlreadySimulated(db_file)
    isWalberlaRun = processes is None
    skipped = 0
    simulated = 0
    if isWalberlaRun:
        processes = waLBerla.mpi.numProcesses()
    for layout in ('fzyx', 'zyxf'):
        for ghost_layers in (1, 2):
            for cartesian_comm in (False, True):
                for stencil in ('D3Q19', 'D3Q27', 'D3Q7'):
                    for cells in BlOCK_SIZES_SQ + BLOCK_SIZES_RECT:
                        for blocksPerProcess in (1, 2, 4, 8, 16):
                            for local_comm in ('start', 'noOptimization', 'buffer'):
                                for single_message in (True, False):
                                    for pdf, pdf_opt, vector, scalar in ([1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1],
                                                                         [0, 0, 0, 3], [0, 0, 0, 19], [2, 0, 0, 0], [0, 2, 0, 0]):
                                        for buffered in (0, 1):
                                            if blocksPerProcess >= 8 and cells[0] >= 64 and cells[1] >= 64 and cells[2] >= 64:
                                                continue

                                            data = (processes, layout, ghost_layers, int(cartesian_comm), stencil, *cells, blocksPerProcess, local_comm,
                                                    int(single_message), pdf, pdf_opt, vector, scalar, buffered)

                                            if simulated_db.in_db(data):
                                                skipped += 1
                                                if skipped % 100 == 0 and isWalberlaRun:
                                                    waLBerla.log_info_on_root("Skipped {} scenarios".format(skipped))
                                                continue
                                            else:
                                                simulated += 1
                                                if not isWalberlaRun:
                                                    continue

                                            cfg = {
                                                'Domain': {
                                                    'cellsPerBlock': cells,
                                                    'domainWeights': (1, 1, 1),
                                                    'blocksPerProcess': blocksPerProcess,
                                                },
                                                'Communication': {
                                                    'buffered': buffered,
                                                    'stencil': stencil,
                                                    'ghostLayers': ghost_layers,
                                                    'cartesianCommunicator': cartesian_comm,
                                                    'singleMessage': single_message,
                                                    'Fields': {
                                                        'pdf': pdf,
                                                        'pdfOpt': pdf_opt,
                                                        'vector': vector,
                                                        'scalar': scalar,
                                                    },
                                                    'layout': layout,
                                                    'localCommunicationMode': local_comm,
                                                },
                                                'Run': {
                                                    'warmupIterations': 3,
                                                    'iterations': 100,
                                                    'outerIterations': outer_iterations,
                                                    'databaseFile': db_file,
                                                    'timeForBenchmark': time_for_benchmark,
                                                    'minIterations': 2,
                                                    'maxIterations': 10000,
                                                },
                                                'Database': {
                                                    'sngNetworkIslands': sng_network['i'],
                                                    'sngNetworkRows': sng_network['r'],
                                                    'sngNetworkCabinets': sng_network['c'],
                                                    'sngNetworkSlots': sng_network['s'],
                                                }
                                            }
                                            yield cfg
    if not isWalberlaRun:
        print("Skipped", skipped, "to simulate", simulated)
        estimated_seconds = simulated * time_for_benchmark * outer_iterations
        print("Estimated time ", timedelta(seconds=estimated_seconds))


if __name__ == '__main__':
    for _ in config(int(sys.argv[1])):
        pass
