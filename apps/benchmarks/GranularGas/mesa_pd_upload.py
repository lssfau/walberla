#! /usr/bin/env python3

import os
import time
import re
from influxdb import InfluxDBClient


class Upload:
    def __init__(self):
        try:
            self.write_user_pw = os.environ["INFLUXDB_MESAPD_PW"]
        except KeyError:
            import sys
            print('Password for the InfluxDB write_user was not set.\n',
                  'See https://docs.gitlab.com/ee/ci/variables/#secret-variables', file=sys.stderr)
            exc_info = sys.exc_info()
            raise exc_info[0].with_traceback(exc_info[1], exc_info[2])

        self.client = InfluxDBClient('i10grafana.informatik.uni-erlangen.de', 8086,
                                     'mesa_pd', self.write_user_pw, 'mesa_pd')

    def process(self, filename):
        tts = dict()

        with open(filename) as f:
            for s in f.readlines():
                m = re.search(r'\[0\]\s*(\w*)\s*\|[\s\d\.\%]*\|\s*([\d\.]*)', s)

                if m is not None:
                    tts[m.group(1)] = float(m.group(2))

        json_body = [
            {
                'measurement': 'mesa_pd_benchmark',
                'tags': {
                    'host': os.uname()[1],
                    'image': os.environ["DOCKER_IMAGE_NAME"],
                },
                'time': int(time.time()),
                'fields': tts
            }
        ]
        print(tts)
        self.client.write_points(json_body, time_precision='s')


if __name__ == "__main__":
    up = Upload()
    up.process("mesa_pd.txt")
