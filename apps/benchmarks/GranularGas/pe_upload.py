import os
import time
import math
import random
import re
from influxdb import InfluxDBClient
from git import Repo


class Upload:
   def __init__(self):
      try:
         self.write_user_pw = os.environ["INFLUXDB_WRITE_USER"]
      except KeyError:
         import sys
         print('Password for the InfluxDB write_user was not set.\n',
               'See https://docs.gitlab.com/ee/ci/variables/#secret-variables', file=sys.stderr)
         exc_info = sys.exc_info()
         raise exc_info[0].with_traceback(exc_info[1], exc_info[2])

      self.client = InfluxDBClient('i10grafana.informatik.uni-erlangen.de', 8086,
                                   'pe', self.write_user_pw, 'pe')

   def process(self, filename, model, friction, sync, parallelization):
      with open(filename) as f:
          s = f.read()
      m = re.search('PUpS: (\S*)', s)

      json_body = [
          {
              'measurement': 'pe_benchmark',
              'tags': {
                  'host'            : os.uname()[1],
                  'image'           : os.environ["DOCKER_IMAGE_NAME"],
                  'model'           : model,
                  'friction'        : friction,
                  'sync'            : sync,
                  'parallelization' : parallelization
              },
              'time': int(time.time()),
              'fields': {'PUpS': float(m.group(1))}
          }
      ]
      print(float(m.group(1)))
      self.client.write_points(json_body, time_precision='s')

if __name__ == "__main__":
   up = Upload()
   up.process("GranularGas_DEM_NN.txt", "DEM", "Coulomb", "next neighbors", "8P1T")
   up.process("GranularGas_DEM_SO.txt", "DEM", "Coulomb", "shadow owners", "8P1T")
   up.process("GranularGas_HCSITS_NN_IFC.txt", "HCSITS", "InelasticFrictionlessContact", "next neighbors", "8P1T")
   up.process("GranularGas_HCSITS_NN_AICCBD.txt", "HCSITS", "ApproximateInelasticCoulombContactByDecoupling", "next neighbors", "8P1T")
   up.process("GranularGas_HCSITS_NN_ICCBD.txt", "HCSITS", "InelasticCoulombContactByDecoupling", "next neighbors", "8P1T")
   up.process("GranularGas_HCSITS_NN_IGMDC.txt", "HCSITS", "InelasticGeneralizedMaximumDissipationContact", "next neighbors", "8P1T")
   up.process("GranularGas_HCSITS_SO_IFC.txt", "HCSITS", "InelasticFrictionlessContact", "shadow owners", "8P1T")

