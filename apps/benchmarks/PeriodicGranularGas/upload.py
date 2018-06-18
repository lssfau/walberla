import os
import time
import math
import random
import re
from influxdb import InfluxDBClient
from git import Repo


def main():
    try:
        write_user_pw = os.environ["INFLUXDB_WRITE_USER"]
    except KeyError:
        import sys
        print('Password for the InfluxDB write_user was not set.\n',
              'See https://docs.gitlab.com/ee/ci/variables/#secret-variables', file=sys.stderr)
        exc_info = sys.exc_info()
        raise exc_info[0].with_traceback(exc_info[1], exc_info[2])

    client = InfluxDBClient('i10grafana.informatik.uni-erlangen.de', 8086,
                            'pe', write_user_pw, 'pe')

    #repo = Repo(search_parent_directories=True)
    #commit = repo.head.commit

    with open("PeriodicGranularGas_DEM_NN.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'DEM',
                'friction': 'Coulomb',
                'sync'    : 'next neighbor'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')

    #*************************************************

    with open("PeriodicGranularGas_DEM_SO.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'DEM',
                'friction': 'Coulomb',
                'sync'    : 'shadow owner'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')

    #*************************************************

    with open("PeriodicGranularGas_HCSITS_NN_IFC.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'HCSITS',
                'friction': 'InelasticFrictionlessContact',
                'sync'    : 'next neighbor'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')

    #*************************************************

    with open("PeriodicGranularGas_HCSITS_NN_AICCBD.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'HCSITS',
                'friction': 'ApproximateInelasticCoulombContactByDecoupling',
                'sync'    : 'next neighbor'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')

    #*************************************************

    with open("PeriodicGranularGas_HCSITS_NN_ICCBD.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'HCSITS',
                'friction': 'InelasticCoulombContactByDecoupling',
                'sync'    : 'next neighbor'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')

    #*************************************************

    with open("PeriodicGranularGas_HCSITS_NN_IGMDC.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'HCSITS',
                'friction': 'InelasticGeneralizedMaximumDissipationContact',
                'sync'    : 'next neighbor'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')

    #*************************************************

    with open("PeriodicGranularGas_HCSITS_SO_IFC.txt") as f:
        s = f.read()
    m = re.search('runtime: (\d*.\d*)', s)

    json_body = [
        {
            'measurement': 'pe_benchmark',
            'tags': {
                'host'    : os.uname()[1],
                'image'   : os.environ["DOCKER_IMAGE_NAME"],
                'model'   : 'HCSITS',
                'friction': 'InelasticFrictionlessContact',
                'sync'    : 'shadow owner'
            },
            'time': int(time.time()),
            'fields': {'runtime': float(m.group(1))}
        }
    ]
    print(float(m.group(1)))
    client.write_points(json_body, time_precision='s')


if __name__ == "__main__":
    main()
