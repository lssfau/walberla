import unittest
import os
import tempfile
import shutil

from waLBerla.tools.sqlitedb import checkAndUpdateSchema, storeSingle, storeMultiple


class SqliteDBTest(unittest.TestCase):

    def testSimpleInsertion(self):
        try:
            d = tempfile.mkdtemp()

            dbFile = os.path.join(d, "database.sqlite")
            myData = {'integerCol': 5, 'stringCol': 'someTestString', 'realCol': 3.141}

            checkAndUpdateSchema(myData, "sometable", dbFile=dbFile)
            runId = storeSingle(myData, 'sometable', dbFile=dbFile)

            valueData = {'x': [1, 2, 3, 4, 5], 'y': [3, 2, 5, 1, 0]}
            checkAndUpdateSchema(valueData, "data", dbFile=dbFile, referenceRuns=True)
            storeMultiple(valueData, "data", dbFile=dbFile, runId=runId)

        finally:
            shutil.rmtree(d)


if __name__ == '__main__':
    unittest.main()
