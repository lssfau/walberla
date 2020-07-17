"""Small wrapper around Pythons sqlite3 module - simplifies table creation and is mainly tailored to storage of
   aggregated simulation results"""

import sqlite3


def sequenceValuesToScalars(data):
    """Transforms sequence values in tht data dictionary into multiple entries

       Example: { 'domainSize' : (4,5,6) } is transformed int
                { 'domainSize_0' : 4, 'domainSize_1' : 5, 'domainSize_2' : 6 }

       This is useful when using the dictionary with storeSingle(), since each entry gets a separate column
       after the sequences have been separated.
     """
    keysToDelete = []
    newValues = {}
    for key, value in data.items():
        if type(value) in [list, tuple]:
            keysToDelete.append(key)
            for i in range(len(value)):
                newValues["%s_%d" % (key, i)] = value[i]
    for k in keysToDelete:
        del data[k]
    data.update(newValues)


def storeSingle(data, tableName, dbFile="database.sqlite", runId=None):
    """Stores results of a single simulation run to a sqlite database.

    Primary key column 'runId' and a timestamp are automatically added.

        :param data:       Dictionary with data to store in sqlite database.
                           The dictionary keys are assumed to be column names.
        :param tableName: name of sqlite table
        :param dbFile:    database file
        :param runId:     override the otherwise manually created runId which serves as primary key
        :return: runId of inserted run
    """
    if runId:
        data['runId'] = runId

    keyString = ",".join(data.keys())
    valueString = ",".join(["?" for e in data.values()])
    query = "INSERT INTO %s ( %s ) VALUES ( %s )" % (tableName, keyString, valueString)

    conn = sqlite3.connect(dbFile)
    c = conn.cursor()
    c.execute(query, list(data.values()))
    lastrowid = c.lastrowid
    conn.commit()
    conn.close()

    return lastrowid


def storeMultiple(data, tableName, dbFile="database.sqlite", runId=None):
    """Stores data from multiple simulation runs into the database.

       :param data:  similar to storeSingle, but each value of the dictionary has to be a list
                     for all keys the values have to be lists of equal length
       :param runId: leave this to none if runId is the primary key of this table
                     if this table references a table with primary key runId, you can reference a
                     specific run by passing here the id of the run to reference.
       :return: primary key of last inserted run, the other runs have runId: returnValue-1, returnValue-2, ...
    """
    # Check if data is a dictionary that only has lists as values, all these lists have to have the
    # same length
    list_length = -1
    for key, value in data.items():
        assert (type(value) is list)
        if list_length < 0:
            list_length = len(value)
        assert (len(value) == list_length)

    if runId is not None:
        data.update({'runId': [runId] * list_length})

    keyString = ",".join(data.keys())
    valueString = ",".join(["?" for e in data.values()])
    query = "INSERT INTO %s ( %s ) VALUES ( %s )" % (tableName, keyString, valueString)

    conn = sqlite3.connect(dbFile)
    c = conn.cursor()

    sql_data = [tuple([v[i] for v in data.values()]) for i in range(list_length)]

    c.executemany(query, sql_data)
    lastrowid = c.lastrowid
    conn.commit()
    conn.close()

    return lastrowid


def checkAndUpdateSchema(data, tableName, dbFile="database.sqlite", referenceRuns=False):
    """Alters a sqlite table in order to match the given data:

        * if table with given name does not exist yet, it is created
        * keys in the data dictionary correspond to columns
        * columns are added if necessary, existing data has NULL in these new columns

    :param data:          see :func:`storeSingle` or :func:`storeMultiple`
    :param referenceRuns: if False the table gets an autoincrementing integer column 'runId'
                          if True, a normal column runId is created that points into
                          another table with runId as primary key
     """

    def pythonToSqlType(python_type):
        if python_type is int:
            return ("INTEGER")
        elif python_type is bool:
            return ("INTEGER")
        elif python_type is float:
            return ("DOUBLE")
        elif python_type is str:
            return ("TEXT")
        elif python_type is list:
            return (pythonToSqlType(type(value[0])))
        elif python_type is tuple:
            return (pythonToSqlType(type(value[0])))

    if referenceRuns:
        names = ["runId"]
        types = ["INTEGER"]
    else:
        names = ["id", "timestamp"]
        types = ["INTEGER PRIMARY KEY", "DATETIME DEFAULT  (datetime('now','localtime'))"]

    for key, value in data.items():
        names.append(key)
        types.append(pythonToSqlType(type(value)))

    columns = [("%s %s") % e for e in zip(names, types)]

    createQuery = "CREATE TABLE IF NOT EXISTS %s ( %s );" % (tableName, ",".join(columns))
    alterQueries = ["ALTER TABLE %s ADD COLUMN %s %s;" % (tableName, key, typ) for key, typ in zip(names, types)]

    conn = sqlite3.connect(dbFile)
    c = conn.cursor()

    c.execute(createQuery)
    for q in alterQueries:
        try:
            c.execute(q)
        except sqlite3.OperationalError as e:
            print(e)
            pass

    conn.commit()
    conn.close()
