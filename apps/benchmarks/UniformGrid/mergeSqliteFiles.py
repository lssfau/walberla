#!/usr/bin/python

import sqlite3
import sys
import shutil


def getColumnNames(db, tableName, dbName):
    cursor = db.cursor()
    cursor.execute("PRAGMA %s.table_info(%s)" % (dbName, tableName))
    columns = cursor.fetchall()

    res = []
    for e in columns:
        res.append((e[1], e[2].upper()))

    return res


def mergeSqliteFiles(targetFile, fileToMerge):
    db = sqlite3.connect(targetFile)
    db.execute('ATTACH "' + fileToMerge + '" AS toMerge')

    targetColumns = getColumnNames(db, "runs", "main")
    toMergeColumns = getColumnNames(db, "runs", "toMerge")

    columnsToCreate = [e for e in toMergeColumns if e not in targetColumns]

    for column in columnsToCreate:
        print
        "Adding Column %s to run table of %s " % (column[0], targetFile)
        db.execute("ALTER TABLE main.runs ADD COLUMN %s %s" % (column[0], column[1]))

    # Fetch all runs from toMerge,
    # check if an entry with same date exists, if not add the run and the timing pool entries
    # to the targetTable
    c = db.cursor()
    assert (toMergeColumns[0][0] == "runId")
    columns = [e[0] for e in toMergeColumns]
    columnString = ",".join(columns)
    columnStringNoRunId = ",".join(columns[1:])

    query = 'SELECT %s FROM toMerge.runs WHERE timestamp || " " || random NOT IN ' % (columnString,)
    query += '( SELECT timestamp || " " || random FROM main.runs )'

    timingPoolColumnsMain = getColumnNames(db, "timingPool", "main")
    timingPoolColumnsToMerge = getColumnNames(db, "timingPool", "toMerge")
    assert (timingPoolColumnsMain == timingPoolColumnsToMerge)
    timingPoolColumnNames = [e[0] for e in timingPoolColumnsMain]
    assert (timingPoolColumnNames[0] == "runId")

    mergedRuns = 0
    for run in c.execute(query):
        # Build up insert statement for 'runs' table
        questionMarkList = ['?'] * (len(run) - 1)
        questionMarkString = ",".join(questionMarkList)
        insertStatement = "INSERT INTO main.runs (%s) VALUES (%s);" % (columnStringNoRunId, questionMarkString)
        # Execute the insert
        insertCursor = db.cursor()
        insertCursor.execute(insertStatement, run[1:])
        # Insert the corresponding timingPool infos
        insertedRunId = insertCursor.lastrowid
        originalRunId = run[0]

        timingPoolQuery = "SELECT %s FROM toMerge.timingPool WHERE runId=?" % (",".join(timingPoolColumnNames[1:]))
        timingPoolInsertCursor = db.cursor()
        timingPoolQueryCursor = db.cursor()

        for tp in timingPoolQueryCursor.execute(timingPoolQuery, (originalRunId,)):
            questionMarkList = ['?'] * len(timingPoolColumnNames)
            questionMarkString = ",".join(questionMarkList)
            insertQuery = "INSERT INTO main.timingPool (%s) VALUES (%s)" % (",".join(timingPoolColumnNames),
                                                                            questionMarkString)
            timingPoolInsertCursor.execute(insertQuery, (insertedRunId,) + tp)

        mergedRuns = mergedRuns + 1

    print("Merged %s runs from %s to %s " % (mergedRuns, fileToMerge, targetFile))
    db.commit()
    db.close()


if len(sys.argv) < 3:
    print("Usage: mergeSqliteFiles resultFile <filesToMerge>")
else:
    print("Copying " + sys.argv[2] + " to " + sys.argv[1])
    shutil.copy(sys.argv[2], sys.argv[1])
    for i in range(3, len(sys.argv)):
        print("Merging " + sys.argv[i])
        mergeSqliteFiles(sys.argv[1], sys.argv[i])
