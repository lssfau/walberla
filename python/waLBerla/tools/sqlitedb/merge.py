#!/usr/bin/env python

import sqlite3


def getColumnNames(db, tableName, dbName):
    """Returns list of columns for a sqlite3 table
    :param db       :  database connection
    :param dbName   : name of database
    :param tableName: name of table
    :return: list of column names for this table"""
    cursor = db.cursor()
    cursor.execute("PRAGMA %s.table_info(%s)" % (dbName, tableName))
    columns = cursor.fetchall()

    res = []
    for e in columns:
        res.append((e[1], e[2].upper()))

    return res


def mergeSqliteFiles(targetFile, fileToMerge):
    """Merges sqlite3 database 'fileToMerge' into 'targetFile' database.
       Works only if both tables have a 'runs' table.  Other tables are ignored!
       If the runs table in one of the databases has more columns, these columns are created in the merged database.
    """
    db = sqlite3.connect(targetFile)
    db.execute('ATTACH "' + fileToMerge + '" AS toMerge')

    targetColumns = getColumnNames(db, "runs", "main")
    toMergeColumns = getColumnNames(db, "runs", "toMerge")

    columnsToCreate = [e for e in toMergeColumns if e not in targetColumns]

    for column in columnsToCreate:
        print("Adding Column {} to run table of {} ".format(column[0], targetFile))
        db.execute("ALTER TABLE main.runs ADD COLUMN %s %s" % (column[0], column[1]))

    # Fetch all runs from toMerge,
    # check if an entry with same date exists, if not add the run and the timing pool entries
    # to the targetTable
    c = db.cursor()
    assert (toMergeColumns[0][0] == "runId")
    columns = [e[0] for e in toMergeColumns]
    columnString = ",".join(columns)
    columnStringNoRunId = ",".join(columns[1:])

    query = 'SELECT {} FROM toMerge.runs WHERE timestamp || " " || uuid NOT IN '.format(columnString, )
    query += '( SELECT timestamp || " " || uuid FROM main.runs )'

    # associated tables are tables that reference the runs table, having a first column of 'runId' which is a
    # foreign key to runs
    associatedTables = []
    associatedTablesColumnNames = []

    c.execute("SELECT name FROM sqlite_master WHERE type='table';")
    tableNames = [e[0] for e in c.fetchall()]
    for tableName in tableNames:
        if tableName == 'runs':
            continue

        mainColumns = getColumnNames(db, tableName, "main")
        toMergeColumns = getColumnNames(db, tableName, "toMerge")

        if mainColumns != toMergeColumns:
            print("Warning: Not merging associated table %s, since they have different columns." % (tableName,))
            continue

        columnNames = [e[0] for e in mainColumns]
        if columnNames[0] != "runId":
            print("Warning: Not merging table %s, since foreign key column 'runId' not found." % (tableName,))
            continue

        associatedTables.append(tableName)
        associatedTablesColumnNames.append(columnNames)

    mergedRuns = 0
    for run in c.execute(query):
        # Build up insert statement for 'runs' table
        questionMarkList = ['?'] * (len(run) - 1)
        questionMarkString = ",".join(questionMarkList)
        insertStatement = "INSERT INTO main.runs (%s) VALUES (%s);" % (columnStringNoRunId, questionMarkString)
        # Execute the insert
        insertCursor = db.cursor()
        insertCursor.execute(insertStatement, run[1:])
        insertedRunId = insertCursor.lastrowid
        originalRunId = run[0]

        # Insert the corresponding entries from associated tables
        for associatedTable, columnNames in zip(associatedTables, associatedTablesColumnNames):
            assocTableQuery = "SELECT %s FROM toMerge.%s WHERE runId=?" % (",".join(columnNames[1:]), associatedTable)
            assocTableInsertCursor = db.cursor()
            assocTableQueryCursor = db.cursor()

            for vals in assocTableQueryCursor.execute(assocTableQuery, (originalRunId,)):
                questionMarkList = ['?'] * len(columnNames)
                questionMarkString = ",".join(questionMarkList)
                insertQuery = "INSERT INTO main.%s (%s) VALUES (%s)" % (associatedTable, ",".join(columnNames),
                                                                        questionMarkString)
                assocTableInsertCursor.execute(insertQuery, (insertedRunId,) + vals)

        mergedRuns = mergedRuns + 1

    c.execute("SELECT COUNT(*) FROM toMerge.runs")
    totalRows = c.fetchall()[0][0]

    print("Merged {}/{} runs from {} to {} ".format(mergedRuns, totalRows, fileToMerge, targetFile))
    db.commit()
    db.close()


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    import re
    import shutil

    parser = argparse.ArgumentParser()
    parser.add_argument('output_file', type=Path)
    parser.add_argument('input_folder', type=Path)
    parser.add_argument('-r', '--recursive', action='store_true', help='Take all subfolders into account.')
    parser.add_argument('-f', '--filter', default='.*', dest='regex',
                        help='Only files matching this regex are accepted.')
    args = parser.parse_args()

    regex = re.compile(args.regex)

    if not args.input_folder.is_dir():
        raise ValueError(f'{args.input_folder} is not a folder')

    search_pattern = '**/*.sqlite' if args.recursive else '*.sqlite'

    files = filter(lambda x: regex.match(x.name), args.input_folder.glob(search_pattern))
    try:
        first_file = next(files)
        print(f'Copying {first_file}')
        shutil.copy(str(first_file.resolve()), str(args.output_file.resolve()))
    except StopIteration:
        raise ValueError(f'{args.input_folder} has no matching sqlite files')
    for file in files:
        mergeSqliteFiles(str(args.output_file.resolve()), str(file.resolve()))
