-- Communication - Layout
SELECT sweep, processes,cells, average*1e3,FZYX,MLUPS
 FROM runs JOIN timingPool ON runs.runId =timingPool.runId WHERE sweep LIKE "Communication%"


SELECT r1.processes, MAX(r1.MLUPS) AS "SRT MLUPS", MAX(r2.MLUPS) AS "TRT MLUPS"
FROM runs as r1, runs as r2 
WHERE    r1.SPLIT IS NOT NULL and r2.SPLIT IS NOT NULL
       AND r1.FZYX=1 AND r2.FZYX= 1
       AND r1.processes = r2.processes
       AND r1.TRT IS NULL AND r2.TRT IS NOT NULL
      GROUP BY r1.processes;
