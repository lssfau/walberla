#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from ConfigGenerator import Config

cfg = Config()
cfg.addParameter("sorting",                "std::string", '"none"')
cfg.addParameter("normal",                 "Vec3",        "Vec3(real_t(1.0), real_t(1.0), real_t(1.0))")
cfg.addParameter("spacing",                "real_t",      "real_t(1.0)")
cfg.addParameter("shift",                  "Vec3",        "Vec3(real_t(0.1), real_t(0.1), real_t(0.1))")
cfg.addParameter("radius",                 "real_t",      "real_t(0.5)")
cfg.addParameter("bBarrier",               "bool",        "false")
cfg.addParameter("storeNodeTimings",       "bool",        "false")
cfg.addParameter("checkSimulation",        "bool",        "false")
cfg.addParameter("numOuterIterations",     "int64_t",     "10")
cfg.addParameter("initialRefinementLevel", "int64_t",     "0")
cfg.addParameter("simulationSteps",        "int64_t",     "10")
cfg.addParameter("dt",                     "real_t",      "real_t(0.01)")
cfg.addParameter("visSpacing",             "int64_t",     "1000")
cfg.addParameter("vtk_out",                "std::string", '"vtk_out"')
cfg.addParameter("sqlFile",                "std::string", '"benchmark.sqlite"')

cfg.addParameter("recalculateBlockLevelsInRefresh",                "bool", "false")
cfg.addParameter("alwaysRebalanceInRefresh",                       "bool", "true")
cfg.addParameter("reevaluateMinTargetLevelsAfterForcedRefinement", "bool", "false")
cfg.addParameter("allowRefreshChangingDepth",                      "bool", "false")

cfg.addParameter("allowMultipleRefreshCycles",                     "bool", "false")
cfg.addParameter("checkForEarlyOutInRefresh",                      "bool", "true")
cfg.addParameter("checkForLateOutInRefresh",                       "bool", "true")

cfg.addParameter("regridMin",              "uint_t",      'uint_c(100)')
cfg.addParameter("regridMax",              "uint_t",      'uint_c(1000)')
cfg.addParameter("maxBlocksPerProcess",    "int",         'int_c(1000)')
cfg.addParameter("baseWeight",             "real_t",      'real_t(10.0)')
cfg.addParameter("metisipc2redist",        "real_t",      'real_t(1000.0)')
cfg.addParameter("LBAlgorithm",            "std::string", '"Hilbert"')
cfg.addParameter("metisAlgorithm",         "std::string", '"PART_GEOM_KWAY"')
cfg.addParameter("metisWeightsToUse",      "std::string", '"BOTH_WEIGHTS"')
cfg.addParameter("metisEdgeSource",        "std::string", '"EDGES_FROM_EDGE_WEIGHTS"')

cfg.generate()
