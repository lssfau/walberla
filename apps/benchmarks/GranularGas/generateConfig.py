#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from ConfigGenerator import Config

cfg = Config()
cfg.addParameter("sorting",                "std::string", '"none"')
cfg.addParameter("spacing",                "real_t",      "real_t(1.0)")
cfg.addParameter("radius",                 "real_t",      "real_t(0.5)")
cfg.addParameter("bBarrier",               "bool",        "false")
cfg.addParameter("storeNodeTimings",       "bool",        "false")
cfg.addParameter("checkSimulation",        "bool",        "false")
cfg.addParameter("numOuterIterations",     "int64_t",     "10")
cfg.addParameter("initialRefinementLevel", "int64_t",     "0")
cfg.addParameter("simulationSteps",        "int64_t",     "10")
cfg.addParameter("dt",                     "real_t",      "real_t(0.01)")
cfg.addParameter("visSpacing",             "int64_t",     "1000")
cfg.addParameter("path",                   "std::string", '"vtk_out"')
cfg.addParameter("sqlFile",                "std::string", '"benchmark.sqlite"')

cfg.generate()
