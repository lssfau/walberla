#! /usr/bin/env python3
# -*- coding: utf-8 -*-

from FilmSpecimenGenerator import Config

cfg = Config()
cfg.add_parameter("sizeX",                    "real_t",     "0",     "Specimen length (x direction)")
cfg.add_parameter("sizeY",                    "real_t",     "0",     "Specimen length (y direction)")
cfg.add_parameter("sizeZ",                    "real_t",     "0",     "Specimen length (z direction)")
cfg.add_parameter("oopp",                     "bool",       "false", "Specimen out-of-plane periodicity (0 - film, 1 - cnt material)")
cfg.add_parameter("numBlocksX",               "uint_t",     "1",     "Number of blocks in x direction")
cfg.add_parameter("numBlocksY",               "uint_t",     "1",     "Number of blocks in y direction")
cfg.add_parameter("numBlocksZ",               "uint_t",     "1",     "Number of blocks in z direction")
cfg.add_parameter("min_OOP",                  "real_t",     "0",     "Out-of-plane angle minimum")
cfg.add_parameter("max_OOP",                  "real_t",     "0",     "Out-of-plane angle maximum")
cfg.add_parameter("numCNTs",                  "int",        "0",     "Number of CNTs")
cfg.add_parameter("numSegs",                  "int",        "0",     "Number of segments in a CNT")
cfg.add_parameter("spacing",                  "real_t",     "0",     "Segment half-spacing")
cfg.add_parameter("localDamping",             "real_t",     "0",     "Local damping coefficient")
cfg.add_parameter("viscousDamping",           "real_t",     "0",     "Viscous damping coefficient")
cfg.add_parameter("seed",                     "uint_t",     "0",     "random generator seed")
cfg.add_parameter("vdW",                      "int",        "0",     "type of vdW interaction model")
cfg.add_parameter("simulationSteps",          "int",        "0",     "Relaxation duration")
cfg.add_parameter("saveVTKEveryNthStep",      "int",        "0",     "timesteps between saving VTK outputs")
cfg.add_parameter("saveEnergyEveryNthStep",   "int",        "0",     "timesteps between saving energies")
cfg.add_parameter("saveConfEveryNthStep",     "int",        "0",     "timesteps between saving confs")
cfg.add_parameter("vtkFolder",                "std::string", '"."',   "Folder for VTK files")
cfg.add_parameter("energyFolder",             "std::string", '"."',   "Folder for energy files")
cfg.add_parameter("confFolder",               "std::string", '"."',   "Folder for conf files")
cfg.add_parameter("useMPIIO",                 "bool",       "false", "Write a single file instead of one file per process")
cfg.add_parameter("sqlFile",                  "std::string", '"cnt.sqlite"',   "database file")
cfg.add_parameter("initialConfigurationFile", "std::string", '""',   "restart from checkpoint")

cfg.generate()


