"""Dynamically generate a job matrix and CMake preset chains
for waLBerla continuous integration and testing.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, asdict, field
import json
import yaml
from pathlib import Path
from typing import Any

IMAGE_PATTERN = "$CI_CONTAINER_REPOSITORY:i10ci-walberla-{compiler_id}"
SCRIPT_DIR = Path(__file__).parent


class Reference:
    """GitLab CI YAML reference Tag recreation."""

    def __init__(self, job_name: str, section_name: str):
        self._job_name: str = job_name
        self._section_name: str = section_name


def reference_representer(dumper: yaml.SafeDumper, reference: Reference) -> yaml.nodes.SequenceNode:
    """Represent an empty GitLab CI YAML reference Tag."""
    return dumper.represent_sequence("!reference", [reference._job_name, reference._section_name], flow_style=True)


def get_dumper() -> yaml.SafeDumper:
    """Get a YAML SafeDumper with custom representers registered."""
    safe_dumper = yaml.SafeDumper
    safe_dumper.default_flow_style = (
        True  # Enable default flow style (inline) representation for collected representers
    )
    safe_dumper.add_representer(Reference, reference_representer)
    return safe_dumper


@dataclass
class CMakePresets:
    version: int
    cmakeMinimumRequired: dict
    include: list[str]

    configurePresets: list = field(default_factory=list)
    buildPresets: list = field(default_factory=list)
    testPresets: list = field(default_factory=list)
    workflowPresets: list = field(default_factory=list)

    @staticmethod
    def create() -> CMakePresets:
        return CMakePresets(
            version=6,
            cmakeMinimumRequired={"major": 3, "minor": 25, "patch": 0},
            include=[str(SCRIPT_DIR / "cmake-fragments.json")],
        )

    def add_preset_chain(
        self,
        cpreset: ConfigurePreset,
        targets: str | list[str] | None = None,
    ):
        if targets is None:
            targets = []

        self.configurePresets.append(cpreset)
        self.buildPresets.append(BuildPreset.for_configure_preset(cpreset, targets=targets))
        self.testPresets.append(TestPreset.for_configure_preset(cpreset))
        self.workflowPresets.append(WorkflowPreset.for_configure_preset(cpreset))

    def export_json(self, fp: Path):
        with fp.open("w", encoding="utf-8") as jsonfile:
            json.dump(asdict(self), jsonfile, indent=2)


@dataclass
class ConfigurePreset:
    name: str
    inherits: list[str] = ()
    displayName: str | None = None
    cacheVariables: dict[str, str] = field(default_factory=dict)
    generator: str = "Ninja"

    @staticmethod
    def from_fragments(*frags: str, name: str | None = None, **kwargs):
        if name is None:
            name = ".ci-" + "-".join(frags)

        inherits = [f".{frag}" for frag in frags] + [".ci-base"]

        if "mac" in frags:  # V8 and SweepGen don't work right with AppleClang
            inherits.remove(".sweepgen")

        return ConfigurePreset(
            name,
            inherits=inherits,
            displayName="-".join(frags),
            **kwargs,
        )


@dataclass
class BuildPreset:
    name: str
    configurePreset: str
    targets: str | list[str] = field(default_factory=list)

    @staticmethod
    def for_configure_preset(preset: ConfigurePreset, **kwargs):
        return BuildPreset(preset.name, preset.name, **kwargs)


@dataclass
class TestPreset:
    name: str
    configurePreset: str
    inherits: list[str] = field(default_factory=lambda: [".ci-test-base"])

    @staticmethod
    def for_configure_preset(preset: ConfigurePreset, **kwargs):
        return TestPreset(preset.name, preset.name, inherits=[".ci-test-base"], **kwargs)


@dataclass
class WorkflowPreset:
    name: str
    steps: list

    @staticmethod
    def for_configure_preset(preset: ConfigurePreset):
        return WorkflowPreset(
            preset.name,
            [{"type": ptype, "name": preset.name} for ptype in ("configure", "build", "test")],
        )


@dataclass
class CompilerSpec:
    id: str
    cc: str
    cxx: str


@dataclass
class CiJobSpec:
    extends: str
    variables: dict[str, str]
    image: str | None

    tags: list[str | None] = field(default_factory=list)
    stage: str = "test"
    script: list[str | Reference] = field(default_factory=list)
    before_script: list[str | Reference] = field(default_factory=list)
    after_script: list[str | Reference] = field(default_factory=list)


@dataclass
class CiConfig:
    jobspecs: dict[str, CiJobSpec] = field(default_factory=dict)

    def add_job_for_preset(
        self,
        preset: ConfigurePreset,
        compiler: CompilerSpec,
    ):
        spec = CiJobSpec(
            extends=".testsuite-base-linux",
            image=IMAGE_PATTERN.format(compiler_id=compiler.id),  # TODO change here the image pattern once necessary
            variables={
                "cmakePresetName": preset.name,
                "CC": compiler.cc,
                "CXX": compiler.cxx,
            },
        )

        # --- Compiler specific adjustments ---

        match compiler.id.split("-"):
            case ["AppleClang", *_]:
                spec.extends = ".testsuite-base-MacOS"
                spec.image = None
            case ["clang", *_]:
                # The clang from spack cannot fine omp without the patch.
                spec.script.append(Reference(".clang-library-path-patch", "before_script"))
                # Since `extends` is overriding lists instead of merging them, we need to re-append the right before_script.
                if "cuda" in preset.name:
                    spec.script.append(Reference(".testsuite-base-linux-cuda", "before_script"))
                else:
                    spec.script.append(Reference(".testsuite-common", "before_script"))
            case ["gcc", int(version), *_] if version < 13 and "cuda" in preset.name:
                # There are problems with older GCC versions as CUDA host compiler,
                # so we set the CUDAHOSTCXX images base compiler wich currently is g++-13.
                # This should be updated if the CUDA host compiler base image is updated.
                spec.variables["CUDAHOSTCXX"] = "g++-13"

        if "cuda" in preset.name:
            spec.extends += "-cuda"

        self.jobspecs[f"{compiler.id} [{preset.displayName}]"] = spec

    def export_yaml(self, fp: Path):
        obj = {
            #   Use a workflow rule to ensure the matrix pipeline always runs, and only runs as
            #   a dynamic child pipeline
            "workflow": {"rules": [{"if": '$CI_PIPELINE_SOURCE == "parent_pipeline"'}]},
            "include": "utilities/workflows/ci-common.yaml",
        }
        obj = obj | {
            name: asdict(
                spec,
                dict_factory=lambda d: {k: v for (k, v) in d if v},  # Remove None or empty lists items
            )
            for name, spec in self.jobspecs.items()
        }
        with fp.open("w", encoding="utf-8") as yamlfile:
            yaml.dump(obj, yamlfile, Dumper=get_dumper())


###################################################################################################
# PRESETS AND COMPILERS
###################################################################################################


MATRIX_CONFIGURE_PRESETS = [
    ConfigurePreset.from_fragments("cuda", "hybrid", "debugOptimized", "singlePrecision"),
    ConfigurePreset.from_fragments("cuda", "omponly", "make"),
    ConfigurePreset.from_fragments("cpuonly", "mpionly", "debugOptimized"),
    ConfigurePreset.from_fragments("cpuonly", "hybrid", "singlePrecision", "make"),
    ConfigurePreset.from_fragments("cpuonly", "omponly", "debugOptimized"),
    ConfigurePreset.from_fragments("cpuonly", "serial"),
]

MATRIX_COMPILERS = [
    CompilerSpec("clang-17", "clang", "clang++"),
    CompilerSpec("clang-21", "clang", "clang++"),
    CompilerSpec("clang-22", "clang", "clang++"),
    CompilerSpec("gcc-bootstrap", "gcc", "g++"),
    CompilerSpec("gcc-14", "gcc", "g++"),
    CompilerSpec("gcc-15", "gcc", "g++"),
    CompilerSpec("icx-2023", "icx", "icpx"),
    CompilerSpec("icx-2025", "icx", "icpx"),
    CompilerSpec("icx-2026", "icx", "icpx"),
]

#   CI Test Matrix.
#   Make sure that each preset and compiler ID used here is defined in the arrays above.
CI_MATRIX = {
    "clang-17": ["cuda-hybrid-debugOptimized-singlePrecision", "cpuonly-serial"],
    "clang-21": ["cpuonly-mpionly-debugOptimized"],
    "icx-2023": ["cpuonly-hybrid-singlePrecision-make"],
    "icx-2025": ["cuda-omponly-make", "cpuonly-mpionly-debugOptimized"],
    "gcc-bootstrap": ["cpuonly-hybrid-singlePrecision-make", "cpuonly-omponly-debugOptimized"],
    "gcc-14": ["cpuonly-mpionly-debugOptimized"],
}


def get_cmake_presets() -> CMakePresets:
    presets = CMakePresets.create()

    #   GCC + Coverage Testsuite preset chain
    presets.add_preset_chain(
        ConfigurePreset.from_fragments("mpionly", "cuda", "gcov", name=".ci-testsuite-gcc-coverage"),
        targets="waLBerlaTestsuite",
    )

    #   Clang Testsuite preset chain
    presets.add_preset_chain(
        ConfigurePreset.from_fragments("hybrid", name=".ci-testsuite-clang"),
        targets="waLBerlaTestsuite",
    )

    #   Intel Testsuite preset chain
    presets.add_preset_chain(
        ConfigurePreset.from_fragments("omponly", name=".ci-testsuite-icx"),
        targets="waLBerlaTestsuite",
    )

    #   clang-tidy configure preset
    presets.configurePresets.append(
        ConfigurePreset.from_fragments(
            "clang",
            "debug",
            "hybrid",
            "cuda",
            name=".ci-clang-tidy",
            generator="Unix Makefiles",
            cacheVariables={
                "CMAKE_EXPORT_COMPILE_COMMANDS": True,
                "WALBERLA_BUFFER_DEBUG": True,
                "WALBERLA_LOGLEVEL": "DETAIL",
            },
        ),
    )

    #   Matrix preset chains
    for cpreset in MATRIX_CONFIGURE_PRESETS:
        presets.add_preset_chain(cpreset)

    return presets


###################################################################################################


def generate_presets(args):
    presets = get_cmake_presets()
    presets_filepath = Path(args.output_file).resolve()
    presets.export_json(presets_filepath)

    if args.patch_main_presets_file:
        presets_file_path = Path(args.patch_main_presets_file).resolve()
        presets_file_path.parent.mkdir(parents=True, exist_ok=True)

        project_dir = presets_file_path.parent
        output_relative_to_proj = presets_filepath.relative_to(project_dir)

        if not presets_file_path.exists():
            #   Create an empty preset file
            data = CMakePresets.create().export_json(presets_file_path)

        with presets_file_path.open("r") as pfile:
            data: dict = json.load(pfile)

        include_file = str(output_relative_to_proj)
        includes = data.setdefault("include", [])

        if include_file not in includes:
            includes.append(str(output_relative_to_proj))

        with presets_file_path.open("w") as pfile:
            json.dump(data, pfile, indent=2)


def generate_ci_test_matrix(args):
    matrix = CI_MATRIX
    ci_config = CiConfig()

    for compiler_id, preset_ids in matrix.items():
        compiler = next(filter(lambda c: c.id == compiler_id, MATRIX_COMPILERS))

        if isinstance(preset_ids, str):
            preset_ids = [preset_ids]

        for preset_id in preset_ids:
            preset = next(filter(lambda p: p.name == f".ci-{preset_id}", MATRIX_CONFIGURE_PRESETS))
            ci_config.add_job_for_preset(preset, compiler)

    config_filepath = Path(args.output_file).resolve()
    config_filepath.parent.mkdir(parents=True, exist_ok=True)
    ci_config.export_yaml(config_filepath)


def main():
    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(required=True)

    p_presets = subparsers.add_parser("cmake-presets")
    p_presets.add_argument("output_file", type=str)
    p_presets.add_argument(
        "-p",
        "--patch-main-presets-file",
        dest="patch_main_presets_file",
        nargs="?",
        type=str,
        help="Patch the named CMakePresets.json to include the generated preset matrix",
    )
    p_presets.set_defaults(func=generate_presets)

    p_ci_matrix = subparsers.add_parser("ci-matrix")
    p_ci_matrix.add_argument("output_file", type=str)
    p_ci_matrix.set_defaults(func=generate_ci_test_matrix)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
