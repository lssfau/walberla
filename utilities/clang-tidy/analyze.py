"""
Script performing fully automated clang-tidy analysis runs
on the waLBerla repository.
"""

import pathlib
import re
import sys
import argparse
import yaml
import json
import shutil
import subprocess
from collections import Counter
from textwrap import indent


def removePrecompiler(x):
    pos = x.find("clang++")
    if pos != -1:
        return x[pos:]
    else:
        return x


def get_directory_filter(include_dirs: list[pathlib.Path]):
    absolute_dirs: list[pathlib.Path] = []
    for src_dir in include_dirs:
        if not src_dir.exists():
            raise ValueError(f"Source directory {src_dir} does not exist")

        absolute_dirs.append(src_dir.absolute())

    def filter(entry: dict):
        fp = pathlib.Path(entry["file"])
        return any(fp.is_relative_to(d) for d in absolute_dirs)

    return filter


def remove_duplicate_commands(db: list):
    seen_files: set[str] = set()

    db_filtered = []
    for x in db:
        if x["file"] not in seen_files:
            seen_files.add(x["file"])
            db_filtered.append(x)
    return db_filtered


WARNING_PATTERN = re.compile(r"\[[a-z-]+,-warnings-as-errors\]\n")
TRAILING = len(",-warnings-as-errors]\n")


def count_diagnostics(output_file: pathlib.Path) -> Counter:
    clang_tidy_log = output_file.read_text()
    matches = WARNING_PATTERN.findall(clang_tidy_log)
    matches = [m[1:-TRAILING] for m in matches]
    counts = Counter(matches)

    return counts


def print_summary(counts: Counter):
    counts = sorted(list(counts.items()), key=lambda it: it[1], reverse=True)
    summary = "\n".join(f"{warning}: {count}" for warning, count in counts)
    return summary


def main():
    parser = argparse.ArgumentParser("analyze.py")

    parser.add_argument(
        "-p",
        "--parameter-file",
        dest="parameter_file",
        type=str,
        default=None,
        help="Parameter file that defines which modules and apps should be analyzed.",
    )
    parser.add_argument(
        "-m",
        "--modules",
        dest="modules",
        type=str,
        nargs="*",
        default=None,
        help="Modules to be analyzed. Overrides modules listed in the parameter file."
    )
    parser.add_argument(
        "-a",
        "--apps",
        dest="apps",
        type=str,
        nargs="*",
        default=None,
        help="Apps to be analyzed. Overrides apps listed in the parameter file."
    )
    parser.add_argument(
        "-r",
        "--project-root",
        dest="project_root",
        type=str,
        default=".",
        help="waLBerla project root directory. If not specified, use the current directory.",
    )
    parser.add_argument(
        "-c",
        "--compile-database",
        required=True,
        dest="compile_database",
        type=str,
        help="Path to the CMake compile command data base.",
    )
    parser.add_argument(
        "-o",
        "--output-directory",
        dest="output_directory",
        type=str,
        default="./clang-tidy-output",
        help="Folder to which the error streams from clang-tidy should be written.",
    )
    parser.add_argument(
        "--checks",
        dest="checks",
        type=str,
        default=None,
        nargs="+",
        help="clang-tidy checks filter. Forwarded to -checks argument of clang-tidy."
    )
    parser.add_argument(
        "--export-fixes",
        dest="export_fixes",
        action="store_true",
        help="Export possible fixes detected by clang-tidy to a YAML file in the output directory"
    )

    args = parser.parse_args()

    walberla_root = pathlib.Path(args.project_root).resolve()
    src_dir = walberla_root / "src"
    tests_dir = walberla_root / "tests"
    apps_dir = walberla_root / "apps"

    output_dir = pathlib.Path(args.output_directory)

    database_fp = pathlib.Path(args.compile_database)
    database_backup = database_fp.parent / f"{database_fp.name}.bak"
    shutil.copy(str(database_fp), str(database_backup))

    success: bool = True
    total_counts = Counter()

    try:
        with database_fp.open() as dbfile:
            orig_db = json.load(dbfile)

        if args.parameter_file is not None:
            parameter_filepath = pathlib.Path(args.parameter_file)
            with parameter_filepath.open() as pfile:
                params = yaml.load(pfile, yaml.SafeLoader)
        else:
            params = dict()

        if args.modules:
            params["modules"] = args.modules

        if args.apps:
            params["apps"] = args.apps

        build_dir = database_fp.parent
        clang_tidy_base_args = ["run-clang-tidy", "-p", str(build_dir)]

        def run_clang_tidy(
            database: list,
            include_dirs: list[pathlib.Path],
            header_filter: str | None,
            output: pathlib.Path,
            info: str,
        ):
            print(f"Analyzing {info}")

            try:
                dir_filter = get_directory_filter(include_dirs)
            except ValueError as e:
                print(e, file=sys.stderr)
                return False

            cc_filtered = list(filter(dir_filter, database))
            for x in cc_filtered:
                x["command"] = removePrecompiler(x["command"])

            cc_filtered = remove_duplicate_commands(cc_filtered)

            print(f"  -- Retained {len(cc_filtered)} compile commands")

            with database_fp.open("w") as db_out:
                json.dump(cc_filtered, db_out)

            clang_tidy_args = clang_tidy_base_args
            if header_filter:
                clang_tidy_args += ["-header-filter", header_filter]
            if args.checks:
                checks_str = ",".join(args.checks)
                clang_tidy_args += [f"-checks={checks_str}"]
            if args.export_fixes:
                fixes_file = output.parent / "fixes.yaml"
                clang_tidy_args += [f"-export-fixes={str(fixes_file)}"]
            output.parent.mkdir(exist_ok=True, parents=True)

            errfile = output.with_name(output.name + ".err")

            print("  -- Running clang-tidy...")
            with output.open("w") as ofile:
                with errfile.open("w") as efile:
                    subprocess.run(clang_tidy_args, stdout=ofile, stderr=efile)

            print(f"  -- clang-tidy output written to {str(output)}")

            counts = count_diagnostics(output)
            total_counts.update(counts)

            if counts.total() == 0:
                print("  -- Success!")
            else:
                summary = print_summary(counts)

                print("  -- Summary:")
                print(indent(summary, "      - "))

            print("\n\n", end="")
            return counts.total() == 0

        for module_spec in params.get("modules", []):
            include_paths: list[pathlib.Path] = []
            module_name: str
            match module_spec:
                case str():
                    module_name = module_spec
                    include_paths = [src_dir / module_name, tests_dir / module_name]
                case dict():
                    module_name, settings = next(iter(module_spec.items()))
                    include_paths = [src_dir / module_name]
                    if not settings.get("exclude-tests", False):
                        include_paths.append(tests_dir / module_name)
            header_filter = rf".*/src/{module_name}/.*"

            succ = run_clang_tidy(
                orig_db,
                include_paths,
                header_filter,
                output_dir / "modules" / f"{module_name}.out",
                f"module {module_name}",
            )
            success = succ and success

        for app_spec in params.get("apps", []):
            include_paths: list[pathlib.Path]
            app_name: str

            match app_spec:
                case str():
                    app_name = app_spec
                    include_paths = [apps_dir / app_name]
                case dict():
                    app_name, settings = next(iter(app_spec.items()))
                    if (only := settings.get("only", None)) is not None:
                        include_paths = [apps_dir / app_name / o for o in only]
                    else:
                        include_paths = [apps_dir / app_name]

            succ = run_clang_tidy(
                orig_db,
                include_paths,
                None,
                output_dir / "apps" / f"{app_name}.out",
                f"application {app_name}",
            )

            success = succ and success

    finally:
        #   Restore the backup
        shutil.move(str(database_backup), str(database_fp))

    print(f"Done. Total number of diagnostics: {total_counts.total()}")
    print()
    print("Summary:")
    print(indent(print_summary(total_counts), "  - "))

    exit(0 if success else 1)


if __name__ == "__main__":
    main()
