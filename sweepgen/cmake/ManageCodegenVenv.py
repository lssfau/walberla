# This file is part of waLBerla. waLBerla is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# waLBerla is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
#
# Author: Frederik Hennig <frederik.hennig@fau.de>
#

"""Virtual environment manager for use within CMake"""

from __future__ import annotations

from typing import Generator
import sys
import subprocess
import json
import shutil
import hashlib
from contextlib import contextmanager
from dataclasses import dataclass, asdict, field
from argparse import ArgumentParser
from pathlib import Path


@dataclass
class VenvState:
    venv_dir: str | None = None
    main_requirements_file: str | None = None
    initialized: bool = False
    populated: bool = False

    user_requirements: list[list[str]] = field(default_factory=list)
    user_requirements_hash: str | None = None

    @property
    def venv_python(self) -> Path:
        assert self.venv_dir is not None
        return Path(self.venv_dir) / "bin" / "python"

    @staticmethod
    @contextmanager
    def lock(statefile: Path) -> Generator[VenvState, None, None]:
        statefile_bak = statefile.with_suffix(".json.bak")
        if statefile.exists():
            statefile.replace(statefile_bak)

            with statefile_bak.open("r") as f:
                state_dict = json.load(f)

            state = VenvState(**state_dict)
        else:
            state = VenvState()

        yield state

        #   If the consumer raises an error, execution terminates here

        state_dict = asdict(state)
        with statefile_bak.open("w") as f:
            json.dump(state_dict, f)
        statefile_bak.replace(statefile)


def reinitialize(args, state: VenvState):
    assert state.venv_dir is not None

    venv_dir = Path(state.venv_dir)
    if venv_dir.exists():
        shutil.rmtree(venv_dir)

    base_py = Path(sys.executable)

    #   Create the virtual environment
    venv_args = [base_py, "-m", "venv", state.venv_dir]
    subprocess.run(venv_args).check_returncode()

    #   Install base requirements
    venv_py = str(state.venv_python)
    install_args = [venv_py, "-m", "pip", "install", "-r", state.main_requirements_file]

    if args.quiet:
        install_args.append("-q")

    subprocess.run(install_args).check_returncode()

    state.initialized = True
    state.user_requirements_hash = None  # Reset installation state


def action_initialize(args):
    statefile = Path(args.statefile)

    with VenvState.lock(statefile) as state:
        if not state.initialized or args.force:
            p_venv_dir = Path(args.venv_dir).resolve()
            reqs_file = Path(args.requirements_file).resolve()

            state.venv_dir = str(p_venv_dir)
            state.main_requirements_file = str(reqs_file)

            reinitialize(args, state)

        #   Reset user requirements
        state.user_requirements = []
        state.populated = False


def action_require(args):
    statefile = Path(args.statefile)

    with VenvState.lock(statefile) as state:
        if not state.initialized:
            raise RuntimeError(
                "venv-require action failed: virtual environment was not initialized"
            )

        if state.populated:
            raise RuntimeError(
                "venv-require action failed: cannot add requirements after venv-populate was run"
            )

        state.user_requirements.append(list(args.requirements))


def action_populate(args):
    statefile = Path(args.statefile)

    with VenvState.lock(statefile) as state:
        if not state.initialized:
            raise RuntimeError(
                "venv-populate action failed: virtual environment was not initialized"
            )

        if state.populated:
            raise RuntimeError(
                "venv-populate action failed: venv-populate action called twice"
            )

        h = hashlib.sha256()
        for req in state.user_requirements:
            h.update(bytes(";".join(str(r) for r in req), encoding="utf8"))
        digest = h.hexdigest()

        if digest != state.user_requirements_hash:
            if state.user_requirements_hash is not None:
                #   Populate was run before -> rebuild entire environment
                print(
                    "User requirements have changed - rebuilding virtual environment.",
                    flush=True,
                )
                reinitialize(args, state)

            pip_args = [str(state.venv_python), "-m", "pip", "install"]

            if args.quiet:
                pip_args.append("-q")

            for req in state.user_requirements:
                install_args = pip_args.copy() + req
                subprocess.run(install_args).check_returncode()

        state.user_requirements_hash = digest
        state.populated = True


def fail(*args: str):
    for arg in args:
        print(arg, file=sys.stderr)
    exit(-1)


def main():
    parser = ArgumentParser("ManageCodegenVenv")
    parser.add_argument(
        "-s",
        "--statefile",
        required=True,
        dest="statefile",
        help="Path to the environment statefile",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Silence output of `pip install` during virtual environment setup",
    )

    subparsers = parser.add_subparsers(required=True)

    parser_initialize = subparsers.add_parser("init")
    parser_initialize.add_argument("-f", "--force", dest="force", action="store_true")
    parser_initialize.add_argument(
        "venv_dir",
        type=str,
        help="Location of the virtual environment in the filesystem",
    )
    parser_initialize.add_argument(
        "requirements_file",
        type=str,
        help="Location of the virtual environment in the filesystem",
    )
    parser_initialize.set_defaults(func=action_initialize)

    parser_require = subparsers.add_parser("require")
    parser_require.add_argument("requirements", nargs="+", type=str)
    parser_require.set_defaults(func=action_require)

    parser_populate = subparsers.add_parser("populate")
    parser_populate.set_defaults(func=action_populate)

    try:
        args = parser.parse_args()
        args.func(args)
    except RuntimeError as e:
        fail(*e.args)
    except subprocess.CalledProcessError as e:
        fail()


if __name__ == "__main__":
    main()
