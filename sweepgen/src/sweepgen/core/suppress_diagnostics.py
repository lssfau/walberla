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


from __future__ import annotations
from abc import ABC, abstractmethod
from warnings import warn

from enum import Enum, auto

from pystencilssfg import SfgComposer

from ..build_config import get_build_config


class Diagnostic(Enum):
    UNUSED_VARIABLE = auto()
    CONVERSION = auto()


class Preprocessor(ABC):
    """Compiler-specific preprocessor handling"""

    @abstractmethod
    def detection(self) -> str:
        """A preprocessor boolean expression that evaluates to true if this preprocessor is active"""

    @abstractmethod
    def enter(self, sfg: SfgComposer, diagnostics: set[Diagnostic], impl: bool = True):
        """Emit preprocessor directives to enter the suppression region."""

    @abstractmethod
    def exit(self, sfg: SfgComposer, diagnostics: set[Diagnostic], impl: bool = True):
        """Emit preprocessor directives to exit the suppression region."""


class GccLike(Preprocessor):

    #   See https://gcc.gnu.org/onlinedocs/gcc/Diagnostic-Pragmas.html

    def detection(self) -> str:
        return "defined(__GNUC__)"

    def enter(self, sfg: SfgComposer, diagnostics: set[Diagnostic], impl: bool = True):
        sfg.code("#pragma GCC diagnostic push", impl=impl)

        for diag in diagnostics:
            diag_text = "-W" + diag.name.lower().replace("_", "-")
            sfg.code(f'#pragma GCC diagnostic ignored "{diag_text}"', impl=impl)

    def exit(self, sfg: SfgComposer, diagnostics: set[Diagnostic], impl: bool = True):
        sfg.code("#pragma GCC diagnostic pop", impl=impl)


class NVCC(Preprocessor):
    def detection(self) -> str:
        return "defined(__NVCC__) && defined(__NVCC_DIAG_PRAGMA_SUPPORT__)"

    def enter(self, sfg: SfgComposer, diagnostics: set[Diagnostic], impl: bool = True):
        sfg.code("#pragma nv_diagnostic push", impl=impl)

        for diag in diagnostics:
            diag_text = self._diag_id(diag)
            sfg.code(f"#pragma nv_diag_suppress {diag_text}", impl=impl)

    def exit(self, sfg: SfgComposer, diagnostics: set[Diagnostic], impl: bool = True):
        sfg.code("#pragma nv_diagnostic pop", impl=impl)

    def _diag_id(self, diag: Diagnostic) -> str:
        match diag:
            case Diagnostic.UNUSED_VARIABLE:
                return "177 // unused variable"
            case _:
                raise ValueError(f"Unknown NVCC diagnostic ID for {diag}")


class SuppressDiagnostics:
    """Context manager that suppresses compiler diagnostics."""

    def __init__(
        self,
        sfg: SfgComposer,
        *diags: str | Diagnostic,
        impl: bool = True,
        preprocessors: list[Preprocessor] | None = None,
    ):
        self._sfg = sfg

        self._diagnostics: set[Diagnostic] = set()
        self._impl = impl

        self._preprocessors: list[Preprocessor]

        if preprocessors is not None:
            self._preprocessors = preprocessors
        else:
            self._preprocessors = []

            cfg = get_build_config(sfg)

            if cfg.cuda_enabled:
                self._preprocessors.append(NVCC())

            match cfg.cxx_compiler_id:
                case "GNU" | "Clang" | "IntelLLVM":
                    self._preprocessors.append(GccLike())
                case other:
                    warn(
                        f"Compiler {other} not supported by SuppressDiagnostics. Diagnostics will not be suppressed.",
                        UserWarning,
                    )

        self.suppress(*diags)

    def suppress(self, *diags: str | Diagnostic) -> SuppressDiagnostics:
        for diag in diags:
            if isinstance(diag, str):
                diag = Diagnostic[diag.upper().replace("-", "_")]
            self._diagnostics.add(diag)
        return self

    def __enter__(self):
        if self._diagnostics:
            for i, prep in enumerate(self._preprocessors):
                if_ = "if" if i == 0 else "elif"
                self._sfg.code(f"#{if_} {prep.detection()}", impl=self._impl)
                prep.enter(self._sfg, self._diagnostics, self._impl)

                if i == len(self._preprocessors) - 1:
                    self._sfg.code("#endif", impl=self._impl)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None and self._diagnostics:
            for i, prep in enumerate(self._preprocessors):
                if_ = "if" if i == 0 else "elif"
                self._sfg.code(f"#{if_} {prep.detection()}", impl=self._impl)
                prep.exit(self._sfg, self._diagnostics, self._impl)

                if i == len(self._preprocessors) - 1:
                    self._sfg.code("#endif", impl=self._impl)
