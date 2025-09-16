from __future__ import annotations

from dataclasses import dataclass, field

from pystencils import CreateKernelConfig, Target
from pystencils.types.quick import Fp
from pystencils.jit import no_jit

from pystencilssfg import SfgContext
from pystencilssfg.composer import SfgIComposer


def cmake_parse_bool(var: str):
    var = var.upper()
    if var in ("ON", "1", "TRUE"):
        return True
    elif var in ("OFF", "0", "FALSE"):
        return False
    else:
        raise ValueError(f"Could not parse cmake value `{var}` as boolean.")


@dataclass
class ConfigOverrides:
    target: Target | None = None
    """Override code generation target"""


@dataclass
class WalberlaBuildConfig:
    """Represents a waLBerla build system configuration.

    Objects of this class represent the build settings which the waLBerla
    CMake system was configured with.
    It provides the appropriate code generator settings for pystencils
    via `get_pystencils_config`.

    The following settings can be overridden for the current code generation script:

     - `target`, to generate code for a different than the default hardware target
    """

    c_compiler_id: str
    """Value of ``CMAKE_C_COMPILER_ID``."""

    cxx_compiler_id: str
    """Value of ``CMAKE_CXX_COMPILER_ID``."""

    use_double_precision: bool
    """Value of ``WALBERLA_DOUBLE_ACCURACY``"""

    optimize_for_localhost: bool
    """Value of ``WALBERLA_OPTIMIZE_FOR_LOCALHOST``."""

    mpi_enabled: bool
    """Value of ``WALBERLA_BUILD_WITH_MPI``."""

    openmp_enabled: bool
    """Value of ``WALBERLA_BUILD_WITH_OPENMP``."""

    cuda_enabled: bool
    """Value of ``WALBERLA_BUILD_WITH_CUDA``."""

    hip_enabled: bool
    """Value of ``WALBERLA_BUILD_WITH_HIP``."""

    likwid_enabled: bool
    """Value of ``WALBERLA_BUILD_WITH_LIKWID_MARKERS``"""

    override: ConfigOverrides = field(default_factory=ConfigOverrides)
    """Override code generator options that would otherwise be inferred from the build system."""

    @staticmethod
    def from_sfg(sfg: SfgContext | SfgIComposer) -> WalberlaBuildConfig:
        """Retrieve the build configuration object from the pystencils-sfg session context."""
        if isinstance(sfg, SfgIComposer):
            ctx = sfg.context
        else:
            ctx = sfg

        if isinstance(ctx.project_info, WalberlaBuildConfig):
            return ctx.project_info
        elif DEBUG.BUILD_CONFIG is not None:
            return DEBUG.BUILD_CONFIG
        else:
            raise ValueError(
                "The given SfgContext does not encapsulate a waLBerla build config object."
            )

    @property
    def target(self) -> Target:
        """Return or override the actual target used for code generation"""

        if self.override.target is not None:
            target = self.override.target
        elif self.cuda_enabled:
            target = Target.CUDA
        elif self.hip_enabled:
            target = Target.HIP
        else:
            #  CPU target
            if self.optimize_for_localhost:
                target = Target.CurrentCPU
            else:
                target = Target.GenericCPU

        if target == Target.CurrentGPU:
            if self.cuda_enabled:
                target = Target.CUDA
            elif self.hip_enabled:
                target = Target.HIP
            else:
                raise ValueError("Target `CurrentGPU` was selected in a non-GPU build.")

        return target

    @target.setter
    def target(self, t: Target | None):
        self.override.target = t

    def get_pystencils_config(self) -> CreateKernelConfig:
        """Get the pystencils code generator configuration from this build configuration."""
        cfg = CreateKernelConfig()
        cfg.default_dtype = Fp(64) if self.use_double_precision else Fp(32)
        cfg.jit = no_jit
        cfg.target = self.target

        if self.openmp_enabled:
            cfg.cpu.openmp.enable = True

        return cfg


def get_build_config(sfg: SfgContext | SfgIComposer):
    """Get the waLBerla build config object for the current generator script."""
    return WalberlaBuildConfig.from_sfg(sfg)


class DEBUG:
    """Settings for debugging generator script outside of the waLBerla build system"""

    BUILD_CONFIG: WalberlaBuildConfig | None = None
    """The mockup build config to be used when running outside of the build system."""

    @staticmethod
    def use_cpu_default():
        """Mimic a default CPU build configuration"""
        DEBUG.BUILD_CONFIG = WalberlaBuildConfig(
            c_compiler_id="GNU",
            cxx_compiler_id="GNU",
            use_double_precision=True,
            optimize_for_localhost=False,
            mpi_enabled=True,
            openmp_enabled=False,
            hip_enabled=False,
            cuda_enabled=False,
            likwid_enabled=False,
        )

    @staticmethod
    def use_cuda_default():
        """Mimic a default build configuration with CUDA enabled"""
        DEBUG.BUILD_CONFIG = WalberlaBuildConfig(
            c_compiler_id="GNU",
            cxx_compiler_id="GNU",
            use_double_precision=True,
            optimize_for_localhost=False,
            mpi_enabled=True,
            openmp_enabled=False,
            hip_enabled=False,
            cuda_enabled=True,
            likwid_enabled=False,
        )

    @staticmethod
    def use_hip_default():
        """Mimic a default build configuration with HIP enabled"""
        DEBUG.BUILD_CONFIG = WalberlaBuildConfig(
            c_compiler_id="GNU",
            cxx_compiler_id="GNU",
            use_double_precision=True,
            optimize_for_localhost=False,
            mpi_enabled=True,
            openmp_enabled=False,
            hip_enabled=True,
            cuda_enabled=False,
            likwid_enabled=False,
        )
