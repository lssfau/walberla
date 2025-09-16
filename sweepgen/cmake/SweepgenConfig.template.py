#   Project and build-system configuration for the code generator

from pystencilssfg import SfgConfig
from sweepgen import WalberlaBuildConfig


def configure_sfg(cfg: SfgConfig):
    cfg.extensions.header = "hpp"
    cfg.extensions.impl = "cpp"


def project_info() -> WalberlaBuildConfig:
    from sweepgen.build_config import cmake_parse_bool

    return WalberlaBuildConfig(
        c_compiler_id="${CMAKE_C_COMPILER_ID}",
        cxx_compiler_id="${CMAKE_CXX_COMPILER_ID}",
        use_double_precision=cmake_parse_bool("${WALBERLA_DOUBLE_ACCURACY}"),
        optimize_for_localhost=cmake_parse_bool("${WALBERLA_OPTIMIZE_FOR_LOCALHOST}"),
        mpi_enabled=cmake_parse_bool("${WALBERLA_BUILD_WITH_MPI}"),
        openmp_enabled=cmake_parse_bool("${WALBERLA_BUILD_WITH_OPENMP}"),
        cuda_enabled=cmake_parse_bool("${WALBERLA_BUILD_WITH_CUDA}"),
        hip_enabled=cmake_parse_bool("${WALBERLA_BUILD_WITH_HIP}"),
        likwid_enabled=cmake_parse_bool("${WALBERLA_BUILD_WITH_LIKWID_MARKERS}"),
    )


def validate():
    _ = project_info()


if __name__ == "__main__":
    validate()
