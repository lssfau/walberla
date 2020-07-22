from .callbacks import callback, ScenarioManager, memberCallback  # noqa:F401
# deprecated, was renamed to memberCallback
from .callbacks import memberCallback as member_callback  # noqa:F401

import sys

try:
    from .walberla_cpp import *  # noqa: F403

    cpp_available = True
except ImportError:
    try:
        from walberla_cpp import *  # noqa: F403

        cpp_available = True
    except ImportError:
        cpp_available = False
        thismodule = sys.modules[__name__]
        thismodule.log_info_on_root = print
        thismodule.log_progress_on_root = print
        thismodule.log_result_on_root = print
        thismodule.log_devel_on_root = print
        thismodule.log_detail_on_root = print
        thismodule.log_warning_on_root = print
        thismodule.log_info = print
        thismodule.log_result = print
        thismodule.log_devel = print
        thismodule.log_progress = print
        thismodule.log_detail = print
        thismodule.log_warning = print

if cpp_available:
    from .core_extension import extend as extend_core

    thismodule = sys.modules[__name__]
    extend_core(thismodule)

    if 'field' in globals():  # check if field was exported
        # Update modules dict to be able to write e.g. from waLBerla import field
        # otherwise "field" would only be a scope not a module
        sys.modules[__name__ + '.field'] = field  # noqa: F405
        # extend the C++ module with some python functions
        from .field_extension import extend as extend_field

        extend_field(field)  # noqa: F405
    if 'cuda' in globals():
        sys.modules[__name__ + '.cuda'] = cuda  # noqa: F405
        from .cuda_extension import extend as extend_cuda

        extend_cuda(cuda)  # noqa: F405
    if 'geometry' in globals():
        sys.modules[__name__ + '.geometry'] = geometry  # noqa: F405
    if 'lbm' in globals():
        sys.modules[__name__ + '.lbm'] = lbm  # noqa: F405
    if 'postprocessing' in globals():
        sys.modules[__name__ + 'postprocessing'] = postprocessing  # noqa: F405
    if 'mpi' in globals():
        sys.modules[__name__ + '.mpi'] = mpi  # noqa: F405
    if 'timeloop' in globals():
        sys.modules[__name__ + '.timeloop'] = timeloop  # noqa: F405
        from .timeloop_extension import extend as extend_timeloop

        extend_timeloop(timeloop)  # noqa: F405
else:
    class Dummy:
        pass

    callbacks = Dummy()
