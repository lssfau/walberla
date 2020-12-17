import numpy

# try:
#     from . import walberla_cpp
# except ImportError:
#     import walberla_cpp


# ----------------------------- Python functions to extend the C++ field module ---------------------------------

def normalize_ghostlayer_info(field, with_ghost_layers):
    """Takes one ghost layer parameter and returns an integer:
        True -> all ghost layers, False->no ghost layers
        Args:
            field: waLberl field object
            with_ghost_layers: see numpy_array_from_walberla_field
    """

    def normalize_component(gl):
        if gl is False:
            return 0
        if gl is True:
            return field.nrOfGhostLayers
        if gl > field.nrOfGhostLayers:
            raise ValueError("Field only has %d ghost layers (requested %d)" % (field.nrOfGhostLayers, gl))
        return gl

    if hasattr(with_ghost_layers, "__len__") and len(with_ghost_layers) == 3:
        ghost_layers = [normalize_component(gl) for gl in with_ghost_layers]
    else:
        ghost_layers = [normalize_component(with_ghost_layers)] * 3
    return ghost_layers


def numpy_array_from_walberla_field(field, with_ghost_layers=False):
    """ Creates a numpy array view on the waLBerla field data
        @field: the waLBerla field
        @with_ghost_layers: Possible values:
                            1. Boolean: False: no ghost layers included
                                        True:  all ghost layers included
                            2. Integer: number of ghost layers to include
                            3. List with three booleans or integers with ghost layer info for x,y,z direction
    """

    if not field:
        return None

    if hasattr(field, 'nrOfGhostLayers'):
        field_gl = field.nrOfGhostLayers
    else:
        field_gl = 0
    ghost_layers = normalize_ghostlayer_info(field, with_ghost_layers)

    if ghost_layers[0] == field_gl and ghost_layers[1] == field_gl and ghost_layers[2] == field_gl:
        return numpy.asarray(field)
    else:
        result = numpy.asarray(field)
        cutoff = [abs(gl - field.nrOfGhostLayers) for gl in ghost_layers]
        if len(result.shape) == 4:
            view = result[cutoff[0]:-cutoff[0] if cutoff[0] > 0 else None,
                          cutoff[1]:-cutoff[1] if cutoff[1] > 0 else None,
                          cutoff[2]:-cutoff[2] if cutoff[2] > 0 else None,
                          :]
        else:
            view = result[cutoff[0]:-cutoff[0] if cutoff[0] > 0 else None,
                          cutoff[1]:-cutoff[1] if cutoff[1] > 0 else None,
                          cutoff[2]:-cutoff[2] if cutoff[2] > 0 else None]
        return view


def copy_array_to_field(dst_field, src_array, idx=None, with_ghost_layers=False):
    """ Copies a numpy array into (part of) a waLBerla field

    Usually no copying has to take place between waLBerla fields and numpy arrays, since an array view can be
    constructed on a field that uses the same memory.
    When running certain numpy operations that cannot be done in-place,however, the data has to be copied back.

    @param dst_field: waLBerla field, where the data is copied to
    @param src_array: numpy array where to copy from
    @param idx:      the numpy array is allowed to be smaller than the field. In this case the target region
                     has to be specified via this 3 dimensional slice
    @param with_ghost_layers: if true the ghost layers of the field are considered as well
    """
    if idx is None:
        idx = [slice(None, None, None)] * 3
    dst_as_array = numpy_array_from_walberla_field(dst_field, with_ghost_layers)
    numpy.copyto(dst_as_array[idx], src_array)


def extend(cpp_field_module):
    def gather_field(blocks, block_data_name, slice_obj, all_gather=False):
        field = cpp_field_module.gather(blocks, block_data_name, slice_obj, targetRank=-1 if all_gather else 0)
        if field is not None:
            field = numpy_array_from_walberla_field(field)
            field.flags.writeable = False
            return field
        else:
            return None

    cpp_field_module.toArray = numpy_array_from_walberla_field
    cpp_field_module.copyArrayToField = copy_array_to_field
    cpp_field_module.gatherField = gather_field
