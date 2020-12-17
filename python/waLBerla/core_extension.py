try:
    from . import walberla_cpp
except ImportError:
    import walberla_cpp


class SliceMaker(object):
    def __getitem__(self, item):
        return item


makeSlice = SliceMaker()


def normalize_slice(slices, sizes):
    """Converts slices with floating point entries to integer slices"""
    assert (len(slices) == len(sizes))

    result = []

    for s, size in zip(slices, sizes):
        if type(s) is int:
            result.append(s)
            continue
        if type(s) is float:
            result.append(int(s * size))
            continue

        assert (type(s) is slice)

        if s.start is None:
            new_start = 0
        elif type(s.start) is float:
            new_start = int(s.start * size)
        else:
            new_start = s.start

        if s.stop is None:
            new_stop = size
        elif type(s.stop) is float:
            new_stop = int(s.stop * size)
        else:
            new_stop = s.stop

        result.append(slice(new_start, new_stop, s.step))

    return tuple(result)


def slice_to_cell_interval(s):
    new_min = [0, 0, 0]
    new_max = [0, 0, 0]
    for i in range(3):
        if type(s[i]) is int:
            new_min[i] = s[i]
            new_max[i] = s[i]
        else:
            new_min[i] = s[i].start
            new_max[i] = s[i].stop - 1
    return walberla_cpp.CellInterval(new_min[0], new_min[1], new_min[2], new_max[0], new_max[1], new_max[2])


def cell_interval_to_slice(cell_interval, collapse_extent_one=True):
    if not hasattr(collapse_extent_one, '__len__'):
        collapse_extent_one = (collapse_extent_one, collapse_extent_one, collapse_extent_one)

    slices = []
    for i, collapseInfo in enumerate(collapse_extent_one):
        if collapseInfo and cell_interval.min[i] == cell_interval.max[i]:
            slices.append(cell_interval.min[i])
        else:
            slices.append(slice(cell_interval.min[i], cell_interval.max[i] + 1, None))
    return tuple(slices)


def extend(core_module):
    core_module.makeSlice = SliceMaker()
    core_module.normalizeSlice = normalize_slice
    core_module.CellInterval.fromSlice = staticmethod(slice_to_cell_interval)
    core_module.CellInterval.toSlice = cell_interval_to_slice
