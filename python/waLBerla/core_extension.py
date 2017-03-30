try:
    from . import walberla_cpp
except ImportError:
    import walberla_cpp

class SliceMaker(object):
    def __getitem__(self, item):
        return item
makeSlice = SliceMaker()


def normalizeSlice( slices, sizes ):
    """Converts slices with floating point entries to integer slices"""
    assert( len(slices) == len(sizes) )

    result = []

    for s, size in zip( slices, sizes ):
        if type(s) is int:
            result.append( s )
            continue
        if type(s) is float:
            result.append( int( s*size ) )
            continue
            
        assert( type(s) is slice )
        
        if s.start is None:
            newStart = 0
        elif type(s.start) is float:
            newStart = int(s.start * size )
        else:
            newStart = s.start
            
        if s.stop is None:
            newStop = size
        elif type(s.stop) is float:
            newStop = int(s.stop * size )
        else:
            newStop = s.stop
        
        result.append( slice(newStart,newStop, s.step) )
    
    return tuple(result)


def sliceToCellInterval( s ):
    newMin = [0,0,0]
    newMax = [0,0,0]
    for i in range(3):
        if type(s[i]) is int:
            newMin[i] = s[i]
            newMax[i] = s[i]
        else:
            newMin[i] = s[i].start
            newMax[i] = s[i].stop-1
    return walberla_cpp.CellInterval( newMin,newMax)
    

def cellIntervalToSlice( cellInterval ):
    slices = []
    for i in range(3):
        if cellInterval.min[i] == cellInterval.max[i]:
            slices.append( cellInterval.min[i] )
        else:
            slices.append( slice(cellInterval.min[i], cellInterval.max[i]+1,None ) )
    return slices
    #return makeSlice[ cellInterval.min[0]:cellInterval.max[0]+1,
    #                  cellInterval.min[1]:cellInterval.max[1]+1,
    #                  cellInterval.min[2]:cellInterval.max[2]+1 ]





def extend(coreModule):
    coreModule.makeSlice = SliceMaker()
    coreModule.normalizeSlice = normalizeSlice
    coreModule.CellInterval.fromSlice = staticmethod(sliceToCellInterval)
    coreModule.CellInterval.toSlice = cellIntervalToSlice