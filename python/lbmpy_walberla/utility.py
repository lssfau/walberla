from lbmpy.advanced_streaming import Timestep


def timestep_suffix(timestep: Timestep):
    """ get the suffix as string for a timestep

    :param timestep: instance of class lbmpy.advanced_streaming.Timestep
    :return: either "even", "odd" or an empty string
    """
    return ("_" + str(timestep)) if timestep != Timestep.BOTH else ''

