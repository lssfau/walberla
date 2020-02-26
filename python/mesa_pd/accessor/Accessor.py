# -*- coding: utf-8 -*-


def create_access(name, type, access):
    """requires that a certain property is accessible

    Parameters
    ----------
    name : str
       name of the property requested
    type : str
       type of the requested property
    access : str
       'g' for getter (getName)
       'r' for reference (getNameRef)
       's' for setter (setName)
       any combination is possible

    Example
    -------
    create_access("position", "walberla::mesa_pd::Vec3", "sg")
    """

    for acc in access:
        if not (acc in ["g", "s", "r"]):
            raise RuntimeError(f"{acc} is not a valid access specifier in {access}")

    return {'name': name, 'type': type, 'access': access}
