import math


def calculate_parameters_taylor_bubble(reference_length=128,
                                       reference_time=16000,
                                       density_heavy=1.0,
                                       d1=0.0254,
                                       d2=0.0127):
    r"""
    Calculate the simulation parameters for a rising Taylor bubble. The calculation can be found in
    10.1016/S0009-2509(97)00210-8 by G. Das

    Args:
        reference_length: chosen reference length
        reference_time: chosen reference time
        density_heavy: density of the heavier fluid
        d1: diameter of the outer tube
        d2: diameter of the inner cylinder
    """

    water_rho = 998  # kg/m3
    air_rho = 1.2047  # kg/m3
    surface_tension = 0.07286  # kg/s2
    water_mu = 1.002e-3  # kg/ms

    water_nu = water_mu / water_rho  # m2/s
    air_mu = 1.8205e-5  # kg/ms
    air_nu = air_mu / air_rho  # m2/s
    gravity = 9.81  # m/s2

    dh = d1 - d2
    dr = d1 / d2
    de = d1 + d2
    # ur = 0.1695  # (0.28913, 0.23882, 0.1695)

    inverse_viscosity_number = math.sqrt((water_rho - air_rho) * water_rho * gravity * dh ** 3) / water_mu
    bond_number = (water_rho - air_rho) * gravity * dh ** 2 / surface_tension
    morton_number = gravity * water_mu ** 4 * (water_rho - air_rho) / (water_rho ** 2 * surface_tension ** 3)

    d = reference_length / dr

    density_light = 1.0 / (water_rho / air_rho)
    dh = reference_length - d
    g = dh / reference_time ** 2

    mu_h = math.sqrt((density_heavy - density_light) * density_heavy * g * dh ** 3) / inverse_viscosity_number
    mu_l = mu_h / (water_mu / air_mu)

    dynamic_viscosity_heavy = mu_h / density_heavy
    dynamic_viscosity_light = mu_l / density_light

    relaxation_time_heavy = 3 * dynamic_viscosity_heavy
    relaxation_time_light = 3 * dynamic_viscosity_light

    sigma = (density_heavy - density_light) * g * dh ** 2 / bond_number

    parameters = {
        "inverse_viscosity_number": inverse_viscosity_number,
        "bond_number": bond_number,
        "morton_number": morton_number,
        "density_light": density_light,
        "dynamic_viscosity_heavy": dynamic_viscosity_heavy,
        "dynamic_viscosity_light": dynamic_viscosity_light,
        "relaxation_time_heavy": relaxation_time_heavy,
        "relaxation_time_light": relaxation_time_light,
        "gravitational_acceleration": -g,
        "surface_tension": sigma
    }
    return parameters
