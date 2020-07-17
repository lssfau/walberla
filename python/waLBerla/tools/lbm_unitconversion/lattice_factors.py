import sympy as sp
from waLBerla.tools.lbm_unitconversion import PintUnitConverter


def computeLatticeFactors(constant=True, **kwargs):
    """
    Computes lattice factors with given physical input values

    :param kwargs: dictionary with symbol name as key and the physical value as value
    :return: lattice factors as dictionary with name, value
    """

    def __create_symbols(d, *args):
        """
        Create sympy symbols and add them to the dictionary d if they are not already inside

        :param d: dictionary where all symbols are stored
        :param args: keys for all symbols to create
        """
        for arg in args:
            if arg not in d.keys():
                d[arg] = sp.symbols(arg, positive=True)

    def __append_equation(d, eq):
        if eq not in d:
            d.append(eq)

    val = {}
    eqn = []

    for k, v in kwargs.items():
        if k in ['time', 'l_time']:
            __create_symbols(val, 'time', 'l_time', 'l_s')
            __append_equation(eqn, sp.Eq(val['time'], val['l_time'] * val['l_s']))
        elif k in ['size', 'l_size']:
            __create_symbols(val, 'size', 'l_size', 'l_m')
            __append_equation(eqn, sp.Eq(val['size'], val['l_size'] * val['l_m']))
        elif k in ['omega', 'tau', 'nu', 'eta', 'l_nu', 'rho']:
            if constant:
                __create_symbols(val, 'omega', 'tau', 'nu', 'eta', 'l_nu', 'rho', 'l_rho', 'l_m', 'l_kg', 'l_s')
                __append_equation(eqn, sp.Eq(val['nu'], val['eta'] / val['rho']))
            else:
                __create_symbols(val, 'omega', 'tau', 'nu', 'l_nu', 'rho', 'l_rho', 'l_m', 'l_kg', 'l_s')
            __append_equation(eqn, sp.Eq(val['omega'], 1 / val['tau']))
            __append_equation(eqn, sp.Eq(val['tau'], 3 * val['l_nu'] + 0.5))
            __append_equation(eqn, sp.Eq(val['nu'], val['l_nu'] * val['l_m'] ** 2 / val['l_s']))
            __append_equation(eqn, sp.Eq(val['l_rho'], 1))
            __append_equation(eqn, sp.Eq(val['rho'], val['l_rho'] * val['l_kg'] / val['l_m'] ** 3))
        elif k in ['omega_sol', 'tau_sol', 'a_sol', 'l_a_sol']:
            __create_symbols(val, 'omega_sol', 'tau_sol', 'a_sol', 'l_a_sol', 'l_m', 'l_s')
            __append_equation(eqn, sp.Eq(val['omega_sol'], 1 / val['tau_sol']))
            __append_equation(eqn, sp.Eq(val['tau_sol'], 3 * val['l_a_sol'] + 0.5))
            __append_equation(eqn, sp.Eq(val['a_sol'], val['l_a_sol'] * val['l_m'] ** 2 / val['l_s']))
        elif k in ['omega_liq', 'tau_liq', 'a_liq', 'l_a_liq']:
            __create_symbols(val, 'omega_liq', 'tau_liq', 'a_liq', 'l_a_liq', 'l_m', 'l_s')
            __append_equation(eqn, sp.Eq(val['omega_liq'], 1 / val['tau_liq']))
            __append_equation(eqn, sp.Eq(val['tau_liq'], 3 * val['l_a_liq'] + 0.5))
            __append_equation(eqn, sp.Eq(val['a_liq'], val['l_a_liq'] * val['l_m'] ** 2 / val['l_s']))
        else:
            __create_symbols(val, k)

        if k in val:
            __append_equation(eqn, sp.Eq(val[k], v))

    solutions = sp.solve(eqn, tuple(val.values()), dict=True, force=True)

    if len(solutions) == 1:
        # print(solutions[0])
        values = dict()
        for var in val.keys():
            if val[var] in solutions[0]:
                values[var] = solutions[0][val[var]]
        return values
    else:
        print(eqn)
        print(val.values())
        print("\n".join(str(item) for item in solutions))


def extractLatticeFactors(config, constant=True):
    """
    takes config dictionary and extracts lattice factors
    """
    puc = PintUnitConverter()

    def __scan_dict(src, dst):
        """
        Scans the whole src dictionary for special keys and store
        the magnitudes of the values after conversion to SI units

        :param src: src dictionary to walk through
        :param dst: dst dictionary, where all converted quantities are stored
        """
        for key, value in src.items():
            if type(value) is dict:
                __scan_dict(value, dst)
            elif key in ['omega', 'tau', 'nu', 'eta', 'l_nu', 'rho', 'time', 'l_time', 'size', 'l_size',
                         'l_m', 'l_s', 'l_kg', 'l_K', 'l_mol', 'l_A', 'l_cd']:
                dst[key] = puc.to_si_units(puc.ureg.Quantity(value)).magnitude

    si_config = {}
    __scan_dict(config, si_config)
    return computeLatticeFactors(constant, **si_config)
