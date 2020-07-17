import pint
from pint import UnitRegistry


class PintUnitConverter:

    def __init__(self, ureg=UnitRegistry(), **kwargs):
        """Possible arguments (all in SI units):
                * l_m  (lattice meter, same as dx)
                * l_s  (lattice second, same as dt)
                * l_kg (lattice kilogram)
                * same for other SI units

            Example:
                conv = PintUnitConverter( l_m = 1e-6, l_s = 1e-7 )

        """
        self.ureg = ureg

        # if not 'l_m'   in kwargs: kwargs['l_m']   = 1
        # if not 'l_s'   in kwargs: kwargs['l_s']   = 1
        # if not 'l_kg'  in kwargs: kwargs['l_kg']  = 1
        # if not 'l_K'   in kwargs: kwargs['l_K']   = 1
        # if not 'l_mol' in kwargs: kwargs['l_mol'] = 1
        # if not 'l_A'   in kwargs: kwargs['l_A']   = 1
        # if not 'l_cd'  in kwargs: kwargs['l_cd']  = 1

        self.ureg.define('l_radian = radian')
        self.define_lattice_units(**kwargs)

    def define_lattice_units(self, **kwargs):
        """Defines lattice units"""
        for key, value in kwargs.items():
            if key == 'l_m':
                self.ureg.define('l_meter    = %g * m            = l_m = l_meter   ' % value)
            elif key == 'l_s':
                self.ureg.define('l_second   = %g * s            = l_s = l_sec     ' % value)
            elif key == 'l_kg':
                self.ureg.define('l_kilogram = %g * kg           = l_kg            ' % value)
            elif key == 'l_K':
                self.ureg.define('l_kelvin   = %g * K  ;offset:0 = l_K = l_degK    ' % value)
            elif key == 'l_mol':
                self.ureg.define('l_mole     = %g * mol          = l_mol           ' % value)
            elif key == 'l_A':
                self.ureg.define('l_ampere   = %g * A            = l_A = l_amp     ' % value)
            elif key == 'l_cd':
                self.ureg.define('l_candela  = %g * cd           = l_cd = l_candle ' % value)

    def to_si_units(self, value):
        """SI= Systeme international d'unites, International System of Units"""
        if not hasattr(value, 'units'):
            return value
        # convert to base units (m,s,g,K,mol,A,cd)
        value.ito_base_units()
        return value

    def to_si_units_mag(self, value):
        if not hasattr(value, 'units'):
            return value
        else:
            return self.to_si_units(value).magnitude

    def to_sl_units(self, value):
        """SL = Lattice System of Units"""

        if not hasattr(value, 'units'):
            return value
        # convert to SI units (m,s,kg,K,mol,A,cd)
        value = self.to_si_units(value)

        latticeUnits = {
            '[length]': 'l_m',
            '[time]': 'l_s',
            '[mass]': 'l_kg',
            '[temperature]': 'l_K',
            '[substance]': 'l_mol',
            '[current]': 'l_A',
            '[luminosity]': 'l_cd',
        }
        units = dict()
        for unit, power in value.dimensionality.items():
            units[latticeUnits[unit]] = power
        return value.to(units)

    def to_sl_units_mag(self, value):
        if not hasattr(value, 'units'):
            return value
        else:
            return self.to_sl_units(value).magnitude

    def conv_config(self, config):
        result = {}

        for key, value in config.items():
            if type(value) is dict:
                result[key] = self.conv_config(value)
            elif type(value) is self.ureg.Quantity:
                result[key] = self.to_sl_units_mag(value)
            elif type(value) is tuple:
                result[key] = tuple([self.to_sl_units_mag(e) for e in value])
            elif type(value) is list:
                result[key] = [self.to_sl_units_mag(e) for e in value]
            else:
                result[key] = value

        return result

    def __conf_config_strings_to_units(self, config):
        result = {}

        for key, value in config.items():
            if type(value) is dict:
                result[key] = self.__conf_config_strings_to_units(value)
            elif key == 'type' or key == 'direction':
                result[key] = value
            elif type(value) is tuple:
                try:
                    result[key] = tuple([self.ureg.Quantity(e) for e in value])
                except (pint.UndefinedUnitError, ValueError, TypeError):
                    result[key] = value
            else:
                try:
                    result[key] = self.ureg.Quantity(value)
                except (pint.UndefinedUnitError, ValueError, TypeError):
                    result[key] = value

        return result

    def conv_config_unit_strings(self, config):
        result = self.__conf_config_strings_to_units(config)
        return self.conv_config(result)
