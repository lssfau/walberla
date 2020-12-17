import unittest


class UnitConversionTest(unittest.TestCase):

    def testExtractLatticeFactors(self):
        try:
            import pint, sympy  # noqa: F401
        except ImportError:
            print("Skipping unit conversion test since pint or sympy module not available")
            return

        from waLBerla.tools.lbm_unitconversion import extractLatticeFactors, computeLatticeFactors
        example = {'rho': '4000kg/m**3',
                   'l_m': '5um',
                   'tau': '0.6',
                   'nu': '1e-6m/s**2',
                   'l_K': '1mK'}

        # from pprint import pprint
        lf1 = extractLatticeFactors(example)
        lf2 = computeLatticeFactors(rho=4000, l_m=5e-6, tau=0.6, nu=1e-6, l_K=1e-3)
        self.assertEqual(lf1, lf2)

        print(computeLatticeFactors(rho=4000, l_m=5e-6, tau=0.6, l_s=1e-8, l_K=1e-3, l_A=1e-2))
        print(computeLatticeFactors(rho=4000, l_m=5e-6, l_s=1e-8, l_K=1e-3))
        print(computeLatticeFactors(rho=4000, l_m=5e-6, time=1e-6, l_time=50))

    def testUnitConverter(self):
        try:
            import pint, sympy  # noqa: F401
        except ImportError:
            print("Skipping unit conversion test since pint or sympy module not available")
            return

        from waLBerla.tools.lbm_unitconversion import PintUnitConverter

        conv = PintUnitConverter(l_m=5e-6, l_s=1e-8, l_kg=5e-13, l_K=1e-3, l_A=1)

        ureg = conv.ureg

        rho = 4000 * ureg.kg / ureg.m ** 3
        vol = 1.0 * ureg.mm ** 3
        mass = rho * vol
        T = ureg.Quantity('1934 K')
        I = ureg.Quantity('1.9 mA')
        print(ureg.Quantity('700.0 J/(kg*K)'))

        self.assertAlmostEqual(rho.to('l_kg/l_m^3').magnitude, 1.0, places=8)
        self.assertAlmostEqual(mass.to('l_kg').magnitude, 8000000, places=8)
        self.assertAlmostEqual(vol.to('l_m**3').magnitude, 8000000, places=8)
        self.assertAlmostEqual(vol.to('l_m**3').magnitude, 8000000, places=8)
        self.assertAlmostEqual(I.to('l_A').magnitude, 0.0019, places=8)

        self.assertAlmostEqual(T.to('l_K').magnitude, 1934000, places=8)
        self.assertAlmostEqual(mass.to_base_units().magnitude, 4e-6, places=8)

        u = conv.ureg
        c1 = {
            'Parameters': {
                'unitlessParam': 2,
                'param2': 4e-6 * u.kg,
            }
        }
        c2 = {
            'Parameters': {
                'unitlessParam': 2,
                'param2': '4e-6 kg',
            }
        }
        c1 = conv.conv_config(c1)
        c2 = conv.conv_config_unit_strings(c2)
        self.assertEqual(c1, c2)


if __name__ == '__main__':
    unittest.main()
