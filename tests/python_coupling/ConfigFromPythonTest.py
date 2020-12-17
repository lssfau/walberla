import waLBerla as wlb


class Scenario:
    def __init__(self, number):
        self.testInt = number
        self.testString = "someString"
        self.testDouble = 42.42
        self.testVector = (0.5, 0.5, 0.7)
        self.testBool = False

    @wlb.member_callback
    def config(self):
        return {
            'DomainSetup': {
                'testInt': self.testInt,
                'testDouble': self.testDouble,
                'testString': self.testString,
                'testVector': self.testVector,
                'testBool': self.testBool
            }
        }


scenarios = wlb.ScenarioManager()
scenarios.add(Scenario(4))
scenarios.add(Scenario(5))
