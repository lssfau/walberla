import waLBerla


@waLBerla.callback("config")
def waLBerlaConfig():
    conf = {
        'testInt': 4,
        'testString': "someString",
        'testDouble': 42.42,
        '44242': 'ohoh_IntegerKey',
        'subBlock': {'subKey1': 'abc',
                     'subKey2': 'def'
                     }
    }
    return conf
