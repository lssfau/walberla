# -*- coding: utf-8 -*-

from ..utility import find, TerminalColor, generate_file


def create_particle_property(name, type, access="grs", defValue="", syncMode="ALWAYS"):
    """
    Parameters
    ----------
    name : str
       name of the property
    type : str
       type of the property
    access : str
       'g' for getter (getName)
       'r' for reference (getNameRef)
       's' for setter (setName)
       any combination is possible
    defValue : str
       default value the property should be initialized with
    syncMode : str
       'NEVER', this property does not have to be synced
       'ON_GHOST_CREATION', this property must be synced on creation
       'ON_OWNERSHIP_CHANGE', this property must be synced when the ownership changes
       'ALWAYS', this property has to be synced in every iteration
    """

    if (type == 'bool'):
        raise RuntimeError("Due to flaws in the implementation of std::vector<bool>, bool is not supported as a "
                           "property type! Please use char instead.")

    # sort access specifier and remove duplicates
    foo = "".join(sorted(access))
    access = ''.join([foo[i] for i in range(len(foo) - 1) if foo[i + 1] != foo[i]] + [foo[-1]])

    for acc in access:
        if not (acc in ["g", "s", "r"]):
            raise RuntimeError(f"{acc} is not a valid access specifier in {access}")

    if (syncMode not in ["NEVER", "ON_GHOST_CREATION", "ON_OWNERSHIP_CHANGE", "ALWAYS"]):
        raise RuntimeError(
            f"{TerminalColor.RED}{syncMode} is no valid sync for property: {name}{TerminalColor.DEFAULT}")

    return {'name': name, 'type': type, 'access': access, 'defValue': defValue, 'syncMode': syncMode}


class ParticleStorage():
    def __init__(self):
        self.context = {'includes': [], 'properties': []}

        self.add_include("mesa_pd/data/Flags.h")

        self.add_property("uid", "walberla::id_t", defValue="UniqueID<data::Particle>::invalidID()", syncMode="ALWAYS")
        self.add_property("position", "walberla::mesa_pd::Vec3", defValue="real_t(0)", syncMode="ALWAYS")
        self.add_property("interactionRadius", "walberla::real_t", defValue="real_t(0)", syncMode="ON_GHOST_CREATION")
        self.add_property("flags", "walberla::mesa_pd::data::particle_flags::FlagT", defValue="",
                          syncMode="ON_GHOST_CREATION")
        self.add_property("owner", "int", defValue="-1", syncMode="ON_GHOST_CREATION")
        self.add_property("ghostOwners", "std::unordered_set<walberla::mpi::MPIRank>", defValue="",
                          syncMode="ON_OWNERSHIP_CHANGE")

    def set_shapes(self, *shapes):
        self.context['shapes'] = list(shapes)

    def add_property(self, name, type, access="grs", defValue="", syncMode="ALWAYS"):
        prop = find(lambda x: x['name'] == name, self.context['properties'])
        if (prop is None):
            # print(f"{TerminalColor.GREEN} creating particle property: {name} {TerminalColor.DEFAULT}")
            self.context['properties'].append(
                create_particle_property(name, type, access=access, defValue=defValue, syncMode=syncMode))
        else:
            if not (prop['type'] == type and prop['name'] == name and prop['defValue'] == defValue):
                new_prop = create_particle_property(name, type, defValue=defValue, syncMode=syncMode)
                raise RuntimeError(
                    f"{TerminalColor.RED} property definition differs from previous one:\n"
                    f"PREVIOUS {prop}\n"
                    f"NEW {new_prop} {TerminalColor.DEFAULT}")
            print(f"{TerminalColor.YELLOW} reusing particle property: {name} {TerminalColor.DEFAULT}")

    def add_include(self, include):
        if (include in self.context['includes']):
            print(f"{TerminalColor.YELLOW} reusing particle include: {include} {TerminalColor.DEFAULT}")
        else:
            # print(f"{TerminalColor.GREEN} creating particle include: {include} {TerminalColor.DEFAULT}")
            self.context['includes'].append(include)

    def print(self):
        print("=" * 90)
        print("Creating Particle Datastructure:")
        print("")
        print("{0: <20}{1: <30}{2: <20}{3: <10}".format("Type", "Name", "Def. Value", "SyncMode"))
        print("=" * 90)
        for prop in self.context['properties']:
            print("{0: <20.19}{1: <30.29}{2: <20.19}{3: <10.9}".format(prop['type'], prop['name'], prop['defValue'],
                                                                       prop['syncMode']))
        print("=" * 90)

    def get_context(self):
        return self.context

    def generate(self, module):
        ctx = {'module': module, **self.context}

        generate_file(module['module_path'], 'data/ParticleStorage.templ.h', ctx)
        generate_file(module['module_path'], 'data/ParticleAccessor.templ.h', ctx)
        generate_file(module['module_path'], 'common/ParticleFunctions.templ.h', ctx)
