# This file is part of waLBerla. waLBerla is free software: you can
# redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.
#
# waLBerla is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
# for more details.
#
# You should have received a copy of the GNU General Public License along
# with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.

from pystencilssfg.lang import CppClass, AugExpr
from pystencilssfg.lang.types import cpptype
from pystencilssfg.ir import SfgClass


class GeneratedClassWrapperBase(CppClass):
    _class: SfgClass

    def __init_subclass__(cls) -> None:
        typename = cls._class.fqname
        cls.template = cpptype(typename)

    def ctor(self, **kwargs) -> AugExpr:
        for candidate_ctor in self._class.constructors():
            ctor_argnames = [param.name for param in candidate_ctor.parameters]
            if set(ctor_argnames) == set(kwargs.keys()):
                break
        else:
            raise Exception(
                f"No constructor of class {self._class.fqname} matches the argument names {kwargs.keys()}"
            )

        ctor_args = [kwargs[name] for name in ctor_argnames]
        return self.ctor_bind(*ctor_args)
