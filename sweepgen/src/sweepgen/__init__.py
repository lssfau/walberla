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

from .api import real_t, Vector3, GhostLayerFieldPtr, glfield, IBlockPtr
from .sweep import Sweep
from .build_config import WalberlaBuildConfig, get_build_config

__all__ = [
    "real_t",
    "Vector3",
    "GhostLayerFieldPtr",
    "glfield",
    "IBlockPtr",
    "Sweep",
    "WalberlaBuildConfig",
    "get_build_config",
]

from . import _version

__version__ = _version.get_versions()["version"]
