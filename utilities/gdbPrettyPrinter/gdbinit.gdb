python
# STL
from gdbPrettyPrinter.STLv6.printers import register_libstdcxx_printers
register_libstdcxx_printers (None)

# Boost
from gdbPrettyPrinter.boost_1_40.printers import register_boost_printers
register_boost_printers(None)

# waLBerla
from gdbPrettyPrinter.walberla.printers import register_walberla_printers
register_walberla_printers(None)

# Qt4
from gdbPrettyPrinter.qt4.printers import register_qt4_printers
register_qt4_printers (None)
end
