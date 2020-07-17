# encoding: utf-8


import gdb
import re


class static:
    "Creates a 'static' method"

    def __init__(self, function):
        self.__call__ = function


walberla_pretty_printers = []


def register_pretty_printer(pretty_printer):
    "Registers a Pretty Printer"
    walberla_pretty_printers.append(pretty_printer)
    return pretty_printer


@register_pretty_printer
class FieldPrinter:
    "Pretty Printer for walberla::Fields"
    regexField = re.compile('walberla::field::Field')
    regexGlField = re.compile('walberla::field::GhostLayerField')

    @static
    def supports(typename):
        return FieldPrinter.regexField.search(typename) or FieldPrinter.regexGlField.search(typename)

    def __init__(self, typename, value):
        self.typename = typename
        self.value = value
        if FieldPrinter.regexGlField.search(typename):
            self.typeStr = "GhostLayerField"
            self.ghostLayers = value.type.template_argument(2)
        else:
            self.typeStr = "Field"

    def to_string(self):
        ret = "%s<%s> of size (x,y,z,f)=(%d,%d,%d,%d) alloc (%d,%d,%d,%d)" % (self.typeStr,
                                                                              self.value.type.template_argument(0),
                                                                              self.value['xSize_'],
                                                                              self.value['ySize_'],
                                                                              self.value['zSize_'],
                                                                              self.value.type.template_argument(1),
                                                                              self.value['xAllocSize_'],
                                                                              self.value['yAllocSize_'],
                                                                              self.value['zAllocSize_'],
                                                                              self.value['fAllocSize_'])
        return ret

    def children(self):
        if self.value['layout_'] == gdb.parse_and_eval("walberla::field::fzyx"):
            arrStr = "field(f,z,y,x)"
            [s1, s2, s3, s4] = [self.value['xSize_'], self.value['ySize_'], self.value['zSize_'],
                                self.value.type.template_argument(1)]
        else:
            arrStr = "field(z,y,x,f)"
            [s1, s2, s3, s4] = [self.value.type.template_argument(1), self.value['xSize_'], self.value['ySize_'],
                                self.value['zSize_']]

        fieldType = self.value.type.template_argument(0)
        yield (arrStr, self.value["values_"].dereference().cast(
            fieldType.array(s1 - 1).array(s2 - 1).array(s3 - 1).array(s4 - 1)))

        memberToDisplay = ['layout_',
                           'xSize_', 'ySize_', 'zSize_',
                           'xAllocSize_', 'yAllocSize_', 'zAllocSize_',
                           'ffact_', 'xfact_', 'yfact_', 'zfact_',
                           'xOff_', 'yOff_', 'zOff_']

        if hasattr(self, 'ghostLayers'):
            yield ("GhostLayers", self.ghostLayers)

        for member in memberToDisplay:
            yield (member, self.value[member])


def find_pretty_printer(value):
    "Find a pretty printer suitable for value"
    type = value.type

    if type.code == gdb.TYPE_CODE_REF:
        type = type.target()

    type = type.unqualified().strip_typedefs()

    typename = type.tag
    if typename is None:
        return None

    for pretty_printer in walberla_pretty_printers:
        if pretty_printer.supports(typename):
            return pretty_printer(typename, value)

    return None


def register_walberla_printers(obj):
    "Register walberla Pretty Printers."
    if obj is None:
        obj = gdb
    obj.pretty_printers.append(find_pretty_printer)
