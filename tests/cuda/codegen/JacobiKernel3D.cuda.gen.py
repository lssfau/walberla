from pystencils_walberla import Sweep

k = Sweep(dim=3)

src = k.field("f1")
dst = k.temporaryField(src)
h = k.constant("h")

rhs = (src[1,0,0] + src[-1,0,0] + src[0,1,0] + src[0, -1, 0] + src[0, 0, 1] + src[0, 0 , -1] ) / (6 * h**2)
k.addEq(dst[0,0,0], rhs)

k.generate()
