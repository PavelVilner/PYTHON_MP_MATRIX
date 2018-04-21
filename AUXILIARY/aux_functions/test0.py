
import gmpy2 as gm

from gmpy2 import mpfr, mpc

from aux_functions import mpr_to_str

from random import random

ctx = gm.get_context()

# print(ctx)

ctx.precision = 500

ctx.subnormalize = True

print(ctx)

gm.set_context(ctx)

qqq = mpc(0.1, 0.1)

print(gm.sqrt(qqq))