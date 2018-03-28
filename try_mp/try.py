from mpmath import *

mp.dps = 300

print("Hello world!")

print(mpf(2) ** mpf('0.5'))

import mpmath.libmp
print(mpmath.libmp.BACKEND)
