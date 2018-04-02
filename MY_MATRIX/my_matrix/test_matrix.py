import cProfile, pstats, io
from my_matrix import my_matrix, _11
import gmpy2
from gmpy2 import mpfr, mpc
from aux_functions import conj, NUMERICAL_TYPES

ctx = gmpy2.get_context()

ctx.precision = 200

gmpy2.set_context(ctx)

N = 5

qqq = my_matrix.crand(N,N)

print(qqq.get_dps())

qqq.set_print_digits(10)

print(qqq)


pr = cProfile.Profile()
pr.enable()

qqq,_ = qqq.get_Hessenberg_form(True)

print(qqq.round_to_tol())

pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print(s.getvalue())





#
# qqq = mpc("1+1j")
#
# print(qqq)
