import cProfile, pstats, io
from my_matrix import my_matrix, _11
#import gmpy2
from gmpy2 import mpfr, mpc 


#ctx = gmpy2.get_context()
#
#ctx.precision = 200
#
#gmpy2.set_context(ctx)


N = 5

qqq = my_matrix.crand(N,N)

qqq.set_print_digits(12)
#
#print(qqq)
#
#www,_ = qqq.get_Hessenberg_form(True)
#
#print(www.round_to_tol())

# print(qqq.get_dps())
#
# qqq.set_print_digits(10)
#
print(qqq.round_to_tol())
#
#
#pr = cProfile.Profile()
#pr.enable()

for _ in range(1):
    www,_ = qqq.get_Hessenberg_form(True)
    print(www.round_to_tol())

#pr.disable()
#s = io.StringIO()
#sortby = 'cumulative'
#ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
#ps.print_stats()
#print(s.getvalue())
#
#
#
#
#
# #
# # qqq = mpc("1+1j")
# #
# # print(qqq)
