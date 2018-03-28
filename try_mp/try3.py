from mpmath import *
import cProfile, pstats, io
from my_matrix import my_matrix



mp.dps = 100

pr = cProfile.Profile()


pr.enable()

N = 30

my_matrix.set_print_digits(10)

qqq = my_matrix.crand(N,N)

qqq = qqq.H() * qqq
  
print(qqq)
   
print(qqq.get_Hessenberg_form().round_to_tol())

pr.disable()


s = io.StringIO()

sortby = 'cumulative'

ps = pstats.Stats(pr, stream=s).sort_stats(sortby)

ps.print_stats()

print(s.getvalue())