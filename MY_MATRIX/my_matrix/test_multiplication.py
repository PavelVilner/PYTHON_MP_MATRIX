import cProfile, pstats, io
from my_matrix import my_matrix

N = 2**6

N_test = 10

A = my_matrix.rand(N,N)
B = my_matrix.rand(N,N)

pr = cProfile.Profile()
pr.enable()

for _ in range(N_test):
    _ = my_matrix.simple_mult(A,B)
    
pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print('SIMPLE')
print(s.getvalue())


pr = cProfile.Profile()
pr.enable()

for _ in range(N_test):
    _ = my_matrix.recursive_mult(A,B)
    
pr.disable()
s = io.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()
print('RECURSIVE')
print(s.getvalue())


