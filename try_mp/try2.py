from mpmath import *
import numpy as np
from my_matrix import my_matrix, _11

from aux_functions import apply_to_iterable, unwind_iterable, test_iterable,\
    S_to_Z, Z_to_S, color_range, mpf_to_str, mpc_to_str

from copy import deepcopy, copy
    
N = 3

mp.pretty = True

mp.dps = 30

# i = mp.mpc("0","1")
# # qqq = -1e-20 * mp.rand() + 1e-20 * mp.rand()*i
# # qqq = mp.mpf("0")
# # qqq = -5.969e19
# # qqq =  10e-2 + i*9e-2
# # qqq_str = mpc_to_str(qqq, 19)
# # 
# # print(qqq_str)
# # 
# # print(len(qqq_str))
# 
my_matrix.set_print_digits(10)

# for _ in range(1):
#    
#     qqq = my_matrix.crand(3,3)
#     
#     www = qqq[1:-1]
#        
#     q,w = qqq.codiagonalize(qqq.to_I())
#          
#     print(qqq)
#        
#     print(q)
#        
#     print(w)
#       
#     print(w * qqq)
#          
#     print((qqq * w).round_to_tol())
# 
# mp.nprint(mp.mpf("1.14815108861"), 10)

# qqq = mp.mpf('0.9999999999999998e0')
#  
# print(mpf_to_str(qqq, 20))
# 
# print(len(mpf_to_str(qqq, 20)))

# qqq = my_matrix([[1,2,3],[4,5,6],[7,8,10]])
# 
# www = qqq.to_np_matrix()
# 
# print(qqq)
# 
# print(www)
#  
# nprint(qqq.det(), 20)
# 
# print(np.linalg.det(www))

# qqq = my_matrix
# 
# www = my_matrix
# 
# print(qqq == www)
# 
# print(floordiv(3,2))

# qqq = my_matrix.crand(10,10)
#   
# print(qqq)
#    
# www = qqq.to_mp_matrix()
#    
# L,V = mp.eig(www)
#    
# print(qqq.eig(1e-17))
#    
# mp.nprint(L)
# www = qqq.renormalize_cols()
#  
# print(www.H()*www)

# qqq = my_matrix([[1,2]])
# 
# www = my_matrix([[5],[7]])
# 
# print(qqq)
# 
# print(www)
# print(qqq * www)

# qqq = my_matrix.rand(2,2)
# 
# www = copy(qqq)
# 
# qqq[0,0] = 2
# 
# print(qqq)
# 
# print(www)

qqq = my_matrix.rand(30,30)

print(qqq)

print(qqq.get_Hessenberg_form().zero_to_tol())

