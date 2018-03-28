import numpy as np
from matplotlib import pyplot as plt
from mpmath import mp
from aux_functions import make_N_dim_list, deg, dB10, apply_to_iterable, color_range




def local_sinc(x):
    if x == 0:
        return 1
    else:
        return mp.sin(x) / x
        


mp.dps = 100


# N_vector = [4, 10]#, 20, 30, 40]

N_vector = [12, 13, 14, 16, 18, 20, 23]

A_length = mp.mpf("11")*mp.mpf("0.5")


i = mp.mpc(0,1)

k = mp.mpf("2") * mp.pi


N_scan = 181

theta_scan_vector = mp.linspace(mp.mpf("0"), mp.pi, N_scan)



Dir_opt = make_N_dim_list([len(N_vector), len(theta_scan_vector)])



figure_handle = plt.figure()

axes_handle = figure_handle.add_subplot(111)


color_range.reset()


for n_N in range(len(N_vector)):
    
    N = N_vector[n_N]
    
    print("  Solving for N=" + str(N))
    
    for n_scan in range(len(theta_scan_vector)):
        
        theta_scan = theta_scan_vector[n_scan]
    
        print("Solving for n_scan=" + str(n_scan))
        
        M_up = mp.zeros(N,N)
         
        M_down = mp.zeros(N,N)
         
        zL_diff = mp.linspace(mp.mpf("0"), A_length, N)
         
        values_up = mp.zeros(1, N)
         
        values_down = mp.zeros(1, N)
         
        for n in range(N):
            
            values_up[n] = mp.exp(i*k*(zL_diff[n] - zL_diff[0])*mp.cos(theta_scan))
             
            values_down[n] = mp.mpf("4")*mp.pi*local_sinc(k*(zL_diff[n] - zL_diff[0]))
             
        for r in range(N):
             
            for c in range(r, N):
                
                M_up[r,c] = values_up[c - r]
                         
                M_up[c,r] = mp.conj(values_up[c - r])
                 
                M_down[r,c] = values_down[c - r]
                 
                M_down[c,r] = mp.conj(values_down[c - r])
                 
        L_down, V_down = mp.eigh(M_down)

        L_aux, V_aux = mp.eigh(((V_down * mp.diag(L_down).apply(mp.sqrt))**-1) * M_up * ((mp.diag(L_down).apply(mp.sqrt) * V_down.H)**-1))

        a_opt = ((mp.diag(L_down).apply(mp.sqrt) * V_down.H)**-1) * V_aux[:, np.argmax(apply_to_iterable(L_aux, mp.re))]
        
        Dir_opt[n_N][n_scan] = dB10((mp.mpf("4") * mp.pi * (a_opt.H * M_up * a_opt) * ((a_opt.H * M_down * a_opt) ** -1))[0,0])
        
    axes_handle.plot(deg(theta_scan_vector), Dir_opt[n_N], color=color_range.make_color(), label="N={}".format(N))
        
plt.xlabel("Theta scanning, Deg")

plt.ylabel("Directivity, dB10")

plt.title("Optimal Directivity")

plt.ylim([10, 30])

plt.xlim([0, 180])

axes_handle.grid(True, linewidth=1)

axes_handle.legend()
         
plt.draw()

plt.show()
        
        
    
    