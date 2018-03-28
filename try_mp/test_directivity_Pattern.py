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

N_pat = 361

N_vector = [4, 10, 20, 30, 40]

A_length = mp.mpf("3")*mp.mpf("0.5")


i = mp.mpc(0,1)

k = mp.mpf("2") * mp.pi


theta_scan_vector = [mp.pi / mp.mpf("2")]

theta_pat = mp.linspace(mp.mpf("0"), mp.pi, N_pat)

# theta_scan_vector = [mp.mpf("0")]

# theta_pat = mp.linspace(-mp.pi/mp.mpf("2"), mp.pi/mp.mpf("2"), N_pat)



Dir_opt = make_N_dim_list([len(theta_scan_vector), len(N_vector), N_pat])



figure_handle = plt.figure()

axes_handle = figure_handle.add_subplot(211)
        
a_axes1 = figure_handle.add_subplot(223)

a_axes2 = figure_handle.add_subplot(224)

color_range.reset()


for n_scan in range(len(theta_scan_vector)):
    
    theta_scan = theta_scan_vector[n_scan]
    
    print("Solving for n_scan=" + str(n_scan))

    for n_N in range(len(N_vector)):
        
        N = N_vector[n_N]
        
        print("  Solving for N=" + str(N))
        
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
        
        a_opt = a_opt / np.amax(a_opt.apply(mp.fabs))
        
        local_color = color_range.make_color()
        
        a_axes1.plot(mp.linspace(-A_length / 2, A_length / 2, N), a_opt.apply(mp.fabs), color=local_color, label="N={}".format(N))
        
        a_axes2.plot(mp.linspace(-A_length / 2, A_length / 2, N), a_opt.apply(lambda x: deg(mp.arg(x))), color=local_color, label="N={}".format(N))
        
        for n_pat in range(N_pat):
            
            print("      Solving for n_pat=" + str(n_pat))
            
            for n in range(N):
                
                values_up[n] = mp.exp(i*k*(zL_diff[n] - zL_diff[0])*mp.cos(theta_pat[n_pat]))
            
            for r in range(N):
                
                for c in range(r, N):
            
                    M_up[r,c] = values_up[c - r]
                         
                    M_up[c,r] = mp.conj(values_up[c - r])
             
            Dir_opt[n_scan][n_N][n_pat] = dB10((mp.mpf("4") * mp.pi * (a_opt.H * M_up * a_opt) * ((a_opt.H * M_down * a_opt) ** -1))[0,0])

        axes_handle.plot(deg(theta_pat), Dir_opt[n_scan][n_N], color=local_color, label="N={}".format(N))
        
axes_handle.set_xlabel("Theta scanning, Deg")

axes_handle.set_ylabel("Directivity, dB10")

axes_handle.set_title("Optimal Directivity")

axes_handle.set_ylim([-20, 40])

axes_handle.set_xlim([0, 180])

# axes_handle.set_xlim([-90, 90])

axes_handle.legend()

axes_handle.grid(True, linewidth=1)

a_axes1.grid(True, linewidth=1)

a_axes1.legend()

a_axes1.set_xlabel("Array Length")

a_axes1.set_ylabel("Abs")

a_axes1.set_title("Optimal Feed, Amplitude")

a_axes2.grid(True, linewidth=1)

a_axes2.legend()
         
a_axes2.set_xlabel("Array Length")

a_axes2.set_ylabel("Phase")

a_axes2.set_title("Optimal Feed, Phase")


plt.show()      
        
    
    