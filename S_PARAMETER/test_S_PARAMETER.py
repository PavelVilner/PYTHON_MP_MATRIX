import numpy as np

from mpmath import *

from S_PARAMETER import S_PAREMETER

import matplotlib.pyplot as plt



plt.close('all')

f_axis = np.linspace(0, 10, 11)

data = np.random.random_sample((5,5,11))

Z0 = 50

# SP = S_PAREMETER(file_name="C:\\Users\\pavelv\\OneDrive - Mellanox\\PERSONAL\\MATLAB\\For_Emmanuel_Cohen\\LNA_S_param\\LNA_2CS_65nm_70GHz.s2p")

SP = S_PAREMETER(file_name="H:\\Users\\Pavel\\HACKATHON17\\FOR_AVNER\\VNA_MEASUREMENTS\\NEW_PADS_11_02_18\\50p_0.5mm.s4p")


# SP.detect_passivity(tolerance=1e-30)

SP.plot_entry([0,2])

SP.cut_freq(min_f=10e9, max_f=20e9)

SP.plot_entry([1,2], [3,0])


# SP.plot_entry([0,0], scale="dB20")

# SP.change_Z0(75)
# 
# SP.plot_entry([0,2], scale="dB20")
# 
plt.show()

# with mp.workdps(1):
#     
#     print(SP.data(1,1))
# 
#     print(SP.is_symmetric(tolerance=1e-30))
#     
#     print(SP.is_passive(tolerance=1e-30))
#     
#     SP.enforce_symmetry()
#     
#     SP.enforce_passivity()
#     
#     print(SP.is_symmetric(tolerance=1e-30))
#     
#     print(SP.is_passive(tolerance=1e-30))

 
# print(SP.is_passive())
# 
# SP.reorder([1, 2, 3, 4, 0])
# 
# # SP.convert_to_Differential()

