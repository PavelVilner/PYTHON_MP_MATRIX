from aux_functions import convert_to_mp_matrix, is_scalar, is_vector, mp_abs, S_to_T, S_to_Z_simple, S_to_Y_simple, S_to_H,\
                          T_to_S, Z_to_S_simple, Y_to_S_simple, H_to_S, rad, is_even, dB20, dB10, deg,\
                          z0_from_args, y0_from_args, get_dps_by_tol
    
import os

import numpy as np

# import mpmath as mp

from mpmath import mp

import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import filedialog
from math import floor
from itertools import compress
from copy import deepcopy



class S_PAREMETER:
    
    default_tolerance = mp.mpf("1e-16")
    
    def __init__(self, **kwargs):
        if len(kwargs) == 0:
            print("Initializing from a file dialog")
            root = tk.Tk()
            root.withdraw()
            file_name = filedialog.askopenfilename(title="Select a Touchstone file", filetypes=(("Touchstone 1 files","*.s*p *.S*P"),("all files","*.*")))
            self.__read_S_parameter_from_file(file_name)
        elif len(kwargs) == 1 and "file_name" in kwargs:
            print("Initializing from a file name \"" + kwargs["file_name"] + "\"")
            self.__read_S_parameter_from_file(kwargs["file_name"])
        elif len(kwargs) >= 2 and "f_axis" in kwargs and "data" in kwargs:
            print("Initializing from a data set")
            self.__f_axis = convert_to_mp_matrix(kwargs["f_axis"])
            self.__N_freq = len(self.__f_axis)
            self.__data = convert_to_mp_matrix(kwargs["data"])
            self.__N_ports = self.__data[0].rows
            if "Z0" in kwargs:
                self.__Z0 = convert_to_mp_matrix(kwargs["Z0"])
            else:
                self.__Z0 = convert_to_mp_matrix(50)
            if "type" in kwargs:
                self.__type = kwargs["type"]
            else:
                self.__type = "S"
            self.__mode = "SE"
        else: 
            raise ValueError("Wrong arguments to S_PARAMETER constructor!")
        self.__check_data_validity_and_reshape()
        self.__throughs = []
        self.detect_throughs()
        self.__symmetric = []
        self.detect_symmetry()
        self.__passive = []
#         self.detect_passivity()
        self.__diff_map = []
        
    def __str__(self):
        if np.size(self.__throughs) == 0:
            throughs_str = "UNKNOWN"
        else:
            throughs_str = str(np.shape(self.__throughs))
            
        if np.size(self.__symmetric) == 0:
            symmetric_str = "UNKNOWN"
        else:
            symmetric_str = str(np.shape(self.__symmetric))
        if np.size(self.__passive) == 0:
            passive_str = "UNKNOWN"
        else:
            passive_str = str(np.shape(self.__symmetric))
        if np.size(self.__diff_map) == 0:
            diff_map_str = "UNKNOWN"
        else:
            diff_map_str = str(np.shape(self.__diff_map))
        return "S_PARAMETER:" +\
               "\n BASIC:" +\
               "\n   f-axis: Array of " + str(len(self.__f_axis)) + " frequency points"\
               "\n   data: List of " + str(len(self.__data)) + " " +str(self.__N_ports) + "x" + str(self.__N_ports) + " matrices"\
               "\n   Z0: Array of " + str(len(self.__Z0)) + " impedance values"\
               "\n   type: " + self.__type +\
               "\n EXTENDED:" +\
               "\n   mode: " + self.__mode +\
               "\n   throughs: " + throughs_str + \
               "\n   symmetric: " + symmetric_str + \
               "\n   passive: " + passive_str + \
               "\n   diff. map: " + diff_map_str +\
               "\n   default tolerance: " + str(self.default_tolerance) +\
               "\n"
                
    def __check_data_validity_and_reshape(self):
        if len(self.__data) != self.__N_freq:
            raise ValueError("Length of f-axis is not equal to the data f-size!")  
        if len(self.__Z0) == 1:
            self.__Z0 = self.__Z0 * mp.matrix(1, self.__N_ports)
        elif len(self.__Z0) != self.__N_ports:
            raise ValueError("Z0 value must be a scalar or a vector of a length equal to a number of ports!")

    def detect_throughs(self, **kwargs):
        if self.__type != "S":
            raise Exception("Must be in S-form for the throuhgs detection to be meaningful!")
        if kwargs and "tolerance" in kwargs:
            tol = mp.mpf(kwargs["tolerance"])
        else:
            tol = mp.mpf("0.9")
        self.__throughs = []
        with mp.workdps(get_dps_by_tol(tol)):
            for n_row in range(self.__data[0].rows):
                for n_col in range(self.__data[0].cols):
                    if mp.fabs(self.__data[0][n_row, n_col]) >= tol:
                        self.__throughs.append([n_row, n_col])
                    
    def detect_symmetry(self, **kwargs):
        if self.__type != "S":
            raise Exception("Must be in S-form for the symmetry detection to be meaningful!")
        self.__symmetric = list([False]) * self.__N_freq
        if kwargs and "tolerance" in kwargs and is_scalar(kwargs["tolerance"]):
            tol = mp.mpf(kwargs["tolerance"])
        else:
            tol = self.default_tolerance
        with mp.workdps(get_dps_by_tol(tol)):
            for n in range(len(self.__data)):
                if max(mp_abs(self.__data[n] - self.__data[n].T)) <= tol:
                    self.__symmetric[n] = True
    
    def symmetric(self):
        if not self.__symmetric:
            self.detect_symmetry()
        return self.__symmetric()
    
    def is_symmetric(self, **kwargs):
        if not self.__symmetric or kwargs:
            self.detect_symmetry(**kwargs)
        if all(self.__symmetric):
            return True
        return False
      
    def enforce_symmetry(self, **kwargs):
        if self.__type != "S":
            raise Exception("Must be in S-form for the symmetry enforcement to be meaningful!")
        if kwargs and "tolerance" in kwargs and is_scalar(kwargs["tolerance"]):
            tol = mp.mpf(kwargs["tolerance"])
        else:
            tol = self.default_tolerance
        with mp.workdps(get_dps_by_tol(tol)):
            for n in range(len(self.__data)):
                self.__data[n] = (self.__data[n] + self.__data[n].T) / mp.mpf("2")
        self.detect_symmetry(**kwargs)
    
    def detect_passivity(self, **kwargs):
        if self.__type != "S":
            raise Exception("Must be in S-form for the passivity detection to be meaningful!")
        self.__passive = list([False]) * self.__N_freq
        if kwargs and "tolerance" in kwargs and is_scalar(kwargs["tolerance"]):
            tol = mp.mpf(kwargs["tolerance"])
        else:
            tol = self.default_tolerance
        with mp.workdps(get_dps_by_tol(tol)):
            for n in range(len(self.__data)):
                if max(mp_abs((mp.eig(self.__data[n].H * self.__data[n]))[0])) <= 1 + tol:
                    self.__passive[n] = True
    
    def passive(self):
        if not self.__passive:
            self.detect_passivity()
        return self.__passive
    
    def is_passive(self, **kwargs):
        if not self.__passive or kwargs:
            self.detect_passivity(**kwargs)
        if all(self.__passive):
            return True
        return False
    
    def enforce_passivity(self, **kwargs):
        if self.__type != "S":
            raise Exception("Must be in S-form for the passivity enforcement to be meaningful!")
        if kwargs and "tolerance" in kwargs:
            tol = kwargs["tolerance"]
        else:
            tol = self.default_tolerance
        with mp.workdps(get_dps_by_tol(tol)):    
            for n in range(len(self.__data)):
                local_max =  mp.sqrt(max(mp.matrix(mp.eig(self.__data[n].H * self.__data[n])[0]).apply(mp.fabs)))
                if local_max > 0:
                    self.__data[n] = self.__data[n] / local_max
    
    def convert_to_Differential(self):
        raise NotImplementedError()
    
        if self.__type != "S":
            raise Exception("Must be in S-form for the conversion to Differential to be meaningful!")
    
    def convert_to_Single_Ended(self):
        raise NotImplementedError()
    
        if self.__type != "S":
            raise Exception("Must be in S-form for the conversion to Single-Ended to be meaningful!")
    
    def set_diff_map(self, *args, **kwargs):
        raise NotImplementedError()
    
    def change_type_to(self, new_type, **kwargs):
        if self.__type == new_type:
            return
        if self.__type != "S":
            self.__return_to_S(**kwargs)
        if new_type == "S":
            return
        elif new_type == "T":
            self.__convert_S_to_T(**kwargs)
        elif new_type == "Z":
            self.__convert_S_to_Z(**kwargs)
        elif new_type == "Y":
            self.__convert_S_to_Y(**kwargs)
        elif new_type == "H":
            self.__convert_S_to_H(**kwargs)
        else:
            raise ValueError("Wrong conversion type!")
        
    def __convert_one_to_another(self, local_function, **kwargs):
        if kwargs and "tolerance" in kwargs:
            tol = kwargs["tolerance"]
        else:
            tol = self.default_tolerance
        with mp.workdps(get_dps_by_tol(tol)):
            for n in range(self.__N_freq):
                self.__data[n] = local_function(self.__data[n])
        
    def __convert_S_to_T(self, **kwargs):
        if self.__type != "S":
            raise RuntimeError("The instance should be of S-type for this to work, we shouldn't be here!")
        self.__convert_one_to_another(S_to_T, **kwargs)
        self.__type = "T"
            
    def __convert_S_to_Z(self, **kwargs):
        if self.__type != "S":
            raise RuntimeError("The instance should be of S-type for this to work, we shouldn't be here!")
        z0 = z0_from_args(self.__N_ports, self.__Z0)
        self.__convert_one_to_another(lambda x: S_to_Z_simple(x, z0), **kwargs)
        self.__type = "Z"
            
    def __convert_S_to_Y(self, **kwargs):
        if self.__type != "S":
            raise RuntimeError("The instance should be of S-type for this to work, we shouldn't be here!")
        y0 = y0_from_args(self.__N_ports, self.__Z0)
        self.__convert_one_to_another(lambda x: S_to_Y_simple(x, y0), **kwargs)
        self.__type = "Y"
            
    def __convert_S_to_H(self, **kwargs):
        if self.__type != "S":
            raise RuntimeError("The instance should be of S-type for this to work, we shouldn't be here!")
        self.__convert_one_to_another(lambda x: S_to_H(x, self.__Z0), **kwargs)
        self.__type = "H"
    
    def __return_to_S(self, *args, **kwargs):
        if args:
            new_Z0 = args[0]
        else:
            new_Z0 = self.__Z0
        if self.__type == "S":
            return
        elif self.__type == "T":
            self.__convert_one_to_another(T_to_S,**kwargs)
        elif self.__type == "Z":
            y0 = y0_from_args(self.__N_ports, new_Z0)
            self.__convert_one_to_another(lambda x: Z_to_S_simple(x, y0),**kwargs)
        elif self.__type == "Y":
            z0 = z0_from_args(self.__N_ports, new_Z0)
            self.__convert_one_to_another(lambda x: Y_to_S_simple(x, z0),**kwargs)
        elif self.__type == "H":
            self.__convert_one_to_another(H_to_S,**kwargs)
        else:
            raise RuntimeError("The instance has un unexpected type, we shouldn't be here!")
        self.__type = "S"
                
    def reorder(self, *args, **kwargs):
        if args and not kwargs:
            local_order = args[0]
        elif kwargs and "new_order" in kwargs and not args:
            local_order = kwargs["new_order"]
        else:
            raise ValueError("Wrong new order vector!")
        if is_vector(local_order, int) and len(args[0]) == self.__N_ports and all(np.sort(local_order) == range(self.__N_ports)):
            tmp = deepcopy(self.__Z0)
            for n in range(self.__N_ports):
                self.__Z0[n] = tmp[local_order[n]]
            for n in range(self.__N_freq):
                tmp = deepcopy(self.__data[n])
                for n_r in range(self.__N_ports):
                    for n_c in range(self.__N_ports):
                        self.__data[n][n_r, n_c] = tmp[local_order[n_r], local_order[n_c]]
            self.detect_throughs()
        else:
            raise ValueError("Wrong new order vector!")
    
    def reorder_by_throughs(self, **kwargs):
        if not self.is_symmetric():
            raise NotImplementedError("Reordering by throughs for non-reciprocal S-parameters is not supported for now!")
        if not self.__throughs:
            self.detect_throughs(**kwargs)
        indexes1 = []
        indexes2 = []
        remaining_indexes = [True]*self.__N_ports
        for through in self.__throughs:
            if remaining_indexes[through[0]]:
                remaining_indexes[through[0]] = False
                indexes1.append(through[0])
            if remaining_indexes[through[1]]:
                remaining_indexes[through[1]] = False
                indexes2.append(through[1])
        self.reorder(indexes1 + indexes2 + list(compress(range(self.__N_ports), remaining_indexes)))
            
            
    def throughs(self):
        if not self.__throughs:
            self.detect_throughs()
        return self.__throughs
            
    def terminate_ports(self, *args, **kwargs):
        raise NotImplementedError()
    
    def change_Z0(self, *args, **kwargs):
        if args:
            local_Z0 = args[0]
        elif kwargs and "new_Z0" in kwargs:
            local_Z0 = kwargs["new_Z0"]
        else:
            raise ValueError("A new Z0 value must be supplied!")
        self.__convert_S_to_Z(**kwargs)
        self.__return_to_S(local_Z0, **kwargs)
        
    def cut_freq(self, **kwargs):
        if kwargs and "min_f" in kwargs:
            min_f = kwargs["min_f"]
        else:
            min_f = self.__f_axis[0,0]
        if kwargs and "max_f" in kwargs:
            max_f = kwargs["max_f"]
        else:
            max_f = self.__f_axis[0,-1]
        local_f = self.__f_axis.tolist()[0]
        indexes_to_remove = []
        for n in range(self.__f_axis.cols):
            if self.__f_axis[0,n] < min_f:
                indexes_to_remove.append(n)
            if self.__f_axis[0,n] > max_f:
                indexes_to_remove.append(n)
        for n in sorted(indexes_to_remove, reverse=True):
            del(local_f[n])
            del(self.__data[n])
            if self.__symmetric:
                del(self.__symmetric[n])
            if self.__passive:
                del(self.__passive[n])
        self.__f_axis = mp.matrix(local_f).T
        self.__N_freq = self.__f_axis.cols
        
        
    
    def interpolate(self, *args, **kwargs):
        raise NotImplementedError()
    
    def TDR_entry(self, *args, **kwargs):
        raise NotImplementedError()
    
    def f_axis(self, *args):
        if args:
            try:
                return self.__f_axis[0, args[0]]
            except:
                raise ValueError("The input value is not a legal index for the f-axis!")
        else:
            return self.__f_axis
    
    def data(self, *args):
        if args:
            aux = []
            try:
                for local_matrix in self.__data:
                    aux.append(local_matrix[args[0], args[1]])
            except:
                raise ValueError("The input values are not legal indexes for the data!")
            return mp.matrix(aux).T
        else:
            return self.__data
    
    def Z0(self):
        return self.__Z0
    
    def plot_entry(self, *args, **kwargs):
        figure_handle = plt.figure()
        axes_handle = figure_handle.add_subplot(111)
        lines_to_plot = []
        if not kwargs:
            scale_converter = dB20
            mag_label = "dB20"
            plot_title = "S-parameter, Magnitude"
        else:
            if not "scale" in kwargs:
                scale_converter = dB20
                mag_label = "dB20"
                plot_title = "S-parameter, Magnitude"
            elif kwargs["scale"] == "dB20":
                scale_converter = dB20
                mag_label = "dB20"
                plot_title = "S-parameter, Magnitude"
            elif kwargs["scale"] == "dB10":
                scale_converter = dB10
                mag_label = "dB10"
                plot_title = "S-parameter, Magnitude"
            elif kwargs["scale"] == "mag" or kwargs["scale"] == "abs":
                scale_converter = mp.fabs
                mag_label = "Abs"
                plot_title = "S-parameter, Magnitude"
            elif kwargs["scale"] == "deg":
                scale_converter = lambda x: deg(mp.phase(x))
                mag_label = "Degrees"
                plot_title = "S-parameter, Phase"
            elif kwargs["scale"] == "rad":
                scale_converter = mp.phase
                mag_label = "Radians"
                plot_title = "S-parameter, Phase"
            else:
                raise ValueError("Wrong scale argument!")
            
        if self.__f_axis[0,self.__N_freq-1] > 1e9:
            f_scale = 1e9
            f_label = "[GHz]"
        elif self.__f_axis[0,self.__N_freq-1] > 1e6:
            f_scale = 1e6
            f_label = "[MHz]"
        elif self.__f_axis[0,self.__N_freq-1] > 1e3:
            f_scale = 1e3
            f_label = "[kHz]"
        else:
            f_scale = 1
            f_label = "[Hz]"
        local_indexes = []
        if args:
            try:
                for indexes in args:
                    local_line = []
                    for n in range(self.__f_axis.cols):
                        local_line.append(scale_converter(self.__data[n][indexes[0], indexes[1]]))
                    lines_to_plot.append(local_line)
                    local_indexes.append(indexes)
            except:
                raise ValueError("The input arguments are not legal indexes for this S-parameter instance!")
        else:
            local_line = []
            for n in range(self.__f_axis.cols):
                local_line.append(scale_converter(self.__data[n][0,0]))
            lines_to_plot.append(local_line)
            local_indexes = [[0,0]]
        line_handles = []
        for n in range(len(lines_to_plot)):
            line_handles.append(axes_handle.plot(self.__f_axis/f_scale, lines_to_plot[n], label="[{},{}]".format(local_indexes[n][0], local_indexes[n][1])))
        axes_handle.legend()
        plt.xlabel(f_label)
        plt.ylabel(mag_label)
        plt.title(plot_title)
        plt.xlim((float(self.__f_axis[0,0]/f_scale), float(self.__f_axis[0, self.__N_freq - 1]/f_scale)))
        axes_handle.grid(True, linewidth=1)
        plt.draw()
    
    def plot_slice(self, *args, **kwargs):
        raise NotImplementedError()
    
    def plot_entry_TDR(self, *args, **kwargs):
        raise NotImplementedError()
    
    def __read_S_parameter_from_file(self, file_name):
        i = mp.mpc("0", "1")
        try:
            _, local_extension = os.path.splitext(file_name)
        except:
            raise ValueError("The given file name couldn't be splitted into the path and the extension!")
        try:
            self.__N_ports = int(local_extension[2:-1])
            if is_even(self.__N_ports):
                lines_per_matrix = floor((self.__N_ports**2)*2/8)
            else:
                lines_per_matrix = floor((self.__N_ports**2)*2/8) + 1
        except:
            raise ValueError("The given file extension is not in correct sNp format!")
        format_flag = True
        file_start_flag = True
        lines_counter = 0
        self.__data = []
        local_freq = []
        local_matrix = mp.matrix(self.__N_ports, self.__N_ports)
        n_row = 0
        n_col = 0
        with open(file_name, "r") as file_handle:
            for line in file_handle:
                if line:
                    if line[0] == "!" and file_start_flag: 
                        continue
                    if line[0] == "!" and not file_start_flag:
                        break
                    local_list = list(filter(None, line.rstrip().split(" ")))
                    if local_list[0] == "#" and format_flag:
                        if len(local_list) < 6:
                            raise RuntimeError("The format string length is wrong!")
                        if local_list[1] == "hz" or local_list[1] == "Hz" or local_list[1] == "HZ":
                            freq_multiplier = mp.mpf("1")
                        elif local_list[1] == "khz" or local_list[1] == "kHz" or local_list[1] == "KHz" or local_list[1] == "KHZ":
                            freq_multiplier = mp.mpf("1e3")
                        elif local_list[1] == "mhz" or local_list[1] == "mHz" or local_list[1] == "MHz" or local_list[1] == "MHZ":
                            freq_multiplier = mp.mpf("1e6")
                        elif local_list[1] == "ghz" or local_list[1] == "gHz" or local_list[1] == "GHz" or local_list[1] == "GHZ":
                            freq_multiplier = mp.mpf("1e9")
                        else:
                            raise RuntimeError("Frequency format is not recognized!")
                        if local_list[2] == "s" or local_list[2] == "S":
                            self.__type = "S"
                        elif local_list[2] == "t" or local_list[2] == "T":
                            self.__type = "T"
                        elif local_list[2] == "z" or local_list[2] == "Z":
                            self.__type = "Z"
                        elif local_list[2] == "y" or local_list[2] == "Y":
                            self.__type = "Y"
                        elif local_list[2] == "h" or local_list[2] == "H":
                            self.__type = "H"
                        else:
                            raise RuntimeError("Matrix type is not recognized!")
                        if local_list[3] == "ma" or local_list[3] == "mA" or local_list[3] == "Ma" or local_list[3] == "MA":
                            local_converter = lambda x,y : mp.mpf(x) * mp.exp(i*rad(mp.mpf(y)))
                        elif local_list[3] == "ri" or local_list[3] == "Ri" or local_list[3] == "rI" or local_list[3] == "RI":
                            local_converter = lambda x,y : mp.mpf(x) + i*mp.mpf(y)
                        elif local_list[3] == "db" or local_list[3] == "Db" or local_list[3] == "dB" or local_list[3] == "DB":
                            local_converter = lambda x,y : (mp.mpf("10")**(mp.mpf(x)/mp.mpf("20"))) * mp.exp(i*rad(mp.mpf(y)))
                        else:
                            raise RuntimeError("Matrix complex numbers convention is not recognized!")
                        self.__Z0 = mp.mpf(local_list[5]) * mp.ones(1,self.__N_ports)
                    elif float(local_list[0]):
                        format_flag = False
                        file_start_flag = False
                        if len(local_list) < 3:
                            raise ValueError("Wrong length of a data line for an sNp file!")
                        if not lines_counter:
                            tmp = mp.mpf(local_list[0]) * freq_multiplier
                            local_freq.append(tmp)
                            if tmp < 1e3:
                                units_str = "[Hz]"
                            elif tmp < 1e6:
                                tmp /= 1e3
                                units_str = "[kHz]"
                            elif tmp < 1e9:
                                tmp /= 1e6
                                units_str = "[MHz]"
                            else:
                                tmp /= 1e9
                                units_str = "[GHz]"
                            print("   Reading " + str(tmp) + units_str + "\n")
                            local_list.pop(0)
                        while local_list:
                            local_matrix[n_row, n_col] = local_converter(local_list.pop(0), local_list.pop(0))
                            n_col += 1
                            if n_col == self.__N_ports:
                                n_col = 0
                                n_row += 1
                        lines_counter += 1
                        if lines_counter == lines_per_matrix:
                            n_row = 0
                            n_col = 0
                            lines_counter = 0
                            self.__data.append(local_matrix)
                            local_matrix = mp.matrix(self.__N_ports, self.__N_ports)
        self.__f_axis = convert_to_mp_matrix(local_freq).T
        self.__N_freq = self.__f_axis.cols
        self.__mode = "SE"
    
    def write(self, **kwargs):
        raise NotImplementedError()
    
    
    
        
    
            