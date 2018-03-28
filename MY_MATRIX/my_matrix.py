
from mpmath import *
import numpy as np
from copy import copy as local_copy
from aux_functions import is_even, NUMERICAL_TYPES, mpc_to_str, mpf_to_str, mpfi_to_str, get_index
from operator import floordiv


class submatrix:
    pass

class _11(submatrix):
    pass

class _12(submatrix):
    pass

class _21(submatrix):
    pass

class _22(submatrix):
    pass

class my_matrix:
    
    __use_tol = True
    
    __tol_overshoot = 5
    
    __print_digits = 5
    
    __max_approximation_steps = 5000
    
    def __init__(self, *args):
        self.__N_rows = 0
        self.__N_cols = 0
        if args:
            if isinstance(args[0], mp.matrix):
                self.__N_rows = args[0].rows
                self.__N_cols = args[0].cols
                for r in self.rows():
                    for c in self.cols():
                        if args[0][r,c] != mp.mpf("0"):
                            self.__data[(r,c)] = args[0][r,c]
            elif isinstance(args[0], list):
                try:
                    self.__N_rows = len(args[0])
                    self.__N_cols = len(args[0][0])
                    self.__data = {}
                    for r in self.rows():
                        for c in self.cols():
                            if args[0][r][c] != mp.mpf("0"):
                                self.__data[(r,c)] = args[0][r][c]
                except:
                    raise ValueError("Wrong list format for my_matrix constructor!")                                
            elif len(args) == 2 and isinstance(args[0], int) and isinstance(args[1], int) and args[0] >=0 and args[1] >= 0:
                self.__N_rows = args[0]
                self.__N_cols = args[1]
                self.__data = {}
            elif isinstance(args[0], mp.matrix):
                self.__N_rows = args[0].rows
                self.__N_cols = args[0].cols
                for r in range(args[0].rows):
                    for c in range(args[0].cols):
                        if args[0][r,c] != mp.mpf("0"):
                            self.__data[(r,c)] = args[0][r,c]
            else:
                raise ValueError("Wrong inputs to my_matrix.__init__!")
        else:
            self.__data = {}
        
    def __str__(self):
        if self.is_Real():
            num_to_str = lambda x: mpf_to_str(x, self.__print_digits)
            n_symbols = self.__print_digits
        elif self.is_Imag():
            num_to_str = lambda x: mpfi_to_str(x, self.__print_digits+1)
            n_symbols = self.__print_digits+1
        else:
            num_to_str = lambda x: mpc_to_str(x, 2*self.__print_digits+5)
            n_symbols = 2*self.__print_digits+5
        sep_str = "-" * (self.__N_cols*(n_symbols+3) + 1)
        local_str = ""
        for r in self.rows():
            local_str = local_str + "\n" + sep_str + "\n|"
            for c in self.cols():
                local_str = local_str + " " + num_to_str(self[r,c]) + " |"
        local_str = local_str + "\n" + sep_str + "\n"
        return local_str
            
    def N_rows(self):
        return self.__N_rows
    
    def rows(self):
        return range(self.__N_rows)
    
    def N_cols(self):
        return self.__N_cols
    
    def cols(self):
        return range(self.__N_cols)
    
    def size(self):
        return self.__N_rows, self.__N_cols
    
    def data(self):
        return local_copy(self.__data)
    
    
    
    def H(self):
        M = my_matrix()
        for keys in self.__data:
            r = keys[0]
            c = keys[1]
            M.__data[(c,r)] = mp.conj(self.__data[keys])
        M.__N_rows = self.__N_cols
        M.__N_cols = self.__N_rows
        return M
    
    def T(self):
        M = my_matrix()
        for keys in self.__data:
            r = keys[0]
            c = keys[1]
            M.__data[(c,r)] = self.__data[keys]
        M.__N_rows = self.__N_cols
        M.__N_cols = self.__N_rows
        return M
        
    def is_Square(self):
        return self.__N_rows == self.__N_cols
    
    def is_Diagonal(self, *args):
        if self.__N_rows != self.__N_cols:
            return False
        tol = self.get_tol(*args)
        for r in self.rows():
            for c in self.cols():
                if r != c and mp.fabs(self[r,c]) > tol:
                    return False
        return True
    
    def is_Symmetric(self, *args):
        if self.__N_rows != self.__N_cols:
            return False
        tol = self.get_tol(*args)
        for r in self.rows():
            for c in range(r, self.__N_cols):
                if mp.fabs(self[r,c] - self[c,r]) > tol:
                    return False
        return True
                    
    def is_Hermitian(self, *args):
        tol = self.get_tol(*args)
        if self.__N_rows != self.__N_cols:
            return False
        for r in self.rows():
            for c in range(r, self.__N_cols):
                if mp.fabs(self[r,c] - mp.conj(self[c,r])) > tol:
                    return False
        return True
    
    def is_Unitary(self, *args):
        tol = self.get_tol(*args)
        if self.is_Diagonal():
            for r in self.rows():
                if mp.fabs(mp.conj(self[r,r])*self[r,r] - mp.mpf("1")) > tol:
                    return False
            return True
        else:
            if (self.H()*self - self.to_I()).max_abs() > tol:
                return False
            return True
        
    def is_Real(self, *args):
        tol = self.get_tol(*args)
        for keys in self.__data:
            if mp.fabs(mp.im(self.__data[keys])) > tol:
                return False
        return True
    
    def is_Imag(self, *args):
        tol = self.get_tol(*args)
        for keys in self.__data:
            if mp.fabs(mp.re(self.__data[keys])) > tol:
                return False
        return True
        
    def is_Positive(self, *args):
        tol = self.get_tol(*args)
        L, _ = (self.H()*self).eig()
        for keys in L.__data:
            if L.__data[keys] < tol:
                return False
        return True
    
    def is_semi_Positive(self, *args):
        tol = self.get_tol(*args)
        L, _ = (self.H()*self).eig()
        for keys in L.__data:
            if L.__data[keys] < -tol:
                return False
        return True
    
    def is_Negative(self, *args):
        tol = self.get_tol(*args)
        L, _ = (self.H()*self).eig()
        for keys in L.__data:
            if L.__data[keys] > -tol:
                return False
        return True
    
    def is_semi_Negative(self, *args):
        tol = self.get_tol(*args)
        L, _ = (self.H()*self).eig()
        for keys in L.__data:
            if L.__data[keys] > tol:
                return False
        return True
    
    def has_orthogonal_columns(self, *args):
        tol = self.get_tol(*args)
        if self.is_Diagonal():
            return True
        else:
            for c1 in range(self.__N_cols-1):
                for c2 in range(c1+1, self.__N_cols):
                    if self["c",c1].H() * self["c",c2] >= tol:
                        return False
            return True
    
    def has_orthogonal_rows(self, *args):
        tol = self.get_tol(*args)
        if self.is_Diagonal():
            return True
        else:
            for r1 in range(self.__N_rows-1):
                for r2 in range(r1+1, self.__N_rows):
                    if self["r", r1] * self["r", r2].H() >= tol:
                        return False
            return True
        
    def rank_Cols(self, *args):
        tol = self.get_tol(*args)
        M = local_copy(self)
        M.renormalize_cols(*args)
        cols_rank = self.__N_cols
        for c in M.cols():
            stop_flag = False
            for r in M.rows():
                if mp.fabs(M[r,c]) > tol:
                    stop_flag = True
                    break
            if not stop_flag:
                cols_rank -= 1
        return cols_rank 
    
    def rank_Rows(self, *args):
        return self.T.rank_Cols(*args)
    
    def rank(self, *args):
        if self.is_Square():
            return self.rank_Cols(*args) == self.__N_cols
        else:
            raise RuntimeError("The matrix is not square, use rank_Rows or rank_Cols instead!")
        
    def diag(self):
        if self.__N_rows == 1:
            M = my_matrix(self.__N_cols, self.__N_cols)
            for c in self.cols():
                if (0,c) in self.__data:
                    M.__data[(c,c)] = self.__data[(0,c)]
            return M
        elif self.__N_cols == 1:
            M = my_matrix(self.__N_rows, self.__N_rows)
            for r in self.rows():
                if (r,0) in self.__data:
                    M.__data[(r,r)] = self.__data[(r,0)]
            return M
        elif self.is_Square():
            M = my_matrix(self.__N_rows,1)
            for r in self.rows():
                if (r,r) in self.__data:
                    M.__data[(r,0)] = self.__data[(r,r)]
            return M            
        else:
            raise RuntimeError("The diag method is applicable only to vector matrices or square matrices!")
    
    def inv(self, *args):
        tol = self.get_tol(*args)
        if self.is_Diagonal():
            M = local_copy(self)
            for r in self.rows():
                if mp.fabs(M[r,r]) < tol:
                    raise RuntimeError("The matrix is singular to the working precision!")
                else:
                    M[r,r] = mp.mpf("1") / M[r,r]
            return M 
        else:
            _, Mi = self.codiagonalize(self.to_I(), *args)
            return Mi
    
    def eig(self, *args):
        tol = self.get_tol(*args)
        Q,R = self.QR_Factorization()
#         prev_guess = R.diag()
        prev_guess = (R*Q).diag()
        for counter in range(self.__max_approximation_steps):
            M = R*Q
#             print(M)
            Q,R = M.QR_Factorization()
            new_guess = (R*Q).diag()
#             print("prev guess: " + prev_guess.__str__())
#             print(str(counter) + ": " + str((prev_guess - new_guess).max_abs()))
            if (prev_guess - new_guess).max_abs() < tol:
                return new_guess.diag()
            prev_guess = local_copy(new_guess)
        raise RuntimeError("QR algorithm didn't converge to eigenvalues in __max_approximation_steps!") 
    
       
    def __add__(self, other):
        if isinstance(other, my_matrix):
            if self.__N_rows == other.__N_rows and self.__N_cols == other.__N_cols:
                M = local_copy(self)
                for keys in other.__data:
                    try:
                        M.__data[keys] += other.__data[keys]
                    except KeyError:
                        M.__data[keys] = other.__data[keys]
                return M
            else:
                raise ValueError("Matrix dimensions must agree!")
        elif isinstance(other, NUMERICAL_TYPES.scalar_types):
            M = local_copy(self)
            for keys in self.__data:
                M.__data[keys] += other 
            return M
        else:
            raise ValueError("Wrong input to __add__!")
        
    def __radd__(self, other):
        return self + other
    
    def __sub__(self, other):
        if isinstance(other, my_matrix):
            if self.__N_rows == other.__N_rows and self.__N_cols == other.__N_cols:
                M = local_copy(self)
                for keys in other.__data:
                    try:
                        M.__data[keys] -= other.__data[keys]
                    except KeyError:
                        M.__data[keys] = -other.__data[keys]
                return M
            else:
                raise ValueError("Matrix dimensions must agree!")
        elif isinstance(other, NUMERICAL_TYPES.scalar_types):
            M = local_copy(self)
            for keys in self.__data:
                M.__data[keys] -= other 
            return M
        else:
            raise ValueError("Wrong input to __sub__!")
        
    def __rsub__(self, other):
        return -self + other
    
    def __neg__(self):
        M = local_copy(self)
        for keys in M.__data:
            M.__data[keys] = -M.__data[keys]
        return M
    
    def __mul__(self, other):
        if isinstance(other, my_matrix):
            if self.__N_cols == other.__N_rows:
                if self.__N_rows == 1 and other.__N_cols == 1:
                    aux = mp.mpf("0")
                    for keys in self.__data:
                        try:
                            aux += self.__data[keys] * other.__data[(keys[1], keys[0])]
                        except KeyError:
                            pass 
                    return aux
                M = my_matrix(self.__N_rows, other.__N_cols)
                for r in self.rows():
                    for c in range(other.__N_cols):
                        aux = mp.mpf("0")
                        for k in self.cols():
                            try:
                                aux += self.__data[(r,k)] * other.__data[(k,c)]
                            except KeyError:
                                pass
                        if aux != mp.mpf("0"):
                            M.__data[(r,c)] = aux
                return M
            else:
                raise ValueError("Matrix dimensions must agree!")
        elif isinstance(other, NUMERICAL_TYPES.scalar_types):
            M = local_copy(self)
            for keys in M.__data:
                M.__data[keys] *= other 
            return M
        else:
            raise ValueError("Wrong input to __add__!")
        
    def __rmul__(self, other):
        if isinstance(other, NUMERICAL_TYPES.scalar_types):
            return self * other
        else:
            raise NotImplementedError
        
    def __truediv__(self, other):
        if isinstance(other, NUMERICAL_TYPES.scalar_types):
            M = local_copy(self)
            if mp.fabs(other) < self.get_tol():
                raise ValueError("Attempting a division by zero!")
            for keys in M.__data:
                M.__data[keys] /= other 
            return M
        else:
            raise NotImplementedError
    
    def __gt__(self, other):
        pass
    
    def __lt__(self, other):
        pass
    
    def __bool__(self):
        if self.__N_rows == 0 and self.__N_cols == 0:
            return False
        return True
    
    def __len__(self):
        return self.__N_rows * self.__N_cols
    
    def __getitem__(self, keys):
        rows_list, cols_list = self.__keys_to_indexes(keys)
        if len(rows_list) == 1 and len(cols_list) == 1:
            if (rows_list[0], cols_list[0]) in self.__data:
                return self.__data[(rows_list[0], cols_list[0])]
            else:
                return mp.mpf("0")
            return self[rows_list[0], cols_list[0]]
        M = my_matrix(len(rows_list), len(cols_list))
        for r in range(len(rows_list)):
            for c in range(len(cols_list)):
                try:
                    M.__data[(r, c)] = self.__data[(rows_list[r], cols_list[c])]
                except KeyError:
                    pass
        return M
    
    def __setitem__(self, keys, value):
        rows_list, cols_list = self.__keys_to_indexes(keys)
        if isinstance(value, NUMERICAL_TYPES.scalar_types):
            if value != mp.mpf("0"):
                for r in rows_list:
                    for c in cols_list:
                        self.__data[(r,c)] = value
            else:
                keys_to_del = []
                for r in rows_list:
                    for c in cols_list:
                        if (r,c) in self.__data:
                            keys_to_del.append((r,c))
                for local_keys in keys_to_del:
                    del self.__data[local_keys]
        elif isinstance(value, my_matrix):
            if len(rows_list) == value.__N_rows and len(cols_list) == value.__N_cols:
                keys_to_del = []
                for r in value.rows():
                    for c in value.cols():
                        if (r,c) in value.__data:
                            self.__data[(rows_list[r], cols_list[c])] = value.__data[(r,c)]
                        else:
                            keys_to_del.append((rows_list[r], cols_list[c]))
                for local_keys in keys_to_del:
                    try:
                        del self.__data[local_keys]
                    except KeyError:
                            pass
            else:
                raise ValueError("Wrong indexes input!")
    
    def __delitem__(self, keys):
        raise NotImplementedError("A matrix entry can't be deleted, replace it by 0 instead!")
    
    def __contains__(self, value):
        raise NotImplementedError
    
    def __copy__(self):
        M = my_matrix()
        M.__N_rows = self.__N_rows
        M.__N_cols = self.__N_cols
        M.__data = self.__data.copy()
        return M
    
    def __pow__(self, degree):
        if isinstance(degree, int) and degree >= 0:
            if self.is_Square():
                if degree == 0:
                    return self.I(self.__N_rows, self.__N_cols)
                if degree == 1:
                    return local_copy(self)
                else:
                    M = local_copy(self)
                    for _ in range(1,degree):
                        M = M * self
                    return M
            else:
                raise ValueError("Power is defiend only for square matrices!")
        else:
            raise NotImplementedError
    
    def __eq__(self, other):
        tol = self.get_tol()
        if isinstance(other, NUMERICAL_TYPES.scalar_types):
            for r in self.rows():
                for c in self.cols():
                    if mp.fabs(self[r,c] - other) > tol:
                        return False
            return True
        if not isinstance(other, my_matrix):
            return False
        if self.__N_cols != other.__N_cols or self.__N_rows != other.__N_rows:
            return False
        for r in self.rows():
            for c in self.cols():
                if mp.fabs(self[r,c] - other[r,c]) > tol:
                    return False
        return True
    
    def __req__(self, other):
        return self == other
    
    def sqrt(self):
        return local_copy(self).apply(mp.sqrt)
    
    def abs(self):
        return local_copy(self).apply(mp.fabs)
    
    def conj(self):
        return local_copy(self).apply(mp.conj)
    
    def angle(self):
        return local_copy(self).apply(mp.phase)
    
    def phase(self):
        return local_copy(self).apply(mp.phase)
    
    def real(self):
        return local_copy(self).apply(mp.re)
    
    def imag(self):
        return local_copy(self).apply(mp.im)
    
    def det(self):
        if not self.is_Square():
            raise RuntimeError("Determinant is not defined for not-square matrices!")
        aux = mp.mpf("0")
        if self.__N_rows == 2 and self.__N_cols == 2:
            return self[0,0]*self[1,1] - self[1,0]*self[0,1]
        for c in self.cols():
            local_c = list(range(self.__N_cols))
            local_c.pop(c)
            if is_even(c):
                aux += self[0,c]*self[1:,local_c].det()
            else:
                aux -= self[0,c]*self[1:,local_c].det()
        return aux                      
    
    def max_abs(self):
        max_val = mp.mpf("0")
        for keys in self.__data:
            if mp.fabs(self.__data[keys]) > max_val:
                max_val = mp.fabs(self.__data[keys])
        return max_val
    
    def min_abs(self):
        min_val = mp.inf
        for r in self.rows():
            for c in self.cols():
                if mp.fabs(self[r,c]) < min_val:
                    min_val = mp.fabs(self[r,c])
        return min_val

    def to_ndarray(self):
        pass
    
    def to_np_matrix(self):
        return np.matrix(self.apply(float).to_list())
    
    def to_mp_matrix(self):
        return mp.matrix(self.to_list())
    
    def to_list(self, **kwargs):
        aux = []
        if kwargs and "cols_form" in kwargs and kwargs["cols_form"]:
            for c in self.cols():
                local_aux = []
                for r in self.rows():
                    local_aux.append(self[r,c])
                aux.append(local_aux)
        else:
            for r in self.rows():
                local_aux = []
                for c in self.cols():
                    local_aux.append(self[r,c])
                aux.append(local_aux)
        return aux
    
    def enforce_Symmetry(self, *args):
        tol = self.get_tol(*args)
        for r in self.rows():
            for c in range(r, self.__N_cols):
                aux1 = self[r,c]
                aux2 = self[c,r]
                val = (aux1 + aux2) / mp.mpf("2")
                if mp.fabs(val) < tol:
                    if (r,c) in self.__data:
                        del self.__data[r,c]
                    if (c,r) in self.__data:
                        del self.__data[c,r]    
                else:
                    self[r,c] = val
                    self[c,r] = val
        return local_copy(self)
    
    def enforce_Hermiticity(self, *args):
        tol = self.get_tol(*args)
        for r in self.rows():
            for c in range(r, self.__N_cols):
                aux1 = self[r,c]
                aux2 = self[c,r]
                val = (aux1 + mp.conj(aux2)) / mp.mpf("2")
                if mp.fabs(val) < tol:
                    if (r,c) in self.__data:
                        del self.__data[r,c]
                    if (c,r) in self.__data:
                        del self.__data[c,r]    
                else:
                    self[r,c] = val
                    self[c,r] = mp.conj(val)
        return local_copy(self)
    
    def enforce_Unitarity(self, *args):
        return self.renormalize_cols()
        if self.has_zero_cols(*args):
            raise RuntimeError("Enforcing unitarity failed, the matrix columns don't constitute a complete basis!")
    
    def get_Hessenberg_form(self):
        i = mp.mpc("0","1")
        if not self.is_Square():
            raise RuntimeError("Can't force a non-square matrix into a Hessenberg form!")
        if self.__N_cols == 2:
            raise RuntimeError("Can't force a 2x2 matrix into a Hessenberg form!")
        M = local_copy(self)
        I = M.to_I()
        for c in range(M.__N_cols-2):
            v = M["c",c]
            v_tar = my_matrix(M.__N_rows,1)
            for loc_c in range(c+1):
                try:
                    v_tar.__data[(loc_c,0)] = v.__data[(loc_c,0)]
                except KeyError:
                    pass
            aux = mp.mpf("0")
            for loc_c in range(c+1,M.__N_rows):
                try:
                    aux += mp.conj(v.__data[(loc_c,0)])*v.__data[(loc_c,0)]
                except KeyError:
                    pass
            v_tar.__data[(c+1,0)] = mp.sqrt(v["r",c+1:].H()*v["r",c+1:])*mp.exp(i*mp.phase(v[c+1,0]))
            u = v - v_tar
            aux = mp.mpf("0")
            for keys in u.__data:
                aux += mp.conj(u.__data[keys]) * u.__data[keys]
            aux = mp.sqrt(aux)
            for keys in u.__data:
                u.__data[keys] = u.__data[keys] / aux
            P = I.to_I()
            for loc_c in range(c+1,P.__N_rows):
                for loc_r in range(loc_c, P.__N_rows):
                    if loc_r == loc_c:
                        P.__data[(loc_r,loc_c)] = mp.mpf("1") - mp.mpf("2") * u.__data[(loc_r,0)] * mp.conj(u.__data[(loc_c,0)])
                    else:
                        P.__data[(loc_r,loc_c)] = - mp.mpf("2") * u.__data[(loc_r,0)] * mp.conj(u.__data[(loc_c,0)])
                        P.__data[(loc_c,loc_r)] = - mp.mpf("2") * u.__data[(loc_c,0)] * mp.conj(u.__data[(loc_r,0)])
            M = M - mp.mpf("2")*u*(u.H()*M)
            M = M - mp.mpf("2")*(M*u)*u.H()
        return M
            
    def zero_to_tol(self, *args):
        tol = self.get_tol(*args)
        keys_to_del = []
        for keys in self.__data:
            if mp.fabs(self.__data[keys]) < tol: 
                    keys_to_del.append(keys)
        for keys in keys_to_del:
            self.__data.pop(keys)
        return local_copy(self)
    
    def round_to_tol(self, *args):
        tol = self.get_tol(*args)
        keys_to_del = []
        for keys in self.__data:
            if mp.fabs(self.__data[keys]) < tol: 
                    keys_to_del.append(keys)
            else:
                aux_re = mp.fabs(mp.re(self.__data[keys]))
                aux_im = mp.fabs(mp.im(self.__data[keys]))
                if mp.fabs(aux_re - mp.nint(aux_re)) < tol:
                    aux_re = mp.nint(aux_re) * mp.sign(mp.re(self.__data[keys]))
                if mp.fabs(aux_im - mp.nint(aux_im)) < tol:
                    aux_im = mp.nint(aux_im) * mp.sign(mp.im(self.__data[keys]))
                if aux_im == mp.mpf("0"):
                    self.__data[keys] = mp.mpf(aux_re)
                else:
                    self.__data[keys] = mp.mpc(aux_re, aux_im)
        for keys in keys_to_del:
            self.__data.pop(keys)
        return local_copy(self)
    
    def apply(self, function_handle):
        try:
            M = local_copy(self)
            for r in range(M.__N_rows):
                for c in range(M.__N_rows):
                    M[r,c] = function_handle(M[r,c])
            return M
        except:
            raise RuntimeError("Application of the given function to matrix elements failed!")
    
    def renormalize_rows(self, *args):
        tol = self.get_tol(*args)
        for r1 in self.rows():
            aux = self["r",r1] * self["r",r1].H()
            if mp.fabs(aux) <= tol:
                continue
            else:
                self["r",r1] = self["r",r1] / mp.sqrt(aux) 
            for r2 in range(r1+1, self.__N_rows):
                self["r",r2] = self["r",r2] - (self["r",r1].conj()*self["r",r2].T())*self["r",r1]
        return local_copy(self)
    
    def renormalize_cols(self, *args):
        tol = self.get_tol(*args)
        for c1 in self.cols():
            aux = self["c",c1].H() * self["c",c1]
            if mp.fabs(aux) <= tol:
                continue
            else:
                self["c",c1] = self["c",c1] / mp.sqrt(aux) 
            for c2 in range(c1+1, self.__N_cols):
                self["c",c2] = self["c",c2] - (self["c",c1].H()*self["c",c2])*self["c",c1]
        return local_copy(self)
    
    def codiagonalize(self, other, *args):
        if not isinstance(other, my_matrix) or not self.is_Square() or not other.is_Square() or self.__N_rows != other.__N_rows:
            raise ValueError("Wrong input to codiagonalize!")
        tol = self.get_tol(*args)
        Ms = local_copy(self)
        Mo = local_copy(other)
        break_flag = False
        for c in range(Ms.__N_cols):
            for r in range(c, Ms.__N_rows):
                if mp.fabs(Ms[r,c]) > tol:
                    break_flag = True
                    break
            if break_flag:
                r_init = r
                aux = Ms[r_init, c]
                Ms["r",r_init] = Ms["r", r_init] / aux
                Mo["r",r_init] = Mo["r", r_init] / aux
                for r in range(r_init+1, self.__N_rows):
                    aux = Ms[r,c]
                    Ms["r",r] = Ms["r",r] - Ms["r", r_init]*aux
                    Mo["r",r] = Mo["r",r] - Mo["r", r_init]*aux
        break_flag = False
        for c in range(Ms.__N_cols-1, -1, -1):
            for r in range(c, -1, -1):
                if mp.fabs(Ms[r,c]) > tol:
                    break_flag = True
                    break
            if break_flag:
                r_init = r
                aux = Ms[r_init, c]
                Ms["r",r_init] = Ms["r", r_init] / aux
                Mo["r",r_init] = Mo["r", r_init] / aux
                for r in range(r_init-1, -1, -1):
                    aux = Ms[r,c]
                    Ms["r",r] = Ms["r",r] - Ms["r", r_init]*aux
                    Mo["r",r] = Mo["r",r] - Mo["r", r_init]*aux
        return Ms, Mo
    
    def AT_Factorization(self):
        if self.is_Symmetric():
            pass
        else:
            raise RuntimeError("Autonne-Takagi factorization is defined only for symmetric matrices!")
        
    def QR_Factorization(self):
        Q = local_copy(self).renormalize_cols()
        R = Q.H() * self
        return Q, R
    
    
        
    def __get_local_row_index(self, ind):
        return get_index(ind, self.__N_rows)
    
    def __get_local_col_index(self, ind):
        return get_index(ind, self.__N_cols)
        
    def __get_local_index(self, *args):
        if len(args) == 2:
            return self.__get_local_row_index(args[0]), self.__get_local_col_index(args[1])
        else:
            raise ValueError("Wrong amount of inputs to __get_local_index!")
        
    def __keys_to_indexes(self, keys):
        if keys == _11:
            if self.is_Square and is_even(self.__N_cols):
                N_M = int(self.__N_cols / 2)
                rows_list = list(range(0, N_M))
                cols_list = list(range(0, N_M))
            else:
                raise RuntimeError("Submatrix keying works only for square matrices with even nomber of rows!")
        elif keys == _12:
            if self.is_Square and is_even(self.__N_cols):
                N_M = int(self.__N_cols / 2)
                rows_list = list(range(0, N_M))
                cols_list = list(range(N_M, self.__N_cols))
            else:
                raise RuntimeError("Submatrix keying works only for square matrices with even nomber of rows!")
        elif keys == _21:
            if self.is_Square and is_even(self.__N_cols):
                N_M = int(self.__N_cols / 2)
                rows_list = list(range(N_M, self.__N_cols))
                cols_list = list(range(0, N_M))
            else:
                raise RuntimeError("Submatrix keying works only for square matrices with even nomber of rows!")
        elif keys == _22:
            if self.is_Square and is_even(self.__N_cols):
                N_M = int(self.__N_cols / 2)
                rows_list = list(range(N_M, self.__N_cols))
                cols_list = list(range(N_M, self.__N_cols))
            else:
                raise RuntimeError("Submatrix keying works only for square matrices with even nomber of rows!")
        elif isinstance(keys, int):
            if keys >= self.__N_rows * self.__N_cols:
                raise StopIteration
            rows_list = [keys % self.__N_rows]
            cols_list = [floordiv(keys, self.__N_rows) % self.__N_cols]
        elif len(keys) == 2:
            if keys[0] == "r" and isinstance(keys[1], (int, range, list, slice)):
                rows_list = self.__get_local_row_index(keys[1])
                cols_list = list(range(self.__N_cols))
            elif keys[0] == "c" and isinstance(keys[1], (int, range, list, slice)):
                rows_list = list(range(self.__N_rows))
                cols_list = self.__get_local_col_index(keys[1])
            else:
                if isinstance(keys[0], (int, range, list, slice)):
                    rows_list = self.__get_local_row_index(keys[0])
                else:
                    raise ValueError("Wrong matrix rows indexes!")
                if isinstance(keys[1], (int, range, list, slice)):
                    cols_list = self.__get_local_col_index(keys[1])
                else:
                    raise ValueError("Wrong matrix columns indexes!")
        else:
            raise ValueError("Wrong matrix indexes!")
        return rows_list, cols_list
            
    @staticmethod
    def get_tol(*args):
        current_dps = mp.dps
        if current_dps <= 5:
            local_tol = mp.mpf("1e-5")
        elif my_matrix.__use_tol:
            local_tol = mp.mpf("1e-" + str(current_dps - 5))
        else:
            local_tol = mp.mpf("0")
        if args:
            if mp.mpf(args[0]) > local_tol:
                    return mp.mpf(args[0])
            else:
                raise RuntimeError("The requested tolerance is finer than the current MPMATH precision!")
        else:
            return local_tol

    @staticmethod    
    def zeros(*args):
        if args:
            try:
                O = my_matrix()
                O.__N_rows = args[0]
                O.__N_cols = args[1]
                return O
            except:
                raise ValueError("Wrong matrix dimensions!")
        else:
            raise ValueError("Wrong matrix dimensions!")
    
    def to_zeros(self):
        return my_matrix(self.__N_rows, self.__N_cols)  
    
    @staticmethod
    def ones(*args):
        if args:
            try:
                O = my_matrix()
                for r in range(args[0]):
                    for c in range(args[1]):
                        O.__data[(r,c)] = mp.mpf("1")
                O.__N_rows = args[0]
                O.__N_cols = args[1]
                return O
            except:
                raise ValueError("Wrong matrix dimensions!")
        else:
            raise ValueError("Wrong matrix dimensions!")
        
    def to_ones(self):
        O = my_matrix(self.__N_rows, self.__N_cols)
        for r in self.rows():
            for c in self.cols():
                O.__data[(r,c)] = mp.mpf("1")
        return O
        
    @staticmethod
    def I(*args):
        if args:
            try:
                if len(args) == 1 or args[0] == args[1]:
                    O = my_matrix()
                    for r in range(args[0]):
                        O.__data[(r,r)] = mp.mpf("1")
                    O.__N_rows = args[0]
                    O.__N_cols = args[0]
                    return O
                else:
                    raise ValueError("Wrong matrix dimensions!")
            except:
                raise ValueError("Wrong matrix dimensions!")
        else:
            raise ValueError("Wrong matrix dimensions!")
    
    def to_I(self):
        if not self.is_Square():
            raise RuntimeError("Can't make a unit matrix out of a not square one!")
        O = my_matrix(self.__N_rows, self.__N_cols)
        for r in self.rows():
            O.__data[(r,r)] = mp.mpf("1")
        return O
        
    @staticmethod
    def rand(*args):
        if args:
            try:
                O = my_matrix()
                for r in range(args[0]):
                    for c in range(args[1]):
                        O.__data[(r,c)] = mp.rand()
                O.__N_rows = args[0]
                O.__N_cols = args[1]
                return O
            except:
                raise ValueError("Wrong matrix dimensions!")
        else:
            raise ValueError("Wrong matrix dimensions!")
        
    def to_rand(self):
        M = my_matrix(self.__N_rows, self.__N_cols)
        for r in self.rows():
            for c in self.cols():
                M[r,c] = mp.rand()
        return M
    
    @staticmethod    
    def crand(*args):
        i = mp.mpc("0", "1")
        if args:
            try:
                O = my_matrix()
                for r in range(args[0]):
                    for c in range(args[1]):
                        O.__data[(r,c)] = mp.rand() * mp.exp(i * mp.mpf("2") * mp.pi * mp.rand())
                O.__N_rows = args[0]
                O.__N_cols = args[1]
                return O
            except:
                raise ValueError("Wrong matrix dimensions!")
        else:
            raise ValueError("Wrong matrix dimensions!")
        
    def to_crand(self):
        i = mp.mpc("0", "1")
        M = my_matrix(self.__N_rows, self.__N_cols)
        for r in self.rows():
            for c in self.cols():
                M[r,c] = mp.rand() * mp.exp(i * mp.mpf("2") * mp.pi * mp.rand())
        return M
    
    @staticmethod
    def set_use_tol(local_in):
        if local_in:
            my_matrix.__use_tol = True
        else:
            my_matrix.__use_tol = False
            
    @staticmethod
    def use_tol():
        return my_matrix.__use_tol
        
    @staticmethod
    def set_tol_overshoot(local_in):
        if isinstance(local_in, int) and local_in > 0: 
            my_matrix.__tol_overshoot = local_in
        else:
            raise ValueError("Wrong input to set_tol_overshot!")
        
    @staticmethod
    def tol_overshoot():
        return my_matrix.__tol_overshoot
    
    @staticmethod    
    def print_digits():
        return my_matrix.__print_digits
        
    @staticmethod    
    def set_print_digits(local_in):
        if isinstance(local_in, int) and local_in > 0: 
            my_matrix.__print_digits = local_in
        else:
            raise ValueError("Wrong input to set_print_digits!")
     
    @staticmethod
    def approximation_steps():
        return my_matrix.__max_approximation_steps
    
    @staticmethod
    def set_approximation_steps(value):
        my_matrix.__max_approximation_steps = value

    
        
        
            
            
        
    
    
    
    
    