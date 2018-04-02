import gmpy2 as gm
from gmpy2 import const_pi as pi, mpfr as mpr, mpc
from copy import deepcopy
from aux_inner_functions import is_lowest_level, is_dict_callable, is_list_callable
from numpy import double


class NUMERICAL_TYPES:
    @staticmethod
    def scalar_types():
        return (int, float, type(mpr(1)), type(mpc(1,1)), double)
    
def is_number(local_in):
    '''
        Must be a numeric scalalr or a numeric iterable.
    '''
    try:
        if is_scalar(local_in) or is_vector(local_in, NUMERICAL_TYPES.scalar_types):
            return True
        else:
            return False
    except:
        return False
    
def is_scalar(local_in):
    '''
        Must be a numeric scalar.
    '''
    if isinstance(local_in, NUMERICAL_TYPES.scalar_types):
        return True
    else:
        return False
    
def is_vector(local_in, *args):
    '''
        Must be an iterable of a given type. 
    '''
    if args:
        local_types = args[0]
    else:
        local_types = NUMERICAL_TYPES.scalar_types
    try:
        if test_iterable(local_in, lambda item: isinstance(item, local_types)):
            return True
        else:
            return False
    except:
        return False

def is_even(local_in):
    if isinstance(local_in, int) and round(local_in/2.0) == local_in/2.0:
        return True
    return False

def is_odd(local_in):
    if isinstance(local_in, int) and round(local_in/2.0) != local_in/2.0:
        return True
    return False

#def mp_matrix_hor_cat(M1, M2):
#    '''
#        Horizontally concatenates two matrices
#    '''
#    if not isinstance(M1, mp.matrix) or not isinstance(M2, mp.matrix):
#        raise ValueError("Both inputs must be mpmath.matrix!")
#    if M1.rows != M2.rows:
#        raise ValueError("Both input matrices must have the same amount of rows!")
#    aux = mp.zeros(M1.rows, M1.cols + M2.cols)
#    for n_rows in range(M1.rows):
#        for n_cols in range(M1.cols + M2.cols):
#            if n_cols < M1.cols:
#                aux[n_rows, n_cols] = M1[n_rows, n_cols]
#            else:
#                aux[n_rows, n_cols] = M2[n_rows, n_cols - M1.cols]
#    return aux   

#def mp_matrix_vert_cat(M1, M2):
#    '''
#        Vertically concatenates two matrices
#    '''
#    if not isinstance(M1, mp.matrix) or not isinstance(M2, mp.matrix):
#        raise ValueError("Both inputs must be mpmath.matrix!")
#    if M1.cols != M2.cols:
#        raise ValueError("Both input matrices must have the same amount of columns!")
#    aux = mp.zeros(M1.rows, M1.cols + M2.cols)
#    for n_rows in range(M1.rows + M2.rows):
#        for n_cols in range(M1.cols):
#            if n_rows < M1.rows:
#                aux[n_rows, n_cols] = M1[n_rows, n_cols]
#            else:
#                aux[n_rows, n_cols] = M2[n_rows - M1.rows, n_cols]
#    return aux

def test_iterable(local_iter, *args):
    '''
        Recursively goes over multi-level iterables, applying a test function to each lowest-level non-iterable object.
        The test function must have boolean return type.
        Returns True only if the test function is True for each lowest-level object    
    '''
    def local_apply(local_iter, local_function, lowest_level):
        if lowest_level(local_iter):
            try:
                if local_function(local_iter):
                    return True
                else:
                    return False
            except:
                raise RuntimeError("The given function can't be applied to the lowest-level variable!")
        else:
            try:
                if is_dict_callable(local_iter):
                    for key in local_iter.keys():
                        if lowest_level(local_iter[key]):
                            if not local_function(local_iter[key]):
                                return False
                        else:
                            if not local_apply(local_iter[key], local_function, lowest_level):
                                return False
                    return True
                elif is_list_callable(local_iter):
                    for key in range(len(local_iter)):
                        if lowest_level(local_iter[key]):
                            if not local_function(local_iter[key]):
                                return False
                        else:
                            if not local_apply(local_iter[key], local_function, lowest_level):
                                return False
                    return True                  
                else:
                    raise RuntimeError("Unspecified error happened, We shouldn't be here!")        
            except Exception as E:
                raise E
    if args:
        local_function = args[0]
        if len(args) > 1:
            string_is_iterable = args[1]
        else:
            string_is_iterable = False
    else:
        local_function = lambda x: x
        string_is_iterable = False
    try:
        return local_apply(local_iter, local_function, lambda x: is_lowest_level(x, string_is_iterable, str))
    except Exception as E:
        raise E
    
def unwind_iterable(local_iter, *args):
    '''
        Recursively goes over multi-level iterables, applying a given function to each lowest-level non-iterable object.
        The function must accept a single argument and return a single output.
        The results are collected into a single-level list
    '''
    def local_apply(local_iter, local_function, lowest_level, local_list):
        if lowest_level(local_iter):
            try:
                return local_function(local_iter)
            except:
                raise RuntimeError("The given function can't be applied to the lowest-level variable!")
        else:
            try:
                if is_dict_callable(local_iter):
                    for key in local_iter.keys():
                        if lowest_level(local_iter[key]):
                            local_list.append(local_function(local_iter[key]))
                        else:
                            local_apply(local_iter[key], local_function, lowest_level, local_list)
                elif is_list_callable(local_iter):
                    for key in range(len(local_iter)):
                        if lowest_level(local_iter[key]):
                            local_list.append(local_function(local_iter[key]))
                        else:
                            local_apply(local_iter[key], local_function, lowest_level, local_list)
                else:
                    raise RuntimeError("Unspecified error happened, We shouldn't be here!")        
            except Exception as E:
                raise E
    if args:
        local_function = args[0]
        if len(args) > 1:
            string_is_iterable = args[1]
        else:
            string_is_iterable = False
    else:
        local_function = lambda x: x
        string_is_iterable = False
    result = []
    try:
        local_apply(local_iter, local_function, lambda x: is_lowest_level(x, string_is_iterable, str), result)
        return result
    except Exception as E:
        raise E

def apply_to_iterable(local_iter, *args):
    '''
        Recursively goes over multi-level iterables, applying a given function to each lowest-level non-iterable object.
        The function must accept a single argument and return a single output of the same type as the input.
        The results are collected into a structure identical to the input one, which is created by a deep copying.
    '''
    def local_apply(local_iter, local_function, lowest_level):
        if lowest_level(local_iter):
            try:
                return local_function(local_iter)
            except:
                raise RuntimeError("The given function can't be applied to the lowest-level variable!")
        else:
            try:
                if is_dict_callable(local_iter):
                    for key in local_iter.keys():
                        if lowest_level(local_iter[key]):
                            local_iter[key] = local_function(local_iter[key])
                        else:
                            local_apply(local_iter[key], local_function, lowest_level)
                elif is_list_callable(local_iter):
                    for key in range(len(local_iter)):
                        if lowest_level(local_iter[key]):
                            local_iter[key] = local_function(local_iter[key])
                        else:
                            local_apply(local_iter[key], local_function, lowest_level)
                else:
                    raise RuntimeError("Unspecified error happened, We shouldn't be here!")        
            except Exception as E:
                raise E
    if args:
        local_function = args[0]
        if len(args) > 1:
            string_is_iterable = args[1]
        else:
            string_is_iterable = False
    else:
        local_function = lambda x: x
        string_is_iterable = False
    result = deepcopy(local_iter)
    try:
        local_apply(result, local_function, lambda x: is_lowest_level(x, string_is_iterable, str))
        return result
    except Exception as E:
        raise E

#def convert_to_mp_matrix(local_in):
#    '''
#        Tries to convert the input into mpmath.matrix.
#        Scalar inputs or inputs which can be converted into numpy.ndarray are accepted.
#        0D, 1D or 2D inputs are converted into mpmath.matrix.
#        3D inputs are converted into a list of mpmath.matrix.
#        4D and larger dimensional iterables are rejected as illegal inputs.  
#    '''
#    if not is_number(local_in):
#        raise ValueError("The input is not a numeric type!")
#    if is_scalar(local_in):
#        return mp.matrix([local_in])
#    if isinstance(local_in, mp.matrix):
#        return local_in
#    try:
#        local_in = np.array(local_in)
#    except:
#        raise ValueError("The input can't be converted to a numeric matrix!") 
#    if np.ndim(local_in) < 3:
#        return mp.matrix(local_in)
#    elif np.ndim(local_in) == 3:
#        result = [mp.matrix(np.shape(local_in)[0], np.shape(local_in)[1])] * np.shape(local_in)[2]
#        for n in range(np.shape(local_in)[2]):
#            result[n] = mp.matrix(local_in[:,:,n])
#        return result
#    else:
#        raise ValueError("Conversion to mp-lists of np-arrays of more than 3 dimensions is not supported for now!")
        
#def mp_abs(local_in):
#    '''
#        Returns an absolute value of input if the input is convertible to mpmath.mpf, mpmath.mpc or mpmath.mpmatrix.
#    '''
#    if is_scalar(local_in):
#        return mp.fabs(local_in)
#    else:
#        if isinstance(local_in, mp.matrix):
#            return local_in.apply(mp.fabs);
#        else:
#            try:
#                return convert_to_mp_matrix(local_in).apply(mp.fabs)
#            except:
#                raise ValueError("The per-element absolute value of the input can't be calculated!")

def z0_from_args(N, *args):
    raise NotImplementedError
#    if not args:
#        return gm.sqrt(mpfr("50"))*mp.diag(mp.ones(1,N))
#    elif is_scalar(args[0]):
#        return mp.sqrt(mp.mpf(args[0]))*mp.diag(mp.ones(1,N))
#    elif isinstance(args[0], mp.matrix) and ((args[0].rows == N and args[0].cols == 1) or (args[0].cols == N and args[0].rows == 1)):
#        aux = mp.matrix(1, N)
#        for n in range(N):
#            aux[0,n] = mp.sqrt(mp.mpf(args[0][n]))
#        return mp.diag(aux)
#    else:
#        raise ValueError("The Z0 argument can't be converted into a proper y0 vector!")
    
def y0_from_args(N, *args):
    raise NotImplementedError
#    if not args:
#        return mp.sqrt(mp.mpf("0.02"))*mp.diag(mp.ones(1,N))
#    elif is_scalar(args[0]):
#        return mp.sqrt(1/mp.mpf(args[0]))*mp.diag(mp.ones(1,N))
#    elif isinstance(args[0], mp.matrix) and ((args[0].rows == N and args[0].cols == 1) or (args[0].cols == N and args[0].rows == 1)):
#        aux = mp.matrix(1, N)
#        for n in range(N):
#            aux[0,n] = mp.sqrt(1/mp.mpf(args[0][n]))
#        return mp.diag(aux)
#    else:
#        raise ValueError("The Z0 argument can't be converted into a proper y0 vector!") 

#def is_square_mp_matrix(local_in):
#    if not isinstance(local_in, mp.matrix):
#        return False
#    if local_in.rows != local_in.cols:
#        return False
#    return True

def S_to_Z(local_S, *args):
    raise NotImplementedError
    '''
        Converts an S-parameter to Z-parameter.
        The matrixes must be in mpmath.matrix.
        The optional argument is a vector of Z0, if not given 50Ohm is assumed.
    '''
#    if not is_square_mp_matrix(local_S):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    z0 = z0_from_args(local_S.rows, *args)
#    return S_to_Z_simple(local_S, z0)

def S_to_Z_simple(local_S, z0):
    raise NotImplementedError
    '''
        Converts an S-parameter to Z-parameter.
        The matrixes must be in mpmath.matrix.
        The second argument is a diagonal matrix of sqrt(Z0).
    '''
#    I = mp.diag(mp.ones(1,local_S.rows))
#    return z0*(I + local_S)*((I - local_S)**-1)*z0 

def Z_to_S(local_Z, *args):
    raise NotImplementedError
    '''
        Converts an Z-parameter to S-parameter.
        The matrixes must be in mpmath.matrix.
        The optional argument is a vector of Z0, if not given 50Ohm is assumed.
    '''
#    if not is_square_mp_matrix(local_Z):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    y0 = y0_from_args(local_Z.rows, *args)
#    return Z_to_S_simple(local_Z, y0)

def Z_to_S_simple(local_Z, y0):
    raise NotImplementedError
    '''
        Converts an Z-parameter to S-parameter.
        The matrixes must be in mpmath.matrix.
        The second argument is a diagonal matrix of 1/sqrt(Z0).
    '''
#    I = mp.diag(mp.ones(1,local_Z.rows))
#    return (y0*local_Z*y0 - I)*((y0*local_Z*y0 + I)**-1)

def S_to_Y(local_S, *args):
    raise NotImplementedError
    '''
        Converts an S-parameter to Y-parameter.
        The matrixes must be in mpmath.matrix.
        The optional argument is a vector of Z0, if not given 50Ohm is assumed.
    '''
#    if not is_square_mp_matrix(local_S):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    y0 = y0_from_args(local_S.rows, *args)
#    return S_to_Y_simple(local_S, y0)
    
def S_to_Y_simple(local_S, y0):
    raise NotImplementedError
    '''
        Converts an S-parameter to Y-parameter.
        The matrixes must be in mpmath.matrix.
        The second argument is a diagonal matrix of 1/sqrt(Z0).
#    '''
#    I = mp.diag(mp.ones(1,local_S.rows))
#    return y0*(I - local_S)*((I + local_S)**-1)*y0

def Y_to_S(local_Y, *args):
    raise NotImplementedError
    '''
        Converts an Y-parameter to S-parameter.
        The matrixes must be in mpmath.matrix.
        The optional argument is a vector of Z0, if not given 50Ohm is assumed.
    '''
#    if not is_square_mp_matrix(local_Y):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    z0 = z0_from_args(local_Y.rows, *args)
#    return Y_to_S_simple(local_Y, z0)

def Y_to_S_simple(local_Y, z0):
    raise NotImplementedError
    '''
        Converts an Y-parameter to S-parameter.
        The matrixes must be in mpmath.matrix.
        The second argument is a diagonal matrix of sqrt(Z0).
    '''
#    I = mp.diag(mp.ones(1,local_Y.rows))
#    return (I - z0*local_Y*z0)*((I + z0*local_Y*z0)**-1)
            
def S_to_T(local_S):
    raise NotImplementedError
    '''
        Converts an S-parameter to T-parameter.
        The matrixes must be in mpmath.matrix.
        The input must have an even amount of ports and
        be organized so that ports 1:N/2 are on one side,
        while ports N+1:2N are on the other:
        1 <--> 3
        2 <--> 4
    '''
#    if not is_square_mp_matrix(local_S):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    if not is_even(local_S.rows):
#        raise ValueError("An S-parameter with odd number of ports can't be converted into a T-parameter")
#    N = int(local_S.rows/2)
#    S11 = local_S[0:N-1, 0:N-1]
#    S12 = local_S[0:N-1, N:-1]
#    S21 = local_S[N:-1, 0:N-1]
#    S22 = local_S[N:-1, N:-1]
#    T11 = -(S12**-1)*S11
#    T12 = S12**-1
#    T21 = S21 - S22*(S12**-1)*S11
#    T22 = S22*(S12**-1)
#    return mp_matrix_vert_cat(mp_matrix_hor_cat(T11, T12), mp_matrix_hor_cat(T21, T22))

def T_to_S(local_T):
    raise NotImplementedError
    '''
        Converts a T-parameter to S-parameter.
        The matrixes must be in mpmath.matrix.
        The input must have an even amount of ports and
        be organized so that in S-form the ports 1:N/2 are on one side,
        while ports N+1:2N are on the other:
        a1 --> b2
        b1 <-- a2
    '''
#    if not is_square_mp_matrix(local_T):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    if not is_even(local_T.rows):
#        raise ValueError("A T-parameter with odd number of ports shouldn't exist!")
#    N = int(local_T.rows/2)
#    T11 = local_T[0:N-1, 0:N-1]
#    T12 = local_T[0:N-1, N:-1]
#    T21 = local_T[N:-1, 0:N-1]
#    T22 = local_T[N:-1, N:-1]
#    S11 = -(T22**-1)*T21
#    S12 = T22**-1
#    S21 = T11-T12*(T22**-1)*T21
#    S22 = T12*(T22**-1)
#    return mp_matrix_vert_cat(mp_matrix_hor_cat(S11, S12), mp_matrix_hor_cat(S21, S22))

def Z_to_H(local_Z):
    raise NotImplementedError
    '''
        Converts a Z-parameter to an H-parameter
        The matrixes must be in mpmath.matrix.
        The input must have an even amount of ports and
        be organized so that ports 1:N/2 are on one side,
        while ports N+1:2N are on the other:
        1 <--> 3
        2 <--> 4
    '''
#    if not is_square_mp_matrix(local_Z):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    if not is_even(local_Z.rows):
#        raise ValueError("A Z-parameter with odd number of ports can't be converted into an H-parameter")
#    N = int(local_Z.rows/2)
#    Z11 = local_Z[0:N-1, 0:N-1]
#    Z12 = local_Z[0:N-1, N:-1]
#    Z21 = local_Z[N:-1, 0:N-1]
#    Z22 = local_Z[N:-1, N:-1]
#    H11 = Z22*(Z12**-1)
#    H12 = Z21 -Z22*(Z12**-1)*Z11
#    H21 = -Z12**-1
#    H22 = (Z12**-1)*Z11
#    return mp_matrix_vert_cat(mp_matrix_hor_cat(H11, H12), mp_matrix_hor_cat(H21, H22))

def H_to_Z(local_H):
    raise NotImplementedError
    '''
        Converts an H-parameter to a Z-parameter
        The matrixes must be in mpmath.matrix.
        The input must have an even amount of ports and
        be organized so that ports 1:N/2 are on one side,
        while ports N+1:2N are on the other:
        1 <--> 3
        2 <--> 4
    '''
#    if not is_square_mp_matrix(local_H):
#        raise ValueError("The first input must be a square mpmath.matrix!")
#    if not is_even(local_H.rows):
#        raise ValueError("An H-parameter with odd number of ports can't be converted into a Z-parameter")
#    N = int(local_H.rows/2)
#    H11 = local_H[0:N-1, 0:N-1]
#    H12 = local_H[0:N-1, N:-1]
#    H21 = local_H[N:-1, 0:N-1]
#    H22 = local_H[N:-1, N:-1]
#    Z11 = -(H21**-1)*H22
#    Z12 = - (H21**-1)
#    Z21 = H12 - H11*(H21**-1)*H22
#    Z22 = H11*(H21**-1)
#    return mp_matrix_vert_cat(mp_matrix_hor_cat(Z11, Z12), mp_matrix_hor_cat(Z21, Z22))

def S_to_H(local_S, *args):
    return Z_to_H(S_to_Z(local_S, args))

def H_to_S(local_H, *args):
    return Z_to_S(H_to_Z(local_H), args)  

def deg(local_in):
    '''
        Converts an input from radianes to degrees
    '''
    if is_scalar(local_in):
        return local_in / pi() * mpr(180)
    else:
        return apply_to_iterable(local_in, lambda x: x / pi() * mpr(180))

def rad(local_in):
    '''
        Converts an input from degrees to radianes
    '''
    if is_scalar(local_in):
        return local_in * pi() / mpr(180)
    else:
        raise NotImplemented()

def dB20(local_in):
    '''
        Converts an input absolute value into dB20 form
    '''
    if is_scalar(local_in):
        return mpr(20) * gm.log10(abs(local_in))
    else:
        raise NotImplemented()

def dB10(local_in):
    '''
        Converts an input absolute value into dB10 form
    '''
    if is_scalar(local_in):
        return mpr(10) * gm.log10(abs(local_in))
    else:
        raise NotImplemented()
    
def get_dps_by_tol(tol):
    raise NotImplementedError
#    current_precision = mp.dps
#    if tol <= mp.mpf("10")**mp.mpf(-current_precision+1):
#        return int(mp.ceil(-mp.log10(tol)) + 1)
#    else:
#        return current_precision
    
def make_N_dim_list(*args):
    if not args:
        return []
    try:
        if len(args[0]) == 1:
            return [None]*args[0][0]
        else: 
            local_list = [make_N_dim_list(args[0][1:len(args[0])])]*args[0][0]
            return local_list
    except Exception as E:
        raise E
    
def mpr_to_str(num, symbols):
    num = num.real
    if symbols <= 3:
        raise ValueError("The value of mp.mpf can't be reproduced with less than 3 symbols!")
    if num >= 0:
        num_sign = " "
    else:
        num_sign = "-"
    if num ==mpr(0):
        return num_sign + "0." + "0"*(symbols-3)
    num_char = gm.floor(gm.log10(abs(num)))
    num_val = abs(num) / (mpr(10)**num_char)
    num_char_str = str(int(num_char))
    symbols_left = symbols - 4 - len(num_char_str)
    val_string = ("{0:-1.%dNf}" % (symbols_left+1)).format(num_val)
    if val_string[1] != ".":
        num_val /= mpr(10)
        num_char += 1
        num_char_str = str(int(num_char))
        val_string = ("{0:-1.%dNf}" % (symbols_left+1)).format(num_val)
    if len(val_string) + 2 + len(num_char_str) < symbols:
        val_string += "0"
    elif len(val_string) + 2 + len(num_char_str) > symbols:
        val_string = val_string[:-1]
    if symbols_left < 1:
        raise ValueError("Too few given symbols for the given value of mp.mpf!")
    return num_sign + val_string + "e" + num_char_str

def mpi_to_str(num, symbols):
    local_str = mpr_to_str(num.real, symbols)
    return local_str[0]+"j"+local_str[1:]

def mpc_to_str(num, symbols):
    if is_even(symbols) or symbols < 7:
        raise ValueError("mpc_to_str requires an odd amount of symbols more than 5 to operate!")
    abs_re_str = mpr_to_str(abs(num.real), int((symbols-5)/2+1))[1:]
    abs_im_str = mpr_to_str(abs(num.imag), int((symbols-5)/2+1))[1:]
    if num.real >= mpr(0):
        num_re_sign = " "
    else:
        num_re_sign = "-"
    if num.imag >= mpr(0):
        num_im_sign = "+"
    else:
        num_im_sign = "-"
    return num_re_sign + abs_re_str + " " + num_im_sign + " j" + abs_im_str

def get_index(ind, length):
    if isinstance(ind, int):
        return [ind % length]
    if isinstance(ind, slice):
        if ind.start != None:
            if length > 1:
                start_val = ind.start % length
            else:
                start_val = ind.start
        else:
            start_val = 0
        if ind.stop != None:
            if length > 1:
                stop_val = ind.stop % length
            else:
                stop_val = ind.stop
        else:
            stop_val = length
        if ind.step != None:
            step_val = ind.step
        else:
            step_val = 1
        local_iter = range(start_val, stop_val, step_val)
    if isinstance(ind, (range, list)):
        local_iter = ind
    aux = []
    for n in local_iter:
        aux.append(n % length)
    return aux

def real(num):
    return num.real()

def imag(num):
    return num.imag()

def conj(num):
    return num.real - mpc(0,1)*num.imag
    
class color_range:
    
    __current_index = 0
    
    @classmethod
    def reset(cls):
        cls.__current_index = 0
    @classmethod
    def make_color(cls):
        if cls.__current_index == 0:
            result = (1, 0, 0)
        elif cls.__current_index == 1:
            result = (0, 0, 1)
        elif cls.__current_index == 2:
            result = (0, 1, 0)
        elif cls.__current_index == 3:
            result = (0, 0, 0)
        elif cls.__current_index == 4:
            result = (1, 0, 1)
        elif cls.__current_index == 5:
            result = (1, 0.8, 0)
        elif cls.__current_index == 6:
            result = (0, 1, 1)
        elif cls.__current_index == 7:
            result = (0.5, 0.5, 0.5)
        elif cls.__current_index == 8:
            result = (1, 0.5, 0)
        elif cls.__current_index == 9:
            result = (0, 0.5, 1)
        elif cls.__current_index == 10:
            result = (0.5, 1, 0)
        elif cls.__current_index == 11:
            result = (0.3, 0, 1)
        elif cls.__current_index == 12:
            result = (0, 1, 0.3)
        elif cls.__current_index == 13:
            result = (1, 0, 0.3)
        else:
            result = (1, 0, 0.3)
        cls.__current_index += 1
        _, cls.__current_index = divmod(cls.__current_index, 14)
        return result
        
         
    
    
