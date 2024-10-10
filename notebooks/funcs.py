import numpy as np
from matplotlib import pyplot as plt
import scipy as sc
from sklearn.preprocessing import PolynomialFeatures
from scipy.sparse import coo_matrix, block_diag


def deriv1(X, dt):
    '''
    Calculates the derivative of column-matrix X with respect to t. It uses second-order accurate derivative
    (central derivative) for internal values and a first order derivative for boundaries, keeping the original dimensions.
    
    Note that this function differentiates for n-vertically-concatenated experiments, since it uses time.
    '''
    
#     t = t.reshape((len(t), 1))             # Avoid the use of (n,)
      
    # Initialize the return array.
    dX = np.zeros_like(X)                 # Array of the shape of X, in this case <--> (N, dims) 
    
     # --> Second-order accurate derivative (central derivative)
    dX[1:-1, :] = (X[2:, :] - X[:-2, :]) / 2  # f(3)-f(1)/(t(3)-t(1)). No 2* needed. Remember i[:-1] doesnt
                                                             # include the last element, includes the second to last
    # --> Treats the boundary points with a first-order derivation
    dX[0, :]    = (X[1, :] - X[0, :])
    dX[-1, :]   = (X[-1, :] - X[-2, :])
    
    return dX/dt

def deriv1_time(X, t):
    '''
    Calculates the derivative for multiple experiments of the same duration n. Segments each experiment and calculates its derivative independently.

    If an error occurs when calculating n, falls back to using deriv4 directly on the data.

    Parameters:
    X (numpy array): 2D numpy array of shape (n_experiments, n_timepoints) containing the data to be differentiated
    t (numpy array): 1D numpy array of time corresponding to the data

    Returns:
    final_dX (numpy array): 2D numpy array of shape (n_experiments, n_timepoints) containing the differentiated data
    '''
    try:
        dt = t[1] - t[0]
        n = np.where(t == t[0])[0][1]
    except Exception as e:
#         print(f"Error: {e}")
#         print("Falling back to deriv4...")
        return deriv1(X, dt)

    final_dX = np.zeros_like(X)
    for i in range(len(np.where(t == t[0])[0])):
        final_dX[i*n:(i+1)*n, :] = deriv1(X[i*n:(i+1)*n, :], dt)

    return final_dX

def deriv4(X, dt):
    # Initialize the return array.
    dX = np.zeros_like(X)                 

    # 4th order forward
    dX[:2, :] = (-(25/12)*X[:2, :] + 4*X[1:3, :] - 3*X[2:4, :] + (4/3)*X[3:5, :]-(1/4)*X[4:6, :])    
    # 4th order central
    dX[2:-2, :] = (-(1/12)*X[4:, :] + (2/3)*X[3:-1, :] - (2/3)*X[1:-3, :] + (1/12)*X[:-4, :])
    # 4th order backward
    dX[-2:, :] = ((1/4)*X[-6:-4, :] - (4/3)*X[-5:-3, :] + 3*X[-4:-2, :] - 4*X[-3:-1, :] + (25/12)*X[-2:, :])

    return dX/dt

def trans(x, n, s, N):
    '''
    Extracts initial data points UP to data "n" of each experiment from a full array,
    where each experimental data splits at s data points, for a total of N experiments
    
    Input:    1d vector
    Returns:  1d vector
    '''
    final = np.zeros((n*N, x.shape[1]))
    for i in range(N):
        final[i*n:(i+1)*n, :] = x[i*s:i*s+n, :]
    return final

def est(x, n, s, N):
    '''
    Extracts final data points starting from data "n" of each experiment from a full array,
    where each experimental data splits at s data points, for a total of N experiments
    
    Input:    1d vector
    Returns:  1d vector
    '''
    final = np.zeros(((s-n)*N, x.shape[1]))
    for i in range(N):
        final[i*(s-n):(i+1)*(s-n), :] = x[n+i*s:(i+1)*s, :]
    return final




def deriv4_time(X, t):
    '''
    Calculates the derivative for multiple experiments of the same duration n. Segments each experiment and calculates its derivative independently.

    If an error occurs when calculating n, falls back to using deriv4 directly on the data.

    Parameters:
    X (numpy array): 2D numpy array of shape (n_experiments, n_timepoints) containing the data to be differentiated
    t (numpy array): 1D numpy array of time corresponding to the data

    Returns:
    final_dX (numpy array): 2D numpy array of shape (n_experiments, n_timepoints) containing the differentiated data
    '''
    try:
        dt = t[1] - t[0]
        n = np.where(t == t[0])[0][1]
    except Exception as e:
#         print(f"Error: {e}")
#         print("Falling back to deriv4...")
        return deriv4(X, dt)

    final_dX = np.zeros_like(X)
    for i in range(len(np.where(t == t[0])[0])):
        final_dX[i*n:(i+1)*n, :] = deriv4(X[i*n:(i+1)*n, :], dt)

    return final_dX

    
def BurstSampling(arr, t, count, step):
    '''
    arr: matriz de datos
    t = vector de tiempos
    count = muestreo consecutivo
    step = distancia entre muestreos
    dt = delta t
    
    returns: dx0, CENTRAL_PTS, CENTRAL_T, x0, t0
    '''
    dt = t[1]-t[0]
    # --------------------------------------------------------------------------
    def separate_matrices(arr, t):
        last_column = arr[:, -1]
        change_indices = list(np.where(np.diff(last_column) != 0)[0] + 1)
        change_indices.insert(0, 0)
        change_indices.append(arr.shape[0])

        x_result = []
        t_result = []
        for i in range(len(change_indices)-1):
            x_result.append(arr[change_indices[i]:change_indices[i+1], :])
            t_result.append(t[change_indices[i]:change_indices[i+1]])
        return x_result, t_result
    # --------------------------------------------------------------------------
    def generate_consecutive_list(count, step, max_length):
        lista = []
        n_iterations = max_length//step
        for i in range(n_iterations):
            for j in range(count):
                lista.append(i*step + j)
        return lista
    # --------------------------------------------------------------------------    
    def dx_maker_4th(X, dt):
        m, n = X.shape
        dx = np.zeros((m//5, n))
        CENTRAL_PTS = []
        CENTRAL_T = []
        idx = 0
        for i in range(0, m, 5):
            dx[idx, :] = ((1/12)*X[i, :] - (2/3)*X[i+1, :] + (2/3)*X[i+3, :] - (1/12)*X[i+4, :])/dt
            CENTRAL_PTS.append(X[i+2, :])
            CENTRAL_T.append(t0[i+2])
            idx += 1
        return dx, CENTRAL_PTS, CENTRAL_T
    # --------------------------------------------------------------------------
    matrices = separate_matrices(arr, t)  # splits matrices by parameter
    indices = []
    
    for i in matrices[0]:
        indices.append(generate_consecutive_list(count, step, len(i)))
        
    x0 = matrices[0][0][indices[0]]
    t0 = matrices[1][0][indices[0]]
    dx0 = dx_maker_4th(x0, dt)[0]
    CENTRAL_PTS = dx_maker_4th(x0, dt)[1]
    CENTRAL_T = dx_maker_4th(x0, dt)[2]
    
    for i in range(1, len(matrices[0])):
        x0 = np.r_[x0, matrices[0][i][indices[i]]]  # consective points
        t0 = np.r_[t0, matrices[1][i][indices[i]]]
        dx0 = np.r_[dx0, dx_maker_4th(matrices[0][i][indices[i]], dt)[0]]
        CENTRAL_PTS = np.append(CENTRAL_PTS, dx_maker_4th(matrices[0][i][indices[i]], dt)[1], axis=0)
        CENTRAL_T = np.append(CENTRAL_T, dx_maker_4th(matrices[0][i][indices[i]], dt)[2], axis=0)
    
    return dx0, np.array(CENTRAL_PTS), np.array(CENTRAL_T), x0, t0

def mysindy2(X, t, lamb, polyorder, normalize=False):
    '''
    Applies a least-squared-sequential-threshholded SINDy algorithm to a matrix of spatial derivatives
    
    Input:
    
    X: Concatenated column-vector-matrix of each spatial data vector. Array type
    t: Column vector with the time of data acquisition. 
    lamb: Sparsity knob
    polyorder: order of the maximum polynomial for the library.
    
    Output:
    
    Column-vector-matrix with each constant of the respective polynomial term
    '''
    
    dims = X.shape[1]  # Number of equations
    vals = X.shape[0]     # Number of elements in each dim
    
    # --> Initialize the return array.
    dX = deriv4_time(X, t)
    
    THETA = PolynomialFeatures(polyorder)
    variable_names = THETA.fit(X).get_feature_names_out(['x', 'y', 'u'][0:dims])
    THETA = THETA.fit_transform(X)
    
    if normalize == True:
    
        #NORMALIZE COLUMNS
        THETAp = np.zeros_like(THETA)
        norms = np.zeros(THETA.shape[1])
        for i in range(THETA.shape[1]):
            norms[i] = np.max(THETA[:, i])
            THETAp[:, i] = THETA[:, i]/np.max(THETA[:, i])
        THETA = THETAp

        # Initial guess
        CHI = np.linalg.lstsq(THETA, dX, rcond=None)[0]

        # Seq. Threshholding    
        for k in range(19):
            for i in range(dims):                        # n is state dimension
                smallinds = np.abs(CHI[:,i]) < np.mean(np.abs(CHI[:,i]))*lamb   # Find small coefficients. This returns a boolean array
                CHI[smallinds, i] = 0                             # add threshold. CHI must be an array for this to work
                biginds = smallinds == 0           # == False
                # Regress dynamics onto remaining terms to find sparse Xi
                CHI[biginds, i] = np.linalg.lstsq(THETA[:, biginds], dX[:, i], rcond=None)[0]


        for i in range(dims):
            CHI[:, i] = CHI[:, i]/norms
    
    else:
        # Initial guess
        CHI = np.linalg.lstsq(THETA, dX, rcond=None)[0]

        # Seq. Threshholding    
        for k in range(19):
            for i in range(dims):                        # n is state dimension
                smallinds = np.abs(CHI[:,i]) < np.max(np.abs(CHI[:,i]))*lamb   # Find small coefficients. This returns a boolean array
                CHI[smallinds, i] = 0                             # add threshold. CHI must be an array for this to work
                biginds = smallinds == 0           # == False
                # Regress dynamics onto remaining terms to find sparse Xi
                CHI[biginds, i] = np.linalg.lstsq(THETA[:, biginds], dX[:, i], rcond=None)[0]
    
    return CHI, variable_names, lamb, THETA