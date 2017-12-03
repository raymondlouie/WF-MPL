import scipy.linalg
from numpy import *
import numpy as np

def estimate_MPL(mu,dt_array,single_mut,double_mut):
 
    L =  single_mut.shape[1]
    T = len(dt_array)
    sum_cov = zeros((L,L))
    sum_b = zeros((1,L))
    
    for ite_T in range(1,T):
        dt = dt_array[ite_T] - dt_array[ite_T-1]
        
        prev_single_mut = single_mut[ite_T-1,:]
        prev_single_mut =prev_single_mut[:,np.newaxis] 

        mut_mat = double_mut[ite_T-1,:,:]
        
        cov_mat = mut_mat - prev_single_mut @ prev_single_mut.T # covariance matrix
        sum_cov = sum_cov + cov_mat*dt
        sum_b = sum_b + dt*(1 - 2*prev_single_mut.T) 


    reg_term = 0.1
    s_MPL = np.linalg.solve(sum_cov + reg_term*np.eye(L),(single_mut[T-1,:]-single_mut[0,:] - mu*sum_b).T)
    return s_MPL