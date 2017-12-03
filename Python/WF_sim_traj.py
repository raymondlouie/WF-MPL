from numpy import *
import numpy as np

# calculate Hamming distance between strings s1 and s2
def hamming_calculate(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def WF_sim_traj(s,mu,L,N,p_init,dt_array):

    
    K = 2**L
    T_tot = max(dt_array)

    # Generate  binary matrix
    U_list = []
    for index in range(K):
        bin_string=format(index, '0'+str(L)+'b') # string list
        bin_int = [int(numeric_string) for numeric_string in bin_string] # integer list
        U_list.append(bin_int) # append to list
        
    U = array(U_list)

    g = U@s # genotype selection coefficients
    w = g + ones((K,1))

    # mutation matrix
    Q= np.empty(shape=[K,K],dtype=float)
    for index1, seq1 in enumerate(U_list):
        for index2, seq2 in enumerate(U_list):
            hamming_dist = hamming_calculate(seq1,seq2)
            Q[index1,index2]=(mu**hamming_dist)*((1-mu)**(L-hamming_dist))
            
            
    Wmat = Q@np.diag(np.squeeze(w))
    
    # initialize first time point
    p  = p_init
    mut_mat = U.T@np.diag(p)@U
    
    double_mut  = zeros((T_tot,L,L))
    single_mut  = zeros((T_tot,L))
    
    double_mut[0,:,:] = mut_mat
    single_mut[0,:] = diag(mut_mat) 
    
    for ite_T in range(1,T_tot):
        
        p=p[:,np.newaxis] 
        f = Wmat@p
        f = f/sum(f)
        p = np.random.multinomial(N, list(f[:,0]))/N
        
        mut_mat = U.T@np.diag(p)@U
        double_mut[ite_T,:,:] = mut_mat
        single_mut[ite_T,:] = diag(mut_mat)
        

    double_mut = double_mut[dt_array-1,:,:]
    single_mut = single_mut[dt_array-1,:]
    
    return single_mut, double_mut